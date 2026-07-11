#!/usr/bin/env python3
"""Summarize reverse-polarity BFBArchitect calls from existing AC outputs."""

from __future__ import annotations

import argparse
import base64
import html
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
CENTROMERE_BEDS = {
    "GRCh37": REPO_ROOT / "bfbarchitect" / "resources" / "GRCh37_centromere.bed",
    "GRCh38": REPO_ROOT / "bfbarchitect" / "resources" / "GRCh38_centromere.bed",
}


@dataclass(frozen=True)
class Cohort:
    name: str
    root: Path
    ref: str | None = None


REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>[0-9]+)-(?P<end>[0-9]+):(?P<score>[^:]+):(?P<multiplicity>.+)$")
SCORE_CHECK_COLUMNS = ["sample_name", "amplicon_number", "mode", "run_region", "multiplicity", "score"]
EXAMPLE_COLUMNS = [
    "rank",
    "cohort",
    "ref",
    "sample_name",
    "amplicon_number",
    "score",
    "mode",
    "native_region",
    "native_arm",
    "BFB_source",
    "aa_plot",
    "bfb_plot",
    "aa_asset",
    "bfb_asset",
]


def clean_bool(value: object) -> str:
    text = str(value)
    if text == "True":
        return "True"
    if text == "False":
        return "False"
    return "NA"


def read_delimited_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
    df.columns = [str(col).strip() for col in df.columns]
    return df


def optional_text(value: object) -> str | None:
    text = str(value).strip()
    if text in {"", "NA", "nan", "None"}:
        return None
    return text


def read_manifest(path: Path) -> list[Cohort]:
    df = read_delimited_table(path)
    df = df.rename(
        columns={
            "name": "cohort",
            "reference": "ref",
            "ref_genome": "ref",
            "genome": "ref",
        }
    )
    missing = {"cohort", "root"} - set(df.columns)
    if missing:
        raise ValueError(f"manifest {path} is missing required columns: {', '.join(sorted(missing))}")

    cohorts = []
    for row_num, row in df.iterrows():
        name = optional_text(row["cohort"])
        root_text = optional_text(row["root"])
        if not name or not root_text:
            raise ValueError(f"manifest {path} has an empty cohort or root on row {row_num + 2}")
        root = Path(root_text).expanduser()
        if not root.exists():
            raise FileNotFoundError(f"cohort root does not exist for {name}: {root}")
        cohorts.append(Cohort(name=name, root=root, ref=optional_text(row.get("ref", ""))))
    if not cohorts:
        raise ValueError(f"manifest {path} did not contain any cohorts")
    return cohorts


def parse_ref(root: Path) -> str:
    for log_path in sorted(root.glob("*.log")):
        text = log_path.read_text(errors="replace")
        match = re.search(r"(?:^|\s)--ref\s+(\S+)", text)
        if match:
            return match.group(1)
    root_text = str(root)
    if "GRCh38" in root_text or "hg38" in root_text:
        return "GRCh38"
    if "GRCh37" in root_text or "hg19" in root_text:
        return "GRCh37"
    return "unknown"


def profile_path(root: Path) -> Path:
    paths = sorted(root.glob("*amplicon_classification_profiles.tsv"))
    if len(paths) != 1:
        raise FileNotFoundError(f"expected one profile table in {root}, found {len(paths)}")
    return paths[0]


def normalize_chrom(chrom: object) -> str:
    text = str(chrom).strip()
    if text.startswith("chr"):
        text = text[3:]
    return text


def load_centromeres(refs: set[str]) -> dict[tuple[str, str], tuple[int, int]]:
    out: dict[tuple[str, str], tuple[int, int]] = {}
    for ref in refs:
        bed = CENTROMERE_BEDS.get(ref)
        if not bed or not bed.exists():
            continue
        rows = []
        with bed.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end, *_ = line.rstrip("\n").split("\t")
                rows.append((normalize_chrom(chrom), int(start), int(end)))
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in rows:
            by_chrom.setdefault(chrom, []).append((start, end))
        for chrom, spans in by_chrom.items():
            out[(ref, chrom)] = (min(s for s, _ in spans), max(e for _, e in spans))
    return out


def parse_regions(value: object) -> list[dict[str, object]]:
    text = str(value)
    if text in {"", "NA", "nan", "None"}:
        return []
    regions = []
    for token in text.split(";"):
        match = REGION_RE.match(token)
        if not match:
            continue
        score_text = match.group("score")
        try:
            score = float(score_text)
        except ValueError:
            score = None
        regions.append(
            {
                "region": token,
                "chrom": normalize_chrom(match.group("chrom")),
                "start": int(match.group("start")),
                "end": int(match.group("end")),
                "score": score,
                "multiplicity": match.group("multiplicity"),
            }
        )
    return regions


def pick_region(row: pd.Series) -> dict[str, object] | None:
    regions = parse_regions(row.get("BFBArchitect_regions", "NA"))
    if not regions:
        return None
    min_score = pd.to_numeric(row.get("BFBArchitect_min_score"), errors="coerce")
    if pd.notna(min_score):
        for region in regions:
            if region["score"] is not None and abs(float(region["score"]) - float(min_score)) < 0.005:
                return region
    scored = [r for r in regions if r["score"] is not None]
    if scored:
        return min(scored, key=lambda r: float(r["score"]))
    return regions[0]


def chromosome_arm(ref: str, region: dict[str, object] | None, centromeres: dict[tuple[str, str], tuple[int, int]]) -> str:
    if not region:
        return "unknown"
    chrom = str(region["chrom"])
    start = int(region["start"])
    end = int(region["end"])
    cent = centromeres.get((ref, chrom))
    if not cent:
        return f"chr{chrom}:unknown-arm"
    cent_start, cent_end = cent
    if end <= cent_start:
        return f"chr{chrom}p"
    if start >= cent_end:
        return f"chr{chrom}q"
    return f"chr{chrom}:centromere-spanning"


def bfb_mode(row: pd.Series) -> str:
    whole = clean_bool(row.get("BFBArchitect_whole_graph_used"))
    reverse = clean_bool(row.get("BFBArchitect_reverse_polarity_used"))
    if whole == "NA" or reverse == "NA":
        return "BFBArchitect not scored"
    if whole == "True" and reverse == "True":
        return "whole graph + reverse_polarity"
    if whole == "True":
        return "whole graph"
    if reverse == "True":
        return "default + reverse_polarity"
    return "default"


def bfb_source_group(row: pd.Series) -> str:
    source = str(row.get("BFB_source", "NA"))
    if row.get("BFB+") != "Positive":
        return "No final BFB"
    if source == "AC":
        return "AC only"
    if source == "BFBArchitect":
        return "BFBArchitect only"
    if source == "AC|BFBArchitect":
        return "AC + BFBArchitect"
    return source


def cohort_ref(cohort: Cohort) -> str:
    return cohort.ref or parse_ref(cohort.root)


def read_calls(cohorts: list[Cohort]) -> tuple[pd.DataFrame, dict[tuple[str, str], tuple[int, int]]]:
    refs = {cohort_ref(cohort) for cohort in cohorts}
    centromeres = load_centromeres(refs)
    frames = []
    for cohort in cohorts:
        ref = cohort_ref(cohort)
        df = pd.read_csv(profile_path(cohort.root), sep="\t", dtype=str).fillna("NA")
        df["cohort"] = cohort.name
        df["cohort_root"] = str(cohort.root)
        df["ref"] = ref
        df["BFBArchitect_min_score_numeric"] = pd.to_numeric(df["BFBArchitect_min_score"].replace("NA", pd.NA), errors="coerce")
        df["BFBArchitect_whole_graph_used"] = df["BFBArchitect_whole_graph_used"].map(clean_bool)
        df["BFBArchitect_reverse_polarity_used"] = df["BFBArchitect_reverse_polarity_used"].map(clean_bool)
        df["bfbarchitect_mode"] = df.apply(bfb_mode, axis=1)
        df["bfb_source_group"] = df.apply(bfb_source_group, axis=1)
        df["bfbarchitect_final_positive"] = (df["BFB+"] == "Positive") & df["BFB_source"].str.contains("BFBArchitect", na=False)
        df["reverse_polarity_final_positive"] = df["bfbarchitect_final_positive"] & (
            df["BFBArchitect_reverse_polarity_used"] == "True"
        )
        df["bfbarchitect_scored"] = df["BFBArchitect_min_score_numeric"].notna()
        regions = df.apply(pick_region, axis=1)
        df["native_region"] = regions.map(lambda r: r["region"] if r else "NA")
        df["native_chrom"] = regions.map(lambda r: f"chr{r['chrom']}" if r else "NA")
        df["native_start"] = regions.map(lambda r: r["start"] if r else pd.NA)
        df["native_end"] = regions.map(lambda r: r["end"] if r else pd.NA)
        df["native_arm"] = regions.map(lambda r: chromosome_arm(ref, r, centromeres))
        frames.append(df)
    return pd.concat(frames, ignore_index=True), centromeres


def find_bfb_plot(row: pd.Series) -> Path | None:
    root = Path(row["cohort_root"])
    out_dir = root / "bfbarchitect_outputs"
    prefix = f"{row['sample_name']}_{row['amplicon_number']}_"
    candidates = sorted(out_dir.glob(prefix + "*_BFB_1.png"))
    if not candidates:
        return None
    mode = row["bfbarchitect_mode"]
    if mode.startswith("whole graph"):
        preferred = [p for p in candidates if "_whole_graph_" in p.name]
    else:
        preferred = [p for p in candidates if "_region" in p.name]
    return (preferred or candidates)[0]


def find_aa_plot(row: pd.Series) -> Path | None:
    root = Path(row["cohort_root"])
    base = root / "files" / f"{row['sample_name']}_{row['amplicon_number']}"
    for suffix in (".png", ".pdf"):
        path = base.with_suffix(suffix)
        if path.exists():
            return path
    return None


def copy_or_convert_image(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix.lower() == ".png":
        shutil.copy2(src, dst.with_suffix(".png"))
        return dst.with_suffix(".png")
    if src.suffix.lower() == ".pdf":
        out_base = dst.with_suffix("")
        subprocess.run(["pdftoppm", "-png", "-singlefile", "-r", "180", str(src), str(out_base)], check=True)
        return out_base.with_suffix(".png")
    raise ValueError(f"unsupported image type: {src}")


def make_summary_tables(calls: pd.DataFrame, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    calls = calls.copy()
    calls["polarity_assessment_category"] = "Not final BFB"
    final_bfb = calls["BFB+"] == "Positive"
    ac_only = final_bfb & (calls["BFB_source"] == "AC")
    assessed = final_bfb & calls["bfbarchitect_final_positive"]
    calls.loc[ac_only, "polarity_assessment_category"] = "AC-native only (polarity not assessed)"
    calls.loc[assessed, "polarity_assessment_category"] = calls.loc[assessed, "bfbarchitect_mode"]

    summary = (
        calls.groupby(["cohort", "ref"], dropna=False)
        .agg(
            amplicons=("sample_name", "size"),
            final_bfb_positive=("BFB+", lambda s: int((s == "Positive").sum())),
            bfbarchitect_assessed_final_bfb=("bfbarchitect_final_positive", "sum"),
            ac_native_only_not_assessed=("polarity_assessment_category", lambda s: int((s == "AC-native only (polarity not assessed)").sum())),
            reverse_polarity_final_positive=("reverse_polarity_final_positive", "sum"),
        )
        .reset_index()
    )
    summary["reverse_fraction_of_assessed_final_bfb"] = (
        summary["reverse_polarity_final_positive"] / summary["bfbarchitect_assessed_final_bfb"].replace(0, pd.NA)
    )
    total_assessed = int(summary["bfbarchitect_assessed_final_bfb"].sum())
    total_reverse = int(summary["reverse_polarity_final_positive"].sum())
    summary.loc["All cohorts"] = [
        "All cohorts",
        "mixed",
        int(summary["amplicons"].sum()),
        int(summary["final_bfb_positive"].sum()),
        total_assessed,
        int(summary["ac_native_only_not_assessed"].sum()),
        total_reverse,
        total_reverse / total_assessed if total_assessed else pd.NA,
    ]

    source = (
        calls[calls["BFB+"] == "Positive"]
        .groupby(["cohort", "ref", "bfb_source_group"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    source_all = source.groupby("bfb_source_group", as_index=False)["count"].sum()
    source_all.insert(0, "ref", "mixed")
    source_all.insert(0, "cohort", "All cohorts")
    source = pd.concat([source, source_all], ignore_index=True)

    mode = (
        calls[final_bfb]
        .groupby(["cohort", "ref", "polarity_assessment_category"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    mode_all = mode.groupby("polarity_assessment_category", as_index=False)["count"].sum()
    mode_all.insert(0, "ref", "mixed")
    mode_all.insert(0, "cohort", "All cohorts")
    mode = pd.concat([mode, mode_all], ignore_index=True)

    assessed_mode = (
        calls[assessed]
        .groupby(["cohort", "ref", "bfbarchitect_mode"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    assessed_mode_all = assessed_mode.groupby("bfbarchitect_mode", as_index=False)["count"].sum()
    assessed_mode_all.insert(0, "ref", "mixed")
    assessed_mode_all.insert(0, "cohort", "All cohorts")
    assessed_mode = pd.concat([assessed_mode, assessed_mode_all], ignore_index=True)

    summary.to_csv(out_dir / "reverse_polarity_usage_summary.tsv", sep="\t", index=False)
    source.to_csv(out_dir / "final_bfb_source_breakdown.tsv", sep="\t", index=False)
    assessed_mode.to_csv(out_dir / "bfbarchitect_assessed_final_mode_breakdown.tsv", sep="\t", index=False)
    assessed_mode.to_csv(out_dir / "bfbarchitect_final_mode_breakdown.tsv", sep="\t", index=False)
    mode.to_csv(out_dir / "final_bfb_polarity_assessment_breakdown.tsv", sep="\t", index=False)
    return summary, source, mode


def plot_mode_breakdown(mode: pd.DataFrame, out: Path) -> None:
    order = [
        "AC-native only (polarity not assessed)",
        "default",
        "default + reverse_polarity",
        "whole graph",
        "whole graph + reverse_polarity",
    ]
    colors = {
        "AC-native only (polarity not assessed)": "#9ca3af",
        "default": "#4e79a7",
        "default + reverse_polarity": "#f28e2b",
        "whole graph": "#59a14f",
        "whole graph + reverse_polarity": "#e15759",
    }
    cur = mode[mode["cohort"] != "All cohorts"]
    pivot = cur.pivot_table(index="cohort", columns="polarity_assessment_category", values="count", fill_value=0, aggfunc="sum")
    pivot = pivot.reindex(columns=order, fill_value=0)
    ax = pivot.plot(kind="barh", stacked=True, figsize=(10, 4.8), color=[colors[o] for o in order])
    ax.set_xlabel("Final BFB-positive amplicons")
    ax.set_ylabel("")
    ax.set_title("Final BFB Calls by Polarity Assessment Category")
    ax.legend(loc="lower right", fontsize=7)
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


def plot_reverse_fraction(summary: pd.DataFrame, out: Path) -> None:
    cur = summary[summary["cohort"] != "All cohorts"].copy()
    cur["reverse_percent"] = 100 * cur["reverse_fraction_of_assessed_final_bfb"]
    ax = cur.plot(kind="bar", x="cohort", y="reverse_percent", legend=False, color="#e15759", figsize=(8, 4.5))
    ax.set_ylabel("Reverse-polarity share of BFBArchitect-assessed final BFB calls (%)")
    ax.set_xlabel("")
    ax.set_title("Reverse-Polarity Use Among Assessed Final BFB Calls")
    max_percent = cur["reverse_percent"].max(skipna=True)
    max_percent = 0 if pd.isna(max_percent) else float(max_percent)
    ax.set_ylim(0, max(100, max_percent * 1.15))
    ax.grid(axis="y", alpha=0.25)
    for idx, row in cur.reset_index(drop=True).iterrows():
        value = row["reverse_percent"]
        if pd.notna(value):
            ax.text(idx, value + 1, f"{value:.1f}%", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


def selected_examples(calls: pd.DataFrame, out_dir: Path, n: int) -> pd.DataFrame:
    examples = calls[
        calls["bfbarchitect_final_positive"] & (calls["BFBArchitect_reverse_polarity_used"] == "True")
    ].copy()
    examples = examples.sort_values(["BFBArchitect_min_score_numeric", "cohort", "sample_name", "amplicon_number"]).head(n)
    rows = []
    assets = out_dir / "assets"
    for idx, (_, row) in enumerate(examples.iterrows(), 1):
        aa_src = find_aa_plot(row)
        bfb_src = find_bfb_plot(row)
        aa_asset = copy_or_convert_image(aa_src, assets / f"case_{idx:02d}_aa.png") if aa_src else None
        bfb_asset = copy_or_convert_image(bfb_src, assets / f"case_{idx:02d}_bfb.png") if bfb_src else None
        rows.append(
            {
                "rank": idx,
                "cohort": row["cohort"],
                "ref": row["ref"],
                "sample_name": row["sample_name"],
                "amplicon_number": row["amplicon_number"],
                "score": row["BFBArchitect_min_score_numeric"],
                "mode": row["bfbarchitect_mode"],
                "native_region": row["native_region"],
                "native_arm": row["native_arm"],
                "BFB_source": row["BFB_source"],
                "aa_plot": str(aa_src) if aa_src else "missing",
                "bfb_plot": str(bfb_src) if bfb_src else "missing",
                "aa_asset": str(aa_asset.relative_to(out_dir)) if aa_asset else "missing",
                "bfb_asset": str(bfb_asset.relative_to(out_dir)) if bfb_asset else "missing",
            }
        )
    out = pd.DataFrame(rows, columns=EXAMPLE_COLUMNS)
    out.to_csv(out_dir / "reverse_polarity_top_examples.tsv", sep="\t", index=False)
    return out


def read_score_checks(path: Path | None, out_dir: Path) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=SCORE_CHECK_COLUMNS)
    out = read_delimited_table(path)
    missing = set(SCORE_CHECK_COLUMNS) - set(out.columns)
    if missing:
        raise ValueError(f"score-check table {path} is missing required columns: {', '.join(sorted(missing))}")
    out.to_csv(out_dir / "four_mode_score_checks.tsv", sep="\t", index=False)
    return out


def format_count(value: object) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def df_to_html_table(df: pd.DataFrame, columns: list[str]) -> str:
    lines = ["<table>", "<thead><tr>"]
    for col in columns:
        lines.append(f"<th>{html.escape(col)}</th>")
    lines.append("</tr></thead><tbody>")
    for _, row in df.iterrows():
        lines.append("<tr>")
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float) and "fraction" in col:
                text = "NA" if pd.isna(value) else f"{100 * value:.1f}%"
            else:
                text = format_count(value)
            lines.append(f"<td>{html.escape(text)}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")
    return "\n".join(lines)


def df_to_markdown_table(df: pd.DataFrame, columns: list[str]) -> str:
    rows = [[col for col in columns]]
    rows.append(["---"] * len(columns))
    for _, row in df.iterrows():
        cur = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float) and "fraction" in col:
                text = "NA" if pd.isna(value) else f"{100 * value:.1f}%"
            elif isinstance(value, float) and col == "score":
                text = f"{value:.2f}"
            else:
                text = format_count(value)
            cur.append(text.replace("|", "\\|"))
        rows.append(cur)
    return "\n".join("| " + " | ".join(row) + " |" for row in rows)


def embedded_image_src(path: Path) -> str:
    suffix = path.suffix.lower()
    mime = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".gif": "image/gif",
        ".webp": "image/webp",
    }.get(suffix, "application/octet-stream")
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def figure_html(out_dir: Path, asset: object, alt: str, caption: str) -> str:
    if pd.isna(asset) or str(asset) == "missing":
        return "\n".join(
            [
                "<figure>",
                f"<div class='missing-image'>{html.escape(caption)} missing</div>",
                f"<figcaption>{html.escape(caption)}</figcaption>",
                "</figure>",
            ]
        )
    path = out_dir / str(asset)
    if not path.exists():
        return "\n".join(
            [
                "<figure>",
                f"<div class='missing-image'>{html.escape(caption)} missing: {html.escape(str(asset))}</div>",
                f"<figcaption>{html.escape(caption)}</figcaption>",
                "</figure>",
            ]
        )
    return "\n".join(
        [
            "<figure>",
            f"<img src='{embedded_image_src(path)}' alt='{html.escape(alt)}'>",
            f"<figcaption>{html.escape(caption)}</figcaption>",
            "</figure>",
        ]
    )


def example_sections(out_dir: Path, rows: pd.DataFrame, score_checks: pd.DataFrame | None = None) -> list[str]:
    if rows.empty:
        return ["<p class='note'>No reverse-polarity BFBArchitect final calls were available for example plots.</p>"]

    html_lines: list[str] = []
    for _, row in rows.iterrows():
        title = f"#{int(row['rank'])}: {row['cohort']} {row['sample_name']} {row['amplicon_number']}"
        meta = (
            f"ref={row['ref']} | score={float(row['score']):.2f} | mode={row['mode']} | "
            f"native arm={row['native_arm']} | region={row['native_region']} | source={row['BFB_source']}"
        )
        html_lines.extend(
            [
                "<section class='case'>",
                f"<h3>{html.escape(title)}</h3>",
                f"<div class='meta'>{html.escape(meta)}</div>",
            ]
        )
        if score_checks is not None and len(score_checks):
            cur_checks = score_checks[
                (score_checks["sample_name"] == row["sample_name"])
                & (score_checks["amplicon_number"] == row["amplicon_number"])
            ]
            if len(cur_checks):
                html_lines.extend(
                    [
                        "<h4>Four-mode score check</h4>",
                        "<p class='note'>Recomputed from the original AA graph with BFBArchitect graph mode, Gurobi, and 3 solver threads per mode.</p>",
                        df_to_html_table(cur_checks, ["mode", "run_region", "multiplicity", "score"]),
                    ]
                )
        html_lines.extend(
            [
                "<div class='plot-row'>",
                figure_html(out_dir, row["aa_asset"], "AA amplicon plot", "AA amplicon plot"),
                figure_html(out_dir, row["bfb_asset"], "BFBArchitect plot", "BFBArchitect reconstruction plot"),
                "</div>",
                "</section>",
            ]
        )
    return html_lines


def write_report(
    summary: pd.DataFrame,
    source: pd.DataFrame,
    mode: pd.DataFrame,
    score_checks: pd.DataFrame,
    examples: pd.DataFrame,
    out_dir: Path,
) -> None:
    style = """
    body { font-family: Arial, Helvetica, sans-serif; font-size: 15px; margin: 28px; color: #1f2933; background: #f7f8fa; }
    h1, h2, h3 { color: #111827; }
    .note { max-width: 980px; line-height: 1.45; }
    table { border-collapse: collapse; margin: 14px 0 28px; font-size: 14px; background: white; }
    th, td { border: 1px solid #d9dee7; padding: 7px 9px; text-align: left; vertical-align: top; }
    th { background: #eef2f7; }
    .plot-row { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; align-items: start; }
    .plot-row figure:nth-child(2) { text-align: center; }
    .plot-row figure:nth-child(2) img { width: 80%; }
    .case { background: white; border: 1px solid #d9dee7; border-radius: 6px; padding: 14px; margin: 18px 0 28px; }
    .case h3 { margin: 0 0 8px; font-size: 19px; }
    .meta { color: #4b5563; font-size: 14px; margin-bottom: 12px; }
    figure { margin: 0; }
    figure img { width: 100%; height: auto; border: 1px solid #d1d5db; background: white; }
    .missing-image { min-height: 180px; display: grid; place-items: center; border: 1px dashed #9ca3af; color: #4b5563; background: #f9fafb; padding: 12px; text-align: center; }
    figcaption { font-size: 13px; color: #4b5563; margin-top: 5px; }
    .summary-plots { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; max-width: 1200px; }
    .summary-plots img { width: 100%; border: 1px solid #d1d5db; background: white; }
    @media (max-width: 900px) { .plot-row, .summary-plots { grid-template-columns: 1fr; } }
    """
    html_lines = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        "<title>BFBArchitect Reverse-Polarity Usage</title>",
        f"<style>{style}</style>",
        "</head><body>",
        "<h1>BFBArchitect reverse-polarity usage</h1>",
        "<p class='note'>This report reads existing AmpliconClassifier profile tables and existing plot outputs only. "
        "It does not rerun BFBArchitect. Polarity is summarized only for final BFB calls that are supported by BFBArchitect; "
        "AC-native-only BFB calls are listed as polarity not assessed. Lower BFBArchitect score is treated as the stronger "
        "BFBArchitect fit when ranking example panels.</p>",
        "<h2>Cohort summary</h2>",
        df_to_html_table(
            summary,
            [
                "cohort",
                "ref",
                "amplicons",
                "final_bfb_positive",
                "bfbarchitect_assessed_final_bfb",
                "ac_native_only_not_assessed",
                "reverse_polarity_final_positive",
                "reverse_fraction_of_assessed_final_bfb",
            ],
        ),
        "<div class='summary-plots'>",
        f"<img src='{embedded_image_src(out_dir / 'mode_breakdown.png')}' alt='BFBArchitect mode breakdown'>",
        f"<img src='{embedded_image_src(out_dir / 'reverse_fraction.png')}' alt='Reverse polarity fraction'>",
        "</div>",
        "<h2>Final BFB source breakdown</h2>",
        df_to_html_table(source, ["cohort", "ref", "bfb_source_group", "count"]),
        "<h2>Final BFB polarity assessment breakdown</h2>",
        df_to_html_table(mode, ["cohort", "ref", "polarity_assessment_category", "count"]),
    ]
    html_lines.extend(
        [
        "<h2>Top reverse-polarity BFBArchitect examples</h2>",
        ]
    )
    html_lines.extend(example_sections(out_dir, examples, score_checks))
    html_lines.extend(["</body></html>"])
    (out_dir / "reverse_polarity_report.html").write_text("\n".join(html_lines))

    md_lines = [
        "# BFBArchitect reverse-polarity usage",
        "",
        "This report reads existing AmpliconClassifier profile tables and existing plot outputs only. It does not rerun BFBArchitect.",
        "",
        "Polarity is summarized only for final BFB calls that are supported by BFBArchitect; AC-native-only BFB calls are listed as polarity not assessed.",
        "",
        "## Cohort summary",
        "",
        df_to_markdown_table(
            summary,
            [
                "cohort",
                "ref",
                "amplicons",
                "final_bfb_positive",
                "bfbarchitect_assessed_final_bfb",
                "ac_native_only_not_assessed",
                "reverse_polarity_final_positive",
                "reverse_fraction_of_assessed_final_bfb",
            ],
        ),
        "",
        "## Final BFB polarity assessment breakdown",
        "",
        df_to_markdown_table(mode, ["cohort", "ref", "polarity_assessment_category", "count"]),
        "",
        "## Top reverse-polarity examples",
        "",
        df_to_markdown_table(
            examples,
            ["rank", "cohort", "ref", "sample_name", "amplicon_number", "score", "mode", "native_arm", "native_region"],
        ),
        "",
    ]
    (out_dir / "reverse_polarity_report.md").write_text("\n".join(md_lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Delimited table with cohort, root, and optional ref columns.",
    )
    parser.add_argument(
        "--score-checks",
        type=Path,
        help="Optional delimited table of externally recomputed four-mode scores to embed in matching examples.",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("reports/bfb_reverse_polarity_usage"))
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cohorts = read_manifest(args.manifest)
    calls, _ = read_calls(cohorts)
    calls.to_csv(args.out_dir / "all_profile_calls_with_bfbarchitect_mode.tsv", sep="\t", index=False)
    summary, source, mode = make_summary_tables(calls, args.out_dir)
    plot_mode_breakdown(mode, args.out_dir / "mode_breakdown.png")
    plot_reverse_fraction(summary, args.out_dir / "reverse_fraction.png")
    examples = selected_examples(calls, args.out_dir, args.top_n)
    score_checks = read_score_checks(args.score_checks, args.out_dir)
    write_report(summary, source, mode, score_checks, examples, args.out_dir)

    print(args.out_dir / "reverse_polarity_report.html")
    print(args.out_dir / "reverse_polarity_usage_summary.tsv")
    if args.score_checks:
        print(args.out_dir / "four_mode_score_checks.tsv")
    print(args.out_dir / "reverse_polarity_top_examples.tsv")


if __name__ == "__main__":
    main()
