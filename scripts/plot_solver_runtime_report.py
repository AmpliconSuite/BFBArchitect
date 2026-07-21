#!/usr/bin/env python3
"""Generate plot-first HTML reports for BFBArchitect solver runtime benchmarks."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


SOLVERS = ["gurobi", "mosek", "cbc"]
SOLVER_COLORS = {"gurobi": "#2f80ed", "mosek": "#27ae60", "cbc": "#eb5757"}
THREAD_COLORS = {1: "#1f2933", 2: "#2f80ed", 3: "#27ae60", 4: "#f2994a", 8: "#9b51e0", 16: "#eb5757"}


def score_cmap():
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#d1d5db")
    return cmap


def read_replicates(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for idx, path in enumerate(paths, 1):
        df = pd.read_csv(path, sep="\t")
        df["replicate"] = idx
        df["replicate_file"] = str(path)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["elapsed_s"] = pd.to_numeric(out["elapsed_s"], errors="coerce")
    out["score_min"] = pd.to_numeric(out["score_min"], errors="coerce")
    out["T"] = pd.to_numeric(out["T"], errors="coerce")
    out["threads"] = pd.to_numeric(out["threads"], errors="coerce")
    if "timeout_s" in out:
        out["timeout_s"] = pd.to_numeric(out["timeout_s"], errors="coerce")
    else:
        out["timeout_s"] = np.nan
    if "max_active_threads" in out:
        out["max_active_threads"] = pd.to_numeric(out["max_active_threads"], errors="coerce")
    else:
        out["max_active_threads"] = np.nan
    out["case_key"] = out["source"].astype(str) + " / " + out["case_id"].astype(str)
    return out


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    ok = df[df["status"] == "ok"].copy()
    grouped = ok.groupby(["source", "case_id", "case_key", "mode", "T", "solver", "threads"], dropna=False)
    summary = grouped.agg(
        count=("elapsed_s", "count"),
        median=("elapsed_s", "median"),
        mean=("elapsed_s", "mean"),
        min=("elapsed_s", "min"),
        max=("elapsed_s", "max"),
        std=("elapsed_s", "std"),
        score_median=("score_min", "median"),
        score_min=("score_min", "min"),
        score_max=("score_min", "max"),
        timeout_s=("timeout_s", "max"),
        max_active_threads=("max_active_threads", "max"),
    ).reset_index()
    summary["cv"] = summary["std"] / summary["mean"]
    return summary


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_solver_distribution(summary: pd.DataFrame, out: Path) -> None:
    solvers = SOLVERS
    data = [summary.loc[summary["solver"] == solver, "median"].dropna().values for solver in solvers]
    plt.figure(figsize=(7.5, 4.8))
    ax = plt.gca()
    ax.boxplot(data, tick_labels=solvers, showfliers=True)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("Runtime per case, median of replicates (s)")
    ax.set_title("Solver Runtime Distribution Across Selected Cases")
    ax.grid(True, axis="y", alpha=0.3)
    savefig(out)


def plot_runtime_ecdf(summary: pd.DataFrame, out: Path) -> None:
    linestyles = {1: "-", 2: "--", 3: "-.", 4: ":", 8: (0, (3, 1, 1, 1)), 16: (0, (5, 2))}
    threads = sorted(int(t) for t in summary["threads"].dropna().unique())
    plt.figure(figsize=(9.5, 6.2))
    ax = plt.gca()
    for solver in SOLVERS:
        for thread in threads:
            vals = summary[(summary["solver"] == solver) & (summary["threads"] == thread)]["median"].dropna().sort_values().values
            if len(vals) == 0:
                continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.step(vals, y, where="post", color=SOLVER_COLORS.get(solver), linestyle=linestyles.get(thread, "-"),
                    alpha=0.85, label=f"{solver} t={thread}")
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("Runtime per case, median of replicates (s)")
    ax.set_ylabel("Fraction of cases solved")
    ax.set_title("Runtime ECDF by Solver and Solver Threads Per Job")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=3)
    savefig(out)


def plot_runtime_ecdf_faceted(summary: pd.DataFrame, out: Path) -> None:
    threads = sorted(int(t) for t in summary["threads"].dropna().unique())
    all_vals = summary["median"].dropna()
    max_runtime = max(1.0, float(all_vals.max())) if len(all_vals) else 1.0
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), sharey=True)
    for ax, solver in zip(axes, SOLVERS):
        for thread in threads:
            vals = summary[(summary["solver"] == solver) & (summary["threads"] == thread)]["median"].dropna().sort_values().values
            if len(vals) == 0:
                continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.step(vals, y, where="post", color=THREAD_COLORS.get(thread), alpha=0.9, label=f"{thread}")
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xlim(0, max_runtime * 1.05)
        ax.set_title(solver)
        ax.set_xlabel("Runtime (s)")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Fraction of cases solved")
    axes[-1].legend(title="Solver threads per job", fontsize=8)
    fig.suptitle("Runtime ECDF Small Multiples")
    savefig(out)


def plot_runtime_heatmap(summary: pd.DataFrame, out: Path) -> None:
    table = summary.pivot_table(index="solver", columns="threads", values="median", aggfunc="median")
    table = table.reindex([s for s in SOLVERS if s in table.index])
    cols = sorted(table.columns)
    table = table[cols]
    values = np.log10(table.replace(0, np.nan).values.astype(float))
    plt.figure(figsize=(7.2, 3.4))
    ax = plt.gca()
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([str(int(c)) for c in cols])
    ax.set_yticks(np.arange(len(table.index)))
    ax.set_yticklabels(table.index)
    ax.set_xlabel("Solver threads per job")
    ax.set_title("Median Runtime Across Cases (log10 seconds)")
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            val = table.iloc[i, j]
            label = "" if pd.isna(val) else f"{val:.1f}s"
            ax.text(j, i, label, ha="center", va="center", color="white" if values[i, j] > np.nanmedian(values) else "black", fontsize=8)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("log10(seconds)")
    savefig(out)


def plot_thread_scaling(summary: pd.DataFrame, out: Path) -> None:
    base = summary[summary["threads"] == 1][["source", "case_id", "solver", "median"]].rename(columns={"median": "median_t1"})
    merged = summary.merge(base, on=["source", "case_id", "solver"], how="inner")
    merged = merged[merged["median"] > 0].copy()
    merged["speedup"] = merged["median_t1"] / merged["median"]
    stats = merged.groupby(["solver", "threads"])["speedup"].agg(
        median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75),
    ).reset_index()
    plt.figure(figsize=(7.8, 4.8))
    ax = plt.gca()
    for solver, cur in stats.groupby("solver"):
        cur = cur.sort_values("threads")
        ax.plot(cur["threads"], cur["median"], marker="o", label=solver, color=SOLVER_COLORS.get(solver))
        ax.fill_between(cur["threads"].to_numpy(dtype=float), cur["q25"].to_numpy(dtype=float),
                        cur["q75"].to_numpy(dtype=float), alpha=0.18, color=SOLVER_COLORS.get(solver))
    ax.axhline(1, color="#111827", lw=1, ls="--")
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted(stats["threads"].unique()))
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("Solver threads per job")
    ax.set_ylabel("Speedup vs 1 thread")
    ax.set_title("Thread Scaling by Solver (median, IQR across cases)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(out)


def plot_absolute_runtime_by_threads(summary: pd.DataFrame, out: Path) -> None:
    stats = summary.groupby(["solver", "threads"])["median"].agg(
        median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75),
    ).reset_index()
    plt.figure(figsize=(7.8, 4.8))
    ax = plt.gca()
    for solver, cur in stats.groupby("solver"):
        cur = cur.sort_values("threads")
        ax.plot(cur["threads"], cur["median"], marker="o", label=solver, color=SOLVER_COLORS.get(solver))
        ax.fill_between(cur["threads"].to_numpy(dtype=float), cur["q25"].to_numpy(dtype=float),
                        cur["q75"].to_numpy(dtype=float), alpha=0.18, color=SOLVER_COLORS.get(solver))
    ax.set_xscale("log", base=2)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(sorted(stats["threads"].unique()))
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("Solver threads per job")
    ax.set_ylabel("Runtime per case, median of replicates (s)")
    ax.set_title("Absolute Runtime by Thread Count")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(out)


def plot_win_fraction(summary: pd.DataFrame, out: Path, tolerance: float = 0.0) -> None:
    pivot = summary.pivot_table(
        index=["source", "case_id", "threads"],
        columns="solver",
        values="median",
        aggfunc="first",
    )
    solvers = [s for s in ["gurobi", "mosek", "cbc"] if s in pivot.columns]
    rows = []
    for (source, case_id, threads), row in pivot.dropna(subset=solvers).iterrows():
        best = row[solvers].min()
        winners = [solver for solver in solvers if row[solver] <= best * (1 + tolerance)]
        for solver in winners:
            rows.append({"threads": threads, "solver": solver, "weight": 1 / len(winners)})
    wins = pd.DataFrame(rows)
    if wins.empty:
        return
    table = wins.pivot_table(index="threads", columns="solver", values="weight", aggfunc="sum", fill_value=0)
    counts = pivot.dropna(subset=solvers).groupby(level="threads").size()
    table = table.div(counts, axis=0)
    table = table[[s for s in solvers if s in table.columns]]
    plt.figure(figsize=(7.8, 4.8))
    ax = plt.gca()
    bottom = np.zeros(len(table))
    x = np.arange(len(table.index))
    for solver in table.columns:
        ax.bar(x, table[solver].values, bottom=bottom, label=solver, color=SOLVER_COLORS.get(solver))
        bottom += table[solver].values
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(t)) for t in table.index])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Solver threads per job")
    ax.set_ylabel("Fraction of cases")
    suffix = "" if tolerance == 0 else f" (within {tolerance:.0%} of fastest)"
    ax.set_title(f"Fastest Solver Fraction by Thread Count{suffix}")
    ax.legend()
    savefig(out)


def plot_pairwise_ratios_by_thread(summary: pd.DataFrame, out: Path) -> None:
    pivot = summary.pivot_table(
        index=["source", "case_id", "threads"],
        columns="solver",
        values="median",
        aggfunc="first",
    )
    pairs = [("cbc", "gurobi"), ("mosek", "gurobi"), ("cbc", "mosek")]
    rows = []
    for pair in pairs:
        if not set(pair).issubset(pivot.columns):
            continue
        numerator = pivot[pair[0]].replace(0, np.nan)
        denominator = pivot[pair[1]].replace(0, np.nan)
        ratio = numerator / denominator
        tmp = ratio.reset_index(name="ratio").dropna()
        tmp["pair"] = f"{pair[0]}/{pair[1]}"
        rows.append(tmp)
    if not rows:
        return
    ratios = pd.concat(rows, ignore_index=True)
    ratios = ratios[(ratios["ratio"] > 0) & np.isfinite(ratios["ratio"])].copy()
    ratios["log2_ratio"] = np.log2(ratios["ratio"])
    threads = sorted(ratios["threads"].unique())
    pairs_order = [f"{a}/{b}" for a, b in pairs]
    plt.figure(figsize=(10.0, 5.6))
    ax = plt.gca()
    positions = []
    data = []
    labels = []
    pos = 1
    for thread in threads:
        for pair in pairs_order:
            vals = ratios[(ratios["threads"] == thread) & (ratios["pair"] == pair)]["log2_ratio"].dropna().values
            if len(vals):
                data.append(vals)
                positions.append(pos)
                labels.append(f"{pair}\nt={int(thread)}")
            pos += 1
        pos += 0.8
    ax.boxplot(data, positions=positions, widths=0.65, showfliers=True)
    ax.axhline(0, color="#111827", lw=1, ls="--")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("log2(runtime ratio)")
    ax.set_title("Pairwise Solver Runtime Ratios by Thread Count")
    ax.grid(True, axis="y", alpha=0.3)
    savefig(out)


def plot_case_solver_bars(summary: pd.DataFrame, out: Path, top_n: int = 20) -> None:
    standard = summary[summary["threads"] == 3].copy()
    pivot = standard.pivot_table(index="case_key", columns="solver", values="median", aggfunc="first")
    order = pivot.max(axis=1).sort_values(ascending=True).tail(top_n).index
    pivot = pivot.loc[order, [c for c in ["gurobi", "mosek", "cbc"] if c in pivot.columns]]
    y = np.arange(len(pivot))
    width = 0.25
    plt.figure(figsize=(10.5, max(6, 0.34 * len(pivot) + 1.5)))
    ax = plt.gca()
    offsets = np.linspace(-width, width, len(pivot.columns))
    for offset, solver in zip(offsets, pivot.columns):
        ax.barh(y + offset, pivot[solver], height=width, label=solver, color=SOLVER_COLORS.get(solver))
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("Runtime at 3 solver threads per job, median of replicates (s)")
    ax.set_yticks(y)
    ax.set_yticklabels([label.split(" / ", 1)[1] for label in pivot.index], fontsize=8)
    ax.set_title("Per-Case Solver Runtime Comparison at 3 Solver Threads")
    ax.legend()
    ax.grid(True, axis="x", alpha=0.3)
    savefig(out)


def plot_cbc_slow(summary: pd.DataFrame, out: Path) -> None:
    cbc = summary[(summary["solver"] == "cbc") & (summary["threads"] == 3) & (summary["median"] >= 30)].copy()
    cbc = cbc.sort_values("median", ascending=True)
    plt.figure(figsize=(9.5, max(3.8, 0.55 * len(cbc) + 1.5)))
    ax = plt.gca()
    labels = cbc["case_id"].tolist()
    y = np.arange(len(cbc))
    med = cbc["median"].values
    lower = med - cbc["min"].values
    upper = cbc["max"].values - med
    ax.barh(y, med, color="#eb5757", alpha=0.85)
    ax.errorbar(med, y, xerr=[lower, upper], fmt="none", ecolor="#1f2933", capsize=3, lw=1)
    ax.axvline(30, color="#111827", lw=1, ls="--")
    ax.set_xlabel("CBC runtime at 3 solver threads per job (s), median with min/max across replicates")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title("CBC Slow Cases at 3 Solver Threads (Median >= 30s)")
    ax.grid(True, axis="x", alpha=0.3)
    savefig(out)


def plot_replicate_scatter(df: pd.DataFrame, out: Path) -> None:
    ok = df[df["status"] == "ok"].copy()
    key_cols = ["source", "case_id", "solver", "threads"]
    wide = ok.pivot_table(index=key_cols, columns="replicate", values="elapsed_s", aggfunc="first").reset_index()
    plt.figure(figsize=(6.8, 6.2))
    ax = plt.gca()
    for solver, cur in wide.groupby("solver"):
        if 1 not in cur or 2 not in cur:
            continue
        ax.scatter(cur[1], cur[2], label=f"{solver}: rep2", alpha=0.75, s=28, marker="o", color=SOLVER_COLORS.get(solver))
        if 3 in cur:
            ax.scatter(cur[1], cur[3], label=f"{solver}: rep3", alpha=0.55, s=28, marker="x", color=SOLVER_COLORS.get(solver))
    vals = ok["elapsed_s"].dropna()
    max_v = max(1.0, vals.max())
    ax.plot([0, max_v], [0, max_v], color="#111827", lw=1, ls="--")
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("Replicate 1 runtime (s)")
    ax.set_ylabel("Replicate 2/3 runtime (s)")
    ax.set_title("Replicate Stability")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    savefig(out)


def plot_score_delta_by_solver_thread(summary: pd.DataFrame, out: Path) -> None:
    score = summary.dropna(subset=["score_median"]).copy()
    if score.empty:
        return
    best = score.groupby(["source", "case_id"])["score_median"].min().rename("best_score").reset_index()
    score = score.merge(best, on=["source", "case_id"], how="left")
    score["score_delta"] = score["score_median"] - score["best_score"]
    stats = score.groupby(["solver", "threads"])["score_delta"].agg(
        median="median",
        q75=lambda x: x.quantile(0.75),
        max="max",
    ).reset_index()
    plt.figure(figsize=(8.2, 4.8))
    ax = plt.gca()
    for solver, cur in stats.groupby("solver"):
        cur = cur.sort_values("threads")
        ax.plot(cur["threads"], cur["median"], marker="o", label=f"{solver} median", color=SOLVER_COLORS.get(solver))
        ax.plot(cur["threads"], cur["q75"], marker=".", linestyle="--", alpha=0.7, color=SOLVER_COLORS.get(solver))
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted(stats["threads"].unique()))
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("Solver threads per job")
    ax.set_ylabel("Score delta from best score for the case")
    ax.set_title("BFB Score Delta by Solver and Thread Count")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    savefig(out)


def plot_score_by_amplicon_3threads(summary: pd.DataFrame, out: Path) -> None:
    cur = summary[summary["threads"] == 3].copy()
    if cur.empty:
        return
    cur["amplicon"] = cur["case_id"]
    pivot = cur.pivot_table(index="amplicon", columns="solver", values="score_median", aggfunc="median", dropna=False)
    cols = [solver for solver in SOLVERS if solver in pivot.columns]
    pivot = pivot[cols]
    sort_key = pivot.max(axis=1).fillna(float("inf")).sort_values(ascending=True)
    pivot = pivot.loc[sort_key.index]
    plt.figure(figsize=(7.2, max(5.0, 0.34 * len(pivot) + 1.4)))
    ax = plt.gca()
    values = pivot.values.astype(float)
    im = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=score_cmap())
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("BFBArchitect Score Per Amplicon at 3 Solver Threads")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i, j]):
                ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", fontsize=7)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("BFBArchitect score (lower is better)")
    ax.text(1.01, -0.06, "grey = no score returned", transform=ax.transAxes, fontsize=8, color="#4b5563")
    savefig(out)


def plot_score_by_amplicon_thread_matrix(summary: pd.DataFrame, out: Path) -> None:
    score = summary.copy()
    if score.empty:
        return
    score["amplicon"] = score["case_id"]
    score["column"] = score["solver"] + " t=" + score["threads"].astype(int).astype(str)
    cols = [f"{solver} t={thread}" for thread in sorted(score["threads"].dropna().unique().astype(int)) for solver in SOLVERS]
    pivot = score.pivot_table(index="amplicon", columns="column", values="score_median", aggfunc="median", dropna=False)
    cols = [col for col in cols if col in pivot.columns]
    pivot = pivot[cols]
    sort_key = pivot.max(axis=1).fillna(float("inf")).sort_values(ascending=True)
    pivot = pivot.loc[sort_key.index]
    plt.figure(figsize=(13.5, max(5.5, 0.34 * len(pivot) + 1.6)))
    ax = plt.gca()
    values = pivot.values.astype(float)
    im = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=score_cmap())
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("BFBArchitect Score Per Amplicon by Solver and Thread Count")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("BFBArchitect score (lower is better)")
    ax.text(1.01, -0.06, "grey = no score returned", transform=ax.transAxes, fontsize=8, color="#4b5563")
    savefig(out)


def plot_score_range_by_amplicon(summary: pd.DataFrame, out: Path) -> None:
    cases = summary.groupby(["source", "case_id"], dropna=False).size().reset_index()[["source", "case_id"]]
    if cases.empty:
        return
    score = summary.dropna(subset=["score_median"]).copy()
    if score.empty:
        ranges = cases.assign(min_score=np.nan, max_score=np.nan, median_score=np.nan)
    else:
        scored = score.groupby(["source", "case_id"])["score_median"].agg(
            min_score="min",
            max_score="max",
            median_score="median",
        ).reset_index()
        ranges = cases.merge(scored, on=["source", "case_id"], how="left")
    ranges["score_range"] = ranges["max_score"] - ranges["min_score"]
    ranges["missing"] = ranges["score_range"].isna()
    ranges = ranges.sort_values(["missing", "score_range"], ascending=[False, True])
    plt.figure(figsize=(9.0, max(5.0, 0.34 * len(ranges) + 1.5)))
    ax = plt.gca()
    y = np.arange(len(ranges))
    colors = np.where(ranges["missing"], "#d1d5db", "#2f80ed")
    ax.barh(y, ranges["score_range"].fillna(0), color=colors, alpha=0.85)
    for idx, row in enumerate(ranges.itertuples()):
        if row.missing:
            ax.text(0.02, idx, "no score", va="center", ha="left", fontsize=7, color="#4b5563")
    ax.set_yticks(y)
    ax.set_yticklabels(ranges["case_id"], fontsize=8)
    ax.set_xlabel("Score range across solver/thread combinations")
    ax.set_title("Per-Amplicon BFBArchitect Score Sensitivity")
    ax.legend(handles=[Patch(facecolor="#d1d5db", label="no score returned")], loc="lower right", fontsize=8)
    ax.grid(True, axis="x", alpha=0.3)
    savefig(out)


def write_score_cutoff_summary(summary: pd.DataFrame, out: Path, cutoff: float) -> pd.DataFrame:
    cases = summary.groupby(["source", "case_id"], dropna=False).size().rename("total_configs").reset_index()
    score = summary.dropna(subset=["score_median"]).copy()
    if score.empty:
        rows = cases.assign(
            n_scored_configs=0,
            min_score=np.nan,
            max_score=np.nan,
            score_range=np.nan,
            bfb_call_fraction=np.nan,
            crosses_cutoff=False,
            missing_all_scores=True,
            cutoff=cutoff,
        )
        rows.to_csv(out, sep="\t", index=False)
        return rows
    score["bfb_call"] = score["score_median"] <= cutoff
    scored = score.groupby(["source", "case_id"], dropna=False).agg(
        n_configs=("score_median", "count"),
        min_score=("score_median", "min"),
        max_score=("score_median", "max"),
        score_range=("score_median", lambda x: x.max() - x.min()),
        bfb_call_fraction=("bfb_call", "mean"),
    ).reset_index()
    scored = scored.rename(columns={"n_configs": "n_scored_configs"})
    rows = cases.merge(scored, on=["source", "case_id"], how="left")
    rows["crosses_cutoff"] = (rows["min_score"] <= cutoff) & (rows["max_score"] > cutoff)
    rows["missing_all_scores"] = rows["n_scored_configs"].isna()
    rows["n_scored_configs"] = rows["n_scored_configs"].fillna(0).astype(int)
    rows["cutoff"] = cutoff
    rows = rows.sort_values(["crosses_cutoff", "missing_all_scores", "score_range"], ascending=[False, True, False])
    rows.to_csv(out, sep="\t", index=False)
    return rows


def plot_score_cutoff_fraction(summary: pd.DataFrame, out: Path, cutoff: float) -> None:
    cases = summary.groupby(["source", "case_id"], dropna=False).size().reset_index()[["source", "case_id"]]
    if cases.empty:
        return
    score = summary.dropna(subset=["score_median"]).copy()
    if score.empty:
        rows = cases.assign(call_fraction=np.nan, min_score=np.nan, max_score=np.nan, crosses=False, missing=True)
    else:
        score["bfb_call"] = score["score_median"] <= cutoff
        scored = score.groupby(["source", "case_id"], dropna=False).agg(
            call_fraction=("bfb_call", "mean"),
            min_score=("score_median", "min"),
            max_score=("score_median", "max"),
        ).reset_index()
        rows = cases.merge(scored, on=["source", "case_id"], how="left")
        rows["crosses"] = (rows["min_score"] <= cutoff) & (rows["max_score"] > cutoff)
        rows["missing"] = rows["call_fraction"].isna()
    rows = rows.sort_values(["crosses", "missing", "call_fraction", "case_id"], ascending=[False, True, True, True])
    colors = np.where(rows["missing"], "#d1d5db", np.where(rows["crosses"], "#eb5757", "#2f80ed"))
    plt.figure(figsize=(9.0, max(5.0, 0.34 * len(rows) + 1.5)))
    ax = plt.gca()
    y = np.arange(len(rows))
    ax.barh(y, rows["call_fraction"].fillna(0), color=colors, alpha=0.85)
    ax.set_xlim(0, 1)
    ax.set_yticks(y)
    ax.set_yticklabels(rows["case_id"], fontsize=8)
    ax.set_xlabel(f"Fraction of solver/thread configurations with score <= {cutoff:g}")
    ax.set_title("BFB Cutoff Stability by Amplicon")
    ax.legend(
        handles=[
            Patch(facecolor="#eb5757", label="crosses cutoff"),
            Patch(facecolor="#2f80ed", label="stable scored case"),
            Patch(facecolor="#d1d5db", label="no score returned"),
        ],
        loc="lower right",
        fontsize=8,
    )
    ax.grid(True, axis="x", alpha=0.3)
    savefig(out)


def plot_score_cutoff_matrix(summary: pd.DataFrame, out: Path, cutoff: float) -> None:
    score = summary.copy()
    if score.empty:
        return
    score["column"] = score["solver"] + " t=" + score["threads"].astype(int).astype(str)
    cols = [f"{solver} t={thread}" for thread in sorted(score["threads"].dropna().unique().astype(int)) for solver in SOLVERS]
    pivot = score.pivot_table(index="case_key", columns="column", values="score_median", aggfunc="median", dropna=False)
    cols = [col for col in cols if col in pivot.columns]
    pivot = pivot[cols]
    ranges = pivot.agg(["min", "max"], axis=1)
    crosses = (ranges["min"] <= cutoff) & (ranges["max"] > cutoff)
    order = (
        pd.DataFrame({
            "crosses": crosses,
            "distance": (ranges.mean(axis=1) - cutoff).abs().fillna(float("inf")),
            "range": ranges["max"] - ranges["min"],
            "missing": pivot.isna().all(axis=1),
        })
        .sort_values(["crosses", "missing", "distance", "range"], ascending=[False, True, True, False])
        .index
    )
    pivot = pivot.loc[order]
    calls = np.full(pivot.shape, 2.0)
    values = pivot.values.astype(float)
    calls[np.isfinite(values) & (values <= cutoff)] = 0.0
    calls[np.isfinite(values) & (values > cutoff)] = 1.0
    cmap = ListedColormap(["#2f80ed", "#f2994a", "#d1d5db"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    plt.figure(figsize=(13.5, max(5.5, 0.35 * len(pivot) + 1.6)))
    ax = plt.gca()
    ax.imshow(calls, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([label.split(" / ", 1)[1] for label in pivot.index], fontsize=8)
    ax.set_title(f"BFB Classification Matrix at Score Cutoff {cutoff:g}")
    ax.legend(
        handles=[
            Patch(facecolor="#2f80ed", label=f"BFB score <= {cutoff:g}"),
            Patch(facecolor="#f2994a", label=f"not BFB score > {cutoff:g}"),
            Patch(facecolor="#d1d5db", label="no score returned"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=8,
    )
    savefig(out)


def plot_score_change_heatmap(summary: pd.DataFrame, out: Path) -> None:
    score = summary.dropna(subset=["score_median"]).copy()
    if score.empty:
        return
    best = score.groupby(["source", "case_id"])["score_median"].min().rename("best_score").reset_index()
    score = score.merge(best, on=["source", "case_id"], how="left")
    score["score_delta"] = score["score_median"] - score["best_score"]
    interesting = (
        score.groupby("case_key")["score_delta"].max()
        .sort_values(ascending=False)
        .head(20)
        .index
    )
    score = score[score["case_key"].isin(interesting)]
    score["column"] = score["solver"] + " t=" + score["threads"].astype(int).astype(str)
    cols = [f"{solver} t={thread}" for thread in sorted(score["threads"].dropna().unique().astype(int)) for solver in SOLVERS]
    pivot = score.pivot_table(index="case_key", columns="column", values="score_delta", aggfunc="median")
    cols = [c for c in cols if c in pivot.columns]
    pivot = pivot[cols]
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=True).index]
    plt.figure(figsize=(12.5, max(5.5, 0.35 * len(pivot) + 1.6)))
    ax = plt.gca()
    values = pivot.fillna(0).values
    im = ax.imshow(values, aspect="auto", cmap="magma")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([label.split(" / ", 1)[1] for label in pivot.index], fontsize=8)
    ax.set_title("BFB Score Delta From Best Observed Score")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Score delta (lower is better)")
    savefig(out)


def plot_speedup_heatmap(summary: pd.DataFrame, out: Path) -> None:
    standard = summary[summary["threads"] == 3].copy()
    pivot = standard.pivot_table(index="case_key", columns="solver", values="median", aggfunc="first")
    if not {"gurobi", "mosek", "cbc"}.issubset(pivot.columns):
        return
    nonzero = pivot.replace(0, np.nan)
    ratios = pd.DataFrame({
        "CBC/Gurobi": nonzero["cbc"] / nonzero["gurobi"],
        "MOSEK/Gurobi": nonzero["mosek"] / nonzero["gurobi"],
        "CBC/MOSEK": nonzero["cbc"] / nonzero["mosek"],
    })
    log_ratios = np.log2(ratios.replace([np.inf, -np.inf], np.nan))
    order_metric = log_ratios.apply(
        lambda row: np.nanmax(np.abs(row)) if np.isfinite(row).any() else -1,
        axis=1,
    )
    order = order_metric.sort_values(ascending=False).index
    ratios = ratios.loc[order]
    values = np.log2(ratios.replace([np.inf, -np.inf], np.nan).values.astype(float))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#d1d5db")
    plt.figure(figsize=(7.5, max(5.5, 0.32 * len(ratios) + 1.5)))
    ax = plt.gca()
    im = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=cmap, vmin=-3, vmax=3)
    ax.set_xticks(np.arange(len(ratios.columns)))
    ax.set_xticklabels(ratios.columns)
    ax.set_yticks(np.arange(len(ratios.index)))
    ax.set_yticklabels([label.split(" / ", 1)[1] for label in ratios.index], fontsize=8)
    ax.set_title("Relative Runtime Ratios for All Selected Cases at 3 Solver Threads")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            label = "" if np.isnan(ratios.iloc[i, j]) else f"{ratios.iloc[i, j]:.1f}x"
            ax.text(j, i, label, ha="center", va="center", fontsize=7)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("log2(runtime ratio)")
    ax.text(1.01, -0.06, "grey = invalid ratio because denominator runtime was 0", transform=ax.transAxes, fontsize=8, color="#4b5563")
    savefig(out)


def write_timeout_followup(df: pd.DataFrame, out: Path, score_cutoff: float) -> pd.DataFrame:
    rows = []
    for threads in [2, 3]:
        sub = df[df["threads"] == threads].copy()
        timeout_cases = (
            sub[sub["status"] == "timeout"]
            .groupby(["replicate", "source", "case_id"], dropna=False)["solver"]
            .nunique()
            .reset_index(name="timeout_solver_count")
        )
        timeout_cases = timeout_cases[timeout_cases["timeout_solver_count"] >= 3]
        for item in timeout_cases.itertuples():
            sixteen = df[
                (df["replicate"] == item.replicate)
                & (df["source"] == item.source)
                & (df["case_id"] == item.case_id)
                & (df["threads"] == 16)
            ]
            for row in sixteen.itertuples():
                rows.append({
                    "replicate": item.replicate,
                    "source": item.source,
                    "case_id": item.case_id,
                    "timed_out_threads": threads,
                    "solver_16_threads": row.solver,
                    "status_16_threads": row.status,
                    "elapsed_16_threads": row.elapsed_s,
                    "score_16_threads": row.score_min,
                    "bfb_at_16_threads": row.score_min <= score_cutoff if pd.notna(row.score_min) else "",
                })
    fields = [
        "replicate", "source", "case_id", "timed_out_threads", "solver_16_threads",
        "status_16_threads", "elapsed_16_threads", "score_16_threads", "bfb_at_16_threads",
    ]
    out_df = pd.DataFrame(rows, columns=fields)
    out_df.to_csv(out, sep="\t", index=False)
    return out_df


def write_html(
    out_dir: Path,
    images: dict[str, str],
    summary: pd.DataFrame,
    replicate_paths: list[Path],
    score_cutoff_summary: pd.DataFrame,
    score_cutoff: float,
    timeout_followup: pd.DataFrame,
) -> None:
    thread_values = sorted(int(t) for t in summary["threads"].dropna().unique())
    case_count = summary[["source", "case_id"]].drop_duplicates().shape[0]
    timeout_values = sorted(
        int(t) for t in pd.to_numeric(summary.get("timeout_s"), errors="coerce").dropna().unique()
    ) if "timeout_s" in summary else []
    timeout_note = ", ".join(f"{t}s" for t in timeout_values) if timeout_values else "not recorded"
    slow = summary[(summary["solver"] == "cbc") & (summary["threads"] == 3) & (summary["median"] >= 30)].sort_values("median", ascending=False)
    solver_bits = []
    for solver, cur in summary[summary["threads"] == 3].groupby("solver"):
        solver_bits.append(
            f"<li><b>{html.escape(solver)}</b>: n={len(cur)}, "
            f"median case median at 3 solver threads={cur['median'].median():.2f}s, "
            f"max case median at 3 solver threads={cur['median'].max():.2f}s</li>"
        )
    slow_bits = "".join(
        f"<li>{html.escape(row.case_id)}: CBC median {row.median:.1f}s "
        f"(min {row.min:.1f}, max {row.max:.1f})</li>"
        for row in slow.itertuples()
    ) or "<li>None at 3 solver threads.</li>"
    crossing = score_cutoff_summary[score_cutoff_summary.get("crosses_cutoff", False) == True]  # noqa: E712
    missing_scores = score_cutoff_summary[score_cutoff_summary.get("missing_all_scores", False) == True]  # noqa: E712
    crossing_bits = "".join(
        f"<li>{html.escape(row.case_id)}: score range {row.min_score:.3g}-{row.max_score:.3g}, "
        f"BFB-call fraction {row.bfb_call_fraction:.2f}</li>"
        for row in crossing.itertuples()
    ) or f"<li>No selected amplicons crossed the {score_cutoff:g} cutoff.</li>"
    missing_bits = "".join(
        f"<li>{html.escape(row.case_id)}</li>"
        for row in missing_scores.itertuples()
    ) or "<li>None.</li>"
    if timeout_followup.empty:
        timeout_bits = "<li>No cases timed out on all solvers at 2 or 3 solver threads in this run.</li>"
    else:
        timeout_bits = "".join(
            f"<li>{html.escape(row.case_id)}, t={row.timed_out_threads}: "
            f"{html.escape(row.solver_16_threads)} at 16 threads score={html.escape(str(row.score_16_threads))}, "
            f"status={html.escape(str(row.status_16_threads))}</li>"
            for row in timeout_followup.itertuples()
        )
    image_sections = "\n".join(
        f"<section><h2>{html.escape(title)}</h2><img src=\"{html.escape(path)}\" alt=\"{html.escape(title)}\"></section>"
        for title, path in images.items()
    )
    reps = "".join(f"<li><code>{html.escape(str(path))}</code></li>" for path in replicate_paths)
    text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>BFBArchitect Solver Runtime Plots</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 32px; color: #1f2933; }}
main {{ max-width: 1180px; }}
section {{ margin: 32px 0; }}
img {{ max-width: 100%; border: 1px solid #d9e2ec; }}
.note {{ color: #52606d; }}
</style>
</head>
<body>
<main>
<h1>BFBArchitect Solver Runtime Plots</h1>
<p class="note">{case_count} selected cases, solver threads per job={', '.join(map(str, thread_values))}, timeout={html.escape(timeout_note)}, replicates={len(replicate_paths)}. "Threads" in plot labels means solver threads assigned to one BFBArchitect solve, not total machine threads.</p>
<section><h2>Inputs</h2><ul>{reps}</ul></section>
<section><h2>3-Thread Summary</h2><ul>{''.join(solver_bits)}</ul></section>
<section><h2>CBC Slow Cases at 3 Solver Threads</h2><ul>{slow_bits}</ul></section>
<section><h2>BFB Score Cutoff Crossings</h2><p class="note">BFB cutoff: score <= {score_cutoff:g}.</p><ul>{crossing_bits}</ul></section>
<section><h2>Missing BFB Scores</h2><p class="note">These selected cases completed but returned no BFBArchitect score, so older score-only plots hid them.</p><ul>{missing_bits}</ul></section>
<section><h2>Timeout Follow-Up</h2><p class="note">This checks whether cases that timed out for all solvers at 2 or 3 solver threads later produced scores at 16 solver threads.</p><ul>{timeout_bits}</ul></section>
{image_sections}
</main>
</body>
</html>
"""
    (out_dir / "solver_runtime_plots.html").write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicate", type=Path, action="append", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--score-cutoff", type=float, default=2.8,
                        help="BFB score cutoff used for cutoff-stability plots.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.out_dir / "plots"
    df = read_replicates(args.replicate)
    summary = summarize(df)
    summary.to_csv(args.out_dir / "plot_summary.tsv", sep="\t", index=False)
    score_cutoff_summary = write_score_cutoff_summary(
        summary,
        args.out_dir / "score_cutoff_summary.tsv",
        args.score_cutoff,
    )
    timeout_followup = write_timeout_followup(df, args.out_dir / "timeout_16thread_followup.tsv", args.score_cutoff)

    images = {
        "Runtime ECDF Small Multiples": "plots/runtime_ecdf_faceted.png",
        "Median Runtime Heatmap": "plots/runtime_heatmap.png",
        "Absolute Runtime by Thread Count": "plots/absolute_runtime_by_threads.png",
        "Thread Scaling by Solver": "plots/thread_scaling.png",
        "Fastest Solver Fraction by Thread Count": "plots/win_fraction.png",
        "Within 10% of Fastest Fraction by Thread Count": "plots/win_fraction_10pct.png",
        "Pairwise Solver Runtime Ratios": "plots/pairwise_ratios_by_thread.png",
        "BFBArchitect Score Per Amplicon at 3 Solver Threads": "plots/score_by_amplicon_3threads.png",
        "BFBArchitect Score Per Amplicon by Solver and Thread Count": "plots/score_by_amplicon_thread_matrix.png",
        "Per-Amplicon BFBArchitect Score Sensitivity": "plots/score_range_by_amplicon.png",
        "BFB Score Cutoff Stability": "plots/score_cutoff_fraction.png",
        "BFB Score Cutoff Matrix": "plots/score_cutoff_matrix.png",
        "Per-Case Solver Runtime Comparison": "plots/case_solver_bars.png",
        "CBC Slow Cases": "plots/cbc_slow_cases.png",
        "Replicate Stability": "plots/replicate_stability.png",
        "Relative Runtime Ratios": "plots/runtime_ratio_heatmap.png",
    }
    plot_runtime_ecdf_faceted(summary, args.out_dir / images["Runtime ECDF Small Multiples"])
    plot_runtime_heatmap(summary, args.out_dir / images["Median Runtime Heatmap"])
    plot_absolute_runtime_by_threads(summary, args.out_dir / images["Absolute Runtime by Thread Count"])
    plot_thread_scaling(summary, args.out_dir / images["Thread Scaling by Solver"])
    plot_win_fraction(summary, args.out_dir / images["Fastest Solver Fraction by Thread Count"], tolerance=0.0)
    plot_win_fraction(summary, args.out_dir / images["Within 10% of Fastest Fraction by Thread Count"], tolerance=0.10)
    plot_pairwise_ratios_by_thread(summary, args.out_dir / images["Pairwise Solver Runtime Ratios"])
    plot_score_by_amplicon_3threads(summary, args.out_dir / images["BFBArchitect Score Per Amplicon at 3 Solver Threads"])
    plot_score_by_amplicon_thread_matrix(summary, args.out_dir / images["BFBArchitect Score Per Amplicon by Solver and Thread Count"])
    plot_score_range_by_amplicon(summary, args.out_dir / images["Per-Amplicon BFBArchitect Score Sensitivity"])
    plot_score_cutoff_fraction(summary, args.out_dir / images["BFB Score Cutoff Stability"], args.score_cutoff)
    plot_score_cutoff_matrix(summary, args.out_dir / images["BFB Score Cutoff Matrix"], args.score_cutoff)
    plot_case_solver_bars(summary, args.out_dir / images["Per-Case Solver Runtime Comparison"])
    plot_cbc_slow(summary, args.out_dir / images["CBC Slow Cases"])
    plot_replicate_scatter(df, args.out_dir / images["Replicate Stability"])
    plot_speedup_heatmap(summary, args.out_dir / images["Relative Runtime Ratios"])
    write_html(args.out_dir, images, summary, args.replicate, score_cutoff_summary, args.score_cutoff, timeout_followup)
    print(args.out_dir / "solver_runtime_plots.html")


if __name__ == "__main__":
    main()
