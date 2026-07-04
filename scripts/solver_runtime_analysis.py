#!/usr/bin/env python3
"""Inventory and benchmark BFBArchitect solver runtimes on graph-mode cases."""

from __future__ import annotations

import argparse
import csv
import html
import importlib.util
import logging
import multiprocessing as mp
import statistics
import time
from collections import deque
from pathlib import Path

from bfbarchitect.BFBArchitect import (
    CHR_CENTRO,
    reconstruct_bfb_from_graph,
    trim_background_segments,
)
from bfbarchitect.graph_input import (
    find_bfb_candidate_regions,
    subsect_graph_for_region,
    whole_graph_as_region,
)

INVENTORY_FIELDS = [
    "source",
    "case_id",
    "mode",
    "region_index",
    "segments",
    "T",
    "max_cn",
    "foldback_cn",
    "graph",
]

BENCHMARK_FIELDS = INVENTORY_FIELDS[:-1] + [
    "solver",
    "threads",
    "max_active_threads",
    "timeout_s",
    "status",
    "elapsed_s",
    "result_count",
    "score_min",
    "error",
    "graph",
]


def case_id_from_graph(graph: Path) -> str:
    suffix = "_BFB_graph.txt"
    return graph.name[:-len(suffix)] if graph.name.endswith(suffix) else graph.stem


def mode_from_case_id(case_id: str) -> str:
    return "whole_graph" if "_whole_graph" in case_id else "region"


def source_label(graph_dir: Path) -> str:
    return graph_dir.parent.name if graph_dir.name == "bfbarchitect_outputs" else graph_dir.name


def resolve_graph_dir(path: Path) -> Path:
    if path.is_dir() and any(path.glob("*_BFB_graph.txt")):
        return path
    nested = path / "bfbarchitect_outputs"
    if nested.is_dir() and any(nested.glob("*_BFB_graph.txt")):
        return nested
    raise FileNotFoundError(
        f"{path} is neither a BFBArchitect graph directory nor an AC output "
        "directory containing bfbarchitect_outputs/*_BFB_graph.txt"
    )


def graph_cases(graph: Path, source: str) -> list[dict[str, object]]:
    case_id = case_id_from_graph(graph)
    mode = mode_from_case_id(case_id)
    cases = []
    if mode == "whole_graph":
        new_segments, cn, lf, rf, *_ = whole_graph_as_region(
            str(graph),
            centromere_dict=CHR_CENTRO,
            max_primary_segments=100,
            report_skips=False,
            deletion=True,
        )
        data = [(1, new_segments, cn, lf, rf)] if new_segments else []
    else:
        regions = find_bfb_candidate_regions(str(graph), deletion=True)
        region_data = subsect_graph_for_region(
            str(graph),
            regions,
            max_segments=100,
            report_skips=False,
            deletion=True,
        )
        data = []
        for idx, item in enumerate(region_data, 1):
            if item is None:
                continue
            new_segments, cn, lf, rf, *_ = item
            new_segments, cn, lf, rf = trim_background_segments(new_segments, cn, lf, rf)
            if new_segments:
                data.append((idx, new_segments, cn, lf, rf))
    for region_index, _segments, cn, lf, rf in data:
        cases.append(
            {
                "case_id": case_id,
                "source": source,
                "mode": mode,
                "region_index": region_index,
                "segments": len(cn),
                "T": round(max(sum(lf) + sum(rf) + 1, max(cn) if cn else 0)),
                "max_cn": max(cn) if cn else 0,
                "foldback_cn": sum(lf) + sum(rf),
                "graph": str(graph),
            }
        )
    return cases


def build_inventory(graph_dirs: list[Path]) -> list[dict[str, object]]:
    rows = []
    for raw_graph_dir in graph_dirs:
        graph_dir = resolve_graph_dir(raw_graph_dir)
        source = source_label(graph_dir)
        for graph in sorted(graph_dir.glob("*_BFB_graph.txt")):
            rows.extend(graph_cases(graph, source))
    rows.sort(
        key=lambda r: (int(r["T"]), int(r["segments"]), int(r["foldback_cn"])),
        reverse=True,
    )
    return rows


def write_tsv(rows: list[dict[str, object]], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def append_tsv_row(row: dict[str, object], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def solver_available(solver: str) -> bool:
    if solver == "cbc":
        return True
    if solver == "gurobi":
        return importlib.util.find_spec("gurobipy") is not None
    if solver == "mosek":
        return importlib.util.find_spec("mosek") is not None
    return False


def _run_case_worker(case: dict[str, str], solver: str, threads: int, queue: mp.Queue) -> None:
    logging.getLogger("BFBArchitect").addHandler(logging.NullHandler())
    start = time.perf_counter()
    try:
        result = reconstruct_bfb_from_graph(
            case["graph"],
            solver=solver,
            silent=True,
            threads=threads,
            max_graph_segments=100,
        )
        elapsed = time.perf_counter() - start
        scores = [score for item in result for score in item.get("scores", [])]
        queue.put(
            {
                "status": "ok",
                "elapsed_s": f"{elapsed:.3f}",
                "result_count": len(result),
                "score_min": f"{min(scores):.6g}" if scores else "",
                "error": "",
            }
        )
    except Exception as exc:
        elapsed = time.perf_counter() - start
        queue.put(
            {
                "status": "error",
                "elapsed_s": f"{elapsed:.3f}",
                "result_count": 0,
                "score_min": "",
                "error": f"{type(exc).__name__}: {exc}",
            }
        )


def run_with_timeout(case: dict[str, str], solver: str, threads: int, timeout_s: int) -> dict[str, object]:
    if not solver_available(solver):
        return {
            **case,
            "solver": solver,
            "threads": threads,
            "timeout_s": timeout_s,
            "status": "unavailable",
            "elapsed_s": "",
            "result_count": 0,
            "score_min": "",
            "error": f"{solver} Python package is not installed",
        }
    queue: mp.Queue = mp.Queue()
    proc = mp.Process(target=_run_case_worker, args=(case, solver, threads, queue))
    start = time.perf_counter()
    proc.start()
    proc.join(timeout_s)
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        return {
            **case,
            "solver": solver,
            "threads": threads,
            "timeout_s": timeout_s,
            "status": "timeout",
            "elapsed_s": f"{time.perf_counter() - start:.3f}",
            "result_count": 0,
            "score_min": "",
            "error": "",
        }
    if queue.empty():
        return {
            **case,
            "solver": solver,
            "threads": threads,
            "timeout_s": timeout_s,
            "status": "error",
            "elapsed_s": f"{time.perf_counter() - start:.3f}",
            "result_count": 0,
            "score_min": "",
            "error": f"worker exited with code {proc.exitcode}",
        }
    return {
        **case,
        "solver": solver,
        "threads": threads,
        "timeout_s": timeout_s,
        **queue.get(),
    }


def make_unavailable_row(case: dict[str, str], solver: str, threads: int, timeout_s: int) -> dict[str, object]:
    return {
        **case,
        "solver": solver,
        "threads": threads,
        "max_active_threads": "",
        "timeout_s": timeout_s,
        "status": "unavailable",
        "elapsed_s": "",
        "result_count": 0,
        "score_min": "",
        "error": f"{solver} Python package is not installed",
    }


def select_inventory_cases(
    inventory: list[dict[str, str]],
    limit: int | None,
    min_t: int,
    case_ids: set[str] | None = None,
    include_case_ids: set[str] | None = None,
) -> list[dict[str, str]]:
    selected = [row for row in inventory if int(row["T"]) >= min_t]
    if case_ids is not None:
        selected = [row for row in selected if row["case_id"] in case_ids]
    if limit is not None:
        selected = selected[:limit]
    if include_case_ids:
        existing = {(row.get("source", ""), row["case_id"]) for row in selected}
        extras = [
            row for row in inventory
            if row["case_id"] in include_case_ids
            and (row.get("source", ""), row["case_id"]) not in existing
        ]
        selected.extend(extras)
    return selected


def launch_job(case: dict[str, str], solver: str, threads: int) -> dict[str, object]:
    queue: mp.Queue = mp.Queue()
    proc = mp.Process(target=_run_case_worker, args=(case, solver, threads, queue))
    proc.start()
    return {
        "case": case,
        "solver": solver,
        "threads": threads,
        "queue": queue,
        "proc": proc,
        "start": time.perf_counter(),
    }


def finish_job(job: dict[str, object], timeout_s: int, timed_out: bool = False) -> dict[str, object]:
    proc: mp.Process = job["proc"]  # type: ignore[assignment]
    queue: mp.Queue = job["queue"]  # type: ignore[assignment]
    case: dict[str, str] = job["case"]  # type: ignore[assignment]
    solver = str(job["solver"])
    threads = int(job["threads"])
    start = float(job["start"])
    if timed_out:
        if proc.is_alive():
            proc.terminate()
            proc.join(5)
        return {
            **case,
            "solver": solver,
            "threads": threads,
            "max_active_threads": "",
            "timeout_s": timeout_s,
            "status": "timeout",
            "elapsed_s": f"{time.perf_counter() - start:.3f}",
            "result_count": 0,
            "score_min": "",
            "error": "",
        }
    proc.join()
    if queue.empty():
        return {
            **case,
            "solver": solver,
            "threads": threads,
            "max_active_threads": "",
            "timeout_s": timeout_s,
            "status": "error",
            "elapsed_s": f"{time.perf_counter() - start:.3f}",
            "result_count": 0,
            "score_min": "",
            "error": f"worker exited with code {proc.exitcode}",
        }
    return {
        **case,
        "solver": solver,
        "threads": threads,
        "max_active_threads": "",
        "timeout_s": timeout_s,
        **queue.get(),
    }


def benchmark_cases(
    inventory: list[dict[str, str]],
    solvers: list[str],
    thread_list: list[int],
    timeout_s: int,
    limit: int | None,
    min_t: int,
    case_ids: set[str] | None = None,
    include_case_ids: set[str] | None = None,
    benchmark_path: Path | None = None,
    resume_rows: list[dict[str, str]] | None = None,
    max_active_threads: int = 0,
) -> list[dict[str, object]]:
    selected = select_inventory_cases(inventory, limit, min_t, case_ids, include_case_ids)
    rows: list[dict[str, object]] = list(resume_rows or [])
    completed = {
        (row.get("source", ""), row["case_id"], row["solver"], str(row["threads"]))
        for row in rows
    }
    pending = deque()
    for case in selected:
        for threads in thread_list:
            for solver in solvers:
                pending.append((case, solver, threads))

    if max_active_threads <= 0:
        max_active_threads = max(thread_list) if thread_list else 1
    active: list[dict[str, object]] = []
    active_threads = 0

    def record(row: dict[str, object]) -> None:
        row["max_active_threads"] = max_active_threads
        rows.append(row)
        completed.add((str(row.get("source", "")), str(row["case_id"]), str(row["solver"]), str(row["threads"])))
        if benchmark_path is not None:
            append_tsv_row(row, benchmark_path, BENCHMARK_FIELDS)
        print(
            f"{row['case_id']}\tthreads={row['threads']}\t{row['solver']}\t"
            f"{row['status']}\t{row['elapsed_s']}",
            flush=True,
        )

    while pending or active:
        launched = False
        rotations = len(pending)
        for _ in range(rotations):
            if not pending:
                break
            case, solver, threads = pending.popleft()
            key = (case.get("source", ""), case["case_id"], solver, str(threads))
            if key in completed:
                print(f"{case['case_id']}\tthreads={threads}\t{solver}\tskipped\t", flush=True)
                continue
            if not solver_available(solver):
                record(make_unavailable_row(case, solver, threads, timeout_s))
                continue
            if active_threads + threads <= max_active_threads or not active:
                active.append(launch_job(case, solver, threads))
                active_threads += threads
                launched = True
            else:
                pending.append((case, solver, threads))

        now = time.perf_counter()
        remaining = []
        for job in active:
            proc: mp.Process = job["proc"]  # type: ignore[assignment]
            threads = int(job["threads"])
            if proc.is_alive() and now - float(job["start"]) <= timeout_s:
                remaining.append(job)
                continue
            timed_out = proc.is_alive() and now - float(job["start"]) > timeout_s
            record(finish_job(job, timeout_s, timed_out=timed_out))
            active_threads -= threads
        active = remaining
        if active and not launched:
            time.sleep(1)
    return rows


def fmt_float(value: str) -> str:
    if value == "":
        return ""
    return f"{float(value):.2f}"


def render_report(inventory: list[dict[str, str]], benchmark: list[dict[str, str]]) -> str:
    graphs = len({row["graph"] for row in inventory})
    cases = len(inventory)
    sources = sorted({row.get("source", "") for row in inventory})
    mode_counts = {}
    for row in inventory:
        mode_counts[row["mode"]] = mode_counts.get(row["mode"], 0) + 1
    lines = [
        "# BFBArchitect solver runtime analysis",
        "",
        "## Dataset and scope",
        "",
        f"- Dataset source(s): " + ", ".join(source for source in sources if source),
        f"- Reconstructable graph-mode cases: {cases} from {graphs} graph files.",
        f"- Modes: " + ", ".join(f"{mode}={count}" for mode, count in sorted(mode_counts.items())),
        "- Solver thread setting used for benchmarks: 3.",
        "",
        "## Largest graph-mode ILPs",
        "",
        "| case | mode | segments | T | max_cn | foldback_cn |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in inventory[:15]:
        lines.append(
            f"| {row['case_id']} | {row['mode']} | {row['segments']} | {row['T']} | "
            f"{row['max_cn']} | {row['foldback_cn']} |"
        )
    lines.extend(["", "## Benchmark results", ""])
    if not benchmark:
        lines.append("No benchmark rows were provided.")
    else:
        lines.extend(
            [
                "| source | case | mode | T | threads | solver | status | elapsed_s | score_min |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in benchmark:
            lines.append(
                f"| {row.get('source', '')} | {row['case_id']} | {row['mode']} | {row['T']} | {row['threads']} | "
                f"{row['solver']} | {row['status']} | {fmt_float(row['elapsed_s'])} | "
                f"{row['score_min']} |"
            )
        by_solver: dict[str, list[float]] = {}
        for row in benchmark:
            if row["status"] == "ok" and row["elapsed_s"]:
                by_solver.setdefault(row["solver"], []).append(float(row["elapsed_s"]))
        lines.extend(["", "## Timing summary", ""])
        for solver, values in sorted(by_solver.items()):
            lines.append(
                f"- {solver}: n={len(values)}, median={statistics.median(values):.2f}s, "
                f"max={max(values):.2f}s."
            )
        slow_cbc = [
            row for row in benchmark
            if row["solver"] == "cbc"
            and row["status"] in {"timeout", "ok"}
            and row["elapsed_s"]
            and float(row["elapsed_s"]) >= 30
        ]
        lines.append(f"- CBC cases at or above 30s wall time in this run: {len(slow_cbc)}.")
    mosek_rows = [row for row in benchmark if row.get("solver") == "mosek"]
    if not mosek_rows:
        mosek_note = "- The MOSEK code path sets `numThreads`; no MOSEK benchmark rows were requested."
    elif all(row.get("status") == "unavailable" for row in mosek_rows):
        mosek_note = "- The MOSEK code path sets `numThreads`, but MOSEK was unavailable in this environment."
    else:
        mosek_note = "- MOSEK benchmark rows were run through the code path that sets `numThreads`."
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Gurobi uses `m.Params.Threads = max_threads`; CBC uses `PULP_CBC_CMD(..., threads=max_threads)`.",
            mosek_note,
            "- Automatic fallback order in `detect_solver()` is Gurobi, then MOSEK, then CBC. Explicit `--solver` bypasses fallback.",
            "- Timings are wall-clock timings of graph-mode reconstruction in a child process, including graph parsing, ILP model construction, and solve.",
        ]
    )
    return "\n".join(lines) + "\n"


def render_html(inventory: list[dict[str, str]], benchmark: list[dict[str, str]]) -> str:
    def esc(value: object) -> str:
        return html.escape(str(value))

    ok_rows = [
        row for row in benchmark
        if row.get("status") == "ok" and row.get("elapsed_s")
    ]
    max_elapsed = max((float(row["elapsed_s"]) for row in ok_rows), default=1.0)
    slow_cbc = [
        row for row in ok_rows
        if row.get("solver") == "cbc" and float(row["elapsed_s"]) >= 30
    ]

    rows_html = []
    for row in benchmark:
        elapsed = row.get("elapsed_s", "")
        width = 0 if not elapsed else max(1, float(elapsed) / max_elapsed * 100)
        rows_html.append(
            "<tr>"
            f"<td>{esc(row.get('source', ''))}</td>"
            f"<td>{esc(row.get('case_id', ''))}</td>"
            f"<td>{esc(row.get('mode', ''))}</td>"
            f"<td class='num'>{esc(row.get('T', ''))}</td>"
            f"<td class='num'>{esc(row.get('threads', ''))}</td>"
            f"<td>{esc(row.get('solver', ''))}</td>"
            f"<td>{esc(row.get('status', ''))}</td>"
            f"<td class='num'>{fmt_float(elapsed)}</td>"
            f"<td><div class='bar'><span style='width:{width:.1f}%'></span></div></td>"
            f"<td>{esc(row.get('error', ''))}</td>"
            "</tr>"
        )

    inv_rows = []
    for row in inventory[:20]:
        inv_rows.append(
            "<tr>"
            f"<td>{esc(row['source'])}</td><td>{esc(row['case_id'])}</td><td>{esc(row['mode'])}</td>"
            f"<td class='num'>{esc(row['segments'])}</td><td class='num'>{esc(row['T'])}</td>"
            f"<td class='num'>{esc(row['max_cn'])}</td><td class='num'>{esc(row['foldback_cn'])}</td>"
            "</tr>"
        )

    by_solver = {}
    for row in ok_rows:
        by_solver.setdefault(row["solver"], []).append(float(row["elapsed_s"]))
    summary = "".join(
        f"<li><b>{esc(solver)}</b>: n={len(values)}, median={statistics.median(values):.2f}s, "
        f"max={max(values):.2f}s</li>"
        for solver, values in sorted(by_solver.items())
    )
    slow = "".join(
        f"<li>{esc(row.get('source', ''))}: {esc(row['case_id'])}, threads={esc(row['threads'])}: {float(row['elapsed_s']):.2f}s</li>"
        for row in slow_cbc
    ) or "<li>None in this benchmark.</li>"

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>BFBArchitect solver runtime analysis</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 32px; color: #1f2933; }}
h1, h2 {{ margin: 0 0 12px; }}
section {{ margin: 28px 0; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th, td {{ border-bottom: 1px solid #d9e2ec; padding: 6px 8px; text-align: left; vertical-align: top; }}
th {{ background: #f0f4f8; position: sticky; top: 0; }}
.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.bar {{ width: 140px; height: 12px; background: #eef2f7; border-radius: 2px; overflow: hidden; }}
.bar span {{ display: block; height: 100%; background: #2f80ed; }}
.note {{ color: #52606d; }}
</style>
</head>
<body>
<h1>BFBArchitect Solver Runtime Analysis</h1>
<p class="note">Generated from BFBArchitect graph outputs. Timings are wall-clock graph-mode reconstruction runs.</p>
<section>
<h2>Summary</h2>
<ul>
<li>Reconstructable graph-mode cases inventoried: {len(inventory)}</li>
<li>CBC rows at or above 30 seconds: {len(slow_cbc)}</li>
{summary}
</ul>
</section>
<section>
<h2>CBC >= 30 Seconds</h2>
<ul>{slow}</ul>
</section>
<section>
<h2>Benchmark Rows</h2>
<table>
<thead><tr><th>source</th><th>case</th><th>mode</th><th>T</th><th>threads</th><th>solver</th><th>status</th><th>elapsed_s</th><th>relative</th><th>error</th></tr></thead>
<tbody>
{''.join(rows_html)}
</tbody>
</table>
</section>
<section>
<h2>Largest Inventoried Cases</h2>
<table>
<thead><tr><th>source</th><th>case</th><th>mode</th><th>segments</th><th>T</th><th>max_cn</th><th>foldback_cn</th></tr></thead>
<tbody>
{''.join(inv_rows)}
</tbody>
</table>
</section>
</body>
</html>
"""


def write_replicate_comparison(rep1: Path, rep2: Path, out_tsv: Path, out_md: Path) -> None:
    rows1 = read_tsv(rep1)
    rows2 = read_tsv(rep2)
    key_fields = ["source", "case_id", "solver", "threads"]
    idx1 = {tuple(row.get(field, "") for field in key_fields): row for row in rows1}
    idx2 = {tuple(row.get(field, "") for field in key_fields): row for row in rows2}
    keys = sorted(set(idx1) & set(idx2))
    out_rows = []
    for key in keys:
        r1, r2 = idx1[key], idx2[key]
        elapsed1 = float(r1["elapsed_s"]) if r1.get("elapsed_s") else float("nan")
        elapsed2 = float(r2["elapsed_s"]) if r2.get("elapsed_s") else float("nan")
        ratio = elapsed2 / elapsed1 if elapsed1 > 0 else float("nan")
        out_rows.append({
            "source": key[0],
            "case_id": key[1],
            "solver": key[2],
            "threads": key[3],
            "status_1": r1.get("status", ""),
            "elapsed_1": f"{elapsed1:.3f}" if elapsed1 == elapsed1 else "",
            "status_2": r2.get("status", ""),
            "elapsed_2": f"{elapsed2:.3f}" if elapsed2 == elapsed2 else "",
            "replicate2_over_1": f"{ratio:.3f}" if ratio == ratio else "",
            "abs_delta_s": f"{abs(elapsed2 - elapsed1):.3f}" if elapsed1 == elapsed1 and elapsed2 == elapsed2 else "",
            "T": r1.get("T", ""),
            "mode": r1.get("mode", ""),
        })
    fields = [
        "source", "case_id", "mode", "T", "solver", "threads", "status_1", "elapsed_1",
        "status_2", "elapsed_2", "replicate2_over_1", "abs_delta_s",
    ]
    write_tsv(out_rows, out_tsv, fields)

    ok = [
        row for row in out_rows
        if row["status_1"] == "ok" and row["status_2"] == "ok" and row["elapsed_1"] and row["elapsed_2"]
    ]
    ratios = [float(row["replicate2_over_1"]) for row in ok if float(row["elapsed_1"]) >= 1.0]
    slow = sorted(ok, key=lambda row: float(row["abs_delta_s"]), reverse=True)[:12]
    lines = [
        "# Solver Runtime Replicate Comparison",
        "",
        f"- Replicate 1: `{rep1}`",
        f"- Replicate 2: `{rep2}`",
        f"- Paired rows: {len(out_rows)}",
        f"- Paired ok rows: {len(ok)}",
    ]
    if ratios:
        lines.append(f"- Median replicate2/replicate1 ratio for rows with replicate 1 >=1s: {statistics.median(ratios):.3f}")
    lines.extend([
        "",
        "## Largest Absolute Differences",
        "",
        "| source | case | solver | threads | rep1_s | rep2_s | ratio | abs_delta_s |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in slow:
        lines.append(
            f"| {row['source']} | {row['case_id']} | {row['solver']} | {row['threads']} | "
            f"{row['elapsed_1']} | {row['elapsed_2']} | {row['replicate2_over_1']} | {row['abs_delta_s']} |"
        )
    out_md.write_text("\n".join(lines) + "\n")


def write_replicate_summary(rep_paths: list[Path], out_tsv: Path, out_md: Path) -> None:
    replicate_rows = [read_tsv(path) for path in rep_paths]
    key_fields = ["source", "case_id", "solver", "threads"]
    indexes = [
        {tuple(row.get(field, "") for field in key_fields): row for row in rows}
        for rows in replicate_rows
    ]
    keys = sorted(set.intersection(*(set(index) for index in indexes)))
    out_rows = []
    for key in keys:
        paired = [index[key] for index in indexes]
        vals = [float(row["elapsed_s"]) for row in paired if row.get("status") == "ok" and row.get("elapsed_s")]
        if len(vals) != len(rep_paths):
            continue
        mean = statistics.mean(vals)
        stdev = statistics.stdev(vals) if len(vals) > 1 else 0.0
        out_rows.append({
            "source": key[0],
            "case_id": key[1],
            "mode": paired[0].get("mode", ""),
            "T": paired[0].get("T", ""),
            "solver": key[2],
            "threads": key[3],
            "n": len(vals),
            "median_s": f"{statistics.median(vals):.3f}",
            "mean_s": f"{mean:.3f}",
            "min_s": f"{min(vals):.3f}",
            "max_s": f"{max(vals):.3f}",
            "stdev_s": f"{stdev:.3f}",
            "cv": f"{(stdev / mean):.3f}" if mean else "",
            "replicate_values_s": ",".join(f"{v:.3f}" for v in vals),
        })
    fields = [
        "source", "case_id", "mode", "T", "solver", "threads", "n",
        "median_s", "mean_s", "min_s", "max_s", "stdev_s", "cv", "replicate_values_s",
    ]
    write_tsv(out_rows, out_tsv, fields)

    slow_cbc = sorted(
        [row for row in out_rows if row["solver"] == "cbc" and float(row["median_s"]) >= 30],
        key=lambda row: float(row["median_s"]),
        reverse=True,
    )
    by_solver = {}
    for row in out_rows:
        by_solver.setdefault(row["solver"], []).append(float(row["median_s"]))
    lines = [
        "# Solver Runtime Three-Replicate Summary",
        "",
        "- Replicate files:",
    ]
    lines.extend(f"  - `{path}`" for path in rep_paths)
    lines.extend([
        f"- Complete paired rows: {len(out_rows)}",
        "",
        "## Median Runtime By Solver",
        "",
    ])
    for solver, vals in sorted(by_solver.items()):
        lines.append(
            f"- {solver}: n={len(vals)}, median of case medians={statistics.median(vals):.3f}s, "
            f"max case median={max(vals):.3f}s"
        )
    lines.extend([
        "",
        "## CBC Cases With Median >= 30s",
        "",
        "| source | case | T | median_s | min_s | max_s | cv | replicates_s |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in slow_cbc:
        lines.append(
            f"| {row['source']} | {row['case_id']} | {row['T']} | {row['median_s']} | "
            f"{row['min_s']} | {row['max_s']} | {row['cv']} | {row['replicate_values_s']} |"
        )
    lines.extend([
        "",
        "## Largest CV Rows",
        "",
        "| source | case | solver | median_s | cv | replicates_s |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for row in sorted(out_rows, key=lambda r: float(r["cv"] or 0), reverse=True)[:15]:
        lines.append(
            f"| {row['source']} | {row['case_id']} | {row['solver']} | "
            f"{row['median_s']} | {row['cv']} | {row['replicate_values_s']} |"
        )
    out_md.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", type=Path, action="append", default=None,
                        help=(
                            "Directory containing *_BFB_graph.txt files, or an AC output directory containing "
                            "bfbarchitect_outputs/. Can be supplied multiple times."
                        ))
    parser.add_argument("--out-dir", type=Path, default=Path("reports/solver_runtime_analysis"))
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--thread-list", type=int, nargs="+", default=None)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--solvers", nargs="+", default=["gurobi", "mosek", "cbc"])
    parser.add_argument("--limit", type=int, default=15,
                        help="Benchmark the top N selected cases after sorting by inferred ILP size.")
    parser.add_argument("--all-cases", action="store_true",
                        help="Benchmark every inventoried case passing --min-t instead of applying --limit.")
    parser.add_argument("--min-t", type=int, default=8,
                        help="Minimum inferred T value to include in the benchmark selection.")
    parser.add_argument("--case-id", action="append", default=None,
                        help="Restrict benchmark to a case id. Can be supplied multiple times.")
    parser.add_argument("--include-case-id", action="append", default=None,
                        help="Always include this case id in addition to the normal limit/min-t selection.")
    parser.add_argument("--include-case-id-file", type=Path, default=None,
                        help="File containing one additional case id per line to include beyond the normal selection.")
    parser.add_argument("--benchmark-file", type=Path, default=None)
    parser.add_argument("--resume", action="store_true",
                        help="Append to the benchmark TSV and skip completed source/case/solver/thread rows.")
    parser.add_argument("--max-active-threads", type=int, default=0,
                        help="Maximum sum of solver thread counts to run concurrently. 0 means sequential.")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--compare-replicates", nargs=2, type=Path, metavar=("REP1", "REP2"),
                        help="Write paired replicate comparison files and exit after inventory.")
    parser.add_argument("--summarize-replicates", nargs="+", type=Path, metavar="REP",
                        help="Write multi-replicate summary files and exit after inventory.")
    parser.add_argument("--inventory-only", action="store_true")
    args = parser.parse_args()

    if not args.graph_dir:
        parser.error("at least one --graph-dir containing *_BFB_graph.txt files is required")
    graph_dirs = args.graph_dir
    inventory = build_inventory(graph_dirs)
    inventory_path = args.out_dir / "case_inventory.tsv"
    write_tsv(inventory, inventory_path, INVENTORY_FIELDS)
    print(f"wrote {inventory_path} ({len(inventory)} cases)")

    if args.compare_replicates:
        write_replicate_comparison(
            args.compare_replicates[0],
            args.compare_replicates[1],
            args.out_dir / "replicate_comparison.tsv",
            args.out_dir / "replicate_comparison.md",
        )
        print(f"wrote {args.out_dir / 'replicate_comparison.tsv'}")
        print(f"wrote {args.out_dir / 'replicate_comparison.md'}")
        return

    if args.summarize_replicates:
        write_replicate_summary(
            args.summarize_replicates,
            args.out_dir / "replicate_summary.tsv",
            args.out_dir / "replicate_summary.md",
        )
        print(f"wrote {args.out_dir / 'replicate_summary.tsv'}")
        print(f"wrote {args.out_dir / 'replicate_summary.md'}")
        return

    benchmark = []
    benchmark_path = args.benchmark_file or (args.out_dir / "benchmark.tsv")
    if args.report_only:
        if benchmark_path.exists():
            benchmark = read_tsv(benchmark_path)
        else:
            print(f"benchmark file not found; generating inventory-only report: {benchmark_path}")
    elif not args.inventory_only:
        thread_list = args.thread_list or [args.threads]
        resume_rows = read_tsv(benchmark_path) if args.resume and benchmark_path.exists() else []
        if benchmark_path.exists() and not args.resume:
            benchmark_path.unlink()
        include_case_ids = set(args.include_case_id or [])
        if args.include_case_id_file:
            include_case_ids.update(
                line.strip() for line in args.include_case_id_file.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            )
        inventory_for_benchmark = [{field: str(row[field]) for field in INVENTORY_FIELDS} for row in inventory]
        selected_cases = select_inventory_cases(
            inventory_for_benchmark,
            None if args.all_cases else args.limit,
            args.min_t,
            set(args.case_id) if args.case_id else None,
            include_case_ids=include_case_ids,
        )
        selection_path = args.out_dir / "benchmark_selection.tsv"
        write_tsv(selected_cases, selection_path, INVENTORY_FIELDS)
        print(f"wrote {selection_path} ({len(selected_cases)} selected cases)")
        benchmark = benchmark_cases(
            inventory_for_benchmark,
            args.solvers,
            thread_list,
            args.timeout,
            None if args.all_cases else args.limit,
            args.min_t,
            set(args.case_id) if args.case_id else None,
            include_case_ids=include_case_ids,
            benchmark_path=benchmark_path,
            resume_rows=resume_rows,
            max_active_threads=args.max_active_threads,
        )
        if benchmark_path is None:
            write_tsv(benchmark, benchmark_path, BENCHMARK_FIELDS)

    report = render_report(
        [{field: str(row[field]) for field in INVENTORY_FIELDS} for row in inventory],
        [{field: str(row.get(field, "")) for field in BENCHMARK_FIELDS} for row in benchmark],
    )
    report_path = args.out_dir / "solver_runtime_analysis.md"
    report_path.write_text(report)
    print(f"wrote {report_path}")
    html_report = render_html(
        [{field: str(row[field]) for field in INVENTORY_FIELDS} for row in inventory],
        [{field: str(row.get(field, "")) for field in BENCHMARK_FIELDS} for row in benchmark],
    )
    html_path = args.out_dir / "solver_runtime_analysis.html"
    html_path.write_text(html_report)
    print(f"wrote {html_path}")


if __name__ == "__main__":
    main()
