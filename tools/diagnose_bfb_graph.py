#!/usr/bin/env python3
"""Diagnose BFB evidence in an AmpliconArchitect graph file.

This is an agent-neutral helper for quickly reproducing the graph inspection
workflow used during BFBArchitect development.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from bfbarchitect.BFBArchitect import (  # noqa: E402
    CHR_CENTRO,
    _reverse_polarity_vectors,
    reconstruct_bfb,
    reconstruct_bfb_from_graph,
    trim_background_segments,
)
from bfbarchitect.graph_input import (  # noqa: E402
    _foldback_flank_cn_change,
    _foldback_passes_recon_cn_filter,
    find_bfb_candidate_regions,
    parse_graph_file,
    subsect_graph_for_region,
    whole_graph_as_region,
    write_tst_report,
)
import bfbarchitect.graph_input as graph_input  # noqa: E402


def _parse_region(region_text: str) -> tuple[str, int, int]:
    chrom, coords = region_text.split(":", 1)
    start, end = coords.replace(",", "").split("-", 1)
    return chrom, int(start), int(end)


def _print_header(title: str) -> None:
    print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")


def _print_segments(chrom_segs: dict, chrom: str | None = None) -> None:
    chroms = [chrom] if chrom else sorted(chrom_segs)
    for cur_chrom in chroms:
        if cur_chrom not in chrom_segs:
            continue
        print(f"\n{cur_chrom} sequence segments ({len(chrom_segs[cur_chrom])})")
        for i, (start, end, cn, cov, rc) in enumerate(chrom_segs[cur_chrom]):
            print(
                f"  {i:>3} {cur_chrom}:{start}-{end}"
                f"  size={end - start:>9}  CN={cn:8.3f}  cov={cov:8.1f}  reads={rc}"
            )


def _print_raw_svs(svs, chrom: str | None, broad_fb_cutoff: int) -> None:
    _print_header("Raw SVs touching selected chromosome" if chrom else "Raw SVs")
    for sv, cn, rc in svs:
        if chrom and sv.chrom1 != chrom and sv.chrom2 != chrom:
            continue
        same_strand = sv.chrom1 == sv.chrom2 and sv.strand1 == sv.strand2
        is_fb_default = sv.is_foldback()
        is_fb_broad = sv.is_foldback(max_distance=broad_fb_cutoff)
        fb_label = ""
        if is_fb_default:
            fb_label = "FOLD-BACK"
        elif same_strand and is_fb_broad:
            fb_label = f"BROAD-FOLD-BACK<={broad_fb_cutoff}"
        print(
            f"  {sv}  CN={cn:.3f}  rc={rc:<5} type={sv.get_SV_type():<3}"
            f" {fb_label}"
        )


def _print_foldback_filters(svs, chrom_segs, chrom: str | None, broad_fb_cutoff: int) -> None:
    _print_header("Foldback candidates and filters")
    for sv, cn, rc in svs:
        if chrom and sv.chrom1 != chrom and sv.chrom2 != chrom:
            continue
        same_strand = sv.chrom1 == sv.chrom2 and sv.strand1 == sv.strand2
        if not same_strand:
            continue
        if not sv.is_foldback(max_distance=broad_fb_cutoff):
            continue
        dist = abs(sv.bp2 - sv.bp1)
        delta = _foldback_flank_cn_change(sv, chrom_segs)
        passes = _foldback_passes_recon_cn_filter(sv, chrom_segs)
        print(
            f"  {sv}  strand={sv.strand1}  dist={dist}  CN={cn:.3f}  rc={rc}"
            f"  flank_delta={'NA' if delta is None else f'{delta:.3f}'}"
            f"  pass_recon_filter={passes}"
            f"  default50={sv.is_foldback()}"
        )


def _run_subsection(graph_file: str, regions, fb_cutoff: int, disable_tst: bool,
                    verbose: bool, deletion: bool):
    original = graph_input.find_tst_foldbacks
    if disable_tst:
        graph_input.find_tst_foldbacks = lambda svs, chrom_segs, **kwargs: svs
    try:
        return subsect_graph_for_region(
            graph_file,
            regions,
            fb_dist_cut=fb_cutoff,
            verbose=verbose,
            max_segments=None,
            deletion=deletion,
        )
    finally:
        graph_input.find_tst_foldbacks = original


def _print_region_result(
    graph_file: str, region, data, fb_cutoff: int, solver: str, reverse_polarity: bool
) -> None:
    print(f"\nRegion {region[0]}:{region[1]}-{region[2]}  fb_dist_cut={fb_cutoff}")
    if data is None:
        print("  No region data")
        return

    new_segments, cn, lf, rf, region_svs, sv_info = data
    for i, seg in enumerate(new_segments):
        chrom, start, end, cn_float, cov, rc = seg
        print(
            f"  seg {i:>2} {chrom}:{start}-{end}"
            f"  CN_float={cn_float:8.3f}  cn={cn[i]:>4}  lf={lf[i]:>4}  rf={rf[i]:>4}"
        )
    print(f"  cn={cn}")
    print(f"  lf={lf}")
    print(f"  rf={rf}")
    print("  foldbacks:")
    for sv in region_svs:
        if sv.is_foldback(max_distance=fb_cutoff):
            cn_rc = sv_info.get(sv)
            print(
                f"    {sv}  info={cn_rc}"
                f"  TST={getattr(sv, 'TST', False)}"
                f"  local_resolved={getattr(sv, 'local_resolved', False)}"
            )

    trim_segments, trim_cn, trim_lf, trim_rf = trim_background_segments(
        new_segments, cn, lf, rf
    )
    if not trim_segments:
        print("  no segments after trimming")
        return

    bfb_strings, scores, multiplicity = reconstruct_bfb(
        trim_segments,
        trim_cn,
        trim_lf,
        trim_rf,
        CHR_CENTRO.get(region[0], 0),
        solver=solver,
        silent=True,
        reverse_polarity=reverse_polarity,
    )
    print(f"  trimmed_cn={trim_cn}")
    print(f"  trimmed_lf={trim_lf}")
    print(f"  trimmed_rf={trim_rf}")
    if reverse_polarity:
        reverse_cn, reverse_lf, reverse_rf = _reverse_polarity_vectors(
            trim_cn, trim_lf, trim_rf
        )
        print(f"  reverse_cn={reverse_cn}")
        print(f"  reverse_lf={reverse_lf}")
        print(f"  reverse_rf={reverse_rf}")
    polarity_label = " reverse_polarity" if reverse_polarity else ""
    print(f"  multiplicity={multiplicity}")
    print(f"  scores{polarity_label}={scores}")
    print(f"  bfb_strings={bfb_strings}")


def _print_tst_report(graph_file: str, regions) -> None:
    _print_header("TST report")
    buf_path = Path("/tmp") / f"{Path(graph_file).stem}.tst_report.txt"
    write_tst_report(graph_file, str(buf_path), bfb_regions=regions)
    print(buf_path.read_text())


def _print_whole_graph(
    graph_file: str, solver: str, verbose: bool, deletion: bool, reverse_polarity: bool
) -> None:
    _print_header("Whole-graph reconstruction")
    new_segments, cn, lf, rf, svs, sv_info, chrom = whole_graph_as_region(
        graph_file,
        verbose=verbose,
        max_primary_segments=None,
        deletion=deletion,
    )
    if not new_segments:
        print("No whole-graph region returned")
        return
    region = (chrom, new_segments[0][1], new_segments[-1][2])
    print(f"primary_region={region}")
    for i, seg in enumerate(new_segments):
        seg_chrom, start, end, cn_float, cov, rc = seg
        print(
            f"  seg {i:>2} {seg_chrom}:{start}-{end}"
            f"  CN_float={cn_float:8.3f}  cn={cn[i]:>4}  lf={lf[i]:>4}  rf={rf[i]:>4}"
        )
    bfb_strings, scores, multiplicity = reconstruct_bfb(
        new_segments,
        cn,
        lf,
        rf,
        CHR_CENTRO.get(chrom, 0),
        solver=solver,
        silent=True,
        reverse_polarity=reverse_polarity,
    )
    print(f"cn={cn}")
    print(f"lf={lf}")
    print(f"rf={rf}")
    if reverse_polarity:
        reverse_cn, reverse_lf, reverse_rf = _reverse_polarity_vectors(cn, lf, rf)
        print(f"reverse_cn={reverse_cn}")
        print(f"reverse_lf={reverse_lf}")
        print(f"reverse_rf={reverse_rf}")
    polarity_label = " reverse_polarity" if reverse_polarity else ""
    print(f"multiplicity={multiplicity}")
    print(f"scores{polarity_label}={scores}")
    print(f"bfb_strings={bfb_strings}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Diagnose BFB evidence in an AA graph file."
    )
    parser.add_argument("graph", help="AA graph file")
    parser.add_argument("--chrom", help="Restrict raw summaries to one chromosome")
    parser.add_argument(
        "--region",
        action="append",
        help="Region to evaluate, e.g. chr7:52926926-55529670. May repeat.",
    )
    parser.add_argument(
        "--fb-cutoff",
        type=int,
        action="append",
        default=None,
        help="Foldback distance cutoff to evaluate. May repeat. Default: 50000.",
    )
    parser.add_argument(
        "--broad-fb-cutoff",
        type=int,
        default=100000,
        help="Distance cutoff used only for raw broad-foldback diagnostics.",
    )
    parser.add_argument("--no-tst", action="store_true", help="Disable TST injection")
    parser.add_argument("--whole-graph", action="store_true", help="Run whole graph mode")
    parser.set_defaults(deletion=True)
    parser.add_argument("--deletion", dest="deletion", action="store_true",
                        help="Apply graph deletion-edge CN correction (default; retained for compatibility)")
    parser.add_argument("--no-deletion", dest="deletion", action="store_false",
                        help="Disable graph deletion-edge CN correction")
    parser.add_argument("--show-tst-report", action="store_true", help="Print TST report")
    parser.add_argument("--show-segments", action="store_true", help="Print raw sequence segments")
    parser.add_argument("--verbose", action="store_true", help="Show verbose segmentation trace")
    parser.add_argument("--solver", default="cbc", choices=["cbc", "gurobi", "mosek"])
    parser.add_argument(
        "--reverse_polarity",
        action="store_true",
        help="Run the opposite of the computed BFB polarity.",
    )
    args = parser.parse_args()

    graph_file = str(Path(args.graph).expanduser())
    svs, chrom_segs = parse_graph_file(graph_file)
    regions = [_parse_region(r) for r in args.region] if args.region else find_bfb_candidate_regions(
        graph_file, deletion=args.deletion
    )
    cutoffs = args.fb_cutoff or [50000]

    _print_header("Graph summary")
    print(f"graph={graph_file}")
    print(f"chrom_segment_counts={{{', '.join(f'{c}: {len(v)}' for c, v in sorted(chrom_segs.items()))}}}")
    print("candidate_regions:")
    for region in find_bfb_candidate_regions(graph_file, deletion=args.deletion):
        print(f"  {region[0]}:{region[1]}-{region[2]}")
    print("selected_regions:")
    for region in regions:
        print(f"  {region[0]}:{region[1]}-{region[2]}")

    if args.show_segments:
        _print_header("Raw sequence segments")
        _print_segments(chrom_segs, args.chrom)

    _print_raw_svs(svs, args.chrom, args.broad_fb_cutoff)
    _print_foldback_filters(svs, chrom_segs, args.chrom, args.broad_fb_cutoff)

    if args.show_tst_report:
        _print_tst_report(graph_file, regions)

    for cutoff in cutoffs:
        label = f"Region reconstruction fb_dist_cut={cutoff}"
        if args.no_tst:
            label += " no_TST"
        if args.reverse_polarity:
            label += " reverse_polarity"
        _print_header(label)
        with contextlib.redirect_stderr(io.StringIO()):
            data = _run_subsection(
                graph_file,
                regions,
                cutoff,
                disable_tst=args.no_tst,
                verbose=args.verbose,
                deletion=args.deletion,
            )
        for region, region_data in zip(regions, data):
            _print_region_result(
                graph_file,
                region,
                region_data,
                cutoff,
                args.solver,
                args.reverse_polarity,
            )

    if args.whole_graph:
        _print_whole_graph(
            graph_file,
            args.solver,
            args.verbose,
            args.deletion,
            args.reverse_polarity,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
