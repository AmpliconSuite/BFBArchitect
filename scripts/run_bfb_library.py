import sys
import os
import argparse

# Ensure we can import bfbarchitect
# If not installed via pip, add the project root to sys.path
try:
    import bfbarchitect
except ImportError:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

from bfbarchitect import (
    reconstruct_bfb_from_graph,
    write_bfb_graph,
    write_bfb_cycles,
    visualize_BFB,
)

def _parse_region(region_str):
    chrom, rest = region_str.split(':')
    start, end = rest.split('-')
    return (chrom, int(start), int(end))

def run_bfb_library(graph_file, output_prefix, whole_graph=False, region=None,
                    multiple=False, solver=None, verbose=False,
                    max_graph_segments=100, max_whole_graph_segments=None,
                    reverse_polarity=False):
    """
    Demonstrate how to use the BFBArchitect library API to reconstruct BFB sequences
    from an AA-format _graph.txt file.
    """
    print(f"--- BFBArchitect Library API ---")
    print(f"Input graph: {graph_file}")
    if max_whole_graph_segments is not None:
        max_graph_segments = max_whole_graph_segments

    if not os.path.exists(graph_file):
        print(f"Error: Graph file {graph_file} not found.")
        return

    results = reconstruct_bfb_from_graph(
        graph_file,
        whole_graph=whole_graph,
        region=region,
        solver=solver,
        multiple=multiple,
        verbose=verbose,
        max_graph_segments=max_graph_segments,
        reverse_polarity=reverse_polarity,
    )

    for i, res in enumerate(results):
        region_prefix = output_prefix if (whole_graph or region is not None) else f"{output_prefix}_region{i+1}"
        graph_out = f"{region_prefix}_graph.txt"
        cycle_out = f"{region_prefix}_cycles.txt"
        
        write_bfb_graph(graph_out, res['new_segments'], res['svs'], res['sv_info'])
        write_bfb_cycles(cycle_out, res['new_segments'], res['bfb_strings'], res['scores'], res['multiplicity'])
        print(f"  Written: {graph_out}, {cycle_out}")

        # Optional visualization
        print(f"  Visualizing {region_prefix}...")
        visualize_BFB(
            cycle_file=cycle_out,
            graph_file=graph_out,
            cnr_file=None,
            output_prefix=f"{region_prefix}_BFB",
            multiple=multiple
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Invoke BFBArchitect library API on a graph file.")
    parser.add_argument("graph", help="Path to AA-format _graph.txt file.")
    parser.add_argument("--output_prefix", help="Prefix for output files.", default="bfb_output")
    parser.add_argument("--whole_graph", action="store_true", help="Treat all segments as one region.")
    parser.add_argument("--region", help="Process a specific region only (chr:start-end).", default=None)
    parser.add_argument("--multiple", action="store_true", help="Reconstruct multiple candidates.")
    parser.add_argument("--reverse_polarity", action="store_true",
                        help="Run the opposite of the computed BFB polarity.")
    parser.add_argument("--solver", help="Solver to use (gurobi or cbc).", default=None)
    parser.add_argument("--verbose", action="store_true", help="Print per-step segment transforms and CN/LF/RF vectors.")
    parser.add_argument("--max-graph-segments", "--max-whole-graph-segments",
                        type=int, default=100, dest="max_graph_segments",
                        help="Maximum number of graph segments allowed per graph-mode region "
                             "(default: 100). Use 0 or a negative value "
                             "to disable this cutoff.")

    args = parser.parse_args()

    if args.whole_graph and args.region:
        parser.error("--whole_graph and --region are mutually exclusive.")

    parsed_region = _parse_region(args.region) if args.region else None
    max_graph_segments = args.max_graph_segments
    if max_graph_segments is not None and max_graph_segments <= 0:
        max_graph_segments = None

    run_bfb_library(args.graph, args.output_prefix, args.whole_graph, parsed_region,
                    args.multiple, args.solver, args.verbose, max_graph_segments,
                    reverse_polarity=args.reverse_polarity)
