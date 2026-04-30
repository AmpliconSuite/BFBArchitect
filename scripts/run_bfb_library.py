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
    find_bfb_candidate_regions, 
    subsect_graph_for_region, 
    reconstruct_bfb, 
    write_bfb_graph, 
    write_bfb_cycles,
    CHR_CENTRO
)

def run_bfb_library(graph_file, output_prefix, multiple=False, solver=None):
    """
    Demonstrate how to use the BFBArchitect library API to reconstruct BFB sequences
    from an AA-format _graph.txt file.
    """
    print(f"--- BFBArchitect Library API ---")
    print(f"Input graph: {graph_file}")
    
    if not os.path.exists(graph_file):
        print(f"Error: Graph file {graph_file} not found.")
        return

    # 1. Detect candidate BFB regions in the graph
    regions = find_bfb_candidate_regions(graph_file)
    print(f"Detected {len(regions)} candidate region(s): {regions}")
    
    # 2. Extract and pre-process segment data for each region
    region_data = subsect_graph_for_region(graph_file, regions)
    
    for i, (region, data) in enumerate(zip(regions, region_data)):
        if data is None:
            print(f"Skipping region {i+1}: {region} (no data extracted)")
            continue
        
        new_segments, cn, lf, rf, region_svs, sv_info = data
        chrom = region[0]
        
        print(f"\nProcessing region {i+1}: {region}")
        
        # 3. Reconstruct BFB strings using the ILP solver
        # Arguments: (segments, cn_vector, lf_vector, rf_vector, centromere_pos, solver, multiple)
        BFB_strings, scores, multiplicity = reconstruct_bfb(
            new_segments, 
            cn, 
            lf, 
            rf, 
            CHR_CENTRO.get(chrom, 0),
            solver=solver,
            multiple=multiple
        )
        
        # 4. Save results
        region_prefix = f"{output_prefix}_region{i+1}"
        write_bfb_graph(f"{region_prefix}_graph.txt", new_segments, region_svs, sv_info)
        write_bfb_cycles(f"{region_prefix}_cycles.txt", new_segments, BFB_strings, scores, multiplicity)
        
        print(f"Found {len(BFB_strings)} BFB candidate(s).")
        print(f"Results written to: {region_prefix}_graph.txt and {region_prefix}_cycles.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Invoke BFBArchitect library API on a graph file.")
    parser.add_argument("graph", help="Path to AA-format _graph.txt file.")
    parser.add_argument("--output_prefix", help="Prefix for output files.", default="bfb_output")
    parser.add_argument("--multiple", action="store_true", help="Reconstruct multiple candidates.")
    parser.add_argument("--solver", help="Solver to use (gurobi or cbc).", default=None)
    
    args = parser.parse_args()
    
    run_bfb_library(args.graph, args.output_prefix, args.multiple, args.solver)
