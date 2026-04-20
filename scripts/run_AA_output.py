import os
import sys
import argparse
import time
import re

try:
    from src.datatypes import SV, CHR_CENTRO
    from src.BFBSolver import reconstruct_BFB_string, reconstruct_BFB_strings, check_BFB_string, print_BFB_string
    from src.utils import create_logger
    from src.BFBVisualizer import visualize_BFB
except ImportError:
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from src.datatypes import SV, CHR_CENTRO
    from src.BFBSolver import reconstruct_BFB_string, reconstruct_BFB_strings, check_BFB_string, print_BFB_string
    from src.utils import create_logger
    from src.BFBVisualizer import visualize_BFB

def compute_BFB_scores(new_segments, cn0, lf0, rf0, BFB_strings, multiplicity, logger):
    scores = []
    for idx, BFB_string in enumerate(BFB_strings):
        logger.info('----------------------------------------')
        if check_BFB_string(BFB_string) == False:
            print(f'{print_BFB_string(BFB_string)} is not a valid BFB sequence.')
            continue
        logger.info(f'BFB string {idx+1}: {print_BFB_string(BFB_string, print_to_console=False)}')
        logger.info(f'Scores for BFB string {idx+1}')
        # Discrepency between observed and predicted CN
        cn = [0 for _ in range(len(cn0))]
        lf = [0 for _ in range(len(lf0))]
        rf = [0 for _ in range(len(rf0))]
        for i in range(len(BFB_string) - 1):
            seg1, seg2 = BFB_string[i], BFB_string[i+1]
            cn[abs(seg1)-1] += 1
            if abs(seg1) == abs(seg2):
                if seg1 > 0 and seg2 < 0:
                    rf[abs(seg1)-1] += 1
                elif seg1 < 0 and seg2 > 0:
                    lf[abs(seg1)-1] += 1
        cn[abs(BFB_string[-1])-1] += 1
        cn = [c * multiplicity for c in cn]
        lf = [c * multiplicity for c in lf]
        rf = [c * multiplicity for c in rf]

        CN_score = sum([abs(cn0[i]-cn[i])/cn[i] for i in range(len(cn))])
        logger.info(f'CN discrepancy: {CN_score}')
        fb_dist = sum([(lf0[i]-lf[i])**2 + (rf0[i]-rf[i])**2 for i in range(len(cn0))])**0.5 / len(cn0)
        logger.info(f'Foldback Euclidean distance: {fb_dist}')
        missing_fb_score = 0
        for i in range(len(cn0)):
            if (lf0[i] == 0 and lf[i] != 0):
                missing_fb_score += 0.5 * lf[i]
            if (rf0[i] == 0 and rf[i] != 0):
                missing_fb_score += 0.5 * rf[i]
        # fb_count = sum([1 for i in range(len(cn0)) if (lf0[i] > 0 or rf0[i] > 0)])
        fb_count = sum([1 for i in range(len(lf0)) if lf0[i] > 0]) + sum([1 for i in range(len(rf0)) if rf0[i] > 0])
        if fb_count < 2 or len(cn0) < 2:
            missing_fb_score += 2
        logger.info(f'Missing foldback score: {missing_fb_score}')
        observed_cn = [seg[3] - 1 for seg in new_segments]
        sequencing_error = sum([abs(observed_cn[i]-cn0[i])/cn0[i] for i in range(len(cn0)) if cn0[i] > 0])
        logger.info(f'CN estimation error of sequencing data: {sequencing_error}')
        score = CN_score + fb_dist + missing_fb_score + sequencing_error
        logger.info(f'Total score: {score}')
        scores.append(score)
    return scores

def generate_graph_file(output_fn, new_segments, SVs, sv_fino):
    # Generate graph.txt
    out_file = open(output_fn, 'w')
    out_file.write('SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberReadsMapped\n')
    for i, seg in enumerate(new_segments):
        size = seg[2] - seg[1] + 1
        entry = f'sequence	{seg[0]}:{seg[1]}-	{seg[0]}:{seg[2]}+	{seg[3]}	{seg[4]}	{size}	{seg[5]}\n'
        out_file.write(entry)
    out_file.write('BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs\n')
    for i in range(1, len(new_segments)):
        seg1, seg2 = new_segments[i-1], new_segments[i]
        if seg1[0] != seg2[0] or seg1[2]+1 != seg2[1]: 
            continue
        cn = min(new_segments[i-1][3], new_segments[i][3])
        read_count = int((new_segments[i-1][5]+new_segments[i][5])/2)
        entry = f'concordant	{seg1[0]}:{seg1[2]}+->{seg2[0]}:{seg2[1]}-	{cn}	{read_count}\n'
        out_file.write(entry)
    for sv in SVs:
        entry = f'discordant	{sv}	{sv_fino[sv][0]}	{sv_fino[sv][1]}\n'
        out_file.write(entry)
    out_file.close()

def generate_cycle_file(output_fn, new_segments, BFB_strings, scores, multiplicity):
    segments = [(seg[0], seg[1], seg[2]) for seg in new_segments]
    # Find all intervals in the BFB amplicon
    intervals = []
    (chr, start, end) = segments[0]
    for i in range(0, len(segments)-1):
        (next_chr, next_start, next_end) = segments[i+1]
        if next_chr == chr and next_start == end + 1:
            end = next_end
        else:
            intervals.append((chr, start, end))
            chr, start, end = next_chr, next_start, next_end
    intervals.append((chr, start, end))
    # Output cycles.txt
    out_file = open(output_fn, 'w')
    for i, interval in enumerate(intervals):
        (chr, start, end) = interval
        out_file.write(f'Interval	{i+1}	{chr}	{start}	{end}\n')
    out_file.write('List of cycle segments\n')
    for i, segment in enumerate(segments):
        (chr, start, end) = segment
        out_file.write(f'Segment	{i+1}	{chr}	{start}	{end}\n')
    out_file.write('List of longest subpath constraints\n')
    for i, BFB_string in enumerate(BFB_strings):
        if check_BFB_string(BFB_string) == False:
            print('Non-BFB string:')
            print_BFB_string(BFB_string)
            continue
        else:
            print(f'BFB string {i+1} saved to cycles.txt file:')
            print_BFB_string(BFB_string)
        path = []
        for seg in BFB_string:
            if seg > 0:
                path.append(f'{seg}+')
            else:
                path.append(f'{-seg}-')
        out_file.write(f'Path={i+1};Copy_count=1;Segments={','.join(path)};Path_constraints_satisfied=;Score={scores[i]};Multiplicity={multiplicity}\n')
    out_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BFBArchitect based on AA output.")
    parser.add_argument("--graph", help="Path to the graph.txt file from AA output.", required=True)
    parser.add_argument("--cycle", help="Path to the cycles.txt file from AA output.", required=True)
    parser.add_argument("--output_prefix", help = "Prefix of output files.", required=True)
    parser.add_argument("--multiple", help="Reconstruct multiple BFB candidates", action='store_true')
    parser.add_argument("--solver", help="ILP solver to use: 'gurobi' or 'cbc' (default: gurobi)", default='gurobi')
    parser.add_argument("-g", "--gene", help="Gene annotation", default=None)
    args = parser.parse_args()

    logger = create_logger('BFBArchitect', f'{args.output_prefix}.log')
    start_time = time.time()
    logger.info(f'Command: python run_AA_output.py {" ".join(sys.argv[1:])}')

    # Extract amplified region from the cycle file
    interval = open(args.cycle, 'r').readline().strip() # Expected format: "Interval	1	chr19	27645416	30219166"
    _, _, chrom, start, end = interval.split()
    region = (chrom, int(start), int(end))
    # Extract segments and foldback SVs from the graph file
    SVs = []
    SV_info = {}
    segments = []
    with open(args.graph, 'r') as f:
        for line in f:
            if line.startswith('discordant'):
                sv_str = line.split()[1] # The SV string, e.g., "chr3:24524687-->chr3:24523614-"
                pattern = r"(\w+):(\d+)([+-])->(\w+):(\d+)([+-])"
                match = re.match(pattern, sv_str)
                if not match:
                    raise ValueError(f"Invalid SV string format: {sv_str}")
                chrom1, bp1, strand1, chrom2, bp2, strand2 = match.groups()
                sv = SV(chrom1, int(bp1), strand1, chrom2, int(bp2), strand2)
                SVs.append(sv)
                SV_info[sv] = (float(line.split()[2]), int(line.split()[3]))
            elif line.startswith('sequence'):
                parts = line.split()
                chrom = parts[1].split(':')[0]
                start = int(parts[1].split(':')[1].split('-')[0])
                end = int(parts[2].split(':')[1].split('+')[0])
                seg_cn = float(parts[3])
                segments.append((chrom, start, end, seg_cn, float(parts[4]), int(parts[6])))
    # Segmentation
    breakpoints = set()
    for sv in SVs:
        flag1, flag2 = sv.is_in_region(region)
        if sv.is_foldback() and flag1 and flag2:
            if sv.strand1 == '+':
                breakpoints.add(sv.bp2)
            else:
                breakpoints.add(sv.bp1 - 1)
    segments.sort(key=lambda x: (x[1]))
    new_segments = []
    chrom = region[0]
    start = region[1]
    total_length, weighted_cn_sum = 0, 0
    total_bp, total_read_count = 0, 0
    for (_, seg_start, seg_end, seg_cn, seg_cov, seg_read_count) in segments:
        if seg_end in breakpoints:
            new_cn = weighted_cn_sum / total_length if total_length > 0 else 0
            new_cov, new_read_count = (total_bp/total_length, total_read_count) if total_length > 0 else (0, 0)
            new_segments.append((chrom, start, seg_end, new_cn, new_cov, new_read_count))
            start = seg_end + 1
            total_length, weighted_cn_sum = 0, 0
            total_bp, total_read_count = 0, 0
        else:
            seg_size = seg_end - seg_start + 1
            total_length += seg_size
            weighted_cn_sum += seg_cn * seg_size
            total_bp += seg_cov * seg_size
            total_read_count += seg_read_count
    # Add the last segment
    if total_length > 0:
        new_cn = weighted_cn_sum / total_length
        new_cov, new_read_count = (total_bp/total_length, total_read_count)
        new_segments.append((chrom, start, region[2], new_cn, new_cov, new_read_count))
    # Determine the chromosome arm
    if CHR_CENTRO[region[0]] is None or region[2] < CHR_CENTRO[region[0]]:
        new_segments.pop(0)
    else:
        new_segments.pop(-1)
    # Build the copy number profile and foldback SV list for BFB reconstruction
    cn = [round(seg[3]) - 1 for seg in new_segments]
    l_bp, r_bp = [seg[1] for seg in new_segments], [seg[2] for seg in new_segments]
    lf, rf = [0 for _ in range(len(cn))], [0 for _ in range(len(cn))]
    for sv in SVs:
        flag1, flag2 = sv.is_in_region(region)
        if sv.is_foldback() and flag1 and flag2:
            if sv.strand1 == '-' and flag1:
                i = l_bp.index(sv.bp1)
                lf[i] += round(SV_info[sv][0])
            elif sv.strand1 == '+' and flag2:
                i = r_bp.index(sv.bp2)
                rf[i] += round(SV_info[sv][0])
    cn0, lf0, rf0 = cn[:], lf[:], rf[:]
    max_pos = max(l_bp + r_bp)
    start_segment = -len(segments) if max_pos < CHR_CENTRO[chrom] else 1
    multiplicity = 1
    cn_bound = 15 if args.solver == 'gurobi' else 12
    while max(cn)/multiplicity > cn_bound or (sum(lf0) + sum(rf0) + 1)/multiplicity > cn_bound:
        multiplicity += 1
    logger.info(f'Start segment: {start_segment}')
    logger.info(f'cn0: {cn0}')
    logger.info(f'lf0: {lf0}')
    logger.info(f'rf0: {rf0}')
    logger.info(f'Multiplicity: {multiplicity}')
    max_cn, total_foldbacks = max(cn), sum(lf) + sum(rf) + 1
    scale_factor = max(1, (multiplicity * min(1, total_foldbacks/max_cn)))
    cn = [c / multiplicity for c in cn]
    lf = [c / multiplicity for c in lf]
    rf = [c / multiplicity for c in rf]
    # Reconstruct BFB strings
    print("Reconstructing BFB sequences using ILP...")
    if args.multiple:
        BFB_strings, obj_val = reconstruct_BFB_strings(cn, lf, rf, start_segment)
    else:
        if args.solver == 'gurobi':
            BFB_strings, obj_val = reconstruct_BFB_strings(cn, lf, rf, start_segment, pool_solutions=1)
        else: # COIN-OR CBC solver does not support solution pool, so we only reconstruct one BFB string
            BFB_string, obj_val = reconstruct_BFB_string(cn, lf, rf, start_segment)
            BFB_strings = [BFB_string]
    print("BFB reconstruction completed.")
    logger.info(f'ILP objective value: {obj_val}')
    scores = compute_BFB_scores(new_segments, cn0, lf0, rf0, BFB_strings, multiplicity, logger)
    generate_graph_file(f'{args.output_prefix}_BFB_graph.txt', new_segments, SVs, SV_info)
    generate_cycle_file(f'{args.output_prefix}_BFB_cycles.txt', new_segments, BFB_strings, scores, multiplicity)
    logger.info(f'Total time: {time.time() - start_time} seconds')
    # Visualize the reconstructed BFB sequences
    visualize_BFB(cycle_file=f'{args.output_prefix}_BFB_cycles.txt',
            graph_file=f'{args.output_prefix}_BFB_graph.txt',
            cnr_file=None,
            output_prefix=f'{args.output_prefix}_BFB',
            gene_annotation=args.gene,
            multiple=args.multiple
        )
