import argparse
import time
import pandas as pd
import pysam
from collections import defaultdict
from pathlib import Path

try:
    from src.SVCaller import call_SVs
    from src.BFBSolver import reconstruct_BFB_string, check_BFB_string, print_BFB_string
    from src.datatypes import CHR_CENTRO
    from src.utils import create_logger
except:
    from SVCaller import call_SVs
    from BFBSolver import reconstruct_BFB_string, check_BFB_string, print_BFB_string
    from datatypes import CHR_CENTRO
    from utils import create_logger

def get_normal_coverage(cns_fn, bam_fn):
    # Get normal genome regions
    cns = pd.read_csv(cns_fn, sep="\t")
    segments = cns.sort_values(by='log2').reset_index(drop=True) # sort all segments by log2
    l = int(len(segments) / 2.4)
    r = l + 1
    total_length = 0
    log2_cn = []
    regions = []
    while total_length < 10_000_000:
        l_chr, l_start, l_end = segments.loc[l].chromosome, int(segments.loc[l].start), int(segments.loc[l].end)
        regions.append((l_chr, l_start, l_end))
        total_length += l_end - l_start
        log2_cn.append(segments.loc[l].log2)
        r_chr, r_start, r_end = segments.loc[r].chromosome, int(segments.loc[r].start), int(segments.loc[r].end)
        regions.append((r_chr, r_start, r_end))
        total_length += r_end - r_start
        log2_cn.append(segments.loc[r].log2)
        l -= 1
        r += 1
    # Get normal coverage
    (normal_cov, _) = get_coverage_and_rc(bam_fn, regions)
    return normal_cov

def get_coverage_and_rc(bam_fn, intervals, qc_threshold=0):
    total_length, total_bases = 0, 0
    bam = pysam.AlignmentFile(bam_fn, "rb")
    read_count = 0
    for (chrom, start, end) in intervals:
        total_length += end - start + 1
        for read in bam.fetch(chrom, start, end):
            if read.mapping_quality < qc_threshold or read.seq == None:
                continue
            read_count += 1
            for block_start, block_end in read.get_blocks():
                if block_end < start or block_start > end:
                    continue
                total_bases += min(block_end, end) - max(block_start, start)
    coverage = total_bases / total_length
    return (coverage, read_count)

def segmentation(cns_fn, bam_fn, regions, SVs, normal_cov, bkp_distance=50000, tolerance=0.1, CNV_segmentation=False):
    # Get SV breakpoints in the BFB regions and extra regions
    SV_breakpoints = defaultdict(list)
    for sv, count in SVs.items():
        flag1, flag2 = sv.is_in_regions(regions)
        if sv.type == 'FBI':
            if sv.strand1 == '-':
                if flag1:
                    SV_breakpoints[sv.chrom1].append(-sv.bp1)
                else:
                    SV_breakpoints[sv.chrom1+'_extra'].append(-sv.bp1)
            else:
                if flag2:
                    SV_breakpoints[sv.chrom2].append(sv.bp2)
                else:
                    SV_breakpoints[sv.chrom2+'_extra'].append(sv.bp2)
        elif flag1 == False:
            SV_breakpoints[sv.chrom1+'_extra'].append(sv.bp1 if sv.strand1 == '+' else -sv.bp1)
        elif flag2 == False:
            SV_breakpoints[sv.chrom2+'_extra'].append(sv.bp2 if sv.strand2 == '+' else -sv.bp2)
    for chrom, bkps in SV_breakpoints.items():
        SV_breakpoints[chrom] = list(set(bkps))
    # Build extra regions
    extra_regions = []
    for key in SV_breakpoints.keys():
        if key.endswith('_extra'):
            chrom = key[:-6]
            bkps = sorted([abs(bkp) for bkp in SV_breakpoints[key]])
            left, right = 0, 0
            while right < len(bkps) - 1:
                if bkps[right+1] - bkps[right] > 1000000:
                    extra_regions.append((chrom, max(0, bkps[left]-100000), bkps[right]+100000))
                    right += 1
                    left = right
                else:
                    right += 1
            if left < len(bkps):
                extra_regions.append((chrom, max(0, bkps[left]-100000), bkps[-1]+100000))
    # Segmentation based on CNV and SV
    segments, extra_segments = [], []
    cns = pd.read_csv(cns_fn, sep="\t")
    for (chrom, start, end) in regions + extra_regions:
        # Find segment boundaries
        is_BFB_region = (chrom, start, end) in regions
        bkps = SV_breakpoints.get(chrom if is_BFB_region else chrom+'_extra', [])
        cns_region = cns[(cns.chromosome == chrom) & (start <= cns.start) & (cns.end <= end)]
        CNV_boundaries = []
        if len(cns_region) > 0:
            CNV_boundaries = cns_region['start'].values.tolist()
            CNV_boundaries.append(int(cns_region.iloc[-1]['end']))
            CNV_boundaries = list(list(filter(lambda b: (b in bkps or -b in bkps) == False, CNV_boundaries)))
        if CNV_segmentation:
            boundaries = list(set(CNV_boundaries + bkps + [start, end]))
        else:
            boundaries = list(set(bkps + [start, end]))
        boundaries.sort(key=abs)
        if is_BFB_region and bkps:
            max_pos = max([abs(bkp) for bkp in bkps])
            p_arm = max_pos < CHR_CENTRO[chrom]
            i = 0
            while p_arm and boundaries[i] not in bkps:
                i += 1
            j = len(boundaries) - 1
            while not p_arm and boundaries[j] not in bkps:
                j -= 1
            boundaries = boundaries[i:j+1]
            left, right = boundaries[0], boundaries[-1]
            if p_arm and left > 0:
                boundaries.insert(0, min(start, abs(left)-100000))
            elif not p_arm and right < 0:
                boundaries.append(max(end, abs(right)+100000))
        else:
            l, r = boundaries.index(start), boundaries.index(end)
            boundaries = boundaries[l:r+1]

        # Merge boundaries based on genomic distance 
        i = 0
        while i < len(boundaries) - 1:
            curr_bd, next_bd = boundaries[i], boundaries[i+1]
            if abs(next_bd) - abs(curr_bd) < bkp_distance:
                if curr_bd not in bkps and next_bd not in bkps:
                    boundaries.remove(next_bd)
                elif curr_bd not in bkps:
                    boundaries.remove(curr_bd)
                elif next_bd not in bkps:
                    boundaries.remove(next_bd)
                else:
                    i += 1
            else:
                i += 1
        # Merge segments based on copy numbers
        segment_list = []
        for bp1, bp2 in zip(boundaries[:-1], boundaries[1:]):
            pos1 = bp1 + 1 if bp1 > 0 else -bp1
            pos2 = bp2 if bp2 > 0 else -bp2 - 1
            if pos1 > pos2:
                pos1, pos2 = pos2, pos1
            segment_list.append((chrom, pos1, pos2))
        coverage_and_rc = [get_coverage_and_rc(bam_fn, [segment]) for segment in segment_list]
        cn = [round(2*c/normal_cov)-1 for (c, _) in coverage_and_rc]
        i = 0
        bkps_abs = [abs(bkp) for bkp in bkps]
        while i < len(segment_list) - 1:
            curr_cn, next_cn = cn[i], cn[i+1]
            # max_difference = curr_cn * tolerance if (segment_list[i][2]-segment_list[i][1]) > \
            #                 (segment_list[i+1][2]-segment_list[i+1][1]) else next_cn * tolerance
            max_difference = next_cn // 10 + 1 if next_cn >= 5 else 0
            if abs(curr_cn - next_cn) <= max_difference or curr_cn == 0 or next_cn == 0:
                curr_seg, next_seg = segment_list[i], segment_list[i+1]
                if curr_seg[2] not in bkps_abs and next_seg[1] not in bkps_abs:
                    new_seg = (curr_seg[0], curr_seg[1], next_seg[2])
                    segment_list[i] = new_seg
                    del segment_list[i+1]
                    (coverage, _) = get_coverage_and_rc(bam_fn, [new_seg])
                    cn[i] = round(2*coverage/normal_cov)-1
                    del cn[i+1]
                    continue
            i += 1
        if is_BFB_region:
            segments += segment_list
        else:
            extra_segments += segment_list
    return segments, extra_segments

def compute_BFB_scores(coverage, normal_cov, cn0, lf0, rf0, BFB_palindromes, multiplicity, logger, is_BFB_string=False):
    scores = []
    for BFB_palindrome in BFB_palindromes:
        if is_BFB_string == False:
            BFB_string = BFB_palindrome[:len(BFB_palindrome)//2]
        else:
            BFB_string = BFB_palindrome
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
        alpha = 1
        fb_dist = alpha * sum([(lf0[i]-lf[i])**2 + (rf0[i]-rf[i])**2 for i in range(len(cn0))])**0.5 / len(cn0)
        logger.info(f'Foldback Euclidean distance: {fb_dist}')
        fb_score = 0
        beta = 0.5
        for i in range(len(cn0)):
            if (lf0[i] == 0 and lf[i] != 0):
                fb_score += beta * lf[i]
            if (rf0[i] == 0 and rf[i] != 0):
                fb_score += beta * rf[i]
        # fb_score /= len(cn0)
        # fb_count = sum([1 for i in range(len(cn0)) if (lf0[i] > 0 or rf0[i] > 0)])
        fb_count = sum([1 for i in range(len(lf0)) if lf0[i] > 0]) + sum([1 for i in range(len(rf0)) if rf0[i] > 0])
        if fb_count < 2 or len(cn0) < 2:
            fb_score += 2
        logger.info(f'Foldback score: {fb_score}')
        gamma = 1 / len(cn0)
        nanopore_score = gamma * sum([abs(2*coverage[i]/normal_cov-1-cn0[i])/cn0[i] for i in range(len(cn0)) if cn0[i] > 0])
        logger.info(f'Nanopore score: {nanopore_score}')
        score = CN_score + fb_dist + fb_score + nanopore_score
        logger.info(f'Final score: {score}')
        scores.append(score)
    return scores

def generate_graph_file(output_fn, bam_fn, segments, extra_segments, coverage_and_rc, SVs, normal_cov):
    coverage_and_rc += [get_coverage_and_rc(bam_fn, [segment]) for segment in extra_segments]
    all_segments = segments + extra_segments
    # Generate graph.txt
    out_file = open(output_fn, 'w')
    out_file.write('SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads\n')
    segment_cn = [max(0, round(2*c/normal_cov)-1) for (c, _) in coverage_and_rc]
    for i, seg in enumerate(all_segments):
        size = seg[2] - seg[1] + 1
        coverage, read_count = coverage_and_rc[i]
        entry = f'sequence	{seg[0]}:{seg[1]}-	{seg[0]}:{seg[2]}+	{segment_cn[i]}	{coverage}	{size}	{read_count}\n'
        out_file.write(entry)
    out_file.write('BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads\n')
    for i in range(1, len(all_segments)):
        seg1, seg2 = all_segments[i-1], all_segments[i]
        if seg1[0] != seg2[0] or seg1[2]+1 != seg2[1]: 
            continue
        cn = min(segment_cn[i-1], segment_cn[i])
        read_count = int((coverage_and_rc[i-1][0]+coverage_and_rc[i][0])/2)
        entry = f'concordant	{seg1[0]}:{seg1[2]}+->{seg2[0]}:{seg2[1]}-	{cn}	{read_count}\n'
        out_file.write(entry)
    for sv, count in SVs.items():
        sv_str = str(sv)
        entry = f'discordant	{sv_str}	{round(2*count/normal_cov)}	{count}\n'
        out_file.write(entry)
    out_file.close()

def generate_cycle_file(output_fn, segments, BFB_palindromes, scores, multiplicity, is_BFB_string=False):
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
    for i, BFB_palindrome in enumerate(BFB_palindromes):
        if check_BFB_string(BFB_palindrome) == False:
            print('Non-BFB sequence:')
            print_BFB_string(BFB_palindrome)
            continue
        else:
            print('BFB sequence saved to cycles.txt file:')
            print_BFB_string(BFB_palindrome)
        path = []
        BFB_string = BFB_palindrome[:len(BFB_palindrome)//2] if is_BFB_string == False else BFB_palindrome
        for seg in BFB_string:
            if seg > 0:
                path.append(f'{seg}+')
            else:
                path.append(f'{-seg}-')
        out_file.write(f'Path={i+1};Copy_count=1;Segments={','.join(path)};Path_constraints_satisfied=;Score={scores[i]};Multiplicity={multiplicity}\n')
    out_file.close()

def main():
    parser = argparse.ArgumentParser(description = "Segmentation based on CNV and SV.")
    parser.add_argument("--bam", help = "Path to a sorted bam file", required = True)
    parser.add_argument("--cns", help = "Path to a sorted cns file generated by CNVKit", required = True)
    parser.add_argument("--regions", help = "Segmentation on specific regions (e.g. chr1:1000000-2000000;chr2:1000000-2000000)", required = True)
    parser.add_argument("--output_prefix", help = "Prefix of output files.", required=True)
    parser.add_argument("--segmentation", help="Consider CNV in segmentation", action='store_true')
    parser.add_argument("--deletion", help="Deletion handling", action='store_true')
    parser.add_argument("--coverage", help="Sequencing coverage (if provided, estimation from cns will be skipped)", type=float, default=None)
    args = parser.parse_args()

    logger = create_logger('BFBArchitect', f'{args.output_prefix}.log')
    start_time = time.time()
    logger.info(f'Command: python {Path(__file__).resolve()} --bam {args.bam} --cns {args.cns} --regions {args.regions} --output_prefix {args.output_prefix}' +
                 (' --segmentation' if args.segmentation else '') + (' --deletion' if args.deletion else '') + (' --coverage ' + str(args.coverage) if args.coverage != None else ''))
    normal_cov = get_normal_coverage(args.cns, args.bam) if args.coverage == None else args.coverage
    logger.info(f'Normal coverage: {normal_cov}')
    regions = [] # regions = [('chr11', 70000000, 75000000), ('chr17', 40000001, 43000001)]
    for info in args.regions.split(';'):
        chrom = info.split(':')[0]
        start = int(info.split(':')[1].split('-')[0])
        end = int(info.split('-')[1])
        regions.append((chrom, start, end))
    # Call SVs from amplicon regions
    print("Calling SVs in the amplicon regions...")
    SVs = {}
    output_read_fn = None if args.output_prefix == None else f'{args.output_prefix}_reads.txt'
    if output_read_fn != None:
        output_file = open(output_read_fn, 'w')
        output_file.close()
    for region in regions:
        SV_dict = call_SVs(args.bam, region, normal_cov=normal_cov, output_fn=output_read_fn)
        SVs.update(SV_dict)
    print(f'Saved structural variants to {args.output_prefix}_reads.txt.')
    foldback_flag = False
    for sv in SVs.keys():
        if sv.type == 'FBI':
            foldback_flag = True
            break
    if foldback_flag == False:
        print('No foldback inversion found in this region. ')
        exit(0)
    # Segmentation 
    print("Segmenting the amplicon region...")
    segments, extra_segments = segmentation(args.cns, args.bam, regions, SVs, normal_cov, tolerance=0.1, CNV_segmentation=args.segmentation)
    # CNV calling 
    coverage_and_rc = [get_coverage_and_rc(args.bam, [segment]) for segment in segments]
    cn = [max(0, round(2*c/normal_cov)-1) for (c, _) in coverage_and_rc]
    # Restimate segment CN based on deletions
    if args.deletion:
        print("Handling deletions...")
        for i, segment in enumerate(segments):
            missing_bases = 0
            deletion_length = 0
            (chrom, start, end) = segment

            for sv, count in SVs.items():
                if sv.TST == False and sv.type == 'DEL' and sv.chrom1 == chrom and abs(sv.bp1 - sv.bp2) <= 10000000:
                    del_start, del_end = sv.bp1, sv.bp2
                    if del_end < start or del_start > end:
                        continue
                    overlap_start = max(start, del_start)
                    overlap_end = min(end, del_end)
                    missing_bases += (overlap_end - overlap_start + 1) * count
                    deletion_length += overlap_end - overlap_start + 1
            if missing_bases > 0:
                segment_length = end - start + 1 # - deletion_length
                coverage_and_rc[i] = ((segment_length * coverage_and_rc[i][0] + missing_bases) / (segment_length), coverage_and_rc[i][1])
                cn[i] = round(coverage_and_rc[i][0] * 2 / normal_cov) - 1
    # Get vectors for CN, left foldbacks, and right foldbacks
    l_bp, r_bp = [bp1 for (_, bp1, _) in segments], [bp2 for (_, _, bp2) in segments]
    lf, rf = [0 for _ in range(len(cn))], [0 for _ in range(len(cn))]
    for sv, count in SVs.items():
        if sv.type == 'FBI':
            flag1, flag2 = sv.is_in_regions(regions)
            if sv.strand1 == '-' and flag1:
                i = l_bp.index(sv.bp1)
                lf[i] += round(2*count/normal_cov)
            elif sv.strand1 == '+' and flag2:
                i = r_bp.index(sv.bp2)
                rf[i] += round(2*count/normal_cov)
    cn0, lf0, rf0 = cn[:], lf[:], rf[:]

    max_pos = max(l_bp + r_bp)
    start_segment = -len(segments) if max_pos < CHR_CENTRO[chrom] else 1
    multiplicity = 1
    while max(cn) / multiplicity > 12:
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
    
    print("Reconstructing BFB sequences using ILP...")
    BFB_string, obj_val = reconstruct_BFB_string(cn, lf, rf, start_segment)
    print("BFB reconstruction completed.")

    logger.info(f'ILP objective value: {obj_val}')
    coverage = [c for (c, _) in coverage_and_rc]
    scores = compute_BFB_scores(coverage, normal_cov, cn0, lf0, rf0, [BFB_string], multiplicity, logger, is_BFB_string=True)
    if args.output_prefix != None:
        generate_graph_file(f'{args.output_prefix}_graph.txt', args.bam, segments, extra_segments, coverage_and_rc, SVs, normal_cov)
        print(f'Generated {args.output_prefix}_graph.txt file.')
        generate_cycle_file(f'{args.output_prefix}_cycles.txt', segments+extra_segments, [BFB_string], scores, multiplicity, is_BFB_string=True)
        print(f'Generated {args.output_prefix}_cycles.txt file.')
    logger.info(f'Total time: {time.time() - start_time} seconds')

if __name__ == "__main__":
    main()