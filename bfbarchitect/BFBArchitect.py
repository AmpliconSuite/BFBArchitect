import argparse
import logging
import os
import time
import pandas as pd

from collections import defaultdict
from pathlib import Path

try:
    from bfbarchitect.SVCaller import call_SVs
    from bfbarchitect.BFBSolver import reconstruct_BFB_cbc, reconstruct_BFB_gurobi, reconstruct_BFB_mosek, check_BFB_string, print_BFB_string
    from bfbarchitect.datatypes import CHR_CENTRO, build_centromere_dict
    from bfbarchitect.utils import create_logger, get_normal_coverage, get_coverage_and_rc, get_chrom_length
except:
    from SVCaller import call_SVs
    from BFBSolver import reconstruct_BFB_cbc, reconstruct_BFB_gurobi, reconstruct_BFB_mosek, check_BFB_string, print_BFB_string
    from datatypes import CHR_CENTRO, build_centromere_dict
    from utils import create_logger, get_coverage_and_rc, get_normal_coverage, get_chrom_length


def _lookup_chrom(d, chrom, label):
    """Dict lookup with an informative error when the chromosome is missing.

    GRCh37/hg19 data uses bare chromosome names (e.g. '7') while the
    built-in reference tables use 'chr'-prefixed names (e.g. 'chr7').
    If a lookup fails, direct the user to supply a matching BED file.
    """
    try:
        return d[chrom]
    except KeyError:
        raise KeyError(
            f"Chromosome '{chrom}' not found in {label}. "
            "If using GRCh37/hg19 chromosome naming (without the 'chr' prefix), "
            "supply a matching centromere BED file via --centromere."
        ) from None

def expand_amplicon_region(bam_fn, normal_cov, region, CN_threshold=2, size=100000, centromere_dict=None):
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    chrom, start, end = region
    centro = _lookup_chrom(centromere_dict, chrom, 'centromere data')
    if end <= centro:
        left_bound, right_bound = 0, centro
    else:
        chr_size = get_chrom_length(bam_fn, chrom)
        if chr_size is None:
            raise KeyError(
                f"Chromosome '{chrom}' not found in BAM header. "
                "Check that --region uses the same chromosome naming as the BAM file."
            )
        left_bound, right_bound = centro, chr_size
    left = start
    while left-size >= left_bound and start - left < 100*size:
        coverage, _ = get_coverage_and_rc(bam_fn, (chrom, left-size, left))
        if round(coverage * 2 / normal_cov - 1) <= CN_threshold:
            break
        left -= size
    right = end
    while right+size <= right_bound and right - end < 100*size:
        coverage, _ = get_coverage_and_rc(bam_fn, (chrom, right, right+size))
        if round(coverage * 2 / normal_cov - 1) <= CN_threshold:
            break
        right += size
    return (chrom, left, right)

def segment_region(cns_fn, bam_fn, region, SVs, normal_cov, bkp_distance=50000, tolerance=0.1, CNV_segmentation=False, centromere_dict=None):
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    # Get SV breakpoints in the BFB region and extra regions
    SV_breakpoints = defaultdict(list)
    for sv, count in SVs.items():
        flag1, flag2 = sv.is_in_region(region)
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
    for (chrom, start, end) in [region] + extra_regions:
        # Find segment boundaries
        is_BFB_region = (chrom, start, end) == region
        bkps = SV_breakpoints.get(chrom if is_BFB_region else chrom+'_extra', [])
        cns_region = cns[(cns.chromosome == chrom) & (start <= cns.start) & (cns.end <= end)]
        CNV_boundaries = []
        if len(cns_region) > 0:
            CNV_boundaries = cns_region['start'].values.tolist()
            CNV_boundaries.append(int(cns_region.iloc[-1]['end']))
            CNV_boundaries = list(filter(lambda b: (b in bkps or -b in bkps) == False, CNV_boundaries))
        if CNV_segmentation:
            boundaries = list(set(CNV_boundaries + bkps + [start, end]))
        else:
            boundaries = list(set(bkps + [start, end]))
        boundaries.sort(key=abs)
        if is_BFB_region and bkps:
            max_pos = max([abs(bkp) for bkp in bkps])
            p_arm = max_pos < _lookup_chrom(centromere_dict, chrom, 'centromere data')
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
        coverage_and_rc = [get_coverage_and_rc(bam_fn, segment) for segment in segment_list]
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
                    (coverage, _) = get_coverage_and_rc(bam_fn, new_seg)
                    cn[i] = round(2*coverage/normal_cov)-1
                    del cn[i+1]
                    continue
            i += 1
        if is_BFB_region:
            segments += segment_list
        else:
            extra_segments += segment_list
    return segments, extra_segments

def compute_bfb_scores(cn0, lf0, rf0, BFB_strings, multiplicity, logger,
                       observed_cn=None, normal_cov=None, silent=False):
    """
    Score BFB candidate strings against segment CN and foldback counts.

    Parameters
    ----------
    cn0 : list of int
        Expected integer copy number per segment.
    lf0, rf0 : list of int
        Expected left/right foldback counts per segment.
    BFB_strings : list of list[int]
        Candidate BFB strings to score.
    multiplicity : int
        The multiplicity factor applied to the ILP solution.
    logger : logging.Logger
        Logger for reporting scores.
    observed_cn : list of float, optional
        Observed float copy number per segment.
        If None, the CN divergence term is omitted.
    normal_cov : float, optional
        Normal (diploid) coverage. If provided, weights are applied to
        CN discrepancy and foldback distance when normal_cov < 7.
    silent : bool, optional
        If True, suppress print statements.

    Returns
    -------
    list of float
        Scores for each BFB string (lower is better).
    """
    scores = []
    for idx, BFB_string in enumerate(BFB_strings):
        logger.info('----------------------------------------')
        if not check_BFB_string(BFB_string):
            logger.warning(f'BFB string {idx+1} is not a valid BFB sequence; skipping.')
            continue

        logger.info(f'Scores for BFB string {idx+1}')
        # Calculate CN and foldback vectors for the BFB string
        cn = [0] * len(cn0)
        lf = [0] * len(lf0)
        rf = [0] * len(rf0)
        for i in range(len(BFB_string) - 1):
            seg1, seg2 = BFB_string[i], BFB_string[i+1]
            cn[abs(seg1)-1] += 1
            if abs(seg1) == abs(seg2):
                if seg1 > 0 and seg2 < 0:
                    rf[abs(seg1)-1] += 1
                elif seg1 < 0 and seg2 > 0:
                    lf[abs(seg1)-1] += 1
        cn[abs(BFB_string[-1])-1] += 1

        # Scale by multiplicity
        cn = [c * multiplicity for c in cn]
        lf = [c * multiplicity for c in lf]
        rf = [c * multiplicity for c in rf]

        # 1. CN discrepancy score (predicted integer CN vs expected integer CN)
        CN_score = sum(abs(cn0[i] - cn[i]) / cn[i] if cn[i] > 0 else abs(cn0[i] - cn[i])
                       for i in range(len(cn)))
        if normal_cov is not None and normal_cov < 7:
            CN_score *= 0.5
            logger.info(f'CN discrepancy (weight = 0.5 due to low coverage): {CN_score}')
        else:
            logger.info(f'CN discrepancy: {CN_score}')

        # 2. Foldback Euclidean distance
        fb_dist = sum((lf0[i] - lf[i])**2 + (rf0[i] - rf[i])**2 for i in range(len(cn0)))**0.5 / len(cn0)
        if normal_cov is not None and normal_cov < 7:
            fb_dist *= 0.3
            logger.info(f'Foldback Euclidean distance (weight = 0.3 due to low coverage): {fb_dist}')
        else:
            logger.info(f'Foldback Euclidean distance: {fb_dist}')

        # 3. Missing foldback penalty
        missing_fb_score = 0
        for i in range(len(cn0)):
            if lf0[i] == 0 and lf[i] != 0:
                missing_fb_score += 0.5 * lf[i]
            if rf0[i] == 0 and rf[i] != 0:
                missing_fb_score += 0.5 * rf[i]
        fb_count = sum(1 for x in lf0 if x > 0) + sum(1 for x in rf0 if x > 0)
        if fb_count < 2 or len(cn0) < 2:
            missing_fb_score += 2
        logger.info(f'Missing foldback score: {missing_fb_score}')

        # 4. CN divergence score (observed float CN vs expected integer CN)
        cn_divergence = 0
        if observed_cn is not None:
            cn_divergence = sum(abs(observed_cn[i] - cn0[i]) / cn0[i]
                               for i in range(len(cn0)) if cn0[i] > 0)
            logger.info(f'CN divergence score: {cn_divergence}')

        total_score = CN_score + fb_dist + missing_fb_score + cn_divergence
        label = print_BFB_string(BFB_string, print_to_console=False)
        logger.info(f'BFB string {idx+1}: {label}')
        logger.info(f'Total score: {total_score}')
        cn_div_str = f", cn_div={cn_divergence:.4f}" if observed_cn is not None else ""
        if not silent:
            print(f"  BFB {idx+1}: score={total_score:.4f}  (CN={CN_score:.4f}, fb={fb_dist:.4f}, miss_fb={missing_fb_score:.4f}{cn_div_str})  {label}")
        scores.append(total_score)

    return scores

def _add_stream_handler(logger):
    """Attach a stdout StreamHandler to logger if none is present."""
    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
               for h in logger.handlers):
        sh = logging.StreamHandler()
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s'))
        logger.addHandler(sh)

def detect_solver() -> str:
    """Return 'gurobi' if a Gurobi license and package are available,
    otherwise 'mosek' if available, otherwise 'cbc'."""
    # 1. Gurobi
    try:
        import gurobipy  # noqa: F401
        if os.path.exists(os.path.expanduser('~/gurobi.lic')):
            return 'gurobi'
    except ImportError:
        pass

    # 2. MOSEK
    try:
        import mosek  # noqa: F401
        if os.path.exists(os.path.expanduser('~/mosek/mosek.lic')):
            return 'mosek'
    except ImportError:
        pass

    return 'cbc'

def write_bfb_graph(output_fn, new_segments, SVs, sv_info):
    """Write a BFB graph file from pre-segmented data."""
    with open(output_fn, 'w') as f:
        f.write('SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberReadsMapped\n')
        for seg in new_segments:
            size = seg[2] - seg[1] + 1
            f.write(f'sequence\t{seg[0]}:{seg[1]}-\t{seg[0]}:{seg[2]}+\t{seg[3] - 1}\t{seg[4]}\t{size}\t{seg[5]}\n')
        f.write('BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs\n')
        for i in range(1, len(new_segments)):
            seg1, seg2 = new_segments[i-1], new_segments[i]
            if seg1[0] != seg2[0] or seg1[2]+1 != seg2[1]:
                continue
            cn = min(seg1[3] - 1, seg2[3] - 1)
            read_count = int((seg1[5] + seg2[5]) / 2)
            f.write(f'concordant\t{seg1[0]}:{seg1[2]}+->{seg2[0]}:{seg2[1]}-\t{cn}\t{read_count}\n')
        for sv in SVs:
            f.write(f'discordant\t{sv}\t{sv_info[sv][0]}\t{sv_info[sv][1]}\n')

def write_bfb_cycles(output_fn, new_segments, BFB_strings, scores, multiplicity):
    """Write a BFB cycles file from pre-segmented data."""
    segments = [(seg[0], seg[1], seg[2]) for seg in new_segments]
    intervals = []
    chrom, start, end = segments[0]
    for next_chrom, next_start, next_end in segments[1:]:
        if next_chrom == chrom and next_start == end + 1:
            end = next_end
        else:
            intervals.append((chrom, start, end))
            chrom, start, end = next_chrom, next_start, next_end
    intervals.append((chrom, start, end))
    with open(output_fn, 'w') as f:
        for i, (chrom, start, end) in enumerate(intervals):
            f.write(f'Interval\t{i+1}\t{chrom}\t{start}\t{end}\n')
        f.write('List of cycle segments\n')
        for i, (chrom, start, end) in enumerate(segments):
            f.write(f'Segment\t{i+1}\t{chrom}\t{start}\t{end}\n')
        f.write('List of longest subpath constraints\n')
        for i, BFB_string in enumerate(BFB_strings):
            if not check_BFB_string(BFB_string):
                continue
            path = [f'{seg}+' if seg > 0 else f'{-seg}-' for seg in BFB_string]
            f.write(f"Path={i+1};Copy_count=1;Segments={','.join(path)};Path_constraints_satisfied=;Score={scores[i]};Multiplicity={multiplicity}\n")

def trim_background_segments(new_segments, cn, lf, rf):
    """Strip leading and trailing segments with cn <= 0 (non-amplified background)."""
    while new_segments and cn[0] <= 0:
        new_segments, cn, lf, rf = new_segments[1:], cn[1:], lf[1:], rf[1:]
    while new_segments and cn[-1] <= 0:
        new_segments, cn, lf, rf = new_segments[:-1], cn[:-1], lf[:-1], rf[:-1]
    return new_segments, cn, lf, rf

def reconstruct_bfb(new_segments, cn, lf, rf, centromere_pos, solver=None, multiple=False,
                    verbose=False, log_file=None, silent=False, threads=8, track_solve=False,
                    min_lp_bound=None):
    """
    Reconstruct BFB sequences from pre-segmented copy-number and foldback data.
    """
    logger = logging.getLogger('BFBArchitect')
    old_level = logger.level
    if silent:
        logger.setLevel(logging.WARNING)
    elif verbose:
        _add_stream_handler(logger)

    try:
        if solver is None:
            solver = detect_solver()
        cn0, lf0, rf0 = cn[:], lf[:], rf[:]
        max_pos = max(seg[2] for seg in new_segments)
        start_segment = -len(new_segments) if max_pos < centromere_pos else 1
        multiplicity = 1
        cn_bound = 15 if solver == 'gurobi' else 12
        while max(cn) / multiplicity > cn_bound or (sum(lf0) + sum(rf0) + 1) / multiplicity > cn_bound:
            multiplicity += 1
        logger.info(f'Solver: {solver}')
        logger.info(f'Start segment: {start_segment}')
        logger.info(f'cn0: {cn0}')
        logger.info(f'lf0: {lf0}')
        logger.info(f'rf0: {rf0}')
        logger.info(f'Multiplicity: {multiplicity}')
        cn_scaled = [c / multiplicity for c in cn]
        lf_scaled = [c / multiplicity for c in lf]
        rf_scaled = [c / multiplicity for c in rf]
        if not silent:
            print(f"Reconstructing BFB sequences using ILP (solver={solver}, multiplicity={multiplicity})...")
        _verbose = verbose if not silent else False
        _score_fn = None
        if (_verbose or track_solve) and solver == 'gurobi':
            _null_log = logging.getLogger('_bfb_score_null')
            _null_log.addHandler(logging.NullHandler())
            _cn0, _lf0, _rf0, _mult = cn0, lf0, rf0, multiplicity
            _score_fn = lambda s: (compute_bfb_scores(_cn0, _lf0, _rf0, [s], _mult,
                                                      _null_log, silent=True) or [float('inf')])[0]
        if multiple:
            BFB_strings, obj_val = reconstruct_BFB_gurobi(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                            max_threads=threads, log_file=log_file,
                                                            verbose=_verbose, score_fn=_score_fn,
                                                            min_lp_bound=min_lp_bound)
        elif solver == 'gurobi':
            BFB_strings, obj_val = reconstruct_BFB_gurobi(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                            pool_solutions=1, max_threads=threads, log_file=log_file,
                                                            verbose=_verbose, score_fn=_score_fn,
                                                            min_lp_bound=min_lp_bound)
        elif solver == 'mosek':
            BFB_strings, obj_val = reconstruct_BFB_mosek(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                        max_threads=threads, log_file=log_file, verbose=_verbose)
        else:
            BFB_string, obj_val = reconstruct_BFB_cbc(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                          max_threads=threads)
            BFB_strings = [BFB_string]
        logger.info(f'ILP objective value: {obj_val}')
        scores = compute_bfb_scores(cn0, lf0, rf0, BFB_strings, multiplicity, logger,
                                    observed_cn=[seg[3] - 1 for seg in new_segments], silent=silent)
        return BFB_strings, scores, multiplicity
    finally:
        if silent:
            logger.setLevel(old_level)

def reconstruct_bfb_from_bam(bam_fn, cns_fn, region, output_prefix, segmentation=False, deletion=False, coverage=None, multiple=False, no_expansion=False, min_sv_cn=0.75, min_mapq=20, solver=None, centromere_dict=None, verbose=False, threads=8, min_lp_bound=None):
    if solver is None:
        solver = detect_solver()
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    log_file = f'{output_prefix}.log'
    logger = create_logger('BFBArchitect', log_file)
    if verbose:
        _add_stream_handler(logger)
    start_time = time.time()
    logger.info(f'Command: python {Path(__file__).resolve()} --bam {bam_fn} --cns {cns_fn} --region {region} --output_prefix {output_prefix}' +
                 (' --segmentation' if segmentation else '') + (' --deletion' if deletion else '') + (' --coverage ' + str(coverage) if coverage != None else '') + 
                 (' --multiple' if multiple else '') + (' --no_expansion' if no_expansion else '') + (f' --min_sv_cn {min_sv_cn}' if min_sv_cn != 0.75 else '') + 
                 (f' --min_mapq {min_mapq}' if min_mapq != 20 else ''))
    normal_cov = get_normal_coverage(cns_fn, bam_fn) if coverage == None else coverage
    logger.info(f'Normal coverage: {normal_cov}')
    if min_sv_cn * normal_cov / 2 < 3:
        min_sv_cn = 6 / normal_cov
        logger.info(f'Adjusted minimum copy number for SV calling to {min_sv_cn} to ensure at least 3 supporting reads based on normal coverage.')
    else:
        logger.info(f'Minimum copy number for SV calling: {min_sv_cn}')
    logger.info(f'Minimum number of supporting reads for SV calling: {min_sv_cn * normal_cov / 2}')
    # Parse the amplified region
    chrom = region.split(':')[0]
    start = int(region.split(':')[1].split('-')[0])
    end = int(region.split('-')[1])
    region = (chrom, start, end)
    if no_expansion == False:
        region = expand_amplicon_region(bam_fn, normal_cov, (chrom, start, end), centromere_dict=centromere_dict)
    print(f'Amplified region: {region[0]}:{region[1]}-{region[2]}')
    # Call SVs from the amplified region
    print("Calling SVs in the amplified region...")
    output_read_fn = None if output_prefix == None else f'{output_prefix}_reads.txt'
    if output_read_fn != None:
        output_file = open(output_read_fn, 'w')
        output_file.close()
    SVs = call_SVs(bam_fn, region, min_mapq=min_mapq, normal_cov=normal_cov, output_fn=output_read_fn, min_cn = min_sv_cn)
    print(f'Saved structural variants to {output_prefix}_reads.txt.')
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
    segments, extra_segments = segment_region(cns_fn, bam_fn, region, SVs, normal_cov, tolerance=0.1, CNV_segmentation=segmentation, centromere_dict=centromere_dict)
    # CNV calling 
    coverage_and_rc = [get_coverage_and_rc(bam_fn, segment) for segment in segments]
    cn = [max(0, round(2*c/normal_cov)-1) for (c, _) in coverage_and_rc]
    # Restimate segment CN based on deletions
    if deletion:
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
                segment_length = end - start + 1
                coverage_and_rc[i] = ((segment_length * coverage_and_rc[i][0] + missing_bases) / (segment_length), coverage_and_rc[i][1])
                cn[i] = round(coverage_and_rc[i][0] * 2 / normal_cov) - 1
    # Get vectors for CN, left foldbacks, and right foldbacks
    l_bp, r_bp = [bp1 for (_, bp1, _) in segments], [bp2 for (_, _, bp2) in segments]
    lf, rf = [0 for _ in range(len(cn))], [0 for _ in range(len(cn))]
    for sv, count in SVs.items():
        if sv.type == 'FBI':
            flag1, flag2 = sv.is_in_region(region)
            if sv.strand1 == '-' and flag1:
                i = l_bp.index(sv.bp1)
                lf[i] += round(2*count/normal_cov)
            elif sv.strand1 == '+' and flag2:
                i = r_bp.index(sv.bp2)
                rf[i] += round(2*count/normal_cov)
    cn0, lf0, rf0 = cn[:], lf[:], rf[:]
    max_pos = max(l_bp + r_bp)
    start_segment = -len(segments) if max_pos < _lookup_chrom(centromere_dict, chrom, 'centromere data') else 1
    multiplicity = 1
    cn_bound = 15 if solver == 'gurobi' else 12
    while max(cn)/multiplicity > cn_bound or (sum(lf0) + sum(rf0) + 1)/multiplicity > cn_bound:
        multiplicity += 1
    logger.info(f'Start segment: {start_segment}')
    logger.info(f'cn0: {cn0}')
    logger.info(f'lf0: {lf0}')
    logger.info(f'rf0: {rf0}')
    logger.info(f'Multiplicity: {multiplicity}')
    cn_scaled = [c / multiplicity for c in cn]
    lf_scaled = [c / multiplicity for c in lf]
    rf_scaled = [c / multiplicity for c in rf]
    print(f"Reconstructing BFB sequences using ILP (solver={solver}, multiplicity={multiplicity})...")
    if multiple:
        BFB_strings, obj_val = reconstruct_BFB_gurobi(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                        max_threads=threads, log_file=log_file, verbose=verbose,
                                                        min_lp_bound=min_lp_bound)
    else:
        if solver == 'gurobi':
            BFB_strings, obj_val = reconstruct_BFB_gurobi(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                            pool_solutions=1, max_threads=threads, log_file=log_file,
                                                            verbose=verbose, min_lp_bound=min_lp_bound)
        elif solver == 'mosek':
            BFB_strings, obj_val = reconstruct_BFB_mosek(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                        max_threads=threads, log_file=log_file, verbose=verbose)
        else:
            BFB_string, obj_val = reconstruct_BFB_cbc(cn_scaled, lf_scaled, rf_scaled, start_segment,
                                                          max_threads=threads)
            BFB_strings = [BFB_string]
    logger.info(f'ILP objective value: {obj_val}')
    coverage_vals = [c for (c, _) in coverage_and_rc]
    scores = compute_bfb_scores(cn0, lf0, rf0, BFB_strings, multiplicity, logger,
                               observed_cn=[2*c/normal_cov - 1 for c in coverage_vals],
                               normal_cov=normal_cov)
    if output_prefix != None:
        all_segments = segments + extra_segments
        full_segments = []
        for i, (chrom, start, end) in enumerate(all_segments):
            cov, rc = get_coverage_and_rc(bam_fn, (chrom, start, end))
            full_segments.append((chrom, start, end, 2*cov/normal_cov, cov, rc))
        sv_info = {sv: (round(2*count/normal_cov), count) for sv, count in SVs.items()}
        write_bfb_graph(f'{output_prefix}_graph.txt', full_segments, SVs, sv_info)
        print(f'Generated {output_prefix}_graph.txt file.')
        write_bfb_cycles(f'{output_prefix}_cycles.txt', full_segments, BFB_strings, scores, multiplicity)
        print(f'Generated {output_prefix}_cycles.txt file.')
    logger.info(f'Total time: {time.time() - start_time} seconds')

def reconstruct_bfb_from_graph(graph_fn, centromere_dict=None, solver=None,
                               multiple=False, whole_graph=False, region=None,
                               verbose=False, log_file=None, silent=False, threads=8,
                               track_solve=False, min_lp_bound=None):
    """
    Reconstruct BFB sequences from an AA-format _graph.txt file.
    """
    try:
        from bfbarchitect.graph_input import (find_bfb_candidate_regions,
                                               subsect_graph_for_region,
                                               whole_graph_as_region)
    except ImportError:
        from graph_input import (find_bfb_candidate_regions,
                                  subsect_graph_for_region,
                                  whole_graph_as_region)
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    results = []
    if whole_graph:
        new_segments, cn, lf, rf, svs_list, sv_info, primary_chrom = whole_graph_as_region(
            graph_fn, centromere_dict=centromere_dict)
        if new_segments:
            region = (primary_chrom, new_segments[0][1], new_segments[-1][2])
            if not silent:
                print(f'Processing whole graph as single region: {region[0]}:{region[1]}-{region[2]}')
            BFB_strings, scores, multiplicity = reconstruct_bfb(
                new_segments, cn, lf, rf,
                centromere_dict.get(primary_chrom, 0),
                solver=solver, multiple=multiple, verbose=verbose, log_file=log_file, silent=silent,
                threads=threads, track_solve=track_solve, min_lp_bound=min_lp_bound)
            results.append({
                'region': region,
                'new_segments': new_segments,
                'bfb_strings': BFB_strings,
                'scores': scores,
                'multiplicity': multiplicity,
                'svs': svs_list,
                'sv_info': sv_info,
            })
    else:
        if region is not None:
            regions = [region]
        else:
            regions = find_bfb_candidate_regions(graph_fn)
            if not regions:
                if not silent:
                    print('No BFB candidate regions found in the graph file.')
                return results
            if not silent:
                print(f'Found {len(regions)} BFB candidate region(s): '
                      + ', '.join(f'{r[0]}:{r[1]}-{r[2]}' for r in regions))
        region_data = subsect_graph_for_region(graph_fn, regions, verbose=verbose if not silent else False)
        for i, (cur_region, data) in enumerate(zip(regions, region_data)):
            if data is None:
                continue
            new_segments, cn, lf, rf, region_svs, sv_info = data
            chrom = cur_region[0]
            new_segments, cn, lf, rf = trim_background_segments(new_segments, cn, lf, rf)
            if not new_segments:
                continue
            if not silent:
                print(f'\nProcessing region {i+1}: {chrom}:{cur_region[1]}-{cur_region[2]}')
            BFB_strings, scores, multiplicity = reconstruct_bfb(
                new_segments, cn, lf, rf,
                centromere_dict.get(chrom, 0),
                solver=solver, multiple=multiple, verbose=verbose, log_file=log_file, silent=silent,
                threads=threads, track_solve=track_solve, min_lp_bound=min_lp_bound)
            results.append({
                'region': cur_region,
                'new_segments': new_segments,
                'bfb_strings': BFB_strings,
                'scores': scores,
                'multiplicity': multiplicity,
                'svs': region_svs,
                'sv_info': sv_info,
            })
    return results

def run_bfb_from_graph(graph_fn, output_prefix, multiple=False, solver=None,
                       whole_graph=False, region=None, gene=None,
                       centromere_dict=None, verbose=False, threads=8, min_lp_bound=None):
    """
    CLI entry point to reconstruct BFB sequences from an AA-format _graph.txt file.
    """
    try:
        from bfbarchitect.BFBVisualizer import visualize_BFB
    except ImportError:
        from BFBVisualizer import visualize_BFB
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    log_file = f'{output_prefix}.log'
    logger = create_logger('BFBArchitect', log_file)
    if verbose:
        _add_stream_handler(logger)
    start_time = time.time()
    results = reconstruct_bfb_from_graph(
        graph_fn, centromere_dict=centromere_dict, solver=solver,
        multiple=multiple, whole_graph=whole_graph, region=region,
        verbose=verbose, log_file=log_file, threads=threads, min_lp_bound=min_lp_bound
    )
    if not results:
        return
    for i, res in enumerate(results):
        region_prefix = output_prefix if (whole_graph or region is not None) else f'{output_prefix}_region{i + 1}'
        write_bfb_graph(f'{region_prefix}_BFB_graph.txt', res['new_segments'], res['svs'], res['sv_info'])
        write_bfb_cycles(f'{region_prefix}_BFB_cycles.txt', res['new_segments'], res['bfb_strings'], res['scores'], res['multiplicity'])
        print(f'Generated {region_prefix}_BFB_graph.txt and {region_prefix}_BFB_cycles.txt')
        visualize_BFB(
            cycle_file=f'{region_prefix}_BFB_cycles.txt',
            graph_file=f'{region_prefix}_BFB_graph.txt',
            cnr_file=None,
            output_prefix=f'{region_prefix}_BFB',
            gene_annotation=gene,
            multiple=multiple
        )
    logger.info(f'Total time: {time.time() - start_time:.1f} seconds')

def _parse_region_string(region_str):
    """Parse 'chrN:start-end' into (chrom, int_start, int_end). Supports M/K suffixes."""
    def _parse_pos(s):
        s = s.strip()
        if s.upper().endswith('M'):
            return int(float(s[:-1]) * 1_000_000)
        if s.upper().endswith('K'):
            return int(float(s[:-1]) * 1_000)
        return int(s)
    chrom, rest = region_str.split(':')
    start_str, end_str = rest.split('-')
    return (chrom, _parse_pos(start_str), _parse_pos(end_str))

def main():
    parser = argparse.ArgumentParser(
        description="BFBArchitect for detecting and reconstructing BFB sequences in an amplicon region.")
    parser.add_argument("--graph", help="Path to an AA-format _graph.txt file.", default=None)
    parser.add_argument("--whole_graph", help="Treat all segments as one region.", action='store_true')
    parser.add_argument("-g", "--gene", help="Gene annotation GTF file for visualization.", default=None)
    parser.add_argument("--bam", help="Path to a sorted bam file.", default=None)
    parser.add_argument("--cns", help="Path to a sorted cns file.", default=None)
    parser.add_argument("--region", help="The amplified region (chr:start-end). For --graph: process this specific region only.", default=None)
    parser.add_argument("--segmentation", help="Consider CNV in segmentation", action='store_true')
    parser.add_argument("--deletion", help="Deletion handling", action='store_true')
    parser.add_argument("--coverage", help="Sequencing coverage.", type=float, default=None)
    parser.add_argument("--no_expansion", help="Keep the specified region without expansion", action='store_true')
    parser.add_argument("--min_sv_cn", type=float, default=0.75, help="Minimum copy number for SV calling.")
    parser.add_argument("--min_mapq", type=int, default=20, help="Minimum mapping quality for SV calling.")
    parser.add_argument("--output_prefix", help="Prefix of output files.", required=True)
    parser.add_argument("--multiple", help="Reconstruct multiple BFB candidates", action='store_true')
    parser.add_argument("--solver", help="ILP solver to use. Options: gurobi (defualt), mosek, and cbc", default=None)
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads for the ILP solver (default: 8).")
    parser.add_argument("--centromere", help="Path to a BED file of centromere regions.", default=None)
    parser.add_argument("--verbose", help="Print all log messages to stdout in addition to the log file.", action='store_true')
    parser.add_argument("--min-lp-bound", type=float, default=None, dest='min_lp_bound',
                        help="Stop early if the root LP bound exceeds this value (Gurobi only). "
                             "Uncapped by default. A value such as 25 can prevent wasting time on clear non-BFBs.")
    args = parser.parse_args()
    centromere_dict = build_centromere_dict(args.centromere)
    if args.graph:
        if args.whole_graph and args.region:
            parser.error("--whole_graph and --region are mutually exclusive.")
        parsed_region = _parse_region_string(args.region) if args.region else None
        run_bfb_from_graph(args.graph, args.output_prefix, args.multiple, args.solver,
                           args.whole_graph, region=parsed_region, gene=args.gene,
                           centromere_dict=centromere_dict, verbose=args.verbose, threads=args.threads,
                           min_lp_bound=args.min_lp_bound)
    elif args.bam and args.cns and args.region:
        reconstruct_bfb_from_bam(args.bam, args.cns, args.region, args.output_prefix, args.segmentation, args.deletion, args.coverage, args.multiple, args.no_expansion, args.min_sv_cn, args.min_mapq, args.solver, centromere_dict, verbose=args.verbose, threads=args.threads, min_lp_bound=args.min_lp_bound)
    else:
        parser.error("Provide either --graph or all of --bam, --cns, --region.")

if __name__ == "__main__":
    main()
