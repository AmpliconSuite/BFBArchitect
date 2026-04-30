import re
from collections import defaultdict

try:
    from bfbarchitect.datatypes import SV, CHR_CENTRO
except ImportError:
    from datatypes import SV, CHR_CENTRO


# ── helpers (ported from ampclasslib/bfb_regions.py) ─────────────────────────

def _weighted_avg_cn(group):
    """Weighted-average CN for a group of (start, end, cn) segments."""
    total = sum(e - s for s, e, _ in group)
    return sum(cn * (e - s) for s, e, cn in group) / total if total > 0 else 0.0


def _merge_similar_segments(segs, gap_cutoff, cn_tol):
    """
    Merge consecutive segments whose gap is <= gap_cutoff bp AND whose CN
    differs from the current group's weighted-average CN by <= cn_tol.
    Returns a new sorted list of (start, end, cn) tuples.
    """
    if not segs:
        return []
    groups = [[segs[0]]]
    for s, e, cn in segs[1:]:
        prev_e = groups[-1][-1][1]
        prev_cn = _weighted_avg_cn(groups[-1])
        if s - prev_e <= gap_cutoff and abs(cn - prev_cn) <= cn_tol:
            groups[-1].append((s, e, cn))
        else:
            groups.append([(s, e, cn)])
    return [(g[0][0], g[-1][1], _weighted_avg_cn(g)) for g in groups]


def _find_gap_index(fb_p1, fb_p2, segs):
    """
    Return index i such that the foldback midpoint falls in the gap
    [segs[i].end, segs[i+1].start]. Returns None if no such gap exists.
    """
    fb_mid = (fb_p1 + fb_p2) / 2.0
    for i in range(len(segs) - 1):
        if segs[i][1] <= fb_mid <= segs[i + 1][0]:
            return i
    return None


# ── graph file parser ─────────────────────────────────────────────────────────

def parse_graph_file(graph_file):
    """
    Parse sequence edges and discordant breakpoint edges from an AA _graph.txt file.

    Returns
    -------
    svs : list of (SV, cn_float, read_count)
    chrom_segs : dict[chrom, list[(start, end, cn_float, coverage, read_count)]]
        Segments sorted by start position per chromosome.
    """
    svs = []
    chrom_segs = defaultdict(list)
    with open(graph_file) as f:
        for line in f:
            if line.startswith('sequence'):
                parts = line.split()
                chrom = parts[1].split(':')[0]
                start = int(parts[1].split(':')[1][:-1])   # strip trailing '-'
                end   = int(parts[2].split(':')[1][:-1])   # strip trailing '+'
                cn       = float(parts[3])
                coverage = float(parts[4])
                read_count = int(parts[6])
                chrom_segs[chrom].append((start, end, cn, coverage, read_count))
            elif line.startswith('discordant'):
                parts = line.split()
                m = re.match(r"(\w+):(\d+)([+-])->(\w+):(\d+)([+-])", parts[1])
                if not m:
                    continue
                c1, p1, s1, c2, p2, s2 = m.groups()
                sv = SV(c1, int(p1), s1, c2, int(p2), s2)
                cn = float(parts[2])
                rc = int(parts[3]) if len(parts) > 3 else 0
                svs.append((sv, cn, rc))
    for segs in chrom_segs.values():
        segs.sort(key=lambda x: x[0])
    return svs, chrom_segs


# ── BFB candidate region detection ───────────────────────────────────────────

def find_bfb_candidate_regions(graph_file, min_seg_size=50000, min_boundary_seg_size=10000,
                                fb_dist_cut=25000, merge_gap=50000, merge_cn_tol=0.5,
                                min_cn_step=0.5, merge_padding=150000):
    """
    Identify BFB-like candidate regions from an AA breakpoint graph file.

    Two detection criteria are applied after segment pre-processing:

    Criterion 1 – Monotonic CN triplet + foldback:
      Three consecutive large segments with copy numbers that increase or decrease
      by >= min_cn_step at each step, with at least one foldback SV whose both
      endpoints lie within the span of those three segments.

    Criterion 2 – Opposite-direction foldback pair:
      A '++' foldback and a '--' foldback on the same chromosome with exactly one
      merged segment between them.

    Parameters
    ----------
    graph_file : str
        Path to an AA-format breakpoint graph file.
    min_seg_size : int
        Minimum bp for non-boundary segments (default 50000).
    min_boundary_seg_size : int
        Minimum bp for first/last segment per chromosome (default 10000).
    fb_dist_cut : int
        Max distance between foldback SV endpoints (default 25000).
    merge_gap : int
        Max gap in bp between segments to merge (default 50000).
    merge_cn_tol : float
        Max CN difference for merging adjacent segments (default 0.5).
    min_cn_step : float
        Min CN change between consecutive segments in a monotonic triplet (default 0.5).
    merge_padding : int
        Padding added to each side of a candidate before merging overlapping
        regions (default 150000).

    Returns
    -------
    list of (chrom, start, end) tuples in half-open [start, end) coordinates.
    """
    svs, chrom_segs = parse_graph_file(graph_file)

    # ── segment pre-processing ────────────────────────────────────────────────
    chrom_merged = {}
    for chrom, raw_segs in chrom_segs.items():
        all_segs = [(s, e, cn) for s, e, cn, *_ in raw_segs]
        n_all = len(all_segs)
        filtered = [
            (s, e, cn)
            for i, (s, e, cn) in enumerate(all_segs)
            if e - s >= (min_boundary_seg_size if i == 0 or i == n_all - 1 else min_seg_size)
        ]
        if not filtered:
            continue
        chrom_merged[chrom] = _merge_similar_segments(filtered, merge_gap, merge_cn_tol)

    # ── foldback SVs grouped by chromosome ───────────────────────────────────
    fb_by_chrom = defaultdict(list)
    for sv, _cn, _rc in svs:
        if sv.is_foldback(max_distance=fb_dist_cut):
            # After sort_breakpoints(): bp1 <= bp2, strand1 == strand2
            fb_by_chrom[sv.chrom1].append((sv.bp1, sv.bp2, sv.strand1))

    candidates = []

    for chrom, segs in chrom_merged.items():
        n = len(segs)
        fb_list = fb_by_chrom.get(chrom, [])

        # ── Criterion 1: monotonic CN triplet with a foldback inside ─────────
        if n >= 3:
            for i in range(n - 2):
                s_a, e_a, cn_a = segs[i]
                s_b, e_b, cn_b = segs[i + 1]
                s_c, e_c, cn_c = segs[i + 2]
                ascending  = cn_b - cn_a >= min_cn_step and cn_c - cn_b >= min_cn_step
                descending = cn_a - cn_b >= min_cn_step and cn_b - cn_c >= min_cn_step
                if not (ascending or descending):
                    continue
                region_start, region_end = s_a, e_c
                for p1, p2, _ in fb_list:
                    if p1 >= region_start and p2 <= region_end:
                        candidates.append((chrom, region_start, region_end))
                        break

        # ── Criterion 2: ++ / -- foldback pair with one segment between ──────
        if n >= 3:
            fb_pp = [(p1, p2) for p1, p2, d in fb_list if d == '+']
            fb_mm = [(p1, p2) for p1, p2, d in fb_list if d == '-']

            for fba in fb_pp:
                for fbb in fb_mm:
                    mid_a = (fba[0] + fba[1]) / 2.0
                    mid_b = (fbb[0] + fbb[1]) / 2.0
                    fb_left  = fba if mid_a <= mid_b else fbb
                    fb_right = fbb if mid_a <= mid_b else fba

                    gap_l = _find_gap_index(fb_left[0],  fb_left[1],  segs)
                    gap_r = _find_gap_index(fb_right[0], fb_right[1], segs)

                    if gap_l is None or gap_r is None:
                        continue
                    if gap_r != gap_l + 1:
                        continue

                    outer_l = gap_l
                    outer_r = gap_r + 1
                    if outer_r >= n:
                        continue
                    candidates.append((chrom, segs[outer_l][0], segs[outer_r][1]))

    if not candidates:
        return []

    # ── pad and merge overlapping candidates per chromosome ───────────────────
    chrom_intervals = defaultdict(list)
    for chrom, start, end in candidates:
        chrom_intervals[chrom].append((max(0, start - merge_padding), end + merge_padding))

    result = []
    for chrom, intervals in chrom_intervals.items():
        intervals.sort()
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                result.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        result.append((chrom, cur_s, cur_e))

    result.sort()
    return result


# ── per-region segmentation ───────────────────────────────────────────────────

def subsect_graph_for_region(graph_file, regions, fb_dist_cut=50000, cn_tol=1.0,
                              small_seg_size=50000):
    """
    Extract and resegment AA graph data for each BFB candidate region, returning
    data in the format expected by BFBArchitect's reconstruct_bfb().

    Processing pipeline (per region):
    1. Extract segments overlapping the region.
    2. Merge adjacent segments within cn_tol CN units, respecting foldback cut points
       as hard merge boundaries.
    3. Absorb segments < small_seg_size bp not touching any foldback endpoint into
       a neighboring segment.
    4. Re-segment at each foldback's outer cut point.
    5. Build the cn/lf/rf vectors required by reconstruct_bfb().

    Note: Caller is responsible for chromosome-arm trimming before invoking
    reconstruct_bfb() (remove the outermost background segment on the non-BFB side).

    Parameters
    ----------
    graph_file : str
        Path to an AA _graph.txt file.
    regions : list of (chrom, start, end)
        As returned by find_bfb_candidate_regions().
    fb_dist_cut : int
        Max bp distance between foldback SV endpoints (default 50000).
    cn_tol : float
        Max CN difference for merging adjacent segments (default 1.0).
    small_seg_size : int
        Segments < this many bp not touching a foldback endpoint are absorbed
        into a neighbor (default 50000).

    Returns
    -------
    list aligned with `regions`. Each entry is None (no segments found) or a
    tuple (new_segments, cn, lf, rf, region_svs, sv_info):
        new_segments  list of (chrom, start, end, cn_float, coverage, read_count)
        cn            list[int]  round(seg_cn) - 1 per segment
        lf            list[int]  left-foldback SV CN per segment
        rf            list[int]  right-foldback SV CN per segment
        region_svs    list[SV]   all discordant SVs in the region
        sv_info       dict[SV, (cn_float, read_count)]
    """
    svs, chrom_segs = parse_graph_file(graph_file)

    # Pre-group foldback SVs by chromosome to avoid a full scan per region.
    sv_cn_map = {sv: cn for sv, cn, _rc in svs}
    foldbacks_by_chrom = defaultdict(list)
    for sv, _cn, _rc in svs:
        if sv.is_foldback(max_distance=fb_dist_cut):
            foldbacks_by_chrom[sv.chrom1].append(sv)

    results = []
    for chrom, region_start, region_end in regions:
        raw_segs = [
            (s, e, cn, cov, rc)
            for s, e, cn, cov, rc in chrom_segs.get(chrom, [])
            if s < region_end and e > region_start
        ]

        if not raw_segs:
            results.append(None)
            continue

        seg_span_lo = raw_segs[0][0]
        seg_span_hi = raw_segs[-1][1]

        # Keep foldbacks with at least one endpoint inside the extracted segment span.
        foldback_svs = [
            sv for sv in foldbacks_by_chrom[chrom]
            if (seg_span_lo <= sv.bp1 <= seg_span_hi
                or seg_span_lo <= sv.bp2 <= seg_span_hi)
        ]

        # All discordant SVs in/near the region for graph output
        region_svs = []
        sv_info = {}
        for sv, cn, rc in svs:
            flag1, flag2 = sv.is_in_region((chrom, region_start, region_end))
            if flag1 or flag2:
                region_svs.append(sv)
                sv_info[sv] = (cn, rc)

        # After sort_breakpoints(): bp1 <= bp2, strand1 == strand2 for foldbacks.
        # Both endpoints guard small-segment smoothing.
        # One cut point per foldback (outer endpoint) drives re-segmentation.
        fb_endpoints = set()
        cut_points = set()
        for sv in foldback_svs:
            fb_endpoints.update((sv.bp1, sv.bp2))
            if sv.strand1 == '+':    # right (++) foldback: outer endpoint is bp2
                cut_points.add(sv.bp2)
            else:                     # left (--) foldback: outer endpoint is bp1 - 1
                cut_points.add(sv.bp1 - 1)

        segs = list(raw_segs)

        # ── Step 1: merge adjacent CN-similar segments ─────────────────────────
        # Foldback cut points are hard merge boundaries.
        changed = True
        while changed:
            changed = False
            new_segs = []
            i = 0
            while i < len(segs):
                if i < len(segs) - 1:
                    s1, e1, cn1, cov1, rc1 = segs[i]
                    s2, e2, cn2, cov2, rc2 = segs[i + 1]
                    if e1 not in cut_points and abs(cn1 - cn2) <= cn_tol:
                        l1, l2 = e1 - s1 + 1, e2 - s2 + 1
                        new_segs.append((s1, e2,
                                         (cn1 * l1 + cn2 * l2) / (l1 + l2),
                                         (cov1 * l1 + cov2 * l2) / (l1 + l2),
                                         rc1 + rc2))
                        i += 2
                        changed = True
                    else:
                        new_segs.append(segs[i])
                        i += 1
                else:
                    new_segs.append(segs[i])
                    i += 1
            segs = new_segs

        # ── Step 2: absorb small non-foldback segments into neighbors ──────────
        changed = True
        while changed:
            changed = False
            new_segs = []
            i = 0
            while i < len(segs):
                s, e, seg_cn, seg_cov, seg_rc = segs[i]
                size = e - s + 1
                touches = s in fb_endpoints or e in fb_endpoints
                if size < small_seg_size and not touches:
                    can_left  = bool(new_segs) and new_segs[-1][1] not in cut_points
                    can_right = i + 1 < len(segs) and e not in cut_points
                    if can_left:
                        ps, pe, pcn, pcov, prc = new_segs.pop()
                        pl = pe - ps + 1
                        new_segs.append((ps, e,
                                         (pcn * pl + seg_cn * size) / (pl + size),
                                         (pcov * pl + seg_cov * size) / (pl + size),
                                         prc + seg_rc))
                        changed = True
                    elif can_right:
                        ns, ne, ncn, ncov, nrc = segs[i + 1]
                        nl = ne - ns + 1
                        new_segs.append((s, ne,
                                         (seg_cn * size + ncn * nl) / (size + nl),
                                         (seg_cov * size + ncov * nl) / (size + nl),
                                         seg_rc + nrc))
                        i += 1
                        changed = True
                    else:
                        new_segs.append((s, e, seg_cn, seg_cov, seg_rc))
                else:
                    new_segs.append((s, e, seg_cn, seg_cov, seg_rc))
                i += 1
            segs = new_segs

        # ── Step 3: re-segment at foldback outer cut points ────────────────────
        new_segments = []
        cur_start = segs[0][0]
        total_length = 0
        weighted_cn_sum = 0.0
        weighted_cov_sum = 0.0
        total_rc = 0

        for seg_start, seg_end, seg_cn, seg_cov, seg_rc in segs:
            seg_size = seg_end - seg_start + 1
            total_length += seg_size
            weighted_cn_sum += seg_cn * seg_size
            weighted_cov_sum += seg_cov * seg_size
            total_rc += seg_rc
            if seg_end in cut_points:
                new_segments.append((chrom, cur_start, seg_end,
                                     weighted_cn_sum / total_length,
                                     weighted_cov_sum / total_length,
                                     total_rc))
                cur_start = seg_end + 1
                total_length = 0
                weighted_cn_sum = 0.0
                weighted_cov_sum = 0.0
                total_rc = 0

        if total_length > 0:
            new_segments.append((chrom, cur_start, segs[-1][1],
                                 weighted_cn_sum / total_length,
                                 weighted_cov_sum / total_length,
                                 total_rc))

        # ── Step 4: build cn/lf/rf vectors ─────────────────────────────────────
        cn_vals = [round(seg[3]) - 1 for seg in new_segments]
        l_bp_idx = {seg[1]: i for i, seg in enumerate(new_segments)}
        r_bp_idx = {seg[2]: i for i, seg in enumerate(new_segments)}
        lf = [0] * len(cn_vals)
        rf = [0] * len(cn_vals)

        for sv in foldback_svs:
            sv_cn = sv_cn_map.get(sv, 0.0)
            if sv.strand1 == '-':    # left (--) foldback: segment starts at bp1
                if sv.bp1 in l_bp_idx:
                    lf[l_bp_idx[sv.bp1]] += round(sv_cn)
            else:                     # right (++) foldback: segment ends at bp2
                if sv.bp2 in r_bp_idx:
                    rf[r_bp_idx[sv.bp2]] += round(sv_cn)

        results.append((new_segments, cn_vals, lf, rf, region_svs, sv_info))

    return results


# ── whole-graph fallback ──────────────────────────────────────────────────────

def whole_graph_as_region(graph_file, centromere_dict=None):
    """
    Treat all segments in the graph as a single BFB region (the --whole_graph
    fallback).

    Returns
    -------
    (new_segments, cn, lf, rf, svs_list, sv_info, primary_chrom)
        new_segments  list of (chrom, start, end, cn_float, coverage, read_count)
        cn            list[int]  round(seg_cn) - 1 per segment
        lf            list[int]  left-foldback SV CN per segment
        rf            list[int]  right-foldback SV CN per segment
        svs_list      list[SV]   all discordant SVs
        sv_info       dict[SV, (cn_float, read_count)]
        primary_chrom str        chromosome of the primary BFB region
    """
    svs_raw, chrom_segs = parse_graph_file(graph_file)

    svs_list = [sv for sv, _cn, _rc in svs_raw]
    sv_info  = {sv: (cn, rc) for sv, cn, rc in svs_raw}

    # Derive the amplified region from the sequence edges
    chrom_ranges = {}
    for chrom, segs in chrom_segs.items():
        starts = [s for s, *_ in segs]
        ends   = [e for _, e, *_ in segs]
        chrom_ranges[chrom] = (min(starts), max(ends))
    primary_chrom = next(iter(chrom_ranges))
    region_start, region_end = chrom_ranges[primary_chrom]
    region = (primary_chrom, region_start, region_end)

    # Foldback cut points
    breakpoints = set()
    for sv in svs_list:
        flag1, flag2 = sv.is_in_region(region)
        if sv.is_foldback() and flag1 and flag2:
            if sv.strand1 == '+':
                breakpoints.add(sv.bp2)
            else:
                breakpoints.add(sv.bp1 - 1)

    # Sort original segments and build new_segments via weighted average
    orig_segs = sorted(chrom_segs[primary_chrom], key=lambda x: x[0])
    new_segments = []
    chrom = primary_chrom
    start = region_start
    total_length, weighted_cn_sum = 0, 0
    total_cov, total_rc = 0.0, 0
    for seg_start, seg_end, seg_cn, seg_cov, seg_rc in orig_segs:
        if seg_end in breakpoints:
            if total_length > 0:
                new_cn  = weighted_cn_sum / total_length
                new_cov = total_cov / total_length
                new_segments.append((chrom, start, seg_end, new_cn, new_cov, total_rc))
            start = seg_end + 1
            total_length, weighted_cn_sum = 0, 0
            total_cov, total_rc = 0.0, 0
        else:
            seg_size = seg_end - seg_start + 1
            total_length += seg_size
            weighted_cn_sum += seg_cn * seg_size
            total_cov += seg_cov * seg_size
            total_rc += seg_rc
    if total_length > 0:
        new_segments.append((chrom, start, region_end,
                             weighted_cn_sum / total_length,
                             total_cov / total_length, total_rc))

    # Arm trimming: remove the segment on the centromere-proximal side
    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    centro = centromere_dict.get(primary_chrom)
    if new_segments:
        if centro is None or new_segments[-1][2] < centro:
            new_segments.pop(0)   # p-arm: remove leftmost (centromere is to the right)
        else:
            new_segments.pop(-1)  # q-arm: remove rightmost (centromere is to the left)

    # Build cn/lf/rf vectors
    cn = [round(seg[3]) - 1 for seg in new_segments]
    l_bp = [seg[1] for seg in new_segments]
    r_bp = [seg[2] for seg in new_segments]
    lf = [0] * len(cn)
    rf = [0] * len(cn)
    for sv in svs_list:
        flag1, flag2 = sv.is_in_region(region)
        if sv.is_foldback() and flag1 and flag2:
            if sv.strand1 == '-' and sv.bp1 in l_bp:
                lf[l_bp.index(sv.bp1)] += round(sv_info[sv][0])
            elif sv.strand1 == '+' and sv.bp2 in r_bp:
                rf[r_bp.index(sv.bp2)] += round(sv_info[sv][0])

    return new_segments, cn, lf, rf, svs_list, sv_info, primary_chrom
