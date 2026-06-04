import re
from collections import defaultdict, deque

try:
    from bfbarchitect.datatypes import SV, CHR_CENTRO
except ImportError:
    from datatypes import SV, CHR_CENTRO

MAX_FOLDBACKS_FOR_GRAPH_RECON = 50


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


def _find_gap_index(point, segs):
    """
    Return index i such that point falls in the gap
    [segs[i].end, segs[i+1].start]. Returns None if no such gap exists.
    """
    for i in range(len(segs) - 1):
        if segs[i][1] <= point <= segs[i + 1][0]:
            return i
    return None


def _foldback_cn_count(sv_cn):
    """Count every observed graph foldback at least once in lf/rf vectors."""
    return max(1, round(sv_cn))


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


# ── TST-jump foldback detection ──────────────────────────────────────────────

def _build_graph_lookups(svs, chrom_segs):
    """Build (chrom, pos, dir) → segment and → SV-entry lookup tables."""
    bp_to_seg = {}
    for chrom, segs in chrom_segs.items():
        for start, end, *_ in segs:
            bp_to_seg[(chrom, start, '-')] = (start, end)
            bp_to_seg[(chrom, end,   '+')] = (start, end)

    sv_index = defaultdict(list)
    for entry in svs:
        sv = entry[0]
        sv_index[(sv.chrom1, sv.bp1, sv.strand1)].append(entry)
        sv_index[(sv.chrom2, sv.bp2, sv.strand2)].append(entry)

    return bp_to_seg, sv_index


def _shard_path_reachable(start, target, bp_to_seg, sv_index, shard_max_bp, max_hops, exclude_svs=None):
    """
    BFS from start to target through shard-sized segments only.

    Movements allowed:
      - Traverse a shard (<= shard_max_bp) from one end to the other (costs 1 hop).
      - Concordant edge: move between adjacent segment endpoints (pos+ <-> (pos+1)-),
        free if it reaches target or leads into another shard.
      - Discordant SV edge: free if it reaches target or leads into a shard endpoint.

    Returns list of (chrom, start, end) shard tuples traversed, or None.
    """
    queue = deque([(start, [], frozenset())])
    seen = set()

    while queue:
        curr, path, visited = queue.popleft()

        if curr == target:
            return path

        key = (curr, visited)
        if key in seen:
            continue
        seen.add(key)

        if len(path) >= max_hops:
            continue

        chrom, pos, d = curr

        # Traverse the shard whose endpoint is curr
        seg = bp_to_seg.get(curr)
        if seg is not None:
            s, e = seg
            seg_id = (chrom, s, e)
            if e - s <= shard_max_bp and seg_id not in visited:
                other = (chrom, e, '+') if d == '-' else (chrom, s, '-')
                queue.append((other, path + [seg_id], visited | {seg_id}))

        # Concordant edge (AA convention: pos+ connects to (pos+1)-)
        adj = (chrom, pos + 1, '-') if d == '+' else (chrom, pos - 1, '+')
        if adj == target:
            queue.append((adj, path, visited))
        else:
            adj_seg = bp_to_seg.get(adj)
            if adj_seg is not None and adj_seg[1] - adj_seg[0] <= shard_max_bp:
                queue.append((adj, path, visited))

        # Discordant SV edges: only follow into shard-sized segments.
        # Do NOT shortcut directly to target — that would accept any far-jumping
        # SV that happens to land on the target, producing false positives.
        for entry in sv_index.get(curr, []):
            sv = entry[0]
            if exclude_svs and sv in exclude_svs:
                continue
            other = (sv.chrom2, sv.bp2, sv.strand2) \
                    if (sv.chrom1, sv.bp1, sv.strand1) == curr \
                    else (sv.chrom1, sv.bp1, sv.strand1)
            other_seg = bp_to_seg.get(other)
            if other_seg is not None and other_seg[1] - other_seg[0] <= shard_max_bp:
                queue.append((other, path, visited))

    return None


def find_tst_foldbacks(svs, chrom_segs, shard_max_bp=5000, max_hops=5,
                       fb_dist=50000, far_min=500000, verbose=False):
    """
    Detect TST-jump (template-switching) foldbacks and inject synthetic FBI SVs.

    A TST foldback is formed when two far-jumping SVs (inter-chromosomal or
    >far_min bp) have local breakends on the same chromosome within fb_dist of
    each other with opposite orientations (+/-), and their far breakends are
    connected through a chain of shard-sized segments (<= shard_max_bp).

    For each confirmed pair a synthetic '++' FBI SV is created so the normal
    lf/rf vector build in subsect_graph_for_region picks it up.  The synthetic
    SV has .TST = True.

    Two local-breakend topologies are handled:
      Topology A  lo='+', hi='-'  breakpoints are outer boundaries of the shard
                  → FBI: (pos_lo, +, pos_hi-1, +)
      Topology B  lo='-', hi='+'  breakpoints are the shard's own endpoints
                  → FBI: (pos_lo-1, +, pos_hi, +)

    Returns the input svs list extended with any (synthetic_SV, cn, 0) tuples.
    """
    bp_to_seg, sv_index = _build_graph_lookups(svs, chrom_segs)

    # Pre-compute cut-point positions of direct (non-TST) foldbacks so we can
    # skip far-side injections that would duplicate an existing foldback's cut.
    # ++ foldbacks cut at bp2; -- foldbacks cut at bp1.
    existing_fb_cuts = set()
    # Also build a per-chromosome list of real foldback (bp1, bp2) spans for the
    # TST duplicate check below.
    real_fb_spans = defaultdict(list)
    for sv, _cn, _rc in svs:
        if sv.is_foldback():
            if sv.strand1 == '+':
                existing_fb_cuts.add((sv.chrom1, sv.bp2))
            else:
                existing_fb_cuts.add((sv.chrom1, sv.bp1))
            real_fb_spans[sv.chrom1].append((sv.bp1, sv.bp2))

    # Collect far-jumping breakends indexed by local chromosome.
    # An SV contributes two entries (once per end) so both ends are candidates
    # for being the "local" side.
    far_by_chrom = defaultdict(list)
    for entry in svs:
        sv, cn, rc = entry
        for local, far in [
            ((sv.chrom1, sv.bp1, sv.strand1), (sv.chrom2, sv.bp2, sv.strand2)),
            ((sv.chrom2, sv.bp2, sv.strand2), (sv.chrom1, sv.bp1, sv.strand1)),
        ]:
            lc, lp, ld = local
            fc, fp, _  = far
            if fc != lc or abs(fp - lp) >= far_min:
                far_by_chrom[lc].append((lp, ld, far, sv, cn, rc))

    seen_pairs = set()
    synthetic = []

    for chrom, breakends in far_by_chrom.items():
        for i, (pos_i, dir_i, far_i, sv_i, cn_i, rc_i) in enumerate(breakends):
            for pos_j, dir_j, far_j, sv_j, cn_j, rc_j in breakends[i + 1:]:
                if sv_i is sv_j:                       # same SV, skip
                    continue
                if abs(pos_i - pos_j) > fb_dist:       # too far apart locally
                    continue
                if dir_i == dir_j:                     # need opposite orientations
                    continue
                if far_i == far_j:                     # degenerate: same far end
                    continue

                pair_key = (min(id(sv_i), id(sv_j)), max(id(sv_i), id(sv_j)))
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)

                fc_i, fp_i, fd_i = far_i
                fc_j, fp_j, fd_j = far_j

                # Skip: the far ends are already the direct endpoints of a real
                # foldback SV.  The BFS would "find" a path through that foldback's
                # own hairpin segment — the fold is already represented and a
                # synthetic FBI would be a duplicate.
                if fc_i == fc_j:
                    p_lo, p_hi = min(fp_i, fp_j), max(fp_i, fp_j)
                    if any(abs(b1 - p_lo) <= 5 and abs(b2 - p_hi) <= 5
                           for b1, b2 in real_fb_spans.get(fc_i, [])):
                        if verbose:
                            cn_ratio = min(cn_i, cn_j) / max(cn_i, cn_j) if max(cn_i, cn_j) > 0 else 0
                            print(f"[TST] Candidate {chrom}:{pos_i}({dir_i}) <-> "
                                  f"{chrom}:{pos_j}({dir_j}) "
                                  f"[{abs(pos_i - pos_j)} bp apart]  "
                                  f"skipped (far ends are direct foldback endpoints)")
                            print(f"       SV1: {sv_i}  CN={cn_i:.2f}  rc={rc_i}")
                            print(f"       SV2: {sv_j}  CN={cn_j:.2f}  rc={rc_j}  "
                                  f"CN_ratio={cn_ratio:.2f}")
                            print(f"       Far ends: {fc_i}:{fp_i}({fd_i})  "
                                  f"and  {fc_j}:{fp_j}({fd_j})")
                        continue

                path = _shard_path_reachable(far_i, far_j, bp_to_seg, sv_index,
                                             shard_max_bp, max_hops,
                                             exclude_svs={sv_i, sv_j})

                if verbose:
                    status = "confirmed" if path is not None else "rejected (no shard path)"
                    print(f"[TST] Candidate {chrom}:{pos_i}({dir_i}) <-> {chrom}:{pos_j}({dir_j}) "
                          f"[{abs(pos_i - pos_j)} bp apart]  {status}")
                    cn_ratio = min(cn_i, cn_j) / max(cn_i, cn_j) if max(cn_i, cn_j) > 0 else 0
                    print(f"       SV1: {sv_i}  CN={cn_i:.2f}  rc={rc_i}")
                    print(f"       SV2: {sv_j}  CN={cn_j:.2f}  rc={rc_j}  CN_ratio={cn_ratio:.2f}")
                    print(f"       Far ends: {fc_i}:{fp_i}({fd_i})  and  {fc_j}:{fp_j}({fd_j})")
                    if path is not None:
                        if path:
                            total_shard_bp = sum(e - s for _, s, e in path)
                            hops = ' -> '.join(f"{c}:{s}-{e} [{e-s} bp]" for c, s, e in path)
                            print(f"       Shard path ({len(path)} hop{'s' if len(path) != 1 else ''}, "
                                  f"{total_shard_bp} bp total): {hops}")
                        else:
                            print(f"       Shard path: concordant connection (0 shard hops)")

                if path is None:
                    continue

                pos_lo, dir_lo = (pos_i, dir_i) if pos_i <= pos_j else (pos_j, dir_j)
                pos_hi          = pos_j            if pos_i <= pos_j else pos_i

                # Topology A (lo='+', hi='-'): shard fills the gap between them
                # Topology B (lo='-', hi='+'): breakpoints are the shard's endpoints
                if dir_lo == '+':
                    fbi_bp1, fbi_bp2 = pos_lo, pos_hi - 1
                else:
                    fbi_bp1, fbi_bp2 = pos_lo - 1, pos_hi

                if fbi_bp2 < fbi_bp1:
                    if verbose:
                        print(f"       WARNING: degenerate coordinates ({fbi_bp1}, {fbi_bp2}), skipping")
                    continue

                cn_tst = min(cn_i, cn_j)
                synth = SV(chrom, fbi_bp1, '+', chrom, fbi_bp2, '+')
                synth.TST = True

                if verbose:
                    print(f"       -> Injecting synthetic FBI: {chrom}:{fbi_bp1}(+) <-> "
                          f"{chrom}:{fbi_bp2}(+)  CN={cn_tst:.2f}")

                synthetic.append((synth, cn_tst, 0))

                # If the far ends share a chromosome and are within foldback
                # distance, the fold also lives there — inject an FBI on that
                # side so subsect_graph_for_region picks it up for that chrom.
                if fc_i == fc_j and abs(fp_i - fp_j) <= fb_dist:
                    synth_far = None
                    if fd_i == fd_j:
                        synth_far = SV(fc_i, min(fp_i, fp_j), fd_i,
                                       fc_i, max(fp_i, fp_j), fd_i)
                    else:
                        fp_lo, fp_hi = (fp_i, fp_j) if fp_i <= fp_j else (fp_j, fp_i)
                        fd_lo = fd_i if fp_i <= fp_j else fd_j
                        if fd_lo == '+':
                            far_bp1, far_bp2 = fp_lo, fp_hi - 1
                        else:
                            far_bp1, far_bp2 = fp_lo - 1, fp_hi
                        if far_bp2 >= far_bp1:
                            synth_far = SV(fc_i, far_bp1, '+', fc_i, far_bp2, '+')
                    if synth_far is not None:
                        # Skip if a direct foldback already covers the same cut point.
                        far_cut = (fc_i, synth_far.bp2 if synth_far.strand1 == '+'
                                   else synth_far.bp1)
                        if far_cut in existing_fb_cuts:
                            if verbose:
                                print(f"       -> Far-side FBI skipped: direct foldback "
                                      f"already at {fc_i} cut {far_cut[1]}")
                        else:
                            synth_far.TST = True
                            synthetic.append((synth_far, cn_tst, 0))
                            if verbose:
                                d = synth_far.strand1
                                print(f"       -> Also injecting far-side FBI: "
                                      f"{fc_i}:{synth_far.bp1}({d}) <-> "
                                      f"{fc_i}:{synth_far.bp2}({d})  CN={cn_tst:.2f}")

    # Deduplicate: multiple SV pairs may produce the same synthetic FBI at the
    # same coordinates.  Keep the entry with the highest CN.
    deduped = {}
    for entry in synthetic:
        sv, cn, rc = entry
        key = str(sv)
        if key not in deduped or cn > deduped[key][1]:
            deduped[key] = entry
    return svs + list(deduped.values())


# ── BFB candidate region detection ───────────────────────────────────────────

def find_bfb_candidate_regions(graph_file, min_seg_size=50000, min_boundary_seg_size=10000,
                                fb_dist_cut=50000, merge_gap=50000, merge_cn_tol=0.5,
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
        Max distance between foldback SV endpoints (default 50000).
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
                if s_b - e_a > merge_gap or s_c - e_b > merge_gap:
                    continue
                ascending  = cn_b - cn_a >= min_cn_step and cn_c - cn_b >= min_cn_step
                descending = cn_a - cn_b >= min_cn_step and cn_b - cn_c >= min_cn_step
                if not (ascending or descending):
                    continue
                region_start, region_end = s_a, e_c
                for p1, p2, _ in fb_list:
                    if p1 >= region_start - merge_gap and p2 <= region_end + merge_gap:
                        candidates.append((chrom, region_start, region_end))
                        break

        # ── Criterion 2: ++ / -- foldback pair with one segment between ──────
        if n >= 3:
            fb_pp = [(p1, p2) for p1, p2, d in fb_list if d == '+']
            fb_mm = [(p1, p2) for p1, p2, d in fb_list if d == '-']

            for fba in fb_pp:
                for fbb in fb_mm:
                    gap_a = _find_gap_index(fba[0], segs)  # ++ foldback: use bp1
                    gap_b = _find_gap_index(fbb[1], segs)  # -- foldback: use bp2

                    if gap_a is None or gap_b is None:
                        continue
                    gap_l, gap_r = min(gap_a, gap_b), max(gap_a, gap_b)
                    if gap_r != gap_l + 1:
                        continue

                    outer_r = gap_r + 1
                    if outer_r >= n:
                        continue
                    candidates.append((chrom, segs[gap_l][0], segs[outer_r][1]))

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
                              small_seg_size=50000, min_seg_size=10000, verbose=False):
    """
    Extract and resegment AA graph data for each BFB candidate region, returning
    data in the format expected by BFBArchitect's reconstruct_bfb().

    Processing pipeline (per region):
    1. Merge the inner foldback hairpin region (between the two endpoints of each
       foldback SV) into its adjacent amplified segment: -- foldbacks merge the
       hairpin right (into the segment starting at bp2); ++ foldbacks merge left
       (into the segment ending at bp1).
    2. Remove segments < small_seg_size bp whose CN is more than 1 copy lower than
       both neighbors (small deletion artifacts), respecting cut points.
    3. Absorb any remaining segments < min_seg_size bp unconditionally into a
       neighbor, respecting cut points.  These are too small to contribute useful
       CN signal and would only add noise to the ILP.
    4. Merge adjacent segments within cn_tol CN units, respecting foldback outer
       cut points as hard merge boundaries.
    5. Apply cut points: pass all segments through unchanged; split any segment
       that a cut point falls strictly inside (rare after prior merging).
    6. Build the cn/lf/rf vectors required by reconstruct_bfb().

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
        Segments < this many bp with CN >1 below both neighbors are treated as
        deletion artifacts and absorbed (default 50000).
    min_seg_size : int
        Segments < this many bp are unconditionally absorbed into a neighbor
        after deletion removal (default 10000).
    verbose : bool
        Print per-step segment transforms and final CN/LF/RF vectors to stdout.

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

    n_foldbacks = sum(1 for sv, _, _ in svs if sv.is_foldback(max_distance=fb_dist_cut))
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        return [None] * len(regions)

    raw_svs = svs
    svs = find_tst_foldbacks(raw_svs, chrom_segs, fb_dist=fb_dist_cut, verbose=verbose)

    n_foldbacks = sum(1 for sv, _, _ in svs if sv.is_foldback(max_distance=fb_dist_cut))
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        svs = raw_svs

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
            if verbose:
                print(f"\n[subsect_graph_for_region] {chrom}:{region_start}-{region_end}")
                print("  No segments found in region.")
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

        if verbose:
            print(f"\n[subsect_graph_for_region] {chrom}:{region_start}-{region_end}")
            print(f"  Raw segments ({len(raw_segs)}):")
            for i, (s, e, cn, cov, rc) in enumerate(raw_segs):
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}  cov={cov:.1f}  rc={rc}")
            print(f"  Foldback SVs ({len(foldback_svs)}):")
            for sv in foldback_svs:
                sv_cn = sv_cn_map.get(sv, 0.0)
                direction = sv.strand1 * 2
                cut = sv.bp2 if sv.strand1 == '+' else sv.bp1 - 1
                print(f"    {direction}  {sv.chrom1}:{sv.bp1}-{sv.bp2}  CN={sv_cn:.2f}  cut_point={cut}")
            if not foldback_svs:
                print("    (none)")
            print(f"  Cut points:   {sorted(cut_points)}")
            print(f"  FB endpoints: {sorted(fb_endpoints)}")

        segs = list(raw_segs)

        # ── Step 1: merge inner foldback loop segments ─────────────────────────
        # For -- foldbacks the hairpin region [bp1, bp2) is merged RIGHT into the
        # segment starting at bp2 (the amplified side).  For ++ foldbacks the
        # hairpin region (bp1, bp2] is merged LEFT into the segment ending at bp1.
        # This keeps bp1 (--) / bp2 (++) as segment boundaries so lf/rf assignment
        # still works after re-segmentation.
        for sv in foldback_svs:
            if sv.strand1 == '-':
                loop_idx = [i for i, seg in enumerate(segs) if seg[0] >= sv.bp1 and seg[1] < sv.bp2]
                tgt_idx  = [i for i, seg in enumerate(segs) if seg[0] == sv.bp2]
            else:
                loop_idx = [i for i, seg in enumerate(segs) if seg[0] > sv.bp1 and seg[1] <= sv.bp2]
                tgt_idx  = [i for i, seg in enumerate(segs) if seg[1] == sv.bp1]
            if not loop_idx or not tgt_idx:
                continue
            all_idx = sorted(loop_idx + tgt_idx)
            to_merge = [segs[i] for i in all_idx]
            total_len = sum(seg[1] - seg[0] + 1 for seg in to_merge)
            merged = (
                to_merge[0][0],
                to_merge[-1][1],
                sum(seg[2] * (seg[1] - seg[0] + 1) for seg in to_merge) / total_len,
                sum(seg[3] * (seg[1] - seg[0] + 1) for seg in to_merge) / total_len,
                sum(seg[4] for seg in to_merge),
            )
            for i in reversed(all_idx):
                segs.pop(i)
            segs.insert(all_idx[0], merged)

        if verbose:
            print(f"  Step 1 – merge inner foldback loops: {len(raw_segs)} → {len(segs)} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

        n_after_step1 = len(segs)

        # ── Step 2: remove small low-CN (deleted) segments ─────────────────────
        # A segment is absorbed if: size < small_seg_size AND its CN is more than
        # 1 copy lower than BOTH neighbors (deletion artifact).
        # Cut points are hard boundaries for absorption.
        changed = True
        while changed:
            changed = False
            new_segs = []
            i = 0
            while i < len(segs):
                s, e, seg_cn, seg_cov, seg_rc = segs[i]
                size = e - s + 1
                left_cn  = new_segs[-1][2] if new_segs else None
                right_cn = segs[i + 1][2] if i + 1 < len(segs) else None
                is_small_del = (
                    size < small_seg_size
                    and left_cn is not None and left_cn > seg_cn + 1
                    and right_cn is not None and right_cn > seg_cn + 1
                )
                if is_small_del:
                    can_left  = bool(new_segs) and new_segs[-1][1] not in cut_points
                    can_right = (i + 1 < len(segs)) and (e not in cut_points)
                    if can_left and can_right:
                        if abs(left_cn - seg_cn) <= abs(right_cn - seg_cn):
                            can_right = False
                        else:
                            can_left = False
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

        n_after_step2 = len(segs)

        if verbose:
            print(f"  Step 2 – remove small deletions (<{small_seg_size} bp, CN >1 below both neighbors): "
                  f"{n_after_step1} → {n_after_step2} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

        # ── Step 3: absorb remaining tiny segments unconditionally ─────────────
        # Segments < min_seg_size bp add noise to the CN vector regardless of CN.
        # They are absorbed into the nearest neighbor within the same cut-point
        # block.  If a foldback cut point will re-split the merged segment in
        # Step 5, the desired boundary is still preserved.
        changed = True
        while changed:
            changed = False
            new_segs = []
            i = 0
            while i < len(segs):
                s, e, seg_cn, seg_cov, seg_rc = segs[i]
                size = e - s + 1
                if size < min_seg_size:
                    can_left  = bool(new_segs) and new_segs[-1][1] not in cut_points
                    can_right = (i + 1 < len(segs)) and (e not in cut_points)
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

        n_after_step3 = len(segs)

        if verbose:
            print(f"  Step 3 – absorb tiny (<{min_seg_size} bp): {n_after_step2} → {n_after_step3} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

        # ── Step 4: merge adjacent CN-similar segments ─────────────────────────
        # Cut points are hard merge boundaries.
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

        if verbose:
            print(f"  Step 4 – merge CN-similar (tol={cn_tol}): {n_after_step3} → {len(segs)} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

        # ── Step 5: apply foldback cut points ─────────────────────────────────
        # Pass segments through unchanged.  If a cut point falls strictly inside
        # a segment (rare after prior merging), split it; otherwise the cut point
        # is already a segment boundary and the segment passes through as-is.
        new_segments = []
        for seg_start, seg_end, seg_cn, seg_cov, seg_rc in segs:
            inner_cuts = sorted(cp for cp in cut_points if seg_start <= cp < seg_end)
            if not inner_cuts:
                new_segments.append((chrom, seg_start, seg_end, seg_cn, seg_cov, seg_rc))
            else:
                cur = seg_start
                for cp in inner_cuts:
                    new_segments.append((chrom, cur, cp, seg_cn, seg_cov, seg_rc))
                    cur = cp + 1
                new_segments.append((chrom, cur, seg_end, seg_cn, seg_cov, seg_rc))

        if verbose:
            print(f"  Step 5 – apply cut points: {len(segs)} input → {len(new_segments)} segments")
            for i, (ch, s, e, cn, cov, rc) in enumerate(new_segments):
                cut_flag = "  [cut here]" if e in cut_points else ""
                print(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}{cut_flag}")

        # ── Step 6: build cn/lf/rf vectors ─────────────────────────────────────
        cn_vals = [round(seg[3]) - 1 for seg in new_segments]
        l_bp_idx = {seg[1]: i for i, seg in enumerate(new_segments)}
        r_bp_idx = {seg[2]: i for i, seg in enumerate(new_segments)}
        lf = [0] * len(cn_vals)
        rf = [0] * len(cn_vals)

        for sv in foldback_svs:
            sv_cn = sv_cn_map.get(sv, 0.0)
            fb_count = _foldback_cn_count(sv_cn)
            if sv.strand1 == '-':    # left (--) foldback: segment starts at bp1
                if sv.bp1 in l_bp_idx:
                    lf[l_bp_idx[sv.bp1]] += fb_count
            else:                     # right (++) foldback: segment ends at bp2
                if sv.bp2 in r_bp_idx:
                    rf[r_bp_idx[sv.bp2]] += fb_count

        if verbose:
            print(f"  Step 6 – CN/LF/RF vectors ({len(new_segments)} segments):")
            print(f"    {'idx':>3}  {'start':>12}  {'end':>12}  {'CN_float':>9}  {'cn':>4}  {'lf':>4}  {'rf':>4}")
            for i, seg in enumerate(new_segments):
                cn_flag = "  *** cn<=0" if cn_vals[i] <= 0 else ""
                print(f"    {i:>3}  {seg[1]:>12}  {seg[2]:>12}  {seg[3]:>9.3f}  {cn_vals[i]:>4}  {lf[i]:>4}  {rf[i]:>4}{cn_flag}")
            unmatched = []
            for sv in foldback_svs:
                if sv.strand1 == '-' and sv.bp1 not in l_bp_idx:
                    unmatched.append(f"-- foldback bp1={sv.bp1} not a segment start → lf contribution lost")
                elif sv.strand1 == '+' and sv.bp2 not in r_bp_idx:
                    unmatched.append(f"++ foldback bp2={sv.bp2} not a segment end → rf contribution lost")
            if unmatched:
                print("  WARNING – unmatched foldbacks (CN not assigned to lf/rf):")
                for msg in unmatched:
                    print(f"    {msg}")

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

    n_foldbacks = sum(1 for sv, _, _ in svs_raw if sv.is_foldback())
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        return [], [], [], [], [], {}, ''

    raw_svs = svs_raw
    svs_raw = find_tst_foldbacks(raw_svs, chrom_segs)

    n_foldbacks = sum(1 for sv, _, _ in svs_raw if sv.is_foldback())
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        svs_raw = raw_svs

    svs_list = [sv for sv, _cn, _rc in svs_raw]
    sv_info  = {sv: (cn, rc) for sv, cn, rc in svs_raw}

    # Derive the amplified region from the sequence edges.
    # Primary chromosome = highest total CN-weighted length (most amplified content).
    chrom_ranges = {}
    chrom_cn_length = {}
    for chrom, segs in chrom_segs.items():
        starts = [s for s, *_ in segs]
        ends   = [e for _, e, *_ in segs]
        chrom_ranges[chrom] = (min(starts), max(ends))
        chrom_cn_length[chrom] = sum((e - s) * cn for s, e, cn, *_ in segs)
    primary_chrom = max(chrom_cn_length, key=chrom_cn_length.__getitem__)
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
        seg_size = seg_end - seg_start + 1
        total_length += seg_size
        weighted_cn_sum += seg_cn * seg_size
        total_cov += seg_cov * seg_size
        total_rc += seg_rc
        if seg_end in breakpoints:
            new_cn  = weighted_cn_sum / total_length
            new_cov = total_cov / total_length
            new_segments.append((chrom, start, seg_end, new_cn, new_cov, total_rc))
            start = seg_end + 1
            total_length, weighted_cn_sum = 0, 0
            total_cov, total_rc = 0.0, 0
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
            fb_count = _foldback_cn_count(sv_info[sv][0])
            if sv.strand1 == '-' and sv.bp1 in l_bp:
                lf[l_bp.index(sv.bp1)] += fb_count
            elif sv.strand1 == '+' and sv.bp2 in r_bp:
                rf[r_bp.index(sv.bp2)] += fb_count

    return new_segments, cn, lf, rf, svs_list, sv_info, primary_chrom


# ── TST chain report ─────────────────────────────────────────────────────────

def write_tst_report(graph_file, output_file, bfb_regions=None,
                     shard_max_bp=5000, max_hops=5, fb_dist=50000,
                     far_min=500000):
    """
    Write a human-readable TST-chain report for graph_file to output_file.

    Each confirmed TST event is described as a chain:
      local breakpoints → SV anchors → shard path → far ends
    with CNs, read counts, and whether each far end lies in a BFB candidate
    region.

    Parameters
    ----------
    graph_file  : str   Path to AA _graph.txt file.
    output_file : str   Path to write the report (plain text).
    bfb_regions : list of (chrom, start, end), optional
        Pre-computed BFB candidate regions.  If None, auto-detected via
        find_bfb_candidate_regions().
    """
    svs, chrom_segs = parse_graph_file(graph_file)

    if bfb_regions is None:
        bfb_regions = find_bfb_candidate_regions(graph_file)

    def in_bfb(chrom, pos):
        for rc, rs, re in bfb_regions:
            if rc == chrom and rs <= pos <= re:
                return True
        return False

    def seg_for(chrom, pos, direction):
        """Return the segment tuple containing this breakend, or None."""
        for s, e, cn, cov, rc in chrom_segs.get(chrom, []):
            ep = s if direction == '-' else e
            if ep == pos:
                return (s, e, cn)
        return None

    bp_to_seg, sv_index = _build_graph_lookups(svs, chrom_segs)
    existing_fb_cuts = set()
    real_fb_spans = defaultdict(list)
    for sv, _cn, _rc in svs:
        if sv.is_foldback():
            existing_fb_cuts.add((sv.chrom1, sv.bp2 if sv.strand1 == '+' else sv.bp1))
            real_fb_spans[sv.chrom1].append((sv.bp1, sv.bp2))

    far_by_chrom = defaultdict(list)
    for entry in svs:
        sv, cn, rc = entry
        for local, far in [
            ((sv.chrom1, sv.bp1, sv.strand1), (sv.chrom2, sv.bp2, sv.strand2)),
            ((sv.chrom2, sv.bp2, sv.strand2), (sv.chrom1, sv.bp1, sv.strand1)),
        ]:
            lc, lp, ld = local
            fc, fp, _ = far
            if fc != lc or abs(fp - lp) >= far_min:
                far_by_chrom[lc].append((lp, ld, far, sv, cn, rc))

    events = []
    seen_pairs = set()

    for chrom, breakends in far_by_chrom.items():
        for i, (pos_i, dir_i, far_i, sv_i, cn_i, rc_i) in enumerate(breakends):
            for pos_j, dir_j, far_j, sv_j, cn_j, rc_j in breakends[i + 1:]:
                if sv_i is sv_j:                    continue
                if abs(pos_i - pos_j) > fb_dist:    continue
                if dir_i == dir_j:                  continue
                if far_i == far_j:                  continue
                pair_key = (min(id(sv_i), id(sv_j)), max(id(sv_i), id(sv_j)))
                if pair_key in seen_pairs:          continue
                seen_pairs.add(pair_key)

                fc_i, fp_i, fd_i = far_i
                fc_j, fp_j, fd_j = far_j

                # Skip if the far ends are the direct endpoints of a real foldback.
                if fc_i == fc_j:
                    p_lo, p_hi = min(fp_i, fp_j), max(fp_i, fp_j)
                    if any(abs(b1 - p_lo) <= 5 and abs(b2 - p_hi) <= 5
                           for b1, b2 in real_fb_spans.get(fc_i, [])):
                        continue

                path = _shard_path_reachable(far_i, far_j, bp_to_seg, sv_index,
                                             shard_max_bp, max_hops,
                                             exclude_svs={sv_i, sv_j})
                if path is None:
                    continue

                # Determine synthetic FBI coordinates (same as find_tst_foldbacks)
                pos_lo, dir_lo = (pos_i, dir_i) if pos_i <= pos_j else (pos_j, dir_j)
                pos_hi          = pos_j           if pos_i <= pos_j else pos_i
                if dir_lo == '+':
                    fbi_bp1, fbi_bp2 = pos_lo, pos_hi - 1
                else:
                    fbi_bp1, fbi_bp2 = pos_lo - 1, pos_hi
                if fbi_bp2 < fbi_bp1:
                    continue

                cn_tst = min(cn_i, cn_j)
                cn_ratio = min(cn_i, cn_j) / max(cn_i, cn_j) if max(cn_i, cn_j) > 0 else 0

                # Far-side FBI (if applicable)
                far_fbi = None
                far_fbi_reason = None
                if fc_i != fc_j or abs(fp_i - fp_j) > fb_dist:
                    far_fbi_reason = "far ends on different chroms or too far apart"
                else:
                    if fd_i == fd_j:
                        far_fbi = SV(fc_i, min(fp_i, fp_j), fd_i,
                                     fc_i, max(fp_i, fp_j), fd_i)
                    else:
                        fp_lo2 = min(fp_i, fp_j); fp_hi2 = max(fp_i, fp_j)
                        fd_lo2 = fd_i if fp_i <= fp_j else fd_j
                        fbp1 = fp_lo2 if fd_lo2 == '+' else fp_lo2 - 1
                        fbp2 = fp_hi2 - 1 if fd_lo2 == '+' else fp_hi2
                        if fbp2 >= fbp1:
                            far_fbi = SV(fc_i, fbp1, '+', fc_i, fbp2, '+')
                    if far_fbi is not None:
                        far_cut = (fc_i, far_fbi.bp2 if far_fbi.strand1 == '+'
                                   else far_fbi.bp1)
                        if far_cut in existing_fb_cuts:
                            far_fbi_reason = (f"direct foldback already covers "
                                              f"{fc_i} cut {far_cut[1]}")
                            far_fbi = None

                events.append({
                    'chrom': chrom,
                    'pos_lo': min(pos_i, pos_j), 'dir_lo': dir_lo,
                    'pos_hi': max(pos_i, pos_j),
                    'local_gap': abs(pos_i - pos_j),
                    'sv_lo': sv_i if pos_i <= pos_j else sv_j,
                    'sv_hi': sv_j if pos_i <= pos_j else sv_i,
                    'cn_lo': cn_i if pos_i <= pos_j else cn_j,
                    'rc_lo': rc_i if pos_i <= pos_j else rc_j,
                    'cn_hi': cn_j if pos_i <= pos_j else cn_i,
                    'rc_hi': rc_j if pos_i <= pos_j else rc_i,
                    'cn_ratio': cn_ratio,
                    'cn_tst': cn_tst,
                    'far_lo': far_i if pos_i <= pos_j else far_j,
                    'far_hi': far_j if pos_i <= pos_j else far_i,
                    'path': path,
                    'fbi_bp1': fbi_bp1, 'fbi_bp2': fbi_bp2,
                    'far_fbi': far_fbi,
                    'far_fbi_reason': far_fbi_reason,
                })

    import os
    with open(output_file, 'w') as out:
        out.write(f"TST-FOLDBACK REPORT\n")
        out.write(f"Graph:   {os.path.basename(graph_file)}\n")
        if bfb_regions:
            out.write(f"BFB regions detected: "
                      + ", ".join(f"{c}:{s}-{e}" for c, s, e in bfb_regions) + "\n")
        else:
            out.write("BFB regions detected: none\n")
        out.write(f"TST events found: {len(events)}\n")

        for idx, ev in enumerate(events, 1):
            out.write(f"\n{'='*60}\n")
            out.write(f"Event {idx}/{len(events)}  [local chromosome: {ev['chrom']}]\n")
            out.write(f"\n  Local breakpoints: {ev['chrom']}:{ev['pos_lo']}({ev['dir_lo']}) "
                      f"<-> {ev['chrom']}:{ev['pos_hi']}(+)   "
                      f"[{ev['local_gap']} bp apart]\n")

            # SV at low-position local end
            sv_lo = ev['sv_lo']
            far_lo_c, far_lo_p, far_lo_d = ev['far_lo']
            far_lo_seg = seg_for(far_lo_c, far_lo_p, far_lo_d)
            bfb_lo = in_bfb(far_lo_c, far_lo_p)
            out.write(f"\n  SV (lo-end)  CN={ev['cn_lo']:.2f}  rc={ev['rc_lo']}\n")
            out.write(f"    {sv_lo}\n")
            out.write(f"    Far end: {far_lo_c}:{far_lo_p}({far_lo_d})")
            if far_lo_seg:
                s, e, cn = far_lo_seg
                out.write(f"  →  seg {far_lo_c}:{s}-{e}  [{e-s} bp, CN={cn:.2f}]")
            out.write(f"  {'◀ in BFB region' if bfb_lo else ''}\n")

            # SV at high-position local end
            sv_hi = ev['sv_hi']
            far_hi_c, far_hi_p, far_hi_d = ev['far_hi']
            far_hi_seg = seg_for(far_hi_c, far_hi_p, far_hi_d)
            bfb_hi = in_bfb(far_hi_c, far_hi_p)
            out.write(f"\n  SV (hi-end)  CN={ev['cn_hi']:.2f}  rc={ev['rc_hi']}\n")
            out.write(f"    {sv_hi}\n")
            out.write(f"    Far end: {far_hi_c}:{far_hi_p}({far_hi_d})")
            if far_hi_seg:
                s, e, cn = far_hi_seg
                out.write(f"  →  seg {far_hi_c}:{s}-{e}  [{e-s} bp, CN={cn:.2f}]")
            out.write(f"  {'◀ in BFB region' if bfb_hi else ''}\n")

            out.write(f"\n  CN_ratio: {ev['cn_ratio']:.2f}   TST CN: {ev['cn_tst']:.2f}\n")

            # Shard path
            path = ev['path']
            if path:
                total = sum(e - s for _, s, e in path)
                out.write(f"\n  Shard path ({len(path)} hop{'s' if len(path)!=1 else ''}, "
                          f"{total} bp total):\n")
                for pc, ps, pe in path:
                    # Find CN of this shard from chrom_segs
                    shard_cn = next((cn for s, e, cn, *_ in chrom_segs.get(pc, [])
                                     if s == ps and e == pe), None)
                    cn_str = f"  CN={shard_cn:.2f}" if shard_cn is not None else ""
                    bfb_str = "  ◀ in BFB region" if in_bfb(pc, (ps + pe) // 2) else ""
                    out.write(f"    {pc}:{ps}-{pe}  [{pe-ps} bp]{cn_str}{bfb_str}\n")
            else:
                out.write(f"\n  Shard path: concordant (0 shard hops)\n")

            # Injected FBIs
            out.write(f"\n  Injected FBI (local):  "
                      f"{ev['chrom']}:{ev['fbi_bp1']}(+) <-> "
                      f"{ev['chrom']}:{ev['fbi_bp2']}(+)  CN={ev['cn_tst']:.2f}\n")
            if ev['far_fbi'] is not None:
                ff = ev['far_fbi']
                out.write(f"  Injected FBI (far):    "
                          f"{ff.chrom1}:{ff.bp1}({ff.strand1}) <-> "
                          f"{ff.chrom2}:{ff.bp2}({ff.strand2})  CN={ev['cn_tst']:.2f}\n")
            else:
                out.write(f"  Injected FBI (far):    none"
                          f"  [{ev['far_fbi_reason']}]\n")

        out.write(f"\n{'='*60}\n")
