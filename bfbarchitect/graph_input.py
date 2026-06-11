import logging
import re
from collections import defaultdict, deque

try:
    from bfbarchitect.datatypes import SV, CHR_CENTRO
except ImportError:
    from datatypes import SV, CHR_CENTRO

MAX_FOLDBACKS_FOR_GRAPH_RECON = 50
MAX_GRAPH_RECON_SEGMENTS = 100
MAX_WHOLE_GRAPH_PRIMARY_SEGMENTS = MAX_GRAPH_RECON_SEGMENTS
LOCAL_FOLDBACK_RESCUE_CN_RATIO = 1.2
MIN_RAW_FOLDBACK_FLANK_CN_CHANGE = 0.2
MIN_TST_FOLDBACK_FLANK_CN_CHANGE = 1.0
FOLDBACK_FLANK_CN_WINDOW_BP = 100
TST_FOLDBACK_CN_METHOD = 'max'
LOGGER = logging.getLogger('BFBArchitect')


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


def _cn_merge_allowed(cn1, cn2, abs_tol=1.0, cv_tol=0.05):
    """Allow larger absolute CN differences at high copy number."""
    mean_cn = (cn1 + cn2) / 2
    return abs(cn1 - cn2) <= max(abs_tol, cv_tol * mean_cn)


def _tst_foldback_cn(cn_i, cn_j):
    """Estimate synthetic TST foldback CN from its two supporting jump SVs."""
    if TST_FOLDBACK_CN_METHOD == 'mean':
        return (cn_i + cn_j) / 2
    if TST_FOLDBACK_CN_METHOD == 'max':
        return max(cn_i, cn_j)
    return min(cn_i, cn_j)


def _regions_by_chrom(local_regions, padding=0):
    if local_regions is None:
        return None
    region_by_chrom = defaultdict(list)
    for chrom, start, end in local_regions:
        region_by_chrom[chrom].append((start - padding, end + padding))
    return region_by_chrom


def _point_in_regions(region_by_chrom, chrom, pos):
    if region_by_chrom is None:
        return True
    return any(start <= pos <= end for start, end in region_by_chrom.get(chrom, []))


def _sv_in_regions(region_by_chrom, sv):
    if region_by_chrom is None:
        return True
    return (_point_in_regions(region_by_chrom, sv.chrom1, sv.bp1)
            or _point_in_regions(region_by_chrom, sv.chrom2, sv.bp2))


def _near_duplicate_direct_foldback(sv, real_fb_spans,
                                    shared_endpoint_slop=5,
                                    other_endpoint_slop=5000,
                                    same_polarity_span_gap=1000):
    """
    True if sv nearly duplicates a direct foldback.

    The shared-endpoint rule catches synthetic foldbacks that land on the same
    cut point as a direct foldback. The span-gap rule catches same-polarity TST
    foldbacks that sit immediately adjacent to a direct foldback without sharing
    an endpoint, which otherwise double-counts one local foldback structure.
    """
    for bp1, bp2, strand in real_fb_spans.get(sv.chrom1, []):
        same_left = abs(sv.bp1 - bp1) <= shared_endpoint_slop
        same_right = abs(sv.bp2 - bp2) <= shared_endpoint_slop
        if same_left and abs(sv.bp2 - bp2) <= other_endpoint_slop:
            return True
        if same_right and abs(sv.bp1 - bp1) <= other_endpoint_slop:
            return True
        if sv.strand1 == strand:
            span_gap = max(bp1 - sv.bp2, sv.bp1 - bp2, 0)
            if span_gap <= same_polarity_span_gap:
                return True
    return False


def _filter_processable_regions(regions, chrom_segs, max_segments):
    if max_segments is None:
        return list(regions)
    processable = []
    for chrom, region_start, region_end in regions:
        raw_count = sum(
            1 for s, e, *_ in chrom_segs.get(chrom, [])
            if s < region_end and e > region_start
        )
        if raw_count and raw_count <= max_segments:
            processable.append((chrom, region_start, region_end))
    return processable


def _flank_cn_candidates(segs, anchor, window_bp):
    cns = []
    for start, end, cn, _cov, _rc in segs:
        if (start <= anchor <= end
                or abs(start - anchor) <= window_bp
                or abs(end - anchor) <= window_bp):
            cns.append(cn)
    return cns


def _nearest_large_flank_cn(segs, pos, side, min_size):
    # Outermost positions in the region — a segment touching these is the
    # terminal one in the amplicon sub-region and gets the size floor waived.
    if side == 'left':
        outer = min((start for start, end, *_ in segs if end <= pos), default=None)
    else:
        outer = max((end for start, end, *_ in segs if start >= pos), default=None)

    best = None
    best_dist = None
    for start, end, cn, *_ in segs:
        is_terminal = (side == 'left' and start == outer) or \
                      (side == 'right' and end == outer)
        if end - start + 1 < min_size and not is_terminal:
            continue
        if side == 'left':
            if end > pos:
                continue
            dist = pos - end
        else:
            if start < pos:
                continue
            dist = start - pos
        if best_dist is None or dist < best_dist:
            best = cn
            best_dist = dist
    return best


def _tst_recipient_large_flank_cn_change(chrom_segs, chrom, pos_i, pos_j,
                                         min_size=5000):
    pos_lo, pos_hi = sorted((pos_i, pos_j))
    segs = sorted(chrom_segs.get(chrom, []), key=lambda seg: seg[0])
    left_cn = _nearest_large_flank_cn(segs, pos_lo, 'left', min_size)
    right_cn = _nearest_large_flank_cn(segs, pos_hi, 'right', min_size)
    if left_cn is None or right_cn is None:
        return None
    return abs(left_cn - right_cn)


def _foldback_flank_cn_change(sv, chrom_segs,
                              window_bp=FOLDBACK_FLANK_CN_WINDOW_BP):
    """
    CN difference between the left and right flanks near a foldback.

    The graph often places tiny TST-derived shards immediately inside or next
    to the synthetic foldback span.  Use a small window around each flank anchor
    so those local CN states contribute to the boundary-change test.

    Returns None if either flank cannot be found.  Callers should treat unknown
    flanks conservatively and keep the foldback.
    """
    if sv.chrom1 != sv.chrom2 or sv.strand1 != sv.strand2:
        return None

    segs = chrom_segs.get(sv.chrom1, [])
    if sv.strand1 == '+':
        left_anchor = sv.bp1
        right_anchor = sv.bp2 + 1
    else:
        left_anchor = sv.bp1 - 1
        right_anchor = sv.bp2

    left_cns = _flank_cn_candidates(segs, left_anchor, window_bp)
    right_cns = _flank_cn_candidates(segs, right_anchor, window_bp)

    if not left_cns or not right_cns:
        return None
    return max(abs(left_cn - right_cn)
               for left_cn in left_cns
               for right_cn in right_cns)


def _foldback_passes_flank_cn_filter(sv, chrom_segs, min_change):
    change = _foldback_flank_cn_change(sv, chrom_segs)
    return change is None or change >= min_change


def _foldback_passes_recon_cn_filter(sv, chrom_segs):
    tst_change = getattr(sv, 'TST_local_cn_change', None)
    if tst_change is not None:
        return tst_change >= MIN_TST_FOLDBACK_FLANK_CN_CHANGE
    min_change = (MIN_TST_FOLDBACK_FLANK_CN_CHANGE
                  if getattr(sv, 'TST', False)
                  else MIN_RAW_FOLDBACK_FLANK_CN_CHANGE)
    return _foldback_passes_flank_cn_filter(sv, chrom_segs, min_change)


def _native_foldbacks_on_chrom(svs, chrom_segs, chrom, fb_dist_cut=50000):
    return _filter_illegal_nested_foldbacks([
        sv for sv, _cn, _rc in svs
        if (not getattr(sv, 'TST', False)
            and sv.chrom1 == chrom
            and sv.chrom2 == chrom
            and sv.is_foldback(max_distance=fb_dist_cut)
            and _foldback_passes_recon_cn_filter(sv, chrom_segs))
    ])


def _filter_illegal_nested_foldbacks(foldbacks):
    if len(foldbacks) < 2:
        return foldbacks

    illegal = set()
    for i, sv_i in enumerate(foldbacks):
        for sv_j in foldbacks[i + 1:]:
            if sv_i.chrom1 != sv_j.chrom1:
                continue
            if sv_i.strand1 == sv_j.strand1:
                continue
            i_contains_j = sv_i.bp1 <= sv_j.bp1 and sv_j.bp2 <= sv_i.bp2
            j_contains_i = sv_j.bp1 <= sv_i.bp1 and sv_i.bp2 <= sv_j.bp2
            if i_contains_j or j_contains_i:
                illegal.add(sv_i)
                illegal.add(sv_j)

    if not illegal:
        return foldbacks
    return [sv for sv in foldbacks if sv not in illegal]


def _legal_foldbacks_by_chrom(svs, chrom_segs, fb_dist_cut=50000,
                              include_tst=True):
    raw_by_chrom = defaultdict(list)
    for sv, _cn, _rc in svs:
        if (sv.is_foldback(max_distance=fb_dist_cut)
                and (include_tst or not getattr(sv, 'TST', False))
                and _foldback_passes_recon_cn_filter(sv, chrom_segs)):
            raw_by_chrom[sv.chrom1].append(sv)
    return {
        chrom: _filter_illegal_nested_foldbacks(foldbacks)
        for chrom, foldbacks in raw_by_chrom.items()
    }


def _segment_index_by_endpoint(chrom_segs):
    endpoint_to_seg = {}
    for chrom, segs in chrom_segs.items():
        for idx, (start, end, cn, cov, rc) in enumerate(segs):
            endpoint_to_seg[(chrom, start, '-')] = (idx, (start, end, cn, cov, rc))
            endpoint_to_seg[(chrom, end, '+')] = (idx, (start, end, cn, cov, rc))
    return endpoint_to_seg


def _sv_other_endpoint(sv, endpoint):
    chrom, pos, strand = endpoint
    if (sv.chrom1, sv.bp1, sv.strand1) == (chrom, pos, strand):
        return (sv.chrom2, sv.bp2, sv.strand2)
    if (sv.chrom2, sv.bp2, sv.strand2) == (chrom, pos, strand):
        return (sv.chrom1, sv.bp1, sv.strand1)
    return None


def _foldback_vector_slot(sv, l_bp_idx, r_bp_idx):
    if sv.strand1 == '-' and sv.bp1 in l_bp_idx:
        return 'lf', l_bp_idx[sv.bp1]
    if sv.strand1 == '+' and sv.bp2 in r_bp_idx:
        return 'rf', r_bp_idx[sv.bp2]
    return None


def _find_local_foldback_rescue_anchors(svs, sv_cn_map, chrom_segs,
                                        primary_chrom, region=None,
                                        verbose=False,
                                        label='local_foldback_rescue'):
    """
    Find displaced foldbacks before simplification hides their target boundary.

    Pattern:
      candidate boundary with CN drop -> non-foldback SV -> one landed sequence
      block -> foldback at the far end of that block, with a CN drop around it.
    """
    endpoint_to_seg = _segment_index_by_endpoint(chrom_segs)

    sv_by_endpoint = defaultdict(list)
    for sv, cn, rc in svs:
        endpoints = [
            (sv.chrom1, sv.bp1, sv.strand1),
            (sv.chrom2, sv.bp2, sv.strand2),
        ]
        for endpoint in endpoints:
            sv_by_endpoint[endpoint].append((sv, cn, rc))
    legal_foldbacks = {
        sv
        for foldbacks in _legal_foldbacks_by_chrom(svs, chrom_segs).values()
        for sv in foldbacks
    }
    foldbacks_by_endpoint = defaultdict(list)
    for sv, cn, rc in svs:
        if sv not in legal_foldbacks:
            continue
        for endpoint in [
            (sv.chrom1, sv.bp1, sv.strand1),
            (sv.chrom2, sv.bp2, sv.strand2),
        ]:
            foldbacks_by_endpoint[endpoint].append((sv, cn, rc))

    raw_segs = chrom_segs.get(primary_chrom, [])
    anchors = []
    moved_foldbacks = set()

    for raw_idx, (start, end, seg_cn, _cov, _rc) in enumerate(raw_segs):
        boundary_specs = []
        if raw_idx + 1 < len(raw_segs):
            next_cn = raw_segs[raw_idx + 1][2]
            if seg_cn >= LOCAL_FOLDBACK_RESCUE_CN_RATIO * next_cn:
                boundary_specs.append(('right', (primary_chrom, end, '+'), next_cn))
        if raw_idx > 0:
            prev_cn = raw_segs[raw_idx - 1][2]
            if seg_cn >= LOCAL_FOLDBACK_RESCUE_CN_RATIO * prev_cn:
                boundary_specs.append(('left', (primary_chrom, start, '-'), prev_cn))

        for side, boundary_endpoint, outside_cn in boundary_specs:
            if region is not None:
                _chrom, region_start, region_end = region
                if not (region_start <= boundary_endpoint[1] <= region_end):
                    continue
            if seg_cn < LOCAL_FOLDBACK_RESCUE_CN_RATIO * outside_cn:
                continue

            for sv, _sv_cn, _sv_rc in sv_by_endpoint.get(boundary_endpoint, []):
                if sv.is_foldback():
                    continue

                landed = _sv_other_endpoint(sv, boundary_endpoint)
                if landed is None:
                    continue
                landed_info = endpoint_to_seg.get(landed)
                if landed_info is None:
                    continue

                landed_idx, landed_seg = landed_info
                l_start, l_end, landed_cn, _l_cov, _l_rc = landed_seg
                if landed[2] == '+':
                    exit_endpoint = (landed[0], l_start, '-')
                    neighbor_idx = landed_idx - 1
                else:
                    exit_endpoint = (landed[0], l_end, '+')
                    neighbor_idx = landed_idx + 1

                landed_chrom_segs = chrom_segs.get(landed[0], [])
                if not (0 <= neighbor_idx < len(landed_chrom_segs)):
                    continue
                neighbor_cn = landed_chrom_segs[neighbor_idx][2]
                if landed_cn < LOCAL_FOLDBACK_RESCUE_CN_RATIO * neighbor_cn:
                    continue

                for fb_sv, fb_cn, _fb_rc in foldbacks_by_endpoint.get(exit_endpoint, []):
                    if fb_sv in moved_foldbacks:
                        continue
                    if getattr(fb_sv, 'TST', False):
                        continue
                    sv_cn = sv_cn_map.get(fb_sv, fb_cn)
                    if isinstance(sv_cn, tuple):
                        sv_cn = sv_cn[0]
                    fb_count = _foldback_cn_count(sv_cn)
                    moved_foldbacks.add(fb_sv)
                    anchor = {
                        'boundary_side': side,
                        'boundary_endpoint': boundary_endpoint,
                        'boundary_sv': sv,
                        'landed_endpoint': landed,
                        'landed_segment': (landed[0], l_start, l_end, landed_cn),
                        'foldback': fb_sv,
                        'count': fb_count,
                        'boundary_cn': seg_cn,
                        'boundary_outside_cn': outside_cn,
                        'landed_cn': landed_cn,
                        'landed_neighbor_cn': neighbor_cn,
                    }
                    anchors.append(anchor)

                    if verbose:
                        LOGGER.info(f"  [{label}] anchor count {fb_count}: "
                              f"{boundary_endpoint[0]}:{boundary_endpoint[1]}({boundary_endpoint[2]})")
                        LOGGER.info(f"       boundary {side}: raw CN={seg_cn:.3f} "
                              f"-> outside CN={outside_cn:.3f}")
                        LOGGER.info(f"       via SV: {sv}")
                        LOGGER.info(f"       landed block: {landed[0]}:{l_start}-{l_end}  "
                              f"CN={landed_cn:.3f} -> neighbor CN={neighbor_cn:.3f}")
                        LOGGER.info(f"       displaced foldback: {fb_sv}")

    return _dedupe_reciprocal_local_rescue_anchors(
        anchors, chrom_segs, verbose=verbose, label=label
    )


def _local_rescue_core_mass(anchor, chrom_segs, baseline_cn=2.0,
                            max_gap_bp=100000):
    chrom, pos, _strand = anchor['boundary_endpoint']
    segs = chrom_segs.get(chrom, [])
    if not segs:
        return 0.0, 0

    if anchor['boundary_side'] == 'right':
        start_idx = next((i for i, seg in enumerate(segs) if seg[1] == pos), None)
        step = -1
    else:
        start_idx = next((i for i, seg in enumerate(segs) if seg[0] == pos), None)
        step = 1
    if start_idx is None:
        return 0.0, 0

    mass = 0.0
    count = 0
    prev = None
    idx = start_idx
    while 0 <= idx < len(segs):
        seg = segs[idx]
        start, end, cn, *_ = seg
        if cn <= baseline_cn:
            break
        if prev is not None:
            if step == -1:
                gap = prev[0] - end - 1
            else:
                gap = start - prev[1] - 1
            if gap > max_gap_bp:
                break
        mass += (end - start + 1) * (cn - baseline_cn)
        count += 1
        prev = seg
        idx += step

    return mass, count


def _local_rescue_anchor_rank(anchor, chrom_segs):
    core_mass, core_segments = _local_rescue_core_mass(anchor, chrom_segs)
    outside_cn = max(anchor['boundary_outside_cn'], 1e-6)
    landed_neighbor_cn = max(anchor['landed_neighbor_cn'], 1e-6)
    return (
        core_mass,
        core_segments,
        anchor['boundary_cn'] / outside_cn,
        anchor['landed_cn'] / landed_neighbor_cn,
        anchor['count'],
    )


def _dedupe_reciprocal_local_rescue_anchors(anchors, chrom_segs=None,
                                            verbose=False,
                                            label='local_foldback_rescue'):
    """Keep only the stronger direction for reciprocal fragment-rescue anchors."""
    if len(anchors) < 2:
        return anchors
    if chrom_segs is None:
        chrom_segs = {}

    discard = set()
    for i, anchor_i in enumerate(anchors):
        if i in discard:
            continue
        for j in range(i + 1, len(anchors)):
            if j in discard:
                continue
            anchor_j = anchors[j]
            same_bridge_sv = anchor_i['boundary_sv'] is anchor_j['boundary_sv']
            reciprocal_endpoints = (
                anchor_i['boundary_endpoint'] == anchor_j['landed_endpoint']
                and anchor_j['boundary_endpoint'] == anchor_i['landed_endpoint']
            )
            if not (same_bridge_sv and reciprocal_endpoints):
                continue

            rank_i = _local_rescue_anchor_rank(anchor_i, chrom_segs)
            rank_j = _local_rescue_anchor_rank(anchor_j, chrom_segs)
            drop_idx, keep_idx = (j, i) if rank_i >= rank_j else (i, j)
            discard.add(drop_idx)

            if verbose:
                dropped = anchors[drop_idx]
                kept = anchors[keep_idx]
                LOGGER.info(
                    f"  [{label}] dropped reciprocal anchor at "
                    f"{dropped['boundary_endpoint'][0]}:"
                    f"{dropped['boundary_endpoint'][1]}"
                    f"({dropped['boundary_endpoint'][2]}); kept "
                    f"{kept['boundary_endpoint'][0]}:"
                    f"{kept['boundary_endpoint'][1]}"
                    f"({kept['boundary_endpoint'][2]}) "
                    f"rank_kept={_local_rescue_anchor_rank(kept, chrom_segs)} "
                    f"rank_dropped={_local_rescue_anchor_rank(dropped, chrom_segs)}"
                )

    if not discard:
        return anchors
    return [anchor for idx, anchor in enumerate(anchors) if idx not in discard]


def _local_rescue_cut_point(anchor):
    pos = anchor['boundary_endpoint'][1]
    if anchor['boundary_side'] == 'right':
        return pos
    return pos - 1


def _apply_local_foldback_rescue_anchors(new_segments, lf, rf, anchors,
                                         verbose=False,
                                         label='local_foldback_rescue'):
    """
    Move displaced foldback counts to pre-detected anchor boundaries.

    Anchors are found on raw graph segments before simplification.  Their cut
    points are added before vector construction, so the target boundary should
    now be a segment start/end.
    """
    if not new_segments or not anchors:
        return []

    l_bp_idx = {seg[1]: i for i, seg in enumerate(new_segments)}
    r_bp_idx = {seg[2]: i for i, seg in enumerate(new_segments)}
    corrections = []
    moved_foldbacks = set()

    for anchor in anchors:
        fb_sv = anchor['foldback']
        if fb_sv in moved_foldbacks:
            continue
        fb_count = anchor['count']
        side = anchor['boundary_side']
        boundary_endpoint = anchor['boundary_endpoint']

        if side == 'right' and boundary_endpoint[1] in r_bp_idx:
            target_idx = r_bp_idx[boundary_endpoint[1]]
            new_slot = ('rf', target_idx)
        elif side == 'left' and boundary_endpoint[1] in l_bp_idx:
            target_idx = l_bp_idx[boundary_endpoint[1]]
            new_slot = ('lf', target_idx)
        else:
            continue

        old_slot = _foldback_vector_slot(fb_sv, l_bp_idx, r_bp_idx)
        if old_slot is not None:
            old_side, old_idx = old_slot
            old_vec = lf if old_side == 'lf' else rf
            old_vec[old_idx] = max(0, old_vec[old_idx] - fb_count)

        if new_slot[0] == 'rf':
            rf[target_idx] += fb_count
        else:
            lf[target_idx] += fb_count

        fb_sv.local_resolved = True
        moved_foldbacks.add(fb_sv)
        correction = dict(anchor)
        correction['old_slot'] = old_slot
        correction['new_slot'] = new_slot
        correction['target_cn'] = new_segments[target_idx][3]
        corrections.append(correction)

        if verbose:
            old_desc = 'unassigned' if old_slot is None else f"{old_slot[0]}[{old_slot[1]}]"
            LOGGER.info(f"  [{label}] moved foldback count {fb_count}: {old_desc} -> "
                  f"{new_slot[0]}[{new_slot[1]}]")
            LOGGER.info(f"       boundary {side}: {boundary_endpoint[0]}:{boundary_endpoint[1]}"
                  f"({boundary_endpoint[2]})  raw CN={anchor['boundary_cn']:.3f}; "
                  f"target CN={correction['target_cn']:.3f} "
                  f"-> outside CN={anchor['boundary_outside_cn']:.3f}")
            LOGGER.info(f"       via SV: {anchor['boundary_sv']}")
            landed_chrom, l_start, l_end, landed_cn = anchor['landed_segment']
            LOGGER.info(f"       landed block: {landed_chrom}:{l_start}-{l_end}  "
                  f"CN={landed_cn:.3f} -> neighbor CN={anchor['landed_neighbor_cn']:.3f}")
            LOGGER.info(f"       displaced foldback: {fb_sv}")

    return corrections


def _trim_local_rescue_outside_flank(new_segments, anchors, verbose=False):
    """Remove a terminal low-CN outside flank created by a rescue-anchor cut."""
    if not new_segments or not anchors:
        return None

    for anchor in anchors:
        side = anchor['boundary_side']
        chrom, pos, _strand = anchor['boundary_endpoint']
        if side == 'right':
            outside_start = pos + 1
            if (new_segments[-1][0] == chrom
                    and new_segments[-1][1] == outside_start):
                removed = new_segments.pop(-1)
                if verbose:
                    LOGGER.info(f"  [local_foldback_rescue] trimmed outside flank: "
                          f"{removed[0]}:{removed[1]}-{removed[2]}")
                return 'right'
        else:
            outside_end = pos - 1
            if (new_segments[0][0] == chrom
                    and new_segments[0][2] == outside_end):
                removed = new_segments.pop(0)
                if verbose:
                    LOGGER.info(f"  [local_foldback_rescue] trimmed outside flank: "
                          f"{removed[0]}:{removed[1]}-{removed[2]}")
                return 'left'
    return None


def _remove_local_rescue_landed_segments(segs, anchors, chrom, verbose=False,
                                         label='local_foldback_rescue'):
    """Consume graph segments reached by an SV-to-foldback-shard rescue path."""
    if not segs or not anchors:
        return segs, []

    consumed = set()
    for anchor in anchors:
        landed_chrom, l_start, l_end, _landed_cn = anchor['landed_segment']
        if landed_chrom == chrom:
            consumed.add((l_start, l_end))
        fb_sv = anchor['foldback']
        if fb_sv.chrom1 != chrom:
            continue
        for s, e, *_ in segs:
            if fb_sv.strand1 == '-':
                if s >= fb_sv.bp1 and e < fb_sv.bp2:
                    consumed.add((s, e))
            else:
                if s > fb_sv.bp1 and e <= fb_sv.bp2:
                    consumed.add((s, e))

    if not consumed:
        return segs, []

    kept = []
    removed = []
    for seg in segs:
        key = (seg[0], seg[1])
        if key in consumed:
            removed.append(seg)
        else:
            kept.append(seg)

    if verbose and removed:
        for s, e, cn, _cov, _rc in removed:
            LOGGER.info(f"  [{label}] consumed landed foldback shard: "
                  f"{chrom}:{s}-{e} CN={cn:.3f}")

    return kept, removed


def _crosses_blocked_gap(left_seg, right_seg, blocked_intervals):
    if not blocked_intervals:
        return False

    gap_start = left_seg[1] + 1
    gap_end = right_seg[0] - 1
    if gap_start > gap_end:
        return False

    return any(not (gap_end < block_start or gap_start > block_end)
               for block_start, block_end in blocked_intervals)


def _contract_hard_deletion_vectors(new_segments, cn, lf, rf,
                                    max_size=250000, max_cn_float=3.0,
                                    max_cn_int=1, min_flank_cn=10,
                                    max_flank_cn_diff=1.0):
    """
    Remove low-CN intervals between comparable high-CN flanks from BFB vectors.

    This differs from ordinary small-segment smoothing: the deleted interval is
    excluded from the merged CN estimate rather than averaged into a neighbor,
    and the corresponding LF/RF vector entries are contracted with the flanks.
    """
    if len(new_segments) < 3:
        return new_segments, cn, lf, rf, []

    out_segments = []
    out_cn = []
    out_lf = []
    out_rf = []
    contractions = []
    i = 0
    while i < len(new_segments):
        if i + 2 < len(new_segments):
            left = new_segments[i]
            gap = new_segments[i + 1]
            right = new_segments[i + 2]
            lchrom, ls, le, left_cn, left_cov, left_rc = left
            gchrom, gs, ge, gap_cn, _gap_cov, _gap_rc = gap
            rchrom, rs, re, right_cn, right_cov, right_rc = right
            gap_size = ge - gs + 1
            is_hard_deletion = (
                lchrom == gchrom == rchrom
                and le + 1 == gs
                and ge + 1 == rs
                and gap_size < max_size
                and gap_cn < max_cn_float
                and cn[i + 1] <= max_cn_int
                and cn[i] >= min_flank_cn
                and cn[i + 2] >= min_flank_cn
                and left_cn >= min_flank_cn
                and right_cn >= min_flank_cn
                and abs(left_cn - right_cn) <= max_flank_cn_diff
            )
            if is_hard_deletion:
                left_len = le - ls + 1
                right_len = re - rs + 1
                flank_len = left_len + right_len
                merged = (
                    lchrom,
                    ls,
                    re,
                    (left_cn * left_len + right_cn * right_len) / flank_len,
                    (left_cov * left_len + right_cov * right_len) / flank_len,
                    left_rc + right_rc,
                )
                merged_cn = round(merged[3]) - 1
                merged_lf = lf[i] + lf[i + 1] + lf[i + 2]
                merged_rf = rf[i] + rf[i + 1] + rf[i + 2]
                contractions.append({
                    'left': left,
                    'gap': gap,
                    'right': right,
                    'merged': merged,
                    'merged_cn': merged_cn,
                    'merged_lf': merged_lf,
                    'merged_rf': merged_rf,
                })
                out_segments.append(merged)
                out_cn.append(merged_cn)
                out_lf.append(merged_lf)
                out_rf.append(merged_rf)
                i += 3
                continue
        out_segments.append(new_segments[i])
        out_cn.append(cn[i])
        out_lf.append(lf[i])
        out_rf.append(rf[i])
        i += 1

    return out_segments, out_cn, out_lf, out_rf, contractions


def _segment_span_overlaps(seg, span):
    """Return true if a reconstruction segment overlaps a genomic span."""
    chrom, start, end = span
    return seg[0] == chrom and seg[1] <= end and seg[2] >= start


def _contract_local_rescue_consumed_gap_vectors(new_segments, cn, lf, rf,
                                                consumed_segments,
                                                cn_tol=1.0):
    if len(new_segments) < 2 or not consumed_segments:
        return new_segments, cn, lf, rf, []

    blocked = [(seg[0], seg[1]) for seg in consumed_segments]
    out_segments = []
    out_cn = []
    out_lf = []
    out_rf = []
    contractions = []
    i = 0
    while i < len(new_segments):
        if i + 1 < len(new_segments):
            left = new_segments[i]
            right = new_segments[i + 1]
            if (_crosses_blocked_gap(
                    (left[1], left[2]), (right[1], right[2]), blocked)
                    and lf[i] == 0 and rf[i] == 0
                    and lf[i + 1] == 0 and rf[i + 1] == 0
                    and _cn_merge_allowed(left[3], right[3], abs_tol=cn_tol)):
                l_len = left[2] - left[1] + 1
                r_len = right[2] - right[1] + 1
                merged = (
                    left[0],
                    left[1],
                    right[2],
                    (left[3] * l_len + right[3] * r_len) / (l_len + r_len),
                    (left[4] * l_len + right[4] * r_len) / (l_len + r_len),
                    left[5] + right[5],
                )
                out_segments.append(merged)
                out_cn.append(round(merged[3]) - 1)
                out_lf.append(0)
                out_rf.append(0)
                contractions.append({
                    'left': left,
                    'right': right,
                    'merged': merged,
                    'merged_cn': out_cn[-1],
                })
                i += 2
                continue

        out_segments.append(new_segments[i])
        out_cn.append(cn[i])
        out_lf.append(lf[i])
        out_rf.append(rf[i])
        i += 1

    return out_segments, out_cn, out_lf, out_rf, contractions


def _trim_foldback_free_terminal_flanks(new_segments, cn, lf, rf,
                                        max_cn=2):
    trimmed = []
    while new_segments and cn and lf and rf:
        if cn[0] <= max_cn and lf[0] == 0 and rf[0] == 0:
            trimmed.append(('left', new_segments.pop(0), cn.pop(0), lf.pop(0), rf.pop(0)))
            continue
        break
    while new_segments and cn and lf and rf:
        if cn[-1] <= max_cn and lf[-1] == 0 and rf[-1] == 0:
            trimmed.append(('right', new_segments.pop(-1), cn.pop(-1), lf.pop(-1), rf.pop(-1)))
            continue
        break
    return trimmed


def _contract_deletion_bridge_plateaus(new_segments, cn, lf, rf, svs,
                                       min_cn=5):
    """
    Merge adjacent equal-CN slices anchored by graph deletions.

    Foldback/TST cut points intentionally split reconstruction segments so LF/RF
    counts can be assigned.  We therefore allow foldback-free spacer slices to
    be absorbed into a neighboring equal-CN state, but we do not merge two
    foldback-bearing slices into one LF/RF vector entry.
    """
    if len(new_segments) < 2 or not svs:
        return new_segments, cn, lf, rf, []

    deletion_spans = []
    for entry in svs:
        sv = entry[0] if isinstance(entry, tuple) else entry
        span = _deletion_span(sv)
        if span is not None:
            deletion_spans.append(span)
    if not deletion_spans:
        return new_segments, cn, lf, rf, []

    out_segments = []
    out_cn = []
    out_lf = []
    out_rf = []
    contractions = []
    i = 0
    while i < len(new_segments):
        group = [new_segments[i]]
        group_lf = lf[i]
        group_rf = rf[i]
        group_cn = cn[i]
        group_has_foldback = (lf[i] != 0 or rf[i] != 0)
        deletion_touched = any(
            _segment_span_overlaps(new_segments[i], span)
            for span in deletion_spans
        )
        j = i + 1
        while j < len(new_segments):
            prev = group[-1]
            cur = new_segments[j]
            cur_touched = any(
                _segment_span_overlaps(cur, span)
                for span in deletion_spans
            )
            same_high_state = (
                prev[0] == cur[0]
                and group_cn == cn[j]
                and group_cn >= min_cn
            )
            cur_has_foldback = (lf[j] != 0 or rf[j] != 0)
            if (not same_high_state
                    or not (deletion_touched or cur_touched)
                    or (group_has_foldback and cur_has_foldback)):
                break
            group.append(cur)
            group_lf += lf[j]
            group_rf += rf[j]
            group_has_foldback = group_has_foldback or cur_has_foldback
            deletion_touched = deletion_touched or cur_touched
            j += 1

        if len(group) == 1:
            out_segments.append(new_segments[i])
            out_cn.append(cn[i])
            out_lf.append(lf[i])
            out_rf.append(rf[i])
            i += 1
            continue

        total_len = sum(seg[2] - seg[1] + 1 for seg in group)
        merged = (
            group[0][0],
            group[0][1],
            group[-1][2],
            sum(seg[3] * (seg[2] - seg[1] + 1) for seg in group) / total_len,
            sum(seg[4] * (seg[2] - seg[1] + 1) for seg in group) / total_len,
            sum(seg[5] for seg in group),
        )
        contractions.append({
            'segments': group,
            'merged': merged,
            'merged_cn': group_cn,
            'merged_lf': group_lf,
            'merged_rf': group_rf,
        })
        out_segments.append(merged)
        out_cn.append(group_cn)
        out_lf.append(group_lf)
        out_rf.append(group_rf)
        i = j

    return out_segments, out_cn, out_lf, out_rf, contractions


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


def _deletion_span(sv):
    """Return the reference interval skipped by an intrachromosomal deletion SV."""
    if sv.chrom1 != sv.chrom2 or sv.get_SV_type() != 'DEL':
        return None
    start = min(sv.bp1, sv.bp2) + 1
    end = max(sv.bp1, sv.bp2) - 1
    if start > end:
        return None
    return sv.chrom1, start, end


def _deletion_crosses_foldback(span, svs):
    """Return true if a DEL span contains a same-chromosome foldback."""
    chrom, del_start, del_end = span
    for sv, _cn, _rc in svs:
        if sv.chrom1 != chrom or sv.chrom2 != chrom:
            continue
        if not sv.is_foldback():
            continue
        if del_start <= sv.bp1 <= del_end and del_start <= sv.bp2 <= del_end:
            return True
    return False


def _segment_overlaps(segs, start, end):
    overlaps = []
    for seg in segs:
        seg_start, seg_end = seg[0], seg[1]
        overlap_start = max(seg_start, start)
        overlap_end = min(seg_end, end)
        if overlap_start <= overlap_end:
            overlaps.append((seg, overlap_end - overlap_start + 1))
    return overlaps


def _deletion_span_is_cn_drop(span, chrom_segs, min_drop=1.0,
                              min_interior_bp=1000):
    """
    Require the skipped interval to look like a hard deletion in sequence CN.

    The DEL edge correction is intended for intervals whose sequence traversal CN
    is consistently below the neighboring amplified sequence.  If the skipped
    interval remains high CN, or only drops after crossing other rearrangement
    structure, adding the DEL edge CN creates artificial spikes.
    """
    chrom, del_start, del_end = span
    segs = chrom_segs.get(chrom, [])
    if not segs:
        return False

    left_flank = None
    right_flank = None
    for seg in segs:
        seg_start, seg_end = seg[0], seg[1]
        if seg_end < del_start:
            left_flank = seg
        elif seg_start > del_end and right_flank is None:
            right_flank = seg
            break

    if left_flank is None or right_flank is None:
        return False

    flank_floor = min(left_flank[2], right_flank[2])
    inner = _segment_overlaps(segs, del_start, del_end)
    if not inner:
        return False

    for seg, overlap_len in inner:
        if overlap_len < min_interior_bp:
            continue
        if seg[2] > flank_floor - min_drop:
            return False

    return True


def apply_deletion_cn_correction(chrom_segs, svs, verbose=False):
    """
    Add deletion-edge CN back to sequence segments skipped by graph DEL edges.

    AA graph sequence-edge CN measures copies that traverse the reference
    sequence.  A deletion breakpoint edge represents additional copies that skip
    the interval between its breakpoints.  For reconstruction, those skipped
    copies should still count toward the segment copy number.
    """
    corrected = {
        chrom: [tuple(seg) for seg in segs]
        for chrom, segs in chrom_segs.items()
    }
    corrections = []

    for sv, del_cn, rc in svs:
        span = _deletion_span(sv)
        if span is None or del_cn <= 0:
            continue
        chrom, del_start, del_end = span
        if chrom not in corrected:
            continue
        if _deletion_crosses_foldback(span, svs):
            if verbose:
                LOGGER.info(f"  Deletion CN correction skipped: {sv} crosses a foldback")
            continue
        if not _deletion_span_is_cn_drop(span, corrected):
            if verbose:
                LOGGER.info(f"  Deletion CN correction skipped: {sv} span is not a CN drop")
            continue

        new_segs = []
        for seg_start, seg_end, seg_cn, seg_cov, seg_rc in corrected[chrom]:
            overlap_start = max(seg_start, del_start)
            overlap_end = min(seg_end, del_end)
            if overlap_start <= overlap_end:
                seg_len = seg_end - seg_start + 1
                overlap_len = overlap_end - overlap_start + 1
                delta_cn = del_cn * overlap_len / seg_len
                seg_cn += delta_cn
                corrections.append({
                    'sv': sv,
                    'sv_cn': del_cn,
                    'read_count': rc,
                    'segment': (chrom, seg_start, seg_end),
                    'overlap': (overlap_start, overlap_end),
                    'delta_cn': delta_cn,
                    'corrected_cn': seg_cn,
                })
            new_segs.append((seg_start, seg_end, seg_cn, seg_cov, seg_rc))
        corrected[chrom] = new_segs

    if verbose and corrections:
        LOGGER.info(f"  Deletion CN correction: {len(corrections)} segment adjustment(s)")
        for c in corrections:
            chrom, seg_start, seg_end = c['segment']
            overlap_start, overlap_end = c['overlap']
            LOGGER.info(f"    {c['sv']}  DEL_CN={c['sv_cn']:.3f}  "
                  f"{chrom}:{seg_start}-{seg_end} overlap "
                  f"{overlap_start}-{overlap_end}  "
                  f"+{c['delta_cn']:.3f} -> CN={c['corrected_cn']:.3f}")

    return corrected, corrections


# ── TST-jump foldback detection ──────────────────────────────────────────────

def _build_graph_lookups(svs, chrom_segs):
    """Build (chrom, pos, dir) -> segment endpoint lookup."""
    bp_to_seg = {}
    for chrom, segs in chrom_segs.items():
        for start, end, *_ in segs:
            bp_to_seg[(chrom, start, '-')] = (start, end)
            bp_to_seg[(chrom, end,   '+')] = (start, end)

    return bp_to_seg


def _shard_sequence_path_reachable(start, target, bp_to_seg, shard_max_bp,
                                   max_hops):
    """Directional shard path using only sequence and concordant graph edges."""
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

        seg = bp_to_seg.get(curr)
        if seg is not None:
            s, e = seg
            seg_id = (chrom, s, e)
            if e - s <= shard_max_bp and seg_id not in visited:
                other = (chrom, e, '+') if d == '-' else (chrom, s, '-')
                queue.append((other, path + [seg_id], visited | {seg_id}))

        adj = _concordant_adjacent(curr)
        if adj == target:
            queue.append((adj, path, visited))
        else:
            adj_seg = bp_to_seg.get(adj)
            if adj_seg is not None and adj_seg[1] - adj_seg[0] <= shard_max_bp:
                queue.append((adj, path, visited))

    return None


def _concordant_adjacent(endpoint):
    chrom, pos, d = endpoint
    return (chrom, pos + 1, '-') if d == '+' else (chrom, pos - 1, '+')


def _sv_endpoints(sv):
    return [
        (sv.chrom1, sv.bp1, sv.strand1),
        (sv.chrom2, sv.bp2, sv.strand2),
    ]


def _other_sv_endpoint(sv, endpoint):
    endpoints = _sv_endpoints(sv)
    if endpoint == endpoints[0]:
        return endpoints[1]
    if endpoint == endpoints[1]:
        return endpoints[0]
    return None


def _recipient_foldback_sv(ep_a, ep_b, fb_dist):
    chrom_a, pos_a, dir_a = ep_a
    chrom_b, pos_b, dir_b = ep_b
    if chrom_a != chrom_b or abs(pos_a - pos_b) > fb_dist:
        return None

    if pos_a <= pos_b:
        lo, hi = ep_a, ep_b
    else:
        lo, hi = ep_b, ep_a
    chrom, pos_lo, dir_lo = lo
    _chrom, pos_hi, dir_hi = hi

    if dir_lo == dir_hi:
        return SV(chrom, pos_lo, dir_lo, chrom, pos_hi, dir_hi)
    if dir_lo == '+' and dir_hi == '-':
        if pos_hi - 1 < pos_lo:
            return None
        return SV(chrom, pos_lo, '+', chrom, pos_hi - 1, '+')
    if dir_lo == '-' and dir_hi == '+':
        if pos_hi < pos_lo - 1:
            return None
        return SV(chrom, pos_lo - 1, '+', chrom, pos_hi, '+')
    return None


def _endpoint_is_far_jump(local_ep, far_ep, far_min):
    return local_ep[0] != far_ep[0] or abs(local_ep[1] - far_ep[1]) >= far_min


def _path_total_bp(path):
    return sum(e - s for _chrom, s, e in path)


def _directional_tst_pair_candidates(entry_i, entry_j, bp_to_seg, chrom_segs,
                                     fb_dist, far_min, shard_max_bp, max_hops,
                                     region_by_chrom=None):
    sv_i, cn_i, rc_i = entry_i
    sv_j, cn_j, rc_j = entry_j
    if sv_i is sv_j:
        return []

    candidates = []
    for rec_i in _sv_endpoints(sv_i):
        shard_i = _other_sv_endpoint(sv_i, rec_i)
        if shard_i is None or not _endpoint_is_far_jump(rec_i, shard_i, far_min):
            continue
        for rec_j in _sv_endpoints(sv_j):
            shard_j = _other_sv_endpoint(sv_j, rec_j)
            if shard_j is None or shard_i == shard_j:
                continue
            if not _endpoint_is_far_jump(rec_j, shard_j, far_min):
                continue

            synth = _recipient_foldback_sv(rec_i, rec_j, fb_dist)
            if synth is None:
                continue
            if region_by_chrom is not None and not _sv_in_regions(region_by_chrom, synth):
                continue

            path = _shard_sequence_path_reachable(
                shard_i, shard_j, bp_to_seg, shard_max_bp, max_hops
            )
            if path is None:
                continue

            cn_change = _tst_recipient_large_flank_cn_change(
                chrom_segs, synth.chrom1, rec_i[1], rec_j[1]
            )
            synth.TST = True
            synth.TST_local_cn_change = cn_change

            candidates.append({
                'synth': synth,
                'recipient_i': rec_i,
                'recipient_j': rec_j,
                'shard_i': shard_i,
                'shard_j': shard_j,
                'path': path,
                'cn_tst': _tst_foldback_cn(cn_i, cn_j),
                'cn_i': cn_i,
                'cn_j': cn_j,
                'rc_i': rc_i,
                'rc_j': rc_j,
                'sv_i': sv_i,
                'sv_j': sv_j,
                'cn_change': cn_change,
                'touches_region': _sv_in_regions(region_by_chrom, synth),
            })

    return candidates


def _select_directional_tst_candidate(candidates):
    if not candidates:
        return None

    def score(candidate):
        cn_change = candidate['cn_change']
        return (
            1 if candidate['touches_region'] else 0,
            round(cn_change, 3) if cn_change is not None else -1.0,
            -_path_total_bp(candidate['path']),
        )

    return max(candidates, key=score)


def find_tst_foldbacks(svs, chrom_segs, shard_max_bp=5000, max_hops=5,
                       fb_dist=50000, far_min=100000, verbose=False,
                       local_regions=None):
    """
    Detect graph-only TST foldbacks and inject synthetic FBI SVs.

    For each pair of far-jumping discordant SVs, both possible recipient sides
    are tested directionally.  A valid candidate starts at one recipient
    endpoint, crosses one candidate SV, traverses only shard-sized sequence
    segments plus concordant adjacency, crosses the other candidate SV, and
    lands at the expected recipient endpoint.  The recipient endpoints are then
    collapsed into a synthetic same-chromosome foldback.

    Synthetic TST foldbacks are admitted only when they touch the active region,
    are not near-duplicates of a direct foldback, and have a non-neutral CN step
    between the nearest substantial recipient flanks.
    """
    bp_to_seg = _build_graph_lookups(svs, chrom_segs)
    region_by_chrom = _regions_by_chrom(local_regions, padding=fb_dist)

    # Pre-compute cut-point positions of direct (non-TST) foldbacks so we can
    # skip far-side injections that would duplicate an existing foldback's cut.
    # ++ foldbacks cut at bp2; -- foldbacks cut at bp1.
    # Also build a per-chromosome list of real foldback (bp1, bp2) spans for the
    # TST duplicate check below.
    real_fb_spans = defaultdict(list)
    for sv, _cn, _rc in svs:
        if (sv.is_foldback()
                and _foldback_passes_recon_cn_filter(sv, chrom_segs)):
            real_fb_spans[sv.chrom1].append((sv.bp1, sv.bp2, sv.strand1))

    synthetic = []
    for i, entry_i in enumerate(svs):
        for entry_j in svs[i + 1:]:
            pair_candidates = _directional_tst_pair_candidates(
                entry_i, entry_j, bp_to_seg, chrom_segs, fb_dist, far_min,
                shard_max_bp, max_hops, region_by_chrom=region_by_chrom
            )
            candidate = _select_directional_tst_candidate(pair_candidates)
            if candidate is None:
                continue

            synth = candidate['synth']
            if _near_duplicate_direct_foldback(synth, real_fb_spans):
                if verbose:
                    LOGGER.info(f"[TST] Candidate {candidate['recipient_i']} -> "
                          f"{candidate['recipient_j']} skipped: near-duplicate "
                          f"direct foldback {synth.chrom1}:{synth.bp1}-{synth.bp2}")
                continue
            if not candidate['touches_region']:
                if verbose:
                    LOGGER.info(f"[TST] Candidate {candidate['recipient_i']} -> "
                          f"{candidate['recipient_j']} skipped: outside local region")
                continue
            if (candidate['cn_change'] is None
                    or candidate['cn_change'] < MIN_TST_FOLDBACK_FLANK_CN_CHANGE):
                if verbose:
                    cn_msg = ("unknown" if candidate['cn_change'] is None
                              else f"{candidate['cn_change']:.3f}")
                    LOGGER.info(f"[TST] Candidate {candidate['recipient_i']} -> "
                          f"{candidate['recipient_j']} skipped: local-side CN "
                          f"change {cn_msg} < "
                          f"{MIN_TST_FOLDBACK_FLANK_CN_CHANGE:.3f}")
                continue

            if verbose:
                total_shard_bp = _path_total_bp(candidate['path'])
                hops = ' -> '.join(
                    f"{c}:{s}-{e} [{e-s} bp]"
                    for c, s, e in candidate['path']
                )
                alt_count = len(pair_candidates) - 1
                alt_msg = f"; {alt_count} alternate side(s)" if alt_count else ""
                LOGGER.info(f"[TST] Injecting directional FBI {synth} "
                      f"CN={candidate['cn_tst']:.2f} "
                      f"local_side_CN_delta={candidate['cn_change']:.3f} "
                      f"path={total_shard_bp} bp{alt_msg}")
                if hops:
                    LOGGER.info(f"       Shard path: {hops}")

            synthetic.append((synth, candidate['cn_tst'], 0))

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
                                min_cn_step=0.5, merge_padding=150000,
                                deletion=True, verbose=False):
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
    deletion : bool
        If true, add deletion-edge CN back to skipped sequence segments before
        detecting BFB candidate regions.

    Returns
    -------
    list of (chrom, start, end) tuples in half-open [start, end) coordinates.
    """
    svs, chrom_segs = parse_graph_file(graph_file)
    if deletion:
        chrom_segs, _deletion_corrections = apply_deletion_cn_correction(
            chrom_segs, svs, verbose=verbose
        )

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
    legal_fb_by_chrom = _legal_foldbacks_by_chrom(
        svs, chrom_segs, fb_dist_cut=fb_dist_cut
    )
    fb_by_chrom = {
        chrom: [
            (sv.bp1, sv.bp2, sv.strand1)
            for sv in foldbacks
        ]
        for chrom, foldbacks in legal_fb_by_chrom.items()
    }

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
                              cn_cv_tol=0.05, small_seg_size=50000,
                              min_seg_size=10000, verbose=False,
                              max_segments=MAX_GRAPH_RECON_SEGMENTS,
                              report_skips=False, deletion=True,
                              disable_tst=False):
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
    4. Merge adjacent segments within max(cn_tol, cn_cv_tol * mean CN),
       respecting foldback outer cut points as hard merge boundaries.
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
        Absolute CN floor for merging adjacent segments (default 1.0).
    cn_cv_tol : float
        Relative CN tolerance for high-copy segments. Adjacent segments merge
        if their absolute CN difference is <= max(cn_tol, cn_cv_tol * mean CN).
    small_seg_size : int
        Segments < this many bp with CN >1 below both neighbors are treated as
        deletion artifacts and absorbed (default 50000).
    min_seg_size : int
        Segments < this many bp are unconditionally absorbed into a neighbor
        after deletion removal (default 10000).
    verbose : bool
        Print per-step segment transforms and final CN/LF/RF vectors to stdout.
    max_segments : int or None
        Maximum graph/reconstruction segments allowed per region. None disables
        the cutoff.
    deletion : bool
        If true, add deletion-edge CN back to skipped sequence segments before
        building reconstruction segments and CN/LF/RF vectors.

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
    if deletion:
        chrom_segs, _deletion_corrections = apply_deletion_cn_correction(
            chrom_segs, svs, verbose=verbose
        )

    n_foldbacks = sum(
        len(foldbacks)
        for foldbacks in _legal_foldbacks_by_chrom(
            svs, chrom_segs, fb_dist_cut=fb_dist_cut
        ).values()
    )
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        return [None] * len(regions)

    region_raw_segs = []
    has_processable_region = False
    processable_regions = []
    for chrom, region_start, region_end in regions:
        raw_segs = [
            (s, e, cn, cov, rc)
            for s, e, cn, cov, rc in chrom_segs.get(chrom, [])
            if s < region_end and e > region_start
        ]
        too_many_segments = (
            max_segments is not None and len(raw_segs) > max_segments
        )
        region_raw_segs.append((raw_segs, too_many_segments))
        if raw_segs and not too_many_segments:
            has_processable_region = True
            processable_regions.append((chrom, region_start, region_end))

    if not has_processable_region:
        results = []
        for (chrom, region_start, region_end), (raw_segs, too_many_segments) in zip(regions, region_raw_segs):
            if not raw_segs:
                if verbose:
                    LOGGER.info(f"\n[subsect_graph_for_region] {chrom}:{region_start}-{region_end}")
                    LOGGER.info("  No segments found in region.")
            elif too_many_segments and (verbose or report_skips):
                LOGGER.info(f"[subsect_graph_for_region] Skipping "
                      f"{chrom}:{region_start}-{region_end}: "
                      f"{len(raw_segs)} graph segments exceeds "
                      f"--max-graph-segments={max_segments}.")
            results.append(None)
        return results

    raw_svs = svs
    if disable_tst:
        svs = raw_svs
    else:
        svs = find_tst_foldbacks(
            raw_svs, chrom_segs, fb_dist=fb_dist_cut, verbose=verbose,
            local_regions=processable_regions
        )

    n_foldbacks = sum(
        len(foldbacks)
        for foldbacks in _legal_foldbacks_by_chrom(
            svs, chrom_segs, fb_dist_cut=fb_dist_cut
        ).values()
    )
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        svs = raw_svs

    # Pre-group foldback SVs by chromosome to avoid a full scan per region.
    sv_cn_map = {sv: cn for sv, cn, _rc in svs}
    foldbacks_by_chrom = _legal_foldbacks_by_chrom(
        svs, chrom_segs, fb_dist_cut=fb_dist_cut
    )

    results = []
    for (chrom, region_start, region_end), (raw_segs, too_many_segments) in zip(regions, region_raw_segs):
        if not raw_segs:
            if verbose:
                LOGGER.info(f"\n[subsect_graph_for_region] {chrom}:{region_start}-{region_end}")
                LOGGER.info("  No segments found in region.")
            results.append(None)
            continue
        if too_many_segments:
            if verbose or report_skips:
                LOGGER.info(f"[subsect_graph_for_region] Skipping "
                      f"{chrom}:{region_start}-{region_end}: "
                      f"{len(raw_segs)} graph segments exceeds "
                      f"--max-graph-segments={max_segments}.")
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

        local_rescue_anchors = _find_local_foldback_rescue_anchors(
            svs, sv_cn_map, chrom_segs, chrom,
            region=(chrom, region_start, region_end),
            verbose=verbose
        )
        for anchor in local_rescue_anchors:
            cut_points.add(_local_rescue_cut_point(anchor))
            fb = anchor['foldback']
            fb_endpoints.update((fb.bp1, fb.bp2))
        local_rescued_foldbacks = {
            anchor['foldback'] for anchor in local_rescue_anchors
        }

        if verbose:
            LOGGER.info(f"\n[subsect_graph_for_region] {chrom}:{region_start}-{region_end}")
            LOGGER.info(f"  Raw segments ({len(raw_segs)}):")
            for i, (s, e, cn, cov, rc) in enumerate(raw_segs):
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}  cov={cov:.1f}  rc={rc}")
            LOGGER.info(f"  Foldback SVs ({len(foldback_svs)}):")
            for sv in foldback_svs:
                sv_cn = sv_cn_map.get(sv, 0.0)
                direction = sv.strand1 * 2
                cut = sv.bp2 if sv.strand1 == '+' else sv.bp1 - 1
                LOGGER.info(f"    {direction}  {sv.chrom1}:{sv.bp1}-{sv.bp2}  CN={sv_cn:.2f}  cut_point={cut}")
            if not foldback_svs:
                LOGGER.info("    (none)")
            LOGGER.info(f"  Cut points:   {sorted(cut_points)}")
            LOGGER.info(f"  FB endpoints: {sorted(fb_endpoints)}")

        segs, local_rescue_consumed_segments = _remove_local_rescue_landed_segments(
            list(raw_segs), local_rescue_anchors, chrom, verbose=verbose
        )
        local_rescue_blocked_gaps = {
            (seg[0], seg[1]) for seg in local_rescue_consumed_segments
        }

        if verbose and local_rescue_consumed_segments:
            LOGGER.info(f"  Local foldback rescue consumed "
                  f"{len(local_rescue_consumed_segments)} landed shard segment(s)")

        # ── Step 1: merge inner foldback loop segments ─────────────────────────
        # For -- foldbacks the hairpin region [bp1, bp2) is merged RIGHT into the
        # segment starting at bp2 (the amplified side).  For ++ foldbacks the
        # hairpin region (bp1, bp2] is merged LEFT into the segment ending at bp1.
        # This keeps bp1 (--) / bp2 (++) as segment boundaries so lf/rf assignment
        # still works after re-segmentation.
        for sv in foldback_svs:
            if sv in local_rescued_foldbacks:
                continue
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
            LOGGER.info(f"  Step 1 – merge inner foldback loops: {len(raw_segs)} → {len(segs)} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

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
                    can_left = (
                        bool(new_segs)
                        and new_segs[-1][1] not in cut_points
                        and not _crosses_blocked_gap(
                            new_segs[-1], segs[i], local_rescue_blocked_gaps
                        )
                    )
                    can_right = (
                        (i + 1 < len(segs))
                        and (e not in cut_points)
                        and not _crosses_blocked_gap(
                            segs[i], segs[i + 1], local_rescue_blocked_gaps
                        )
                    )
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
            LOGGER.info(f"  Step 2 – remove small deletions (<{small_seg_size} bp, CN >1 below both neighbors): "
                  f"{n_after_step1} → {n_after_step2} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

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
                    can_left = (
                        bool(new_segs)
                        and new_segs[-1][1] not in cut_points
                        and not _crosses_blocked_gap(
                            new_segs[-1], segs[i], local_rescue_blocked_gaps
                        )
                    )
                    can_right = (
                        (i + 1 < len(segs))
                        and (e not in cut_points)
                        and not _crosses_blocked_gap(
                            segs[i], segs[i + 1], local_rescue_blocked_gaps
                        )
                    )
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
            LOGGER.info(f"  Step 3 – absorb tiny (<{min_seg_size} bp): {n_after_step2} → {n_after_step3} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

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
                    if (e1 not in cut_points
                            and not _crosses_blocked_gap(
                                segs[i], segs[i + 1],
                                local_rescue_blocked_gaps
                            )
                            and _cn_merge_allowed(cn1, cn2, abs_tol=cn_tol,
                                                  cv_tol=cn_cv_tol)):
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
            LOGGER.info(f"  Step 4 – merge CN-similar (abs_tol={cn_tol}, cv_tol={cn_cv_tol}): "
                  f"{n_after_step3} → {len(segs)} segments")
            for i, (s, e, cn, cov, rc) in enumerate(segs):
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}")

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
            LOGGER.info(f"  Step 5 – apply cut points: {len(segs)} input → {len(new_segments)} segments")
            for i, (ch, s, e, cn, cov, rc) in enumerate(new_segments):
                cut_flag = "  [cut here]" if e in cut_points else ""
                LOGGER.info(f"    {i:>3}: {s:>12}-{e:<12}  size={e-s:>9}  CN={cn:.3f}{cut_flag}")

        if max_segments is not None and len(new_segments) > max_segments:
            if verbose or report_skips:
                LOGGER.info(f"[subsect_graph_for_region] Skipping "
                      f"{chrom}:{region_start}-{region_end}: "
                      f"{len(new_segments)} reconstruction segments exceeds "
                      f"--max-graph-segments={max_segments}.")
            results.append(None)
            continue

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

        _apply_local_foldback_rescue_anchors(
            new_segments, lf, rf, local_rescue_anchors, verbose=verbose
        )

        new_segments, cn_vals, lf, rf, local_rescue_gap_contractions = (
            _contract_local_rescue_consumed_gap_vectors(
                new_segments, cn_vals, lf, rf,
                local_rescue_consumed_segments, cn_tol=cn_tol
            )
        )
        new_segments, cn_vals, lf, rf, hard_deletion_contractions = (
            _contract_hard_deletion_vectors(new_segments, cn_vals, lf, rf)
        )
        new_segments, cn_vals, lf, rf, deletion_bridge_contractions = (
            _contract_deletion_bridge_plateaus(
                new_segments, cn_vals, lf, rf, svs if deletion else []
            )
        )
        terminal_flank_trims = _trim_foldback_free_terminal_flanks(
            new_segments, cn_vals, lf, rf
        )

        if verbose:
            if local_rescue_gap_contractions:
                LOGGER.info(f"  Step 6a – contract consumed-shard flanks: "
                      f"{len(local_rescue_gap_contractions)} contraction(s)")
                for contraction in local_rescue_gap_contractions:
                    left = contraction['left']
                    right = contraction['right']
                    merged = contraction['merged']
                    LOGGER.info(f"    merged {left[1]}-{left[2]} and "
                          f"{right[1]}-{right[2]} across consumed shard -> "
                          f"{merged[1]}-{merged[2]} "
                          f"CN_float={merged[3]:.3f} "
                          f"cn={contraction['merged_cn']}")
            if hard_deletion_contractions:
                LOGGER.info(f"  Step 6b – contract hard deletion vectors: "
                      f"{len(hard_deletion_contractions)} contraction(s)")
                for contraction in hard_deletion_contractions:
                    _lch, ls, le, left_cn, *_ = contraction['left']
                    _gch, gs, ge, gap_cn, *_ = contraction['gap']
                    _rch, rs, re, right_cn, *_ = contraction['right']
                    _mch, ms, me, merged_cn_float, *_ = contraction['merged']
                    LOGGER.info(f"    contracted {gs}-{ge} CN={gap_cn:.3f} "
                          f"between {ls}-{le} CN={left_cn:.3f} and "
                          f"{rs}-{re} CN={right_cn:.3f} -> "
                          f"{ms}-{me} CN_float={merged_cn_float:.3f} "
                          f"cn={contraction['merged_cn']} "
                          f"lf={contraction['merged_lf']} "
                          f"rf={contraction['merged_rf']}")
            if deletion_bridge_contractions:
                LOGGER.info(f"  Step 6c – contract deletion-bridge plateaus: "
                      f"{len(deletion_bridge_contractions)} contraction(s)")
                for contraction in deletion_bridge_contractions:
                    first = contraction['segments'][0]
                    last = contraction['segments'][-1]
                    _mch, ms, me, merged_cn_float, *_ = contraction['merged']
                    LOGGER.info(f"    merged {len(contraction['segments'])} segment(s) "
                          f"{first[1]}-{last[2]} -> "
                          f"{ms}-{me} CN_float={merged_cn_float:.3f} "
                          f"cn={contraction['merged_cn']} "
                          f"lf={contraction['merged_lf']} "
                          f"rf={contraction['merged_rf']}")
            if terminal_flank_trims:
                LOGGER.info(f"  Step 6d – trim foldback-free terminal flanks: "
                      f"{len(terminal_flank_trims)} trim(s)")
                for side, seg, seg_cn, _seg_lf, _seg_rf in terminal_flank_trims:
                    LOGGER.info(f"    trimmed {side} {seg[0]}:{seg[1]}-{seg[2]} "
                          f"CN_float={seg[3]:.3f} cn={seg_cn}")
            LOGGER.info(f"  Step 6 – CN/LF/RF vectors ({len(new_segments)} segments):")
            LOGGER.info(f"    {'idx':>3}  {'start':>12}  {'end':>12}  {'CN_float':>9}  {'cn':>4}  {'lf':>4}  {'rf':>4}")
            for i, seg in enumerate(new_segments):
                cn_flag = "  *** cn<=0" if cn_vals[i] <= 0 else ""
                LOGGER.info(f"    {i:>3}  {seg[1]:>12}  {seg[2]:>12}  {seg[3]:>9.3f}  {cn_vals[i]:>4}  {lf[i]:>4}  {rf[i]:>4}{cn_flag}")
            unmatched = []
            for sv in foldback_svs:
                if sv.strand1 == '-' and sv.bp1 not in l_bp_idx:
                    unmatched.append(f"-- foldback bp1={sv.bp1} not a segment start → lf contribution lost")
                elif sv.strand1 == '+' and sv.bp2 not in r_bp_idx:
                    unmatched.append(f"++ foldback bp2={sv.bp2} not a segment end → rf contribution lost")
            if unmatched:
                LOGGER.info("  WARNING – unmatched foldbacks (CN not assigned to lf/rf):")
                for msg in unmatched:
                    LOGGER.info(f"    {msg}")

        results.append((new_segments, cn_vals, lf, rf, region_svs, sv_info))

    return results


# ── whole-graph fallback ──────────────────────────────────────────────────────

def whole_graph_as_region(graph_file, centromere_dict=None, verbose=False,
                          max_primary_segments=MAX_WHOLE_GRAPH_PRIMARY_SEGMENTS,
                          report_skips=False, deletion=True,
                          disable_tst=False,
                          require_native_foldback=True):
    """
    Treat all segments in the graph as a single BFB region (the --whole_graph
    fallback).

    Parameters
    ----------
    graph_file : str
        Path to an AA _graph.txt file.
    deletion : bool
        If true, add deletion-edge CN back to skipped sequence segments before
        deriving the whole-graph region and vectors.

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
    if deletion:
        chrom_segs, _deletion_corrections = apply_deletion_cn_correction(
            chrom_segs, svs_raw, verbose=verbose
        )

    # Derive the amplified region from the sequence edges.
    # Primary chromosome = highest total CN-weighted length (most amplified content).
    chrom_ranges = {}
    chrom_cn_length = {}
    for chrom, segs in chrom_segs.items():
        starts = [s for s, *_ in segs]
        ends   = [e for _, e, *_ in segs]
        chrom_ranges[chrom] = (min(starts), max(ends))
        chrom_cn_length[chrom] = sum((e - s) * cn for s, e, cn, *_ in segs)
    native_foldbacks_by_chrom = {
        chrom: _native_foldbacks_on_chrom(svs_raw, chrom_segs, chrom)
        for chrom in chrom_cn_length
    }
    candidate_chroms = [
        chrom for chrom, native_foldbacks in native_foldbacks_by_chrom.items()
        if native_foldbacks
    ]
    if require_native_foldback:
        if not candidate_chroms:
            if verbose or report_skips:
                LOGGER.info("[whole_graph_as_region] Skipping whole graph: "
                      "no chromosomes have native foldbacks passing "
                      "reconstruction filters.")
            return [], [], [], [], [], {}, ''
        primary_chrom = max(candidate_chroms, key=chrom_cn_length.__getitem__)
    else:
        primary_chrom = max(chrom_cn_length, key=chrom_cn_length.__getitem__)
    primary_segment_count = len(chrom_segs[primary_chrom])
    if (max_primary_segments is not None
            and primary_segment_count > max_primary_segments):
        if verbose or report_skips:
            LOGGER.info(f"[whole_graph_as_region] Skipping {primary_chrom}: "
                  f"{primary_segment_count} graph segments exceeds "
                  f"--max-graph-segments={max_primary_segments}.")
        return [], [], [], [], [], {}, ''

    n_foldbacks = sum(1 for sv, _, _ in svs_raw
                      if sv.is_foldback()
                      and _foldback_passes_recon_cn_filter(sv, chrom_segs))
    if n_foldbacks > MAX_FOLDBACKS_FOR_GRAPH_RECON:
        return [], [], [], [], [], {}, ''

    region_start, region_end = chrom_ranges[primary_chrom]
    region = (primary_chrom, region_start, region_end)

    region_data = subsect_graph_for_region(
        graph_file, [region], verbose=verbose, max_segments=max_primary_segments,
        report_skips=report_skips, deletion=deletion, disable_tst=disable_tst
    )
    if not region_data or region_data[0] is None:
        return [], [], [], [], [], {}, ''

    new_segments, cn, lf, rf, svs_list, sv_info = region_data[0]

    if centromere_dict is None:
        centromere_dict = CHR_CENTRO
    centro = centromere_dict.get(primary_chrom)
    if new_segments:
        arm_trim_side = 'left' if centro is None or new_segments[-1][2] < centro else 'right'
        if arm_trim_side == 'left' and cn[0] <= 1 and lf[0] == 0 and rf[0] == 0:
            new_segments.pop(0)
            cn.pop(0)
            lf.pop(0)
            rf.pop(0)
        elif arm_trim_side == 'right' and cn[-1] <= 1 and lf[-1] == 0 and rf[-1] == 0:
            new_segments.pop(-1)
            cn.pop(-1)
            lf.pop(-1)
            rf.pop(-1)

    return new_segments, cn, lf, rf, svs_list, sv_info, primary_chrom

# ── TST chain report ─────────────────────────────────────────────────────────

def write_tst_report(graph_file, output_file, bfb_regions=None,
                     shard_max_bp=5000, max_hops=5, fb_dist=50000,
                     far_min=100000,
                     max_segments=MAX_GRAPH_RECON_SEGMENTS):
    """
    Write a human-readable TST-chain report for graph_file to output_file.

    Each candidate is described as a directional graph walk:
      recipient endpoint -> candidate SV -> shard path -> candidate SV ->
      recipient endpoint.

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
    tst_regions = _filter_processable_regions(
        bfb_regions, chrom_segs, max_segments
    )
    region_by_chrom = _regions_by_chrom(tst_regions, padding=fb_dist)

    bp_to_seg = _build_graph_lookups(svs, chrom_segs)
    real_fb_spans = defaultdict(list)
    for sv, _cn, _rc in svs:
        if sv.is_foldback():
            real_fb_spans[sv.chrom1].append((sv.bp1, sv.bp2, sv.strand1))

    report_region_by_chrom = region_by_chrom if bfb_regions else None
    events = []
    for i, entry_i in enumerate(svs):
        for entry_j in svs[i + 1:]:
            candidates = _directional_tst_pair_candidates(
                entry_i, entry_j, bp_to_seg, chrom_segs, fb_dist, far_min,
                shard_max_bp, max_hops, region_by_chrom=report_region_by_chrom
            )
            selected = _select_directional_tst_candidate(candidates)
            if selected is None:
                continue

            reason = None
            synth = selected['synth']
            if _near_duplicate_direct_foldback(synth, real_fb_spans):
                reason = "near-duplicate of direct foldback"
            elif not selected['touches_region']:
                reason = "outside local region"
            elif (selected['cn_change'] is None
                  or selected['cn_change'] < MIN_TST_FOLDBACK_FLANK_CN_CHANGE):
                cn_msg = ("unknown" if selected['cn_change'] is None
                          else f"{selected['cn_change']:.3f}")
                reason = (f"local-side CN change {cn_msg} < "
                          f"{MIN_TST_FOLDBACK_FLANK_CN_CHANGE:.3f}")

            events.append({
                'selected': selected,
                'candidates': candidates,
                'reason': reason,
            })

    import os
    with open(output_file, 'w') as out:
        out.write("TST-FOLDBACK REPORT\n")
        out.write(f"Graph:   {os.path.basename(graph_file)}\n")
        if bfb_regions:
            out.write("BFB regions detected: "
                      + ", ".join(f"{c}:{s}-{e}" for c, s, e in bfb_regions)
                      + "\n")
        else:
            out.write("BFB regions detected: none\n")
        out.write(f"TST events found: {len(events)}\n")

        for idx, ev in enumerate(events, 1):
            selected = ev['selected']
            synth = selected['synth']
            path = selected['path']
            cn_msg = ("unknown" if selected['cn_change'] is None
                      else f"{selected['cn_change']:.3f}")
            out.write(f"\n{'='*60}\n")
            out.write(f"Event {idx}/{len(events)}  [directional graph traversal]\n")
            out.write(f"\n  Recipient endpoints: "
                      f"{selected['recipient_i'][0]}:{selected['recipient_i'][1]}"
                      f"({selected['recipient_i'][2]}) <-> "
                      f"{selected['recipient_j'][0]}:{selected['recipient_j'][1]}"
                      f"({selected['recipient_j'][2]})\n")
            out.write(f"  Synthetic FBI:      {synth.chrom1}:{synth.bp1}"
                      f"({synth.strand1}) <-> {synth.chrom2}:{synth.bp2}"
                      f"({synth.strand2})  CN={selected['cn_tst']:.2f}  "
                      f"local_CN_delta={cn_msg}\n")
            out.write(f"  SV1: {selected['sv_i']}  CN={selected['cn_i']:.2f}  "
                      f"rc={selected['rc_i']}\n")
            out.write(f"  SV2: {selected['sv_j']}  CN={selected['cn_j']:.2f}  "
                      f"rc={selected['rc_j']}\n")
            if path:
                out.write(f"\n  Shard path ({len(path)} hop"
                          f"{'s' if len(path) != 1 else ''}, "
                          f"{_path_total_bp(path)} bp total):\n")
                for pc, ps, pe in path:
                    shard_cn = next((cn for s, e, cn, *_ in chrom_segs.get(pc, [])
                                     if s == ps and e == pe), None)
                    cn_text = f"  CN={shard_cn:.2f}" if shard_cn is not None else ""
                    out.write(f"    {pc}:{ps}-{pe}  [{pe-ps} bp]{cn_text}\n")
            else:
                out.write("\n  Shard path: concordant (0 shard hops)\n")
            out.write(f"\n  Alternate directional side(s): "
                      f"{max(0, len(ev['candidates']) - 1)}\n")
            if ev['reason'] is None:
                out.write("  Injection status: injected\n")
            else:
                out.write(f"  Injection status: none [{ev['reason']}]\n")

        out.write(f"\n{'='*60}\n")
    return
