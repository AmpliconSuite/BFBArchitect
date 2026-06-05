# TST-Foldback Identification in BFBArchitect

## Overview

Template-switching transposition (TST) events at BFB foldback junctions produce a characteristic pattern in AmpliconArchitect graph files: instead of a single clean foldback inversion (FBI) SV, the fold is represented as **two far-jumping SVs** whose local breakends sit close together on the BFB chromosome, while their far ends connect through one or more tiny "shard" segments on a remote locus. The BFBArchitect graph-input path detects this pattern, injects synthetic FBI SVs to represent the hidden fold, and then proceeds with normal BFB region detection and ILP reconstruction.

The graph-input path also has a second, related rescue pass for cases where a BFB boundary SV jumps out to a short segment whose far end contains a real foldback. This is not marked as a `.TST` synthetic SV. Instead, the observed foldback count is moved back onto the BFB boundary before segment simplification, so the BFB vector reflects the fold that is locally displaced by the jump.

TST detection runs before BFB vector construction. The local displaced-foldback rescue runs after TST reconciliation and before simplification/vector construction; in whole-graph mode it runs after primary chromosome identification.

---

## Where it fits in the graph-input pipeline

```
parse_graph_file()
       │
       ▼
find_tst_foldbacks()        ← TST detection; injects synthetic FBI SVs
       │
       ▼
find_bfb_candidate_regions() or whole_graph_as_region()
       │
       ▼
local foldback rescue       ← follows boundary jump to one segment with a foldback
       │
       ▼
subsect_graph_for_region()  ← builds cn / lf / rf vectors; uses injected SVs
       │
       ▼
BFBSolver (ILP)
```

`write_tst_report()` runs the TST pair-detection logic independently and outputs a human-readable per-event description for inspection. The report describes synthetic TST events; it does not currently report the later local displaced-foldback rescue.

---

## The core structural signature

A TST-masked fold leaves two discordant SVs in the graph:

```
BFB chromosome:
  ──────[pos_lo(+)] ··· [pos_hi(-)]──────
            │                  │
            └──── far SV A ────┘  (jump to remote locus)
                                   remote locus: shard₁ ──── shard₂ ──── …
            ┌──── far SV B ────────────────────────────────────────────┘
            │
  [pos_lo still in BFB region]

```

The two far ends connect through a chain of small segments (shards). The synthetic FBI that BFBArchitect injects spans `pos_lo` to `pos_hi` on the BFB chromosome, standing in for the fold that those two SVs collectively represent.

Two local breakend topologies arise depending on whether the shard lies between the breakpoints or at them:

| Topology | dir_lo | dir_hi | Synthetic FBI coordinates |
|----------|--------|--------|--------------------------|
| A — shard fills the gap | `+` | `−` | `(pos_lo, pos_hi − 1)` |
| B — breakpoints are shard endpoints | `−` | `+` | `(pos_lo − 1, pos_hi)` |

---

## Algorithm: `find_tst_foldbacks`

**Inputs:**
- `svs` — list of `(SV, cn_float, read_count)` from the parsed graph
- `chrom_segs` — dict of chromosome → sorted segment list
- `shard_max_bp` = 5,000 bp — maximum size of a shard segment
- `max_hops` = 5 — maximum shard traversals in the BFS path
- `fb_dist` = 50,000 bp — maximum distance between the two local breakends
- `far_min` = 500,000 bp — minimum jump distance to qualify as "far-jumping"

**Output:** `svs` extended with `(synthetic_SV, cn_tst, 0)` tuples for each confirmed TST event.

```
PROCEDURE find_tst_foldbacks(svs, chrom_segs):

  BUILD lookup tables:
    bp_to_seg : (chrom, pos, strand) → (start, end) of the segment at that breakend
    sv_index  : (chrom, pos, strand) → list of SV entries touching that breakend

  BUILD exclusion sets:
    existing_fb_cuts : set of (chrom, cut_pos) for every direct (non-TST) foldback
      # ++ foldback cuts at bp2; -- foldback cuts at bp1
    real_fb_spans[chrom] : list of (bp1, bp2) for every direct foldback
      # only direct foldbacks that pass the raw flank-CN filter are used

  COLLECT far-jumping breakends, indexed by local chromosome:
    FOR each SV in svs:
      FOR each end (local, far) of the SV:
        IF far.chrom ≠ local.chrom  OR  |far.pos − local.pos| ≥ far_min:
          far_by_chrom[local.chrom].append(
              (local.pos, local.strand, far_endpoint, SV, cn, rc))

  seen_pairs = {}
  synthetic  = []

  FOR chrom, breakends in far_by_chrom:
    FOR every pair (i, j) of breakends on chrom  [i < j]:

      FILTER — skip if:
        sv_i is sv_j                          (same SV object)
        |pos_i − pos_j| > fb_dist             (local ends too far apart)
        dir_i == dir_j                        (need opposite orientations)
        far_i == far_j                        (degenerate: same far end)
        pair (sv_i, sv_j) already processed

      GUARD — skip if far ends are direct foldback endpoints:
        IF far_i.chrom == far_j.chrom:
          p_lo, p_hi = min/max(far_i.pos, far_j.pos)
          IF any real_fb_span on that chrom matches (p_lo ± 5, p_hi ± 5):
            SKIP  # far-end BFS would traverse the foldback's own hairpin

      RUN shard-path BFS:
        path = _shard_path_reachable(
                   start  = far_i endpoint,
                   target = far_j endpoint,
                   exclude_svs = {sv_i, sv_j})   # prevent circular traversal

      IF path is None: SKIP  # no valid shard connection

      COMPUTE synthetic FBI coordinates from (pos_lo, dir_lo, pos_hi):
        IF dir_lo == '+':  fbi_bp1, fbi_bp2 = pos_lo,     pos_hi − 1   # Topology A
        ELSE:              fbi_bp1, fbi_bp2 = pos_lo − 1, pos_hi        # Topology B

      cn_tst = min(cn_i, cn_j)
      synth = SV(chrom, fbi_bp1, '+', chrom, fbi_bp2, '+', TST=True)
      IF synth flank CN change is known AND < 1.0 copy:
        SKIP local synthetic FBI
      ELSE:
        synthetic.append(synth, cn_tst, rc=0)

      IF far_i.chrom == far_j.chrom  AND  |far_i.pos − far_j.pos| ≤ fb_dist:
        # The fold also manifests at the remote locus; inject an FBI there too
        COMPUTE far-side FBI from (far_i, far_j) using the same topology logic
        IF far-side cut point is NOT in existing_fb_cuts:
          IF far-side flank CN change is known AND < 1.0 copy:
            SKIP far-side synthetic FBI
          ELSE:
            synthetic.append(far-side FBI, cn_tst, rc=0)

  DEDUPLICATE synthetic list:
    For SVs with the same (chrom, bp1, bp2) string representation,
    keep the entry with the highest CN.

  RETURN svs + synthetic
```

Synthetic foldbacks are deduplicated by coordinate before returning. If adding TST foldbacks would push the graph above the reconstruction foldback cap, the graph-input path falls back to the raw SV set and does not use the injected TST foldbacks.

---

## Subroutine: `_shard_path_reachable`

BFS over the graph topology from `start` to `target`, restricted to shard-sized segments. Movements allowed at each node `(chrom, pos, strand)`:

1. **Segment traversal** (costs 1 hop): if the segment at `curr` is ≤ `shard_max_bp`, cross it to the other end.
2. **Concordant edge** (free): move from `pos(+)` to `(pos+1)(−)` or vice versa — the AA convention for adjacent sequence. Accepted only if it leads directly to `target` or into another shard endpoint.
3. **Discordant SV edge** (free): follow any SV in `sv_index[curr]` to its other breakend, but **only** if that breakend belongs to a shard-sized segment. SVs in `exclude_svs` are never followed, preventing circular paths where the candidate SVs themselves form the "path."

```
PROCEDURE _shard_path_reachable(start, target, bp_to_seg, sv_index,
                                 shard_max_bp, max_hops, exclude_svs=∅):

  queue = [(curr=start, path=[], visited={})]
  seen  = {}

  WHILE queue not empty:
    curr, path, visited = dequeue

    IF curr == target: RETURN path   # success

    IF (curr, visited) already seen: SKIP
    IF len(path) ≥ max_hops: SKIP

    # 1. Segment traversal
    seg = bp_to_seg[curr]
    IF seg exists AND seg.size ≤ shard_max_bp AND seg not in visited:
      other_end = opposite end of seg
      enqueue(other_end, path + [seg], visited ∪ {seg})

    # 2. Concordant edge
    adj = (chrom, pos+1, '−') if strand=='+' else (chrom, pos−1, '+')
    IF adj == target:
      enqueue(adj, path, visited)
    ELSE IF bp_to_seg[adj] exists AND its size ≤ shard_max_bp:
      enqueue(adj, path, visited)

    # 3. Discordant SV edges
    FOR sv in sv_index[curr]:
      IF sv in exclude_svs: SKIP
      other = other breakend of sv
      IF bp_to_seg[other] exists AND its size ≤ shard_max_bp:
        enqueue(other, path, visited)

  RETURN None   # no path found
```

The `exclude_svs` parameter is the key guard against false positives: without it, the BFS could use sv_i to reach a shard, then use sv_j to reach the target — a circular path that trivially "confirms" any pair of SVs sharing a shard on their far ends, regardless of whether those SVs jointly represent a fold.

---

## Duplicate-foldback guard

A pair is rejected before BFS if the far ends of sv_i and sv_j are within ±5 bp of the endpoints of an existing direct FBI SV on the same chromosome. In this case the BFS would traverse the direct foldback's own tiny hairpin segment and return a spurious "path." The fold is already represented by the direct FBI; injecting a synthetic one would duplicate it.

Direct foldbacks used for this duplicate guard must pass the raw foldback flank-CN filter. This prevents a zero-CN-change hairpin from blocking a better-supported TST interpretation.

---

## Foldback CN-change filters

Foldbacks are used for reconstruction only if they show a copy-number change across their outside flanks, with separate thresholds for raw and synthetic evidence:

| Foldback source | Required flank CN change | Behavior if flanks are missing |
|-----------------|--------------------------|--------------------------------|
| Raw graph FBI | `>= 0.2` copies | Keep conservatively |
| Synthetic TST FBI | `>= 1.0` copy | Keep conservatively |

For a `++` foldback, the flanks are the segment ending at `bp1` and the segment starting at `bp2 + 1`. For a `--` foldback, the flanks are the segment ending at `bp1 - 1` and the segment starting at `bp2`. If either flank is not found exactly in the graph segments, the foldback is retained because the CN context is ambiguous.

Raw graph foldbacks that pass this filter are counted at least once in the `lf`/`rf` vectors, even when their SV CN would otherwise round to 0.

---

## Local displaced-foldback rescue

Some BFB-like graphs do not show the foldback directly at the BFB boundary after TST reconciliation. Instead, the boundary has a non-foldback SV that jumps away to a short local block, and the foldback sits at the far end of that block. Segment simplification can erase the boundary relationship before the graph-input vector is built, so BFBArchitect detects this pattern on the raw graph segments first.

This pass is intentionally narrow:

1. Work on the primary BFB chromosome after TST reconciliation. In default region mode, anchors must lie inside the candidate region. In whole-graph mode, the primary chromosome is chosen first, then anchors are searched.
2. Look for a boundary where a high-CN primary-chromosome segment drops to an outside neighbor by at least `1.25x`.
3. From that boundary endpoint, follow one non-foldback SV to its landed endpoint.
4. Follow exactly one landed sequence block in the direction implied by the landed endpoint.
5. Require the landed block to be at least `1.25x` higher CN than the neighbor on the far side.
6. Require a raw, non-TST foldback at the far end of that landed block.
7. Require that foldback to pass the raw `>= 0.2` flank-CN-change filter.

When the pattern is found, BFBArchitect adds a cut point at the original BFB boundary before simplification. During vector construction, it removes the foldback count from the displaced foldback's original slot if that slot exists, and adds the count to the boundary slot instead. The foldback itself is marked `local_resolved`, not `.TST`.

This implements the "jump out, get one segment with a foldback, and move the fold back to the boundary" case without trusting the small fragment's absolute CN. The 1.25x CN checks are used only as structural support for the boundary and landed-block pattern.

---

## Injected SV properties

Each synthetic SV is a standard `SV` object on the BFB chromosome with:
- `strand1 = strand2 = '+'` (always a `++` foldback)
- `.TST = True` (flag distinguishing it from real graph edges)
- `cn` = `min(cn_i, cn_j)` (conservative: the weaker of the two anchoring SVs)
- `rc = 0` (no direct read count; it is inferred from the SV pair)

The injected SV enters the `svs` list seen by `subsect_graph_for_region`, which then counts it toward the `lf`/`rf` foldback vectors exactly as it would a real FBI, enabling the ILP to reconstruct the correct BFB structure.

---

## Parameters and their role

| Parameter | Default | Effect |
|-----------|---------|--------|
| `shard_max_bp` | 5,000 bp | Upper size limit for a shard segment; larger values allow noisier TST junctions |
| `max_hops` | 5 | Maximum number of shard segments in the connecting path |
| `fb_dist` | 50,000 bp | Maximum separation of the two local breakends (mirrors the `SV.is_foldback` threshold) |
| `far_min` | 500,000 bp | Minimum jump distance; filters out local SVs that are not far-jumping |
| raw foldback flank CN change | 0.2 copies | Minimum CN change across raw foldback flanks to use it in reconstruction |
| synthetic TST flank CN change | 1.0 copy | Minimum CN change across synthetic TST foldback flanks to inject it |
| local rescue CN ratio | 1.25x | Required CN support at the BFB boundary and landed block |
| graph foldback cap | 50 foldbacks | TST/local rescue additions are suppressed if they would exceed this cap |

---

## `write_tst_report`

Runs the candidate-pair enumeration and BFS independently (it does not call `find_tst_foldbacks`), and emits a structured plain-text report per graph file. Each confirmed pair/path event is annotated with:
- Local breakend positions, gap size, and topology
- SV identity, CN, and read count for each anchor SV
- Far-end coordinates, the segment they land in (size, CN), and whether the far end falls inside a detected BFB region
- The shard path (number of hops, total shard bp, per-segment details)
- The synthetic FBI coordinates implied by the pair/path topology
- Whether a far-side FBI was also computed, and why the direct-foldback duplicate guard would block it

This report is the primary tool for manually inspecting TST candidates before running the full reconstruction. The reconstruction path applies additional CN-change and foldback-cap gates, and the report does not currently describe the later local displaced-foldback rescue.
