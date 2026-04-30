# BFBArchitect — Design TODO

## 1. Bidirectional BFB assessment

**Status:** Future

BFB reconstruction currently assumes a fixed arm orientation. However, if a chromosomal arm undergoes translocation and inversion prior to BFB formation, the BFB sequence would run in the opposite direction relative to the reference. The solver should optionally assess both orientations and return the better-scoring result.

Design thoughts:
- Flip the segment order + strand orientation before passing to the ILP, run both, take the lower score.
- Could be a flag (`--bidirectional`) or always-on with the better result reported.
- Need to decide how to handle cases where both orientations score similarly (ties).

---

## 2. Manual thread count control

**Status:** Future

Certain steps (BAM parsing, Gurobi, batch runs) are either single-threaded or use their own internal thread pool. A `--threads` / `-t` flag would let users cap or expand resource usage.

Design thoughts:
- Propagate a `threads` argument through `BFBArchitect.py` → `SVCaller`, Gurobi `Params.Threads`, `batch_run.py` pool size.
- Default: autodetect (e.g. `os.cpu_count()`).

---

## 3. Simplified library invocation (next)

**Status:** In progress — Jens working on this

The current graph-file library API requires callers to chain four separate calls (`find_bfb_candidate_regions` → `subsect_graph_for_region` → `trim_background_segments` → `reconstruct_bfb`) plus manual None-checking and centromere lookup. This is too much boilerplate for external tools (e.g. AmpliconClassifier).

Proposed: a single `reconstruct_bfb_from_graph_file(graph_file, ...)` convenience wrapper that:
1. Detects regions.
2. Subsets + trims each region.
3. Calls the ILP.
4. Returns a list of `(region, BFB_strings, scores, multiplicity)` — skipping None/empty regions internally.

The lower-level functions remain public for callers that need fine-grained control.
