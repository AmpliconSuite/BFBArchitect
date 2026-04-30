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

## 3. Path constraints from graph files

**Status:** Future — format TBD

Some AA-format graph files include a path constraints section; others do not. The current `parse_graph_file` ignores it entirely, and the ILP has no mechanism to enforce such constraints.

Design thoughts: deferred until format documentation is available.
