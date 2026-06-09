# BFBArchitect — Design TODO

## 1. Polarity reversal

**Status:** Done

BFB reconstruction currently assumes a fixed arm orientation. However, if a chromosomal arm undergoes translocation and inversion prior to BFB formation, the BFB sequence would run in the opposite direction relative to the reference. The solver should optionally assess both orientations and return the better-scoring result.

Implemented:
- `--reverse_polarity` flips segment order and swaps left/right foldback vectors before ILP reconstruction.
- Reverse-polarity solver paths are mapped back to the original segment numbering before scoring/output.
- Graph-mode result dictionaries record whether `reverse_polarity` was used.

Possible future extension:
- Add automatic bidirectional assessment that runs both orientations and reports the better-scoring result.
- Decide how to handle cases where both orientations score similarly (ties).

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

---

## 4. Distinguish missing versus undercounted foldbacks in scoring

**Status:** Future

The current BFB score does not clearly separate two different cases:

- A foldback required by the reconstructed BFB string is absent from the graph.
- A foldback exists at the right boundary, but its graph-estimated SV copy number is lower than the BFB string implies.

These should not carry the same evidence weight. An absent foldback is stronger negative evidence than an observed foldback with noisy or underestimated read support. This matters in AA graph mode, where short junctions, repetitive sequence, mapping artifacts, and local graph complexity can make breakpoint copy number noisier than segment copy number.

Proposed scoring model:

- Treat foldback boundary presence as the primary signal.
- Penalize expected foldbacks with no matching observed foldback using a higher `absent_foldback_weight`.
- Penalize observed-but-undercounted foldbacks with a lower `undercounted_foldback_weight`.
- Let read-pair support and boundary proximity soften the undercount penalty.
- Keep overcount penalties separate from undercount penalties; extra observed foldback CN may indicate subclonality or graph noise rather than a missing BFB operation.

Sketch:

```text
if expected_fb > 0 and observed_fb == 0:
    penalty += absent_foldback_weight * expected_fb
elif expected_fb > observed_fb:
    support_factor = f(read_pairs, breakpoint_quality)
    penalty += undercounted_foldback_weight * (expected_fb - observed_fb) * support_factor
```

Initial defaults to evaluate:

- `absent_foldback_weight`: current missing-foldback penalty scale.
- `undercounted_foldback_weight`: 25-50% of absent penalty.
- `support_factor`: decreases with read support, capped so strong read evidence cannot erase all copy-number mismatch.

Validation cases:

- Graphs with low-CN but real foldbacks, such as `LP6005409-DNA_B02_amplicon1_graph.txt`.
- Graphs with no foldback evidence at required boundaries, where the high absent-foldback penalty should remain.
