# BFB Graph Diagnosis Workflow

Use this workflow to evaluate whether an AmpliconArchitect graph supports a BFB call and to diagnose why the score changes under different assumptions.

## Quick Command

```bash
python tools/diagnose_bfb_graph.py path/to/*_graph.txt --show-tst-report --whole-graph
```

For a suspected wider region or broader foldback:

```bash
python tools/diagnose_bfb_graph.py path/to/*_graph.txt \
  --region chr7:52926926-55529670 \
  --fb-cutoff 50000 \
  --fb-cutoff 100000 \
  --show-tst-report
```

For a no-TST control:

```bash
python tools/diagnose_bfb_graph.py path/to/*_graph.txt --no-tst
```

For a no-deletion control:

```bash
python tools/diagnose_bfb_graph.py path/to/*_graph.txt --no-deletion --whole-graph
```

## Checklist

1. **Find candidate regions.**
   Compare `find_bfb_candidate_regions()` output to the visual amplicon span. If the region starts downstream of an obvious foldback, rerun with a manual wider `--region`.

2. **List raw foldbacks and broad foldbacks.**
   Inspect all same-strand intrachromosomal SVs. The current default foldback cutoff is 50 kb; real broad foldbacks may be wider. Test 100 kb when a same-strand SV has strong read support and a compatible CN change.

3. **Run TST report.**
   Check whether paired far-jumping SVs imply hidden foldbacks. Confirm which synthetic foldbacks survive region and flank-CN filters. Always compare with a no-TST control when TST calls are low-CN or numerous.

4. **Inspect CN segmentation.**
   Look for over-splitting of high-CN plateaus or tails. Current reconstruction merging uses an absolute CN floor plus a relative tolerance. Foldback cut points remain hard boundaries.

5. **Check deletion-edge correction when applicable.**
   Same-chromosome `DEL` edges represent copies that skip the interval between breakpoints. Deletion correction is enabled by default. If the skipped sequence segment has depressed CN relative to both flanks, run a `--no-deletion` control to confirm that default correction removes the artifact before segmentation and vector construction.

6. **Inspect `cn / lf / rf` vectors.**
   Verify foldbacks land on the expected segment boundary:
   - `--` foldbacks contribute to `lf` at `bp1`.
   - `++` foldbacks contribute to `rf` at `bp2`.
   If an expected foldback is absent, check foldback distance cutoff, flank-CN filter, ROI truncation, and local displaced-foldback rescue.

7. **Compare score modes.**
   Compare at least:
   - default candidate region
   - manual wider region when visually warranted
   - no-TST
   - no-deletion control when same-chromosome DEL edges are present
   - alternate foldback cutoff, often 100 kb
   - whole-graph mode for compact graphs

8. **Interpret score components.**
   A good CN fit with poor foldback score usually means the graph is BFB-like but the foldback evidence is missing, undercounted, filtered, displaced, or outside the default ROI.

## Common Failure Modes

### Candidate ROI truncates the BFB

Default candidate detection may start too late and exclude an upstream foldback. If whole-graph mode or a manually widened region improves the score and restores expected foldbacks, the default ROI is too narrow.

Example pattern:

```text
default region: sees one left foldback
wider region: sees two left foldbacks and one right foldback
score improves substantially
```

### Foldback cutoff is too strict

`SV.is_foldback()` defaults to 50 kb. Some real foldbacks are broader, especially in noisy or complex AA graphs.

If a same-strand SV is 50-100 kb wide, has read support, and matches the CN structure, rerun:

```bash
python tools/diagnose_bfb_graph.py graph.txt --region chr:start-end --fb-cutoff 100000
```

### TST foldbacks are hidden by far jumps

TST-masked foldbacks can appear as two far-jumping SVs whose far ends connect through small shard segments. The report should show:

```text
local breakpoints
far ends
shard path
injected FBI
```

If default mode should "reach outside" the ROI, enumerate TST candidates globally but only retain synthetic foldbacks that overlap the target region.

### Local displaced-foldback rescue can move evidence

The local rescue moves a foldback count from a landed block back to a boundary when the graph shows:

```text
CN drop at boundary -> non-foldback SV -> high-CN landed block -> foldback
```

Check verbose output for `local_foldback_rescue`. If a foldback disappears, verify that the rescue target slot exists before the original count is removed.

### Deletion edges depress sequence CN

AA graph sequence-edge CN counts copies that traverse the reference sequence. A same-chromosome `DEL` edge counts copies that skip the interval. In graph mode, deletion handling is enabled by default and adds the DEL-edge CN back to skipped sequence segments before merging and BFB vector construction. Use `--no-deletion` for a control run.

Example pattern:

```text
sequence segment:  chr:start-end  CN=1.6
DEL edge:          chr:start-1+->chr:end+1-  CN=5.9
corrected CN:      7.5
```

If correction removes a tiny low-CN shard from the model and improves the score, the uncorrected segmentation was likely an implementation artifact of using sequence-edge CN alone.

### Observed but undercounted foldbacks

Do not treat all foldback deficits as equivalent. There is a biological difference between:

```text
required foldback absent
required foldback present but AA SV CN is low
```

An observed foldback with low CN or modest reads should reduce confidence less severely than a completely absent foldback. See `TODO.md` for the proposed scoring improvement.

## How To Report A Diagnosis

Use this concise structure:

```text
Candidate region(s):
  chr:start-end

Raw foldbacks:
  -- chr:start-end CN=... rc=...
  ++ chr:start-end CN=... rc=...

TST-derived foldbacks:
  none / list

Vectors:
  cn = [...]
  lf = [...]
  rf = [...]

Scores:
  default = ...
  no-deletion = ...
  no-TST = ...
  wider region / alternate cutoff = ...

Interpretation:
  one paragraph explaining whether it is clean BFB, BFB-like complex, or weak/non-BFB,
  and which assumptions drive the score.
```
