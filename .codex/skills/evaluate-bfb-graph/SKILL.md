---
name: evaluate-bfb-graph
description: Diagnose BFB calls from AmpliconArchitect graph files. Use when asked to evaluate an AA *_graph.txt file for BFB evidence, explain why a BFB call scores well or poorly, inspect candidate BFB regions, compare default versus wider context, test foldback cutoff sensitivity, inspect TST-derived foldbacks, or produce cn/lf/rf vectors and BFB scores.
---

# Evaluate BFB Graph

Use this skill to diagnose BFB evidence in AmpliconArchitect graph files in this repository.

## Primary Tool

Run the repository diagnostic CLI:

```bash
python tools/diagnose_bfb_graph.py /path/to/*_graph.txt --show-tst-report --whole-graph
```

For a suspected manual region or foldback cutoff issue:

```bash
python tools/diagnose_bfb_graph.py /path/to/*_graph.txt \
  --region chr7:52926926-55529670 \
  --fb-cutoff 50000 \
  --fb-cutoff 100000 \
  --show-tst-report
```

For a no-TST control:

```bash
python tools/diagnose_bfb_graph.py /path/to/*_graph.txt --no-tst
```

Read `docs/bfb_graph_diagnosis.md` for the full workflow and interpretation checklist.

## Required Analysis Pattern

1. Identify candidate BFB regions.
2. List raw foldbacks, including same-strand broad foldbacks that fail the default 50 kb cutoff.
3. Check TST report and note which synthetic foldbacks survive filters.
4. Inspect segmentation and `cn / lf / rf` vectors.
5. Compare default mode to relevant controls:
   - no-TST
   - wider manual region if ROI seems truncated
   - 50 kb versus 100 kb foldback cutoff when broad same-strand SVs exist
   - whole-graph mode for compact graphs
6. Report score components when available.
7. Distinguish absent foldbacks from observed but undercounted foldbacks.

## Reporting Style

Keep the final diagnosis concise:

```text
Candidate region(s):
Raw foldbacks:
TST-derived foldbacks:
Vectors:
Scores:
Interpretation:
```

Call out likely implementation artifacts separately from biological interpretation, especially:

- ROI truncation.
- Overly strict foldback distance cutoff.
- TST chains outside the ROI.
- CN segment over-splitting.
- Local displaced-foldback rescue moving evidence.
- Foldback CN undercounting versus true missing foldbacks.
