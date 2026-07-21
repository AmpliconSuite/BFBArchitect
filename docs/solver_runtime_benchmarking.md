# Solver Runtime Benchmarking

This document describes a reusable BFBArchitect solver runtime benchmark workflow
for comparing Gurobi, MOSEK, and CBC on graph-mode BFBArchitect cases produced
by AmpliconClassifier or another pipeline.

## Purpose

The benchmark is intended to answer four runtime questions:

1. How does solver choice generally affect runtime across solver thread counts?
2. How much does increasing the solver thread count help each solver type?
3. Is one solver consistently faster across both short and long-running cases?
4. What fraction of cases is won by each solver, optionally by thread count?

The benchmark uses BFBArchitect graph-mode reconstruction rather than invoking
the full AmpliconClassifier workflow. This keeps the benchmark focused on graph
parsing, ILP model construction, and solver runtime.

## Inputs

Provide one or more input directories. Each input may be either:

- an AmpliconClassifier output directory containing `bfbarchitect_outputs/`, or
- a direct BFBArchitect graph directory containing `*_BFB_graph.txt` files.

Set input paths explicitly in your shell:

```bash
export AC_RUN_1=/path/to/ac_run_1
export AC_RUN_2=/path/to/ac_run_2
```

The harness inventories graph files, reconstructs candidate graph-mode cases
using the same BFBArchitect graph-region logic, ranks cases by inferred ILP size
(`T`, segment count, foldback CN), and benchmarks the selected cases. The default
selection is the top `--limit` cases after applying `--min-t`; use `--all-cases`
to benchmark every inventoried case passing `--min-t`.

## Solver And Thread Handling

BFBArchitect solver selection behavior:

- Automatic fallback order in `detect_solver()` is Gurobi, then MOSEK, then CBC.
- Explicit `--solver` bypasses fallback.
- Gurobi receives the thread count via `m.Params.Threads`.
- MOSEK receives the thread count via `m.setSolverParam("numThreads", ...)`.
- CBC receives the thread count via `PULP_CBC_CMD(..., threads=...)`.

The benchmark invokes explicit solvers, so fallback behavior does not affect
which solver is measured.

## Environment

Run from an environment where BFBArchitect and the requested solver bindings are
installed. If MOSEK is not available in that environment, install it with:

```bash
conda install -c mosek mosek
```

Gurobi may need license-server network access. If running through a sandboxed
agent, launch commands with network access enabled or run them directly in a
normal shell.

For non-interactive runs, set `PYTHON_BIN` if the desired Python is not first on
`PATH`:

```bash
export PYTHON_BIN=/path/to/env/bin/python
```

## Quick Inventory

Build only the case inventory:

```bash
python scripts/solver_runtime_analysis.py \
  --inventory-only \
  --graph-dir "$AC_RUN_1" \
  --graph-dir "$AC_RUN_2" \
  --out-dir reports/solver_runtime_analysis
```

Output:

```bash
reports/solver_runtime_analysis/case_inventory.tsv
```

Generated `reports/solver_runtime*` directories are ignored by git.

## Broad 3-Thread Benchmark

This smaller benchmark provides a stable 3-thread baseline. It runs the top 20
combined cases, all three solvers, 3 solver threads, a 900 second timeout, and a
total active solver-thread budget of 15.

```bash
python scripts/solver_runtime_analysis.py \
  --graph-dir "$AC_RUN_1" \
  --graph-dir "$AC_RUN_2" \
  --out-dir reports/solver_runtime_analysis \
  --benchmark-file reports/solver_runtime_analysis/broad_3thread_900s.tsv \
  --threads 3 \
  --timeout 900 \
  --limit 20 \
  --min-t 0 \
  --max-active-threads 15 \
  --resume
```

For replicates, change only the `--benchmark-file`, for example:

```bash
--benchmark-file reports/solver_runtime_analysis/broad_3thread_900s_replicate2.tsv
--benchmark-file reports/solver_runtime_analysis/broad_3thread_900s_replicate3.tsv
```

Summarize three replicates:

```bash
python scripts/solver_runtime_analysis.py \
  --graph-dir "$AC_RUN_1" \
  --graph-dir "$AC_RUN_2" \
  --out-dir reports/solver_runtime_analysis \
  --summarize-replicates \
  reports/solver_runtime_analysis/broad_3thread_900s.tsv \
  reports/solver_runtime_analysis/broad_3thread_900s_replicate2.tsv \
  reports/solver_runtime_analysis/broad_3thread_900s_replicate3.tsv
```

Outputs:

```bash
reports/solver_runtime_analysis/replicate_summary.tsv
reports/solver_runtime_analysis/replicate_summary.md
```

## Full Thread-Matrix Benchmark

The thread-matrix benchmark is designed to answer the thread-scaling and
solver-fraction questions. By default the wrapper runs:

- top 20 combined inventoried cases
- solvers: `gurobi`, `mosek`, `cbc`
- solver threads: `1, 2, 3, 4, 8, 16`
- timeout: 900 seconds per job
- replicates: 3
- max active declared solver threads: 16

Run it with:

```bash
bash scripts/run_solver_thread_matrix.sh \
  reports/solver_runtime_thread_matrix \
  "$AC_RUN_1" \
  "$AC_RUN_2"
```

The script is resumable because each replicate uses
`scripts/solver_runtime_analysis.py --resume`. If interrupted, rerun the same
command and completed `(source, case, solver, threads)` rows will be skipped.

Wrapper settings are controlled with environment variables:

```bash
PYTHON_BIN=/path/to/env/bin/python \
THREADS="1 2 3 4 8 16" \
REPLICATES="1 2 3" \
TIMEOUT=900 \
LIMIT=20 \
MIN_T=0 \
MAX_ACTIVE_THREADS=16 \
bash scripts/run_solver_thread_matrix.sh reports/solver_runtime_thread_matrix "$AC_RUN_1" "$AC_RUN_2"
```

Primary outputs:

```bash
reports/solver_runtime_thread_matrix/thread_matrix_replicate1.tsv
reports/solver_runtime_thread_matrix/thread_matrix_replicate2.tsv
reports/solver_runtime_thread_matrix/thread_matrix_replicate3.tsv
reports/solver_runtime_thread_matrix/solver_runtime_plots.html
reports/solver_runtime_thread_matrix/plots/
reports/solver_runtime_thread_matrix/plot_summary.tsv
```

The plotted HTML report includes:

- runtime ECDF by solver and thread count, including a less-overlapped
  small-multiple ECDF view with shared x limits
- median runtime heatmap by solver and solver-thread count
- absolute runtime by thread count
- speedup vs 1 thread by solver
- fastest-solver fraction by thread count
- within-10%-of-fastest fraction by thread count
- pairwise runtime-ratio distributions
- replicate stability scatter
- BFB score by amplicon, solver, and thread count
- BFB score range by amplicon

## Adding Extra Long-Running Samples

The matrix benchmark uses the top cases ranked by inferred ILP size. Future runs
can include additional known long-running cases without replacing the top
selection.

Create a manifest with one `case_id` per line:

```bash
cat > reports/solver_runtime_thread_matrix/extra_long_cases.txt <<'EOF'
sample_a_amplicon1_whole_graph
sample_b_amplicon3_region1
sample_c_amplicon2_region2
EOF
```

Then pass it to `scripts/solver_runtime_analysis.py`:

```bash
python scripts/solver_runtime_analysis.py \
  --graph-dir "$AC_RUN_1" \
  --graph-dir "$AC_RUN_2" \
  --out-dir reports/solver_runtime_thread_matrix_augmented \
  --benchmark-file reports/solver_runtime_thread_matrix_augmented/thread_matrix_replicate1.tsv \
  --thread-list 1 2 3 4 8 16 \
  --timeout 900 \
  --limit 20 \
  --min-t 0 \
  --max-active-threads 16 \
  --include-case-id-file reports/solver_runtime_thread_matrix/extra_long_cases.txt \
  --resume
```

Use `--include-case-id <case_id>` for one-off additions. Use `--case-id`
instead only when you want to restrict the benchmark to a specific case or set
of cases.

## Sampling Strategy For Future Runs

The top-`N` complexity strategy is useful for stressing long solver runtimes,
but it is not sufficient by itself for score-stability analysis. Future expanded
runs should use a stratified selection:

- hard cases: top 50 or more cases by inferred ILP size
- easy cases: 25-50 low-complexity cases, sampled from the lower end of the
  inventory, but not no-op cases. Require enough structure to exercise
  BFBArchitect, for example at least 3 segments and nonzero foldback signal.
- cutoff-near cases: any cases with prior BFBArchitect scores near the BFB
  cutoff, currently score <= 2.8
- disagreement cases: any cases where solver or thread choice changes whether
  the score is above or below 2.8

This keeps runtime conclusions anchored on difficult examples while also testing
whether solver/thread choices perturb classification for routine cases.

If reusing a completed hard-case run, seed the new output directory with the
existing replicate TSVs by setting `SEED_FROM_DIR`, and pass extra cases through
`INCLUDE_CASE_ID_FILE`. The runner uses `--resume`, so completed
`(source, case, solver, threads)` rows are skipped and only the additional cases
run.

```bash
SEED_FROM_DIR=reports/solver_runtime_thread_matrix_expanded_50 \
INCLUDE_CASE_ID_FILE=reports/solver_runtime_thread_matrix_hard50_easy25_include_cases.txt \
LIMIT=50 \
bash scripts/run_solver_thread_matrix.sh reports/solver_runtime_thread_matrix_hard50_easy25 "$AC_RUN_1" "$AC_RUN_2"
```

## Running Detached

For a long unattended run, use a terminal multiplexer or a user systemd
transient service. A direct shell or multiplexer is simplest:

```bash
PYTHON_BIN=/path/to/env/bin/python \
bash scripts/run_solver_thread_matrix.sh reports/solver_runtime_thread_matrix "$AC_RUN_1" "$AC_RUN_2" \
  > reports/solver_runtime_thread_matrix.log 2>&1
```

If using `systemd-run --user`, make sure `PYTHON_BIN`, `AC_RUN_1`, and
`AC_RUN_2` are set inside the command because user services may not inherit
the interactive environment:

```bash
systemd-run --user \
  --unit=bfb-solver-thread-matrix \
  --collect \
  --same-dir \
  /bin/bash -lc 'PYTHON_BIN=/path/to/env/bin/python AC_RUN_1=/path/to/ac_run_1 AC_RUN_2=/path/to/ac_run_2 bash scripts/run_solver_thread_matrix.sh reports/solver_runtime_thread_matrix "$AC_RUN_1" "$AC_RUN_2" > reports/solver_runtime_thread_matrix.log 2>&1'
```

Check progress:

```bash
tail -f reports/solver_runtime_thread_matrix.log
wc -l reports/solver_runtime_thread_matrix/thread_matrix_replicate*.tsv
```

## Notes On Interpretation

The benchmark supports multiplexing with `--max-active-threads`. This improves
throughput but means individual runtimes reflect a controlled shared-load
setting rather than a fully isolated machine. For final publication-quality
numbers, prefer one of these approaches:

- keep multiplexing fixed and report it explicitly, or
- rerun the final subset sequentially with `--max-active-threads 0`.
