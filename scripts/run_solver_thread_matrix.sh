#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/run_solver_thread_matrix.sh OUT_DIR INPUT_DIR [INPUT_DIR ...]

Runs the BFBArchitect solver/thread benchmark matrix and then renders the plot report.
Each INPUT_DIR may be either an AC output directory containing bfbarchitect_outputs/
or a direct directory of *_BFB_graph.txt files.

Environment overrides:
  PYTHON_BIN          Python executable to use (default: python)
  THREADS            Space-separated solver thread counts (default: "1 2 3 4 8 16")
  REPLICATES         Space-separated replicate ids (default: "1 2 3")
  TIMEOUT            Per-job timeout in seconds (default: 900)
  LIMIT              Number of top inventory cases to benchmark, or "all" (default: 20)
  MIN_T              Minimum inferred T value to include (default: 0)
  MAX_ACTIVE_THREADS Maximum concurrent declared solver threads (default: 16)
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -lt 2 ]]; then
    usage >&2
    exit 2
fi

OUT_DIR="$1"
shift
PYTHON_BIN="${PYTHON_BIN:-python}"
read -r -a THREAD_VALUES <<< "${THREADS:-1 2 3 4 8 16}"
read -r -a REPLICATE_VALUES <<< "${REPLICATES:-1 2 3}"
TIMEOUT="${TIMEOUT:-900}"
LIMIT="${LIMIT:-20}"
MIN_T="${MIN_T:-0}"
MAX_ACTIVE_THREADS="${MAX_ACTIVE_THREADS:-16}"

mkdir -p "$OUT_DIR"

GRAPH_ARGS=()
for input_dir in "$@"; do
    GRAPH_ARGS+=(--graph-dir "$input_dir")
done

LIMIT_ARGS=(--limit "$LIMIT")
if [[ "$LIMIT" == "all" ]]; then
    LIMIT_ARGS=(--all-cases)
fi

for rep in "${REPLICATE_VALUES[@]}"; do
    "$PYTHON_BIN" scripts/solver_runtime_analysis.py \
        "${GRAPH_ARGS[@]}" \
        --out-dir "$OUT_DIR" \
        --benchmark-file "$OUT_DIR/thread_matrix_replicate${rep}.tsv" \
        --thread-list "${THREAD_VALUES[@]}" \
        --timeout "$TIMEOUT" \
        "${LIMIT_ARGS[@]}" \
        --min-t "$MIN_T" \
        --max-active-threads "$MAX_ACTIVE_THREADS" \
        --resume
done

PLOT_ARGS=()
for rep in "${REPLICATE_VALUES[@]}"; do
    PLOT_ARGS+=(--replicate "$OUT_DIR/thread_matrix_replicate${rep}.tsv")
done

"$PYTHON_BIN" scripts/plot_solver_runtime_report.py \
    "${PLOT_ARGS[@]}" \
    --out-dir "$OUT_DIR"
