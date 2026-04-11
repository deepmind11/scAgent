#!/bin/bash
# Run scAgent (via pi) against all 7 Chromium scBench evals.
#
# This launches scAgent as a real pi agent with the full system prompt
# and skills loaded. Each eval gets a fresh session. The agent must
# interpret the task, choose the right tools, and produce an answer.
#
# Usage:
#   cd ~/Projects/scAgent
#   bash eval/run_llm_benchmark.sh
#
# Results saved to eval/results/llm_benchmark_TIMESTAMP/

set -euo pipefail

SCAGENT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCAGENT_DIR"

# Locate pi binary (same as scagent.sh)
FEYNMAN_ROOT="$(ls -d "$HOME/.local/share/feynman/feynman-"*"-darwin-arm64" 2>/dev/null | sort -V | tail -1)"
if [ -z "$FEYNMAN_ROOT" ]; then
  echo "Error: Feynman installation not found" >&2
  exit 1
fi
PI_BIN="$FEYNMAN_ROOT/node/bin/node $FEYNMAN_ROOT/app/.feynman/npm/node_modules/@mariozechner/pi-coding-agent/dist/cli.js"

export PI_CODING_AGENT_DIR="$SCAGENT_DIR/.pi"

# Output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="$SCAGENT_DIR/eval/results/llm_benchmark_$TIMESTAMP"
mkdir -p "$RESULTS_DIR"

# Activate venv for grading
source "$SCAGENT_DIR/.venv/bin/activate"

# Eval definitions: eval_id | data_path | answer_file | grader_file | task_prompt
declare -a EVALS=(
  "chromium_qc_4T1_filter_cells"
  "chromium_4t1_normalization"
  "chromium_4t1_hvg_gene_sets"
  "chromium_celltyping_01_4t1_compartment_fractions"
  "chromium_clustering_01_4t1_pericyte_adjacent_to_caf"
  "chromium_differential_expression_01_contractile_caf_marker_recovery"
  "chromium_trajectory_01_caf_terminal_marker_recovery"
)

EVAL_JSON_DIR="$HOME/Projects/scbench-eval/evals_canonical/chromium"
CACHE_DIR="$SCAGENT_DIR/.eval_cache/cache"

echo "========================================================================"
echo "scAgent LLM Benchmark"
echo "Time: $(date)"
echo "Results: $RESULTS_DIR"
echo "========================================================================"

PASSED=0
FAILED=0
ERRORS=0

for EVAL_ID in "${EVALS[@]}"; do
  echo ""
  echo "──────────────────────────────────────────────────────────"
  echo "EVAL: $EVAL_ID"
  echo "──────────────────────────────────────────────────────────"

  # Load eval JSON
  EVAL_FILE="$EVAL_JSON_DIR/${EVAL_ID}.json"
  if [ ! -f "$EVAL_FILE" ]; then
    echo "  ✗ Eval file not found: $EVAL_FILE"
    ERRORS=$((ERRORS + 1))
    continue
  fi

  # Extract task prompt and data node
  TASK_PROMPT=$(python3 -c "import json; print(json.load(open('$EVAL_FILE'))['task'])")
  DATA_NODE=$(python3 -c "import json; print(json.load(open('$EVAL_FILE'))['data_node'])")

  # Resolve cached data path
  DATA_FILE=$(python3 -c "
import sys
sys.path.insert(0, '$HOME/Projects/scbench-eval')
from latch_eval_tools.harness import download_single_dataset
print(download_single_dataset('$DATA_NODE', cache_name='.eval_cache'))
" 2>/dev/null)

  if [ ! -f "$DATA_FILE" ]; then
    echo "  ✗ Data file not found: $DATA_FILE"
    ERRORS=$((ERRORS + 1))
    continue
  fi

  echo "  Data: $(basename $DATA_FILE)"

  # Build the prompt for scAgent
  ANSWER_FILE="$RESULTS_DIR/${EVAL_ID}_answer.json"

  AGENT_PROMPT="You are running a scBench evaluation task. You must solve this single-cell RNA-seq analysis task.

DATA FILE: $DATA_FILE
Load it with: import anndata; adata = anndata.read_h5ad('$DATA_FILE')

Always inspect the data first: check adata.shape, adata.obs.columns, adata.layers.keys(), adata.X.max() to understand what you're working with.

IMPORTANT: The Python venv is already activated. Use our scagent.tools package when applicable (import from scagent.tools.*). You can also use scanpy directly.

TASK:
$TASK_PROMPT

After completing the analysis, write ONLY the answer JSON to this file:
$ANSWER_FILE

Do NOT include any text or explanation in the answer file — just the raw JSON object."

  # Run scAgent via pi in print mode
  START_TIME=$(date +%s)

  LOG_FILE="$RESULTS_DIR/${EVAL_ID}_agent.log"

  echo "$AGENT_PROMPT" | timeout 300 $PI_BIN \
    --print \
    --no-prompt-templates \
    --no-session \
    --model "claude-opus-4-6" \
    > "$LOG_FILE" 2>&1 || true

  END_TIME=$(date +%s)
  DURATION=$((END_TIME - START_TIME))

  echo "  Agent finished in ${DURATION}s"

  # Check if answer file was created
  if [ ! -f "$ANSWER_FILE" ]; then
    echo "  ✗ FAILED — no answer file produced"
    echo "  Log: $LOG_FILE"
    FAILED=$((FAILED + 1))
    continue
  fi

  # Grade the answer
  GRADER_FILE="$SCAGENT_DIR/eval/results/grader_${EVAL_ID}.json"
  if [ ! -f "$GRADER_FILE" ]; then
    echo "  ✗ Grader config not found: $GRADER_FILE"
    ERRORS=$((ERRORS + 1))
    continue
  fi

  GRADE_OUTPUT=$(python3 -c "
import json, sys
sys.path.insert(0, '$HOME/Projects/scbench-eval')
from latch_eval_tools.graders import GRADER_REGISTRY
try:
    answer = json.load(open('$ANSWER_FILE'))
    grader_config = json.load(open('$GRADER_FILE'))
    grader = GRADER_REGISTRY[grader_config['type']]()
    result = grader.evaluate_answer(answer, grader_config['config'])
    print(f'PASSED:{result.passed}')
    print(result.reasoning[:200])
except Exception as e:
    print(f'PASSED:False')
    print(f'Grading error: {e}')
" 2>&1)

  IS_PASSED=$(echo "$GRADE_OUTPUT" | head -1 | cut -d: -f2)
  REASONING=$(echo "$GRADE_OUTPUT" | tail -n +2)

  if [ "$IS_PASSED" = "True" ]; then
    echo "  ✓ PASSED (${DURATION}s)"
    echo "  $REASONING"
    PASSED=$((PASSED + 1))
  else
    echo "  ✗ FAILED (${DURATION}s)"
    echo "  $REASONING"
    FAILED=$((FAILED + 1))
  fi

  # Save result metadata
  python3 -c "
import json
result = {
    'eval_id': '$EVAL_ID',
    'passed': $( [ \"$IS_PASSED\" = \"True\" ] && echo 'True' || echo 'False' ),
    'duration_s': $DURATION,
    'reasoning': '''$REASONING'''[:500]
}
json.dump(result, open('$RESULTS_DIR/${EVAL_ID}_result.json', 'w'), indent=2)
"
done

# Summary
TOTAL=$((PASSED + FAILED + ERRORS))
echo ""
echo "========================================================================"
echo "RESULTS SUMMARY"
echo "========================================================================"
echo "Passed: $PASSED / $TOTAL"
echo "Failed: $FAILED"
echo "Errors: $ERRORS"
echo "Score:  $PASSED/$TOTAL ($(python3 -c "print(f'{$PASSED/$TOTAL*100:.0f}%')"))"
echo ""
echo "Results directory: $RESULTS_DIR"
echo ""
echo "Compare with direct adapter baseline: 4/7 (57%)"

# Write summary JSON
python3 -c "
import json
summary = {
    'timestamp': '$(date -u +%Y-%m-%dT%H:%M:%SZ)',
    'adapter': 'llm_scagent_pi',
    'model': 'claude-opus-4-6',
    'passed': $PASSED,
    'failed': $FAILED,
    'errors': $ERRORS,
    'total': $TOTAL,
    'score_pct': round($PASSED / max($TOTAL, 1) * 100, 1),
}
json.dump(summary, open('$RESULTS_DIR/summary.json', 'w'), indent=2)
print(f'Summary saved to: $RESULTS_DIR/summary.json')
"
