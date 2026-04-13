#!/bin/bash
# Launch scAgent — single-cell RNA-seq analysis agent
#
# Prerequisites:
#   - Node.js ≥ 20.19 (for the agent runtime)
#   - Python ≥ 3.11 with: pip install -e .
#
# Usage:
#   ./scagent.sh                              # default
#   ./scagent.sh --model opus                 # Opus model
#   ./scagent.sh --model opus --thinking max  # max thinking
#   ./scagent.sh -c                           # continue previous session
#   ./scagent.sh -r                           # resume/pick a session

set -euo pipefail

SCAGENT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# Point the agent at scAgent's config (system prompt, skills, settings)
export PI_CODING_AGENT_DIR="$SCAGENT_DIR/.pi"

# --- Find the pi CLI ---

# 1. pi in PATH (e.g. npm install -g @mariozechner/pi-coding-agent)
if command -v pi >/dev/null 2>&1; then
  exec pi "$@"
fi

# 2. npx fallback (no global install needed — just Node.js)
if command -v npx >/dev/null 2>&1; then
  echo "Note: Running pi-coding-agent via npx (first run may be slow)..." >&2
  exec npx --yes @mariozechner/pi-coding-agent "$@"
fi

echo "Error: Could not find pi or npx." >&2
echo "" >&2
echo "Install Node.js (≥20.19) from https://nodejs.org, then either:" >&2
echo "  npm install -g @mariozechner/pi-coding-agent   # global install" >&2
echo "  npx @mariozechner/pi-coding-agent              # no install needed" >&2
exit 1
