#!/bin/bash
# Launch scAgent — single-cell RNA-seq analysis agent
#
# scAgent runs on pi-coding-agent. On first launch, type /login in the
# REPL to authenticate with your Claude subscription (no API key needed).
#
# Prerequisites:
#   - Node.js ≥ 20.19 (for the agent runtime)
#   - Python ≥ 3.11 with: pip install -e .
#
# Usage:
#   ./scagent.sh                              # default
#   ./scagent.sh --model opus                 # specific model
#   ./scagent.sh --thinking max               # max thinking
#   ./scagent.sh -c                           # continue previous session
#   ./scagent.sh -r                           # resume/pick a session

set -euo pipefail

SCAGENT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# --- Find the pi CLI and launch from the project directory ---

# 1. pi in PATH (e.g. npm install -g @mariozechner/pi-coding-agent)
if command -v pi >/dev/null 2>&1; then
  cd "$SCAGENT_DIR"
  exec pi "$@"
fi

# 2. npx fallback (no global install needed — just Node.js)
if command -v npx >/dev/null 2>&1; then
  echo "Note: Running pi-coding-agent via npx (first run may be slow)..." >&2
  cd "$SCAGENT_DIR"
  exec npx --yes @mariozechner/pi-coding-agent "$@"
fi

echo "Error: Could not find pi or npx." >&2
echo "" >&2
echo "Install Node.js (≥20.19) from https://nodejs.org, then either:" >&2
echo "  npm install -g @mariozechner/pi-coding-agent   # global install" >&2
echo "  npx @mariozechner/pi-coding-agent              # no install needed" >&2
exit 1
