#!/bin/bash
# Launch scAgent — single-cell RNA-seq analysis agent
#
# Prerequisites:
#   - Feynman CLI (https://github.com/mariozechner/feynman)
#     Install: curl -fsSL https://feynman.sh/install | bash
#   - Python ≥3.11 with: pip install -e .
#
# Usage:
#   ./scagent.sh                              # default
#   ./scagent.sh --model opus                 # Opus model
#   ./scagent.sh --model opus --thinking max  # max thinking
#   ./scagent.sh -c                           # continue previous session
#   ./scagent.sh -r                           # resume/pick a session

set -euo pipefail

SCAGENT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# Point the agent at scAgent's config (system prompt, skills, tools)
export PI_CODING_AGENT_DIR="$SCAGENT_DIR/.pi"

# --- Find the pi CLI ---

# 1. Feynman bundles pi — prefer the bundled version
FEYNMAN_ROOT="$(ls -d "$HOME/.local/share/feynman/feynman-"*"-$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)" 2>/dev/null | sort -V | tail -1)"
if [ -n "$FEYNMAN_ROOT" ] && [ -x "$FEYNMAN_ROOT/node/bin/node" ]; then
  PI_CLI="$FEYNMAN_ROOT/app/.feynman/npm/node_modules/@mariozechner/pi-coding-agent/dist/cli.js"
  if [ -f "$PI_CLI" ]; then
    exec "$FEYNMAN_ROOT/node/bin/node" "$PI_CLI" "$@"
  fi
fi

# 2. pi in PATH (standalone install)
if command -v pi >/dev/null 2>&1; then
  exec pi "$@"
fi

# 3. npx fallback
if command -v npx >/dev/null 2>&1; then
  echo "Note: Installing pi-coding-agent via npx (first run may be slow)..." >&2
  exec npx --yes @mariozechner/pi-coding-agent "$@"
fi

echo "Error: Could not find Feynman, pi, or npx." >&2
echo "Install Feynman: curl -fsSL https://feynman.sh/install | bash" >&2
echo "Or install Node.js and run: npx @mariozechner/pi-coding-agent" >&2
exit 1
