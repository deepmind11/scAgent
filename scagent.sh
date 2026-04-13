#!/bin/bash
# Launch scAgent — single-cell RNA-seq analysis agent
#
# Usage:
#   scagent                              # default
#   scagent --model opus                 # Opus 4.6
#   scagent --model opus --thinking max  # max thinking
#   scagent -c                           # continue previous session
#   scagent -r                           # resume/pick a session

SCAGENT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# Don't cd — run from the user's current directory so .scagent/ state
# is created per-project.  PI_CODING_AGENT_DIR tells the agent where
# to find its system prompt, skills, and settings.
export PI_CODING_AGENT_DIR="$SCAGENT_DIR/.pi"

# Use Feynman's bundled pi (has working auth)
FEYNMAN_ROOT="$(ls -d "$HOME/.local/share/feynman/feynman-"*"-darwin-arm64" 2>/dev/null | sort -V | tail -1)"
if [ -z "$FEYNMAN_ROOT" ]; then
  echo "Error: Feynman installation not found" >&2
  exit 1
fi

exec "$FEYNMAN_ROOT/node/bin/node" \
  "$FEYNMAN_ROOT/app/.feynman/npm/node_modules/@mariozechner/pi-coding-agent/dist/cli.js" \
  "$@"
