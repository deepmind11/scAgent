#!/bin/bash
# Launch scAgent — single-cell RNA-seq analysis agent
#
# Requires Feynman: curl -fsSL https://feynman.is/install | bash
#
# Usage:
#   ./scagent.sh                   # launch
#   ./scagent.sh --model opus      # specific model
#   ./scagent.sh --thinking max    # max thinking
#   ./scagent.sh -c                # continue previous session
#   ./scagent.sh -r                # resume/pick a session

set -euo pipefail

# Where the shell script (and package) lives — fallback only
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# Resolve project root: prefer cwd (or ancestor), fall back to script dir
find_root() {
  # Walk up from cwd
  local dir="$PWD"
  while true; do
    if [ -f "$dir/.pi/SYSTEM.md" ]; then
      echo "$dir"
      return
    fi
    local parent="$(dirname "$dir")"
    [ "$parent" = "$dir" ] && break
    dir="$parent"
  done
  # Fallback: script's own directory
  if [ -f "$SCRIPT_DIR/.pi/SYSTEM.md" ]; then
    echo "$SCRIPT_DIR"
    return
  fi
  echo "Error: Could not find scAgent project root (.pi/SYSTEM.md)." >&2
  exit 1
}

SCAGENT_DIR="$(find_root)"

if ! command -v feynman >/dev/null 2>&1; then
  echo "Error: Feynman is not installed." >&2
  echo "Install: curl -fsSL https://feynman.is/install | bash" >&2
  echo "Then run: feynman setup" >&2
  exit 1
fi

export FEYNMAN_CODING_AGENT_DIR="$SCAGENT_DIR/.pi"
cd "$SCAGENT_DIR"
exec feynman "$@"
