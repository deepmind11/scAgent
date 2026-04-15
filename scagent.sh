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

# Where the script (and package) lives — used for config fallback
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")" && pwd)"

# Resolve config root + working directory
# Priority: cwd project > script location (config only)
find_local_project() {
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
  return 1
}

WORK_DIR="$PWD"

if LOCAL="$(find_local_project)"; then
  # Running inside a scAgent project — use it for config + working dir
  CONFIG_ROOT="$LOCAL"
  WORK_DIR="$LOCAL"
elif [ -f "$SCRIPT_DIR/.pi/SYSTEM.md" ]; then
  # No local project — use script's config, stay in user's directory
  CONFIG_ROOT="$SCRIPT_DIR"
else
  echo "Error: Could not find scAgent project root (.pi/SYSTEM.md)." >&2
  exit 1
fi

if ! command -v feynman >/dev/null 2>&1; then
  echo "Error: Feynman is not installed." >&2
  echo "Install: curl -fsSL https://feynman.is/install | bash" >&2
  echo "Then run: feynman setup" >&2
  exit 1
fi

export FEYNMAN_CODING_AGENT_DIR="$CONFIG_ROOT/.pi"
cd "$WORK_DIR"
exec feynman "$@"
