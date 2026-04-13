#!/bin/bash
# scAgent PreCompact Hook — emergency save before context compaction
#
# Fires RIGHT BEFORE conversation gets compressed.  Always blocks.
# Tells the agent to save EVERYTHING before losing detailed context.
#
# === INSTALL ===
# Add to .claude/settings.local.json (or equivalent):
#
#   "hooks": {
#     "PreCompact": [{
#       "hooks": [{
#         "type": "command",
#         "command": "/absolute/path/to/scAgent/hooks/scagent_precompact_hook.sh",
#         "timeout": 30
#       }]
#     }]
#   }

STATE_DIR="$HOME/.mempalace/hook_state"
mkdir -p "$STATE_DIR"

INPUT=$(cat)
SESSION_ID=$(echo "$INPUT" | python3 -c "import sys,json; print(json.load(sys.stdin).get('session_id','unknown'))" 2>/dev/null)

echo "[$(date '+%H:%M:%S')] scAgent PRE-COMPACT triggered for session $SESSION_ID" >> "$STATE_DIR/hook.log"

cat << 'HOOKJSON'
{
  "decision": "block",
  "reason": "COMPACTION IMMINENT — scAgent emergency save. Save ALL analysis context to MemPalace before it is lost:\n1. Every analysis decision and its rationale\n2. Every tool run and its key results\n3. Current DAG position and what steps remain\n4. Current branch name and experiment context\n5. Any rejected alternatives and why\n\nTag every entry with the current analysis branch. Be thorough — after compaction, detailed context will be lost. Then allow compaction to proceed."
}
HOOKJSON
