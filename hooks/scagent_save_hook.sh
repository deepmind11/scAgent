#!/bin/bash
# scAgent Save Hook — auto-save analysis context every N exchanges
#
# Adapted from MemPalace's save hook for scAgent.
# Fires on the "Stop" event.  Every SAVE_INTERVAL human messages it
# blocks the AI and tells it to save analysis decisions, tool results,
# and rationale to the MemPalace palace — tagged with the current branch.
#
# === INSTALL ===
# Add to .claude/settings.local.json (or equivalent):
#
#   "hooks": {
#     "Stop": [{
#       "matcher": "*",
#       "hooks": [{
#         "type": "command",
#         "command": "/absolute/path/to/scAgent/hooks/scagent_save_hook.sh",
#         "timeout": 30
#       }]
#     }]
#   }

SAVE_INTERVAL=15
STATE_DIR="$HOME/.mempalace/hook_state"
mkdir -p "$STATE_DIR"

INPUT=$(cat)

eval $(echo "$INPUT" | python3 -c "
import sys, json, re
data = json.load(sys.stdin)
safe = lambda s: re.sub(r'[^a-zA-Z0-9_/.\-~]', '', str(s))
print(f'SESSION_ID=\"{safe(data.get(\"session_id\", \"unknown\"))}\"')
print(f'STOP_HOOK_ACTIVE=\"{data.get(\"stop_hook_active\", False)}\"')
print(f'TRANSCRIPT_PATH=\"{safe(data.get(\"transcript_path\", \"\"))}\"')
" 2>/dev/null)

TRANSCRIPT_PATH="${TRANSCRIPT_PATH/#\~/$HOME}"

# Prevent infinite loop
if [ "$STOP_HOOK_ACTIVE" = "True" ] || [ "$STOP_HOOK_ACTIVE" = "true" ]; then
    echo "{}"
    exit 0
fi

# Count human messages
if [ -f "$TRANSCRIPT_PATH" ]; then
    EXCHANGE_COUNT=$(python3 - "$TRANSCRIPT_PATH" <<'PYEOF'
import json, sys
count = 0
with open(sys.argv[1]) as f:
    for line in f:
        try:
            entry = json.loads(line)
            msg = entry.get('message', {})
            if isinstance(msg, dict) and msg.get('role') == 'user':
                content = msg.get('content', '')
                if isinstance(content, str) and '<command-message>' in content:
                    continue
                count += 1
        except:
            pass
print(count)
PYEOF
2>/dev/null)
else
    EXCHANGE_COUNT=0
fi

LAST_SAVE_FILE="$STATE_DIR/${SESSION_ID}_last_save"
LAST_SAVE=0
[ -f "$LAST_SAVE_FILE" ] && LAST_SAVE=$(cat "$LAST_SAVE_FILE")
SINCE_LAST=$((EXCHANGE_COUNT - LAST_SAVE))

echo "[$(date '+%H:%M:%S')] scAgent session $SESSION_ID: $EXCHANGE_COUNT exchanges, $SINCE_LAST since last save" >> "$STATE_DIR/hook.log"

if [ "$SINCE_LAST" -ge "$SAVE_INTERVAL" ] && [ "$EXCHANGE_COUNT" -gt 0 ]; then
    echo "$EXCHANGE_COUNT" > "$LAST_SAVE_FILE"
    echo "[$(date '+%H:%M:%S')] TRIGGERING scAgent SAVE at exchange $EXCHANGE_COUNT" >> "$STATE_DIR/hook.log"

    cat << 'HOOKJSON'
{
  "decision": "block",
  "reason": "scAgent AUTO-SAVE checkpoint. Save to MemPalace using the memory skill:\n1. Key analysis decisions and their rationale\n2. Tool results and parameter choices\n3. Important discussion points and rejected alternatives\n\nTag every entry with the current analysis branch. Use store_decision() for decisions, store_step() for tool results, store_exchange() for discussion. Then continue."
}
HOOKJSON
else
    echo "{}"
fi
