# scAgent Hooks

Auto-save analysis context to MemPalace.  Adapted from
[MemPalace hooks](https://github.com/MemPalace/mempalace/tree/main/hooks).

## Install

Add to `.claude/settings.local.json`:

```json
{
  "hooks": {
    "Stop": [{
      "matcher": "*",
      "hooks": [{
        "type": "command",
        "command": "/absolute/path/to/scAgent/hooks/scagent_save_hook.sh",
        "timeout": 30
      }]
    }],
    "PreCompact": [{
      "hooks": [{
        "type": "command",
        "command": "/absolute/path/to/scAgent/hooks/scagent_precompact_hook.sh",
        "timeout": 30
      }]
    }]
  }
}
```

Make executable:

```bash
chmod +x hooks/scagent_save_hook.sh hooks/scagent_precompact_hook.sh
```

## What they do

| Hook | Fires | Action |
|------|-------|--------|
| **Save** | Every 15 human messages | Blocks agent → tells it to save decisions, steps, exchanges to MemPalace tagged with current branch |
| **PreCompact** | Before context compaction | Emergency save → save everything before context is lost |

The agent does the actual saving — the hooks just tell it **when**.
