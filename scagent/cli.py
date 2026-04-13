"""scAgent CLI — launch the single-cell RNA-seq analysis agent.

This finds the scAgent project root (where .pi/SYSTEM.md lives),
then launches pi-coding-agent from that directory so it auto-discovers
the system prompt, skills, and settings.

Auth is stored in ~/.pi/agent/auth.json (managed by /login in the REPL).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def find_scagent_root() -> Path:
    """Find the scAgent project root by locating .pi/SYSTEM.md."""
    # 1. Explicit env var
    env_root = os.environ.get("SCAGENT_ROOT")
    if env_root:
        root = Path(env_root)
        if (root / ".pi" / "SYSTEM.md").exists():
            return root

    # 2. Installed package location (scagent/ is inside the project)
    package_dir = Path(__file__).resolve().parent        # scagent/
    project_dir = package_dir.parent                     # scAgent/
    if (project_dir / ".pi" / "SYSTEM.md").exists():
        return project_dir

    # 3. Walk up from cwd
    cwd = Path.cwd()
    for parent in [cwd, *cwd.parents]:
        if (parent / ".pi" / "SYSTEM.md").exists():
            return parent

    print("Error: Could not find scAgent project root (.pi/SYSTEM.md).", file=sys.stderr)
    print("Either:", file=sys.stderr)
    print("  - Run from inside the scAgent directory", file=sys.stderr)
    print("  - Set SCAGENT_ROOT=/path/to/scAgent", file=sys.stderr)
    sys.exit(1)


def find_pi_cli() -> list[str]:
    """Find the pi CLI. Returns the command as a list."""
    # 1. pi in PATH
    pi_path = shutil.which("pi")
    if pi_path:
        return [pi_path]

    # 2. npx fallback
    npx_path = shutil.which("npx")
    if npx_path:
        return [npx_path, "--yes", "@mariozechner/pi-coding-agent"]

    print("Error: Could not find pi or npx.", file=sys.stderr)
    print("", file=sys.stderr)
    print("Install Node.js (≥20.19) from https://nodejs.org, then either:", file=sys.stderr)
    print("  npm install -g @mariozechner/pi-coding-agent   # global install", file=sys.stderr)
    print("  npx @mariozechner/pi-coding-agent              # no install needed", file=sys.stderr)
    sys.exit(1)


def main():
    """Entry point for the scagent command."""
    root = find_scagent_root()
    cmd = find_pi_cli()

    # Pass through all CLI args to pi
    args = cmd + sys.argv[1:]

    # Launch pi from the scAgent project root so it auto-discovers
    # .pi/SYSTEM.md, .pi/skills/, and .pi/settings.json
    os.chdir(root)
    os.execvp(args[0], args)


if __name__ == "__main__":
    main()
