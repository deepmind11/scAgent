"""scAgent CLI — launch the single-cell RNA-seq analysis agent.

This finds the scAgent project root (where .pi/SYSTEM.md lives),
ensures the agent runtime (pi-coding-agent) is available, then
launches it from the project directory so it auto-discovers the
system prompt, skills, and settings.

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


def check_node() -> bool:
    """Check if Node.js ≥20.19 is available."""
    node = shutil.which("node")
    if not node:
        return False
    try:
        result = subprocess.run(
            [node, "--version"], capture_output=True, text=True, timeout=5
        )
        version = result.stdout.strip().lstrip("v")
        parts = [int(x) for x in version.split(".")[:2]]
        return parts[0] > 20 or (parts[0] == 20 and parts[1] >= 19)
    except Exception:
        return False


def install_pi_coding_agent() -> bool:
    """Attempt to install pi-coding-agent globally via npm."""
    npm = shutil.which("npm")
    if not npm:
        return False
    print("Installing pi-coding-agent...", file=sys.stderr)
    try:
        result = subprocess.run(
            [npm, "install", "-g", "@mariozechner/pi-coding-agent"],
            timeout=120,
        )
        return result.returncode == 0
    except Exception:
        return False


def ensure_auth():
    """Ensure pi-coding-agent has auth credentials.

    If ~/.pi/agent/auth.json is empty/missing but Feynman has working
    OAuth tokens at ~/.feynman/agent/auth.json, symlink them so scAgent
    works immediately without a separate /login step.
    """
    home = Path.home()
    pi_auth_dir = home / ".pi" / "agent"
    pi_auth = pi_auth_dir / "auth.json"
    feynman_auth = home / ".feynman" / "agent" / "auth.json"

    # Already has credentials?
    if pi_auth.exists() and not pi_auth.is_symlink():
        try:
            import json
            data = json.loads(pi_auth.read_text())
            if data:  # non-empty dict
                return
        except Exception:
            pass

    # Feynman has working auth? Symlink it.
    if feynman_auth.exists():
        try:
            import json
            data = json.loads(feynman_auth.read_text())
            if data:
                pi_auth_dir.mkdir(parents=True, exist_ok=True)
                if pi_auth.exists() or pi_auth.is_symlink():
                    pi_auth.unlink()
                pi_auth.symlink_to(feynman_auth)
                print(f"Linked auth from Feynman ({feynman_auth})", file=sys.stderr)
                return
        except Exception:
            pass


def find_pi_cli() -> list[str]:
    """Find the pi CLI, installing it if necessary. Returns the command."""
    # 1. pi in PATH
    pi_path = shutil.which("pi")
    if pi_path:
        return [pi_path]

    # 2. Node.js available? Try to install pi-coding-agent
    if check_node():
        npx_path = shutil.which("npx")

        # Try global install
        if install_pi_coding_agent():
            pi_path = shutil.which("pi")
            if pi_path:
                return [pi_path]

        # Fall back to npx
        if npx_path:
            print("Using npx to run pi-coding-agent...", file=sys.stderr)
            return [npx_path, "--yes", "@mariozechner/pi-coding-agent"]

    # 3. Nothing works
    print("Error: Node.js (≥20.19) is required but not found.", file=sys.stderr)
    print("", file=sys.stderr)
    print("Install Node.js from https://nodejs.org, then rerun scagent.", file=sys.stderr)
    print("The agent runtime will be installed automatically.", file=sys.stderr)
    sys.exit(1)


def main():
    """Entry point for the scagent command."""
    root = find_scagent_root()
    ensure_auth()
    cmd = find_pi_cli()

    # Pass through all CLI args to pi
    args = cmd + sys.argv[1:]

    # Launch pi from the scAgent project root so it auto-discovers
    # .pi/SYSTEM.md, .pi/skills/, and .pi/settings.json
    os.chdir(root)
    os.execvp(args[0], args)


if __name__ == "__main__":
    main()
