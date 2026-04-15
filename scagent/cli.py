"""scAgent CLI — launch the single-cell RNA-seq analysis agent.

Requires Feynman (https://github.com/getcompanion-ai/feynman).
Install: curl -fsSL https://feynman.is/install | bash
"""

import os
import shutil
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

    # 2. Walk up from cwd (preferred — lets you run from any scAgent project)
    cwd = Path.cwd()
    for parent in [cwd, *cwd.parents]:
        if (parent / ".pi" / "SYSTEM.md").exists():
            return parent

    # 3. Installed package location (fallback — editable installs point here)
    package_dir = Path(__file__).resolve().parent
    project_dir = package_dir.parent
    if (project_dir / ".pi" / "SYSTEM.md").exists():
        return project_dir

    print("Error: Could not find scAgent project root (.pi/SYSTEM.md).", file=sys.stderr)
    print("Either:", file=sys.stderr)
    print("  - Run from inside the scAgent directory", file=sys.stderr)
    print("  - Set SCAGENT_ROOT=/path/to/scAgent", file=sys.stderr)
    sys.exit(1)


def main():
    """Entry point for the scagent command."""
    root = find_scagent_root()

    # Find feynman
    feynman = shutil.which("feynman")
    if not feynman:
        print("Error: Feynman is not installed.", file=sys.stderr)
        print("", file=sys.stderr)
        print("Install it with:", file=sys.stderr)
        print("  curl -fsSL https://feynman.is/install | bash", file=sys.stderr)
        print("", file=sys.stderr)
        print("Then run: feynman setup", file=sys.stderr)
        sys.exit(1)

    # Set the agent config to scAgent's .pi/ directory
    os.environ["FEYNMAN_CODING_AGENT_DIR"] = str(root / ".pi")

    # Launch feynman from the project root
    args = [feynman] + sys.argv[1:]
    os.chdir(root)
    os.execvp(args[0], args)


if __name__ == "__main__":
    main()
