"""scAgent CLI — launch the single-cell RNA-seq analysis agent.

Requires Feynman (https://github.com/getcompanion-ai/feynman).
Install: curl -fsSL https://feynman.is/install | bash
"""

import os
import shutil
import sys
from pathlib import Path


def _package_root() -> Path | None:
    """Return the scAgent project root next to the installed package, if valid."""
    package_dir = Path(__file__).resolve().parent  # scagent/
    project_dir = package_dir.parent               # scAgent/
    if (project_dir / ".pi" / "SYSTEM.md").exists():
        return project_dir
    return None


def _find_local_project() -> Path | None:
    """Walk up from cwd to find a directory with .pi/SYSTEM.md."""
    cwd = Path.cwd()
    for parent in [cwd, *cwd.parents]:
        if (parent / ".pi" / "SYSTEM.md").exists():
            return parent
    return None


def resolve_roots() -> tuple[Path, Path]:
    """Return (config_root, work_dir).

    - config_root: where .pi/ lives (system prompt, skills, tool schemas)
    - work_dir:    where feynman should run (where the user's data is)

    Priority:
      1. SCAGENT_ROOT env var         → config + work dir
      2. Local project (cwd ancestor) → config + work dir
      3. Package location             → config only, work dir stays as cwd
    """
    # 1. Explicit env var — full override
    env_root = os.environ.get("SCAGENT_ROOT")
    if env_root:
        root = Path(env_root)
        if (root / ".pi" / "SYSTEM.md").exists():
            return root, root

    # 2. Local project in cwd tree — use it for both
    local = _find_local_project()
    if local:
        return local, local

    # 3. No local project — use package config but stay in cwd
    pkg = _package_root()
    if pkg:
        return pkg, Path.cwd()

    print("Error: Could not find scAgent project root (.pi/SYSTEM.md).", file=sys.stderr)
    print("Either:", file=sys.stderr)
    print("  - Run from inside the scAgent directory", file=sys.stderr)
    print("  - Set SCAGENT_ROOT=/path/to/scAgent", file=sys.stderr)
    sys.exit(1)


def main():
    """Entry point for the scagent command."""
    config_root, work_dir = resolve_roots()

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

    # Point feynman at the scAgent config
    os.environ["FEYNMAN_CODING_AGENT_DIR"] = str(config_root / ".pi")

    # Launch feynman in the working directory
    args = [feynman] + sys.argv[1:]
    os.chdir(work_dir)
    os.execvp(args[0], args)


if __name__ == "__main__":
    main()
