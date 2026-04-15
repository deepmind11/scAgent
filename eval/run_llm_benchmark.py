#!/usr/bin/env python3
"""Run scAgent as a full LLM agent against SC-Bench canonical evals.

Unlike run_benchmark.py (which calls tools directly), this launches the
actual agent via Feynman in non-interactive mode. The agent must interpret
each task prompt, choose the right tools, and produce a correct answer.

Usage:
    python eval/run_llm_benchmark.py [--model MODEL] [--keep-workspace]

Requires:
    - Feynman installed and authenticated (feynman setup)
    - pip install -e ".[eval]"
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from latch_eval_tools import EvalRunner

EVAL_DIR = Path(__file__).parent / "evals_canonical_chromium"
SCAGENT_ROOT = Path(__file__).resolve().parent.parent


def build_agent_prompt(task_prompt: str, data_file: Path, answer_file: Path) -> str:
    """Build the prompt sent to the agent."""
    return f"""You are running a SC-Bench evaluation task. Solve this single-cell RNA-seq analysis task.

DATA FILE: {data_file}
Load it with: import anndata; adata = anndata.read_h5ad('{data_file}')

Always inspect the data first: check adata.shape, adata.obs.columns, adata.layers.keys(), adata.X.max() to understand what you're working with.

IMPORTANT: The Python venv is already activated. You can use scanpy, scagent.tools, or any standard analysis library.

TASK:
{task_prompt}

After completing the analysis, write ONLY the answer JSON to this file:
{answer_file}

Do NOT include any text or explanation in the answer file — just the raw JSON object."""


def _load_auth_env() -> dict[str, str]:
    """Load Anthropic OAuth token from .pi/auth.json into env vars.

    pi -p in subprocess mode doesn't inherit Feynman's OAuth session,
    so we read the access token and pass it as ANTHROPIC_API_KEY.
    Falls back to existing env var if auth.json is missing.
    """
    env = os.environ.copy()
    if env.get("ANTHROPIC_API_KEY"):
        return env  # already set

    auth_file = SCAGENT_ROOT / ".pi" / "auth.json"
    if auth_file.exists():
        try:
            import json as _json
            auth = _json.loads(auth_file.read_text())
            token = auth.get("anthropic", {}).get("access")
            if token:
                env["ANTHROPIC_API_KEY"] = token
        except Exception:
            pass

    if not env.get("ANTHROPIC_API_KEY"):
        print(
            "Error: No Anthropic API key found.\n"
            "Either set ANTHROPIC_API_KEY or run 'feynman setup' to authenticate.",
            file=sys.stderr,
        )
        sys.exit(1)

    return env


def run_agent(prompt: str, model: str) -> tuple[bool, str]:
    """Run the agent via pi in non-interactive mode."""
    pi_bin = shutil.which("pi")
    if not pi_bin:
        print("Error: pi not found. Install Feynman: curl -fsSL https://feynman.is/install | bash", file=sys.stderr)
        sys.exit(1)

    env = _load_auth_env()

    # Write prompt to temp file to avoid shell escaping issues
    prompt_file = Path(tempfile.mktemp(suffix=".md", prefix="eval_prompt_"))
    prompt_file.write_text(prompt, encoding="utf-8")

    try:
        result = subprocess.run(
            [pi_bin, "-p", "--model", model, f"@{prompt_file}"],
            cwd=str(SCAGENT_ROOT),
            capture_output=True,
            text=True,
            timeout=300,
            env=env,
        )
    finally:
        prompt_file.unlink(missing_ok=True)

    log = result.stdout + "\n" + result.stderr
    return result.returncode == 0, log


def llm_agent_function(task_prompt: str, work_dir: Path, model: str = "claude-opus-4-6") -> dict:
    """Agent function for EvalRunner that uses the full LLM agent."""
    # Find data file
    data_files = list(work_dir.glob("**/*.h5ad")) + list(work_dir.glob("**/*.node"))
    if not data_files:
        raise FileNotFoundError(f"No data files in {work_dir}")

    answer_file = work_dir / "eval_answer.json"
    prompt = build_agent_prompt(task_prompt, data_files[0], answer_file)

    print(f"  Running agent ({model})...")
    success, log = run_agent(prompt, model)

    # Save agent log
    log_file = work_dir / "agent.log"
    log_file.write_text(log)

    if not answer_file.exists():
        raise RuntimeError(f"Agent did not produce answer file. Log: {log_file}")

    answer = json.loads(answer_file.read_text())
    return {"answer": answer}


def main():
    model = "claude-opus-4-6"
    keep_workspace = False

    args = sys.argv[1:]
    while args:
        if args[0] == "--model" and len(args) > 1:
            model = args[1]
            args = args[2:]
        elif args[0] == "--keep-workspace":
            keep_workspace = True
            args = args[1:]
        else:
            args = args[1:]

    eval_files = sorted(EVAL_DIR.glob("*.json"))
    if not eval_files:
        print(f"No eval files found in {EVAL_DIR}")
        sys.exit(1)

    print("=" * 80)
    print(f"scAgent LLM Benchmark — {len(eval_files)} Chromium evals")
    print(f"Model: {model}")
    print(f"Time: {datetime.now().isoformat()}")
    print("=" * 80)

    results = []
    results_dir = Path(__file__).parent / "results" / f"llm_benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    results_dir.mkdir(parents=True, exist_ok=True)

    for i, eval_file in enumerate(eval_files, 1):
        eval_name = eval_file.stem
        print(f"\n[{i}/{len(eval_files)}] {eval_name}")
        print("-" * 60)

        start = time.time()
        try:
            runner = EvalRunner(
                eval_file,
                keep_workspace=keep_workspace,
                cache_name=".eval_cache",
            )
            result = runner.run(
                agent_function=lambda prompt, work_dir: llm_agent_function(prompt, work_dir, model)
            )
            duration = time.time() - start

            passed = result.get("passed", False)
            status = "✓ PASSED" if passed else "✗ FAILED"
            print(f"\n  {status} ({duration:.0f}s)")

            results.append({
                "eval": eval_name,
                "passed": passed,
                "duration_s": round(duration, 1),
                "answer": result.get("agent_answer"),
            })

        except Exception as e:
            duration = time.time() - start
            print(f"\n  ✗ ERROR: {e}")
            results.append({
                "eval": eval_name,
                "passed": False,
                "duration_s": round(duration, 1),
                "error": str(e),
            })

    # Summary
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)

    n_passed = sum(1 for r in results if r["passed"])
    n_total = len(results)

    for r in results:
        status = "✓" if r["passed"] else "✗"
        extra = f" — ERROR: {r['error'][:80]}" if "error" in r else ""
        print(f"  {status} {r['eval']} ({r['duration_s']:.0f}s){extra}")

    print(f"\nScore: {n_passed}/{n_total} ({n_passed/n_total*100:.0f}%)")
    total_time = sum(r["duration_s"] for r in results)
    print(f"Total time: {total_time:.0f}s")

    # Save results
    output = {
        "timestamp": datetime.now().isoformat(),
        "adapter": "llm_agent",
        "model": model,
        "n_passed": n_passed,
        "n_total": n_total,
        "score_pct": round(n_passed / n_total * 100, 1),
        "total_duration_s": round(total_time, 1),
        "results": results,
    }

    output_file = results_dir / "summary.json"
    output_file.write_text(json.dumps(output, indent=2, default=str))
    print(f"\nResults saved to: {results_dir}")


if __name__ == "__main__":
    main()
