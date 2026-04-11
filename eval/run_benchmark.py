#!/usr/bin/env python3
"""Run scAgent against all 7 Chromium canonical scBench evals.

Usage:
    python eval/run_benchmark.py [--keep-workspace]
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scbench import EvalRunner
from eval.adapter import scagent_direct_agent

EVAL_DIR = Path(__file__).parent.parent.parent / "scbench-eval" / "evals_canonical" / "chromium"


def main():
    keep_workspace = "--keep-workspace" in sys.argv

    eval_files = sorted(EVAL_DIR.glob("*.json"))
    if not eval_files:
        print(f"No eval files found in {EVAL_DIR}")
        sys.exit(1)

    print("=" * 80)
    print(f"scAgent Benchmark — {len(eval_files)} Chromium evals")
    print(f"Adapter: direct (no LLM)")
    print(f"Time: {datetime.now().isoformat()}")
    print("=" * 80)

    results = []

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
            result = runner.run(agent_function=scagent_direct_agent)
            duration = time.time() - start

            passed = result.get("passed", False)
            grader_result = result.get("grader_result")

            status = "✓ PASSED" if passed else "✗ FAILED"
            print(f"\n  {status} ({duration:.1f}s)")

            results.append({
                "eval": eval_name,
                "passed": passed,
                "duration_s": round(duration, 1),
                "answer": result.get("agent_answer"),
                "grader": _extract_grader_summary(grader_result),
            })

        except Exception as e:
            duration = time.time() - start
            print(f"\n  ✗ ERROR: {e}")
            traceback.print_exc()
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
        extra = ""
        if "grader" in r and r["grader"]:
            extra = f" — {r['grader']}"
        if "error" in r:
            extra = f" — ERROR: {r['error'][:80]}"
        print(f"  {status} {r['eval']}{extra}")

    print(f"\nScore: {n_passed}/{n_total} ({n_passed/n_total*100:.0f}%)")
    total_time = sum(r["duration_s"] for r in results)
    print(f"Total time: {total_time:.0f}s")

    # Save results
    output = {
        "timestamp": datetime.now().isoformat(),
        "adapter": "direct",
        "n_passed": n_passed,
        "n_total": n_total,
        "score_pct": round(n_passed / n_total * 100, 1),
        "total_duration_s": round(total_time, 1),
        "results": results,
    }

    output_dir = Path(__file__).parent.parent / "eval" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    output_file.write_text(json.dumps(output, indent=2, default=str))
    print(f"\nResults saved to: {output_file}")


def _extract_grader_summary(grader_result) -> str | None:
    if grader_result is None:
        return None
    if hasattr(grader_result, "reasoning"):
        return grader_result.reasoning[:150]
    if isinstance(grader_result, dict):
        return str(grader_result.get("reasoning", ""))[:150]
    return str(grader_result)[:150]


if __name__ == "__main__":
    main()
