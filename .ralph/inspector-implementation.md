# Implement AnnData Inspector + Dependencies + DAG Refactor

Following the plan in `notes/inspector-plan.md`.

## Goals
1. New `scagent/inspector.py` with `AnnDataState`, `inspect_adata()`, `find_raw_counts()`, `summarize_state()`
2. New `scagent/dependencies.py` with `DEPENDENCIES`, `check_prerequisites()`, `plan_steps()`, `ensure_ready_for()`
3. Refactor `scagent/dag.py` — add `add_step()`, `mark_precomputed()`, make advisory
4. Refactor `scagent/context.py` — make paradigm optional
5. Fix eval handlers to use inspector
6. All tests passing
7. Update architecture doc
8. Update skills
9. Update `__init__.py` exports
10. Run scBench and verify HVG no longer crashes

## Checklist
- [x] Phase 1: `scagent/inspector.py` — AnnDataState, inspect_adata(), find_raw_counts(), summarize_state()
- [x] Phase 1: `tests/test_inspector.py` — 39 unit tests, all pass
- [x] Phase 1: `tests/test_inspector_integration.py` — 5 integration tests against real scBench data, all pass
- [x] Phase 2: `scagent/dependencies.py` — DEPENDENCIES, check_prerequisites(), plan_steps(), ensure_ready_for()
- [x] Phase 2: `tests/test_dependencies.py` — 21 unit tests, all pass
- [x] Phase 3: `scagent/dag.py` refactor — add_step(), mark_precomputed(), mark_precomputed_from_state()
- [x] Phase 3: `scagent/context.py` refactor — paradigm now optional
- [x] Phase 3: `tests/test_dag.py` updates — 8 new tests for add_step/mark_precomputed
- [x] Phase 3: `.pi/skills/init/SKILL.md` update — Mode B onboarding
- [x] Phase 3: `.pi/skills/dag/SKILL.md` update — advisory behavior, dependency graph
- [x] Phase 4: `eval/handlers/*.py` fixes — all 7 handlers use inspector
- [x] Phase 4: Run scBench — **6/7 (86%)**, up from 4/7 (57%)
- [x] Phase 5: `outputs/architecture.md` update — new pillars, inspector section, scBench table
- [x] Phase 5: `scagent/__init__.py` exports

## Results
- scBench: 4/7 (57%) → **6/7 (86%)**
- HVG eval: CRASH → PASS (100% recall)
- Cell typing: FAIL (18pt off) → PASS (1.8pt off)
- Total tests: 302 passed across all modules
