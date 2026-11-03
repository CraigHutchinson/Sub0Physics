# Sub0Physics Project Rules

Read [AGENTS.md](AGENTS.md) first.

## Build & Test

```bash
cmake --preset default          # Configure (Release + tests + benchmarks)
cmake --build --preset default  # Build
ctest --preset default          # Run tests
```

`build/tests/Sub0Physics_Bench` runs the benchmarks; they are not part of ctest.

The pre-push hook is the performance gate that configuring installs into `.git/hooks` (sub0perf, docs/PERFORMANCE.md); no hook is committed and `core.hooksPath` stays unset. Tests run in local CI (`local_ci.py --level quick`), and `hook_check.py --repo .` from the .github repository's `tools/ci` proves the gate works.

## Commit Rules

- A commit that changes the public API in `include/sub0physics/` says so in its message.
- Tests pass locally before committing.
- Ported code and its `upstream/bikelab.json` record change together.

## Branch Strategy

`main` is the integration branch. Use focused topic branches and pull requests.

## Repository Hygiene

Do not commit handoff notes, scratch documents, probes or build output.
