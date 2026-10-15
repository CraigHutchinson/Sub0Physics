# Sub0Physics Project Rules

Read [AGENTS.md](AGENTS.md) first.

## Build & Test

```bash
cmake --preset default          # Configure (Release + tests + benchmarks)
cmake --build --preset default  # Build
ctest --preset default          # Run tests
```

`build/tests/Sub0Physics_Bench` runs the benchmarks; they are not part of ctest.

A git pre-push hook runs the tests: `git config core.hooksPath .githooks`.

## Commit Rules

- A commit that changes the public API in `include/sub0physics/` says so in its message.
- Tests pass locally before committing.
- Ported code and its `upstream/bikelab.json` record change together.

## Branch Strategy

`main` is the integration branch. Use focused topic branches and pull requests.

## Repository Hygiene

Do not commit handoff notes, scratch documents, probes or build output.
