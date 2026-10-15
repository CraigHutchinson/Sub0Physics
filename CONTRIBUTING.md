# Contributing

Read [AGENTS.md](AGENTS.md). Sub0Physics stays product-agnostic: generic primitives only.

1. **State the need.** Name the real consumer and the observable problem. New primitives arrive with
   that consumer, not ahead of it.
2. **Verify independently.** Each primitive is tested against a closed form, a convergence study, a
   brute-force check or a conservation law, not only against its own output.
3. **Keep the hot path clean.** No allocation, exceptions or I/O while stepping. Faults are values.
4. **Record provenance.** Code ported from elsewhere updates `upstream/*.json`
   ([docs/UPSTREAM.md](docs/UPSTREAM.md)).

By contributing you agree that your contribution is licensed under the [MIT licence](LICENSE.md).
