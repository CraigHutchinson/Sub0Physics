# Repository guidance

Sub0Physics is a product-agnostic C++23 physics library in the Sub0 family. Use generic examples
(oscillators, wheels, contacts, profiles) in source, tests, documentation, issues and pull requests.
Do not introduce consumer product names, consumer-specific models or consumer dependencies.

- Read [STYLE_GUIDE.md](STYLE_GUIDE.md) (profile `sub0`) before writing C++.
- Every primitive is verified against an independent reference in `tests/`.
- Ported code updates [upstream/bikelab.json](upstream/bikelab.json) in the same commit
  ([docs/UPSTREAM.md](docs/UPSTREAM.md)).
- Performance claims come from `Sub0Physics_Bench` with the host and compiler stated.

## C++ Skills

Three skills manage the C++ authoring and review workflow. Full reference material lives in `~/.claude/skills/cpp/references/`.

| When | Skill | Invoke condition |
|---|---|---|
| Before writing a new class, public API, or module | **cpp-write** | Proactively before authoring |
| After substantive C++ changes | **cpp-review** | Proactively after authoring |
| Before implementing -- to review a plan or design doc | **cpp-review** (plan mode) | Explicitly, with plan doc as input |
| At a phase boundary, or after parallel work across modules or repositories | **cpp-architecture-review** | Proactively: approach and algorithms, algorithmic redundancy, one solution per problem across repositories, DRY, contracts the next milestone builds on; consolidated decisions and plan go in `docs/reviews/` of the most private repository whose internals it describes |
| After L0/L1 MUST findings | **cpp-simplify** | With Rewrite Brief from cpp-review as input |

The `.cursor/rules/cpp-standards.mdc` project rule auto-attaches the most important constraints to Agent sessions when C++ files are in context. The skills carry the full reference material and worked examples.
