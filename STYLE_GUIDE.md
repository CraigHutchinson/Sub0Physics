# Sub0Physics Code Style Guide

style-profile: sub0

Sub0Physics follows the Sub0 family style profile. This guide summarises how the rules apply here;
where it and the profile disagree, the profile wins.

## Naming

| Element | Convention | Example |
|---|---|---|
| Namespace | lowercase | `sub0physics` |
| Classes, structs, enums, aliases, concepts | PascalCase | `PolylineProfile`, `Rk4`, `OdeSystem` |
| Functions and methods | camelCase | `measureCircleClearance`, `respond` |
| Members | camelCase with `_` suffix | `cacheValid_` (public aggregate fields have none) |
| Constants and enumerators | PascalCase, no prefix | `StandardGravity`, `FaultKind::OutsideDomain` |
| Macros and CMake options | `SUB0PHYSICS_` prefix | `SUB0PHYSICS_BUILD_TESTING` |
| Files | snake_case `.hpp` / `.cpp` | `circle_clearance.hpp`, `test_integrators.cpp` |

## Layout

4 spaces, no tabs; Allman braces for namespaces, types, functions and control flow. About 120
columns. `template<typename T>` with no inner spaces.

## Includes

System and third-party headers first (`<doctest.h>` included), then own headers quoted and
library-rooted (`"sub0physics/fault.hpp"`); each group alphabetical. Every folder has an umbrella
header (`element/` → `elements.hpp`); `sub0physics.hpp` includes everything.
`Sub0Physics_HeaderCheck` compiles every public header on its own.

## Performance notes in code

Where performance has been measured and iterated, leave a `// PERF:` comment at the code, ten words or
fewer, saying what was found faster than what. Numbers and evidence go in `docs/PERFORMANCE.md` under
Findings. PERF notes are hints that record the area has been worked on; later work may override one
when its own measurement says so, and then updates the note. The family convention is in Sub0Math's
STYLE_GUIDE.md.

## Documentation

Doxygen `/** */` on every public declaration; the first sentence is the brief (no `@brief`);
`@param`, `@return` and `@tparam` throughout. `///<` for data members. State units (SI) and sign
conventions on every physical quantity.

## Physics and numerics

- SI units internally. Signs, axes and frames are stated where a quantity is declared.
- Errors are values: `std::expected<T, Fault>`. No exceptions, no `assert` on physical input.
- Hot-path results are status aggregates (`{value, known}` or a struct with a `status` field), not
  `std::optional` or `std::expected`: measured slower, because the optional keeps the result in memory.
  `std::expected` is for setup and anywhere else off the hot path. A query that returns a variable number
  of results (witnesses, cells) writes them to a caller-owned fixed-capacity buffer and reports the count
  and an `overflowed` flag; it never allocates.
- Stepping paths never allocate; scratch belongs to the integrator or the caller.
- Every model is checked against something independent: a closed form, a finer step, a brute-force
  search or a conserved quantity.
- Fault details are string literals so faults are free to create and copy.

## Portability across the CI compilers

- **Don't iterate a view of a temporary in a range-for** (`for (x : make().view())`). C++23 (P2718) extends the
  temporary's lifetime, but GCC 14 in CI does not implement that yet, so the view dangles there. Name the
  result first: `const auto tiles = make(); for (x : tiles.view())`. This was caught by ASan
  (stack-use-after-scope) on 2026-10-21; GCC 16 passed locally.
