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

## Documentation

Doxygen `/** */` on every public declaration; the first sentence is the brief (no `@brief`);
`@param`, `@return` and `@tparam` throughout. `///<` for data members. State units (SI) and sign
conventions on every physical quantity.

## Physics and numerics

- SI units internally. Signs, axes and frames are stated where a quantity is declared.
- Errors are values: `std::expected<T, Fault>`. No exceptions, no `assert` on physical input.
- Stepping paths never allocate; scratch belongs to the integrator or the caller.
- Every model is checked against something independent: a closed form, a finer step, a brute-force
  search or a conserved quantity.
- Fault details are string literals so faults are free to create and copy.
