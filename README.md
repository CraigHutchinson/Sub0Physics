# Sub0Physics

Product-agnostic physics primitives and their test framework, part of the
[Sub0 family](https://github.com/CraigHutchinson/.github/blob/main/docs/sub0-family.md).

**Status:** new; no code yet. The first primitives will be extracted from
[BikeLab](https://github.com/CraigHutchinson/BikeLab), where they are being proven.

## Scope

In:

- **Primitives:** bodies, constraints, springs and dampers, contact geometry, integrators and fixed-step stepping.
- **Verification:** analytical references, error and energy-balance measurement, convergence and
  step-refinement checks, deterministic replay.
- **Solver comparison:** run the same experiment through different solvers and measure their agreement.

Out:

- **Domain models:** the bike, tyre, suspension and rider solvers live in their products (for example
  PhantomFlow). Sub0Physics knows nothing of bikes.
- **Rendering, UI and product policy.**

## Rules

Follows the [forge rule](https://github.com/CraigHutchinson/.github/blob/main/docs/forge.md): a primitive is
extracted with its tests once a real consumer has proven it, product-neutral names only, and consumers pin
exact versions.

## Licence

[MIT](LICENSE.md). Copyright Craig Hutchinson.
