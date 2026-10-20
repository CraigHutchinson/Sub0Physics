# Sub0Physics

Product-agnostic physics primitives and the tools to verify them, part of the
[Sub0 family](https://github.com/CraigHutchinson/.github/blob/main/docs/sub0-family.md). C++23, no
allocation on stepping paths. Built on [Sub0Math](https://github.com/CraigHutchinson/Sub0Math): element,
reference and geometry interfaces take dimensioned quantities (`f64::Stiffness`, `f64::Length`), so unit
mistakes do not compile.

Domain models (vehicles, tyres, suspensions, riders) live in their products and build on these
primitives. Sub0Physics knows nothing of them.

## What's in it

| Area | Header | Provides |
|---|---|---|
| Integrators | `sub0physics/integrators.hpp` | `OdeSystem` concept; transactional fixed-step `Rk4` and `Midpoint` that reuse the endpoint derivative and never commit a refused or non-finite state |
| Elements | `sub0physics/elements.hpp` | `AxialSpringDamper` (preload, compression/rebound damping), `UnilateralSpringDamper` (push-only penalty contact) and `TravelBuffer` (progressive end-stops inside a finite stroke, built by the validating `TravelBuffer::create`, with exact Jacobians), with stored energy and dissipated power |
| References | `sub0physics/references.hpp` | Closed-form damped linear oscillator for all damping regimes, numerically guarded |
| Verification | `sub0physics/verification.hpp` | `EnergyLedger` (work/energy balance residual), `ErrorStatistics` (peak/RMS against a limit), `observedOrder` (convergence), `ChannelTable` (preallocated column recording) |
| Geometry | `sub0physics/geometry.hpp` | `PolylineProfile` of `ProfilePoint` (`sub0math::Point2<f64::Length>`; a validated height field, searched with `sub0math::closestPoint` with vertical faces) and `sub0math::Circle` clearance with vertex/segment witnesses and approach normals |
| Multibody | `sub0physics/multibody.hpp` | Articulated rigid-body trees in reduced coordinates: a typed description compiled by `ModelBuilder` into an immutable `Model`; `State` and a per-model `Workspace`; kinematics, site poses and Jacobians, energy; CRBA mass matrix, RNEA bias and inverse dynamics, tree-sparse LTDL factorisation; forward dynamics and a semi-implicit Euler step. No allocation once the workspace exists. Design: [docs/MULTIBODY.md](docs/MULTIBODY.md) |
| Terrain | `sub0physics/terrain.hpp` | `GridSpec` (paging of the plane into levels of pages, built by the validating `GridSpec::create`), `TileKey` and `PageFrame` (double-precision terrain frame, float page-local frame), and `HeightPage`: a quantised 257 x 257 height page with an authored-node mask, a min/max pyramid and a fixed-triangulation `heightAt` |
| Profiling | `sub0physics/profiling.hpp` | `ProfileTask`: ITT task markers for VTune, compiled out unless `SUB0PHYSICS_ENABLE_ITT` |
| Level of detail | `sub0physics/lod.hpp` | `ModelDescriptor` and `FidelitySelector`: choose between models of one mechanism from a physical demand, with hysteresis, dwell and commit-after-transfer. See [docs/LOD.md](docs/LOD.md) |

Everything is reachable from `sub0physics/sub0physics.hpp`.

```cpp
struct Wheel                                  // a model owned by the product
{
    sub0physics::UnilateralSpringDamper tyre{newtonsPerMetre(90000.0), newtonSecondsPerMetre(400.0)};
    std::expected<void, sub0physics::Fault> evaluate(double t, const sub0physics::StateVector<2>& x,
                                                     sub0physics::StateVector<2>& dxdt);
};

sub0physics::Rk4<2> rk4;
if (auto stepped = rk4.step(wheel, t, h, state); !stepped)
{
    report(stepped.error());                  // state is unchanged
}
```

## Build and test

```bash
cmake --preset default
cmake --build --preset default
ctest --preset default
```

`build/tests/Sub0Physics_Bench` runs the benchmarks (not part of ctest). Allocation and
evaluation-count gates run in ctest; timing budgets, VTune, Advisor, `perf` and vectorization reports are
one command each, using Sub0Math's tooling: see [docs/PERFORMANCE.md](docs/PERFORMANCE.md).

Develop against local checkouts by setting `CPM_Sub0Math_SOURCE` (and similar) in your environment:
CPM 0.43+ uses them, with a warning, and CI builds the pinned commits. CPM remembers a package's source in
the build cache, so configure a fresh build directory after setting or clearing one. CI builds with GCC, Clang,
MSVC and Apple Clang, runs ASan/UBSan, and compiles every public header on its own.

Consume it with CPM at an exact commit:

```cmake
CPMAddPackage(NAME Sub0Physics GITHUB_REPOSITORY CraigHutchinson/Sub0Physics
    GIT_TAG <full-sha> EXCLUDE_FROM_ALL YES OPTIONS "SUB0PHYSICS_BUILD_TESTING OFF")
target_link_libraries(app PRIVATE Sub0Physics::Sub0Physics)
```

## Where the code came from

The first primitives were ported from [BikeLab](https://github.com/CraigHutchinson/BikeLab), the
JavaScript laboratory where they were proven. [upstream/bikelab.json](upstream/bikelab.json) records the
exact BikeLab commit, each ported file's blob hash, and the files here that carry it. CI checks that
record, and a weekly workflow opens an issue listing BikeLab changes since then. See
[docs/UPSTREAM.md](docs/UPSTREAM.md).

## Rules

Follows the [forge rule](https://github.com/CraigHutchinson/.github/blob/main/docs/forge.md): a primitive
arrives with its tests once a real consumer has proven it, names stay product-neutral, and consumers
pin exact versions. Style: [STYLE_GUIDE.md](STYLE_GUIDE.md).

## Licence

[MIT](LICENSE.md). Copyright Craig Hutchinson.
