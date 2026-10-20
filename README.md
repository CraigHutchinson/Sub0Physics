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
| Integrators | `sub0physics/integrators.hpp` | `OdeSystem` concept; transactional fixed-step `Rk4` and `Midpoint` (step size as `f64::Time`) that reuse the endpoint derivative and never commit a refused or non-finite state; a fixed-size form for small systems and a `std::dynamic_extent` span form for run-time-sized ones |
| Elements | `sub0physics/elements.hpp` | Laws of one coordinate, each built by a validating `create` and answering a `LawResponse` (generalized force on its coordinate, exact Jacobians, `potential`, `dissipationPower`): `AxialSpringDamper` (preload, compression/rebound damping), `UnilateralSpringDamper` (push-only penalty contact), `TravelBuffer` (progressive end-stops inside a finite stroke) and `LinearSpringDamper<Coordinate>` (on a length or an angle) |
| References | `sub0physics/references.hpp` | Closed-form damped linear oscillator for all damping regimes, numerically guarded |
| Verification | `sub0physics/verification.hpp` | `EnergyAccount` and `EnergyLedger` (work/energy balance residual), `ErrorStatistics` (peak/RMS against a limit), `observedOrder` (convergence), `ChannelTable` (preallocated column recording) |
| Geometry | `sub0physics/geometry.hpp` | `PolylineProfile` of `ProfilePoint` (`sub0math::Point2<f64::Length>`; a validated height field, searched with `sub0math::closestPoint` with vertical faces) and `sub0math::Circle` clearance. Tied nearest features are `Witness` values (`point`, approach `normal`, `WitnessKind` vertex/edge/face, `index`) written to a caller-owned buffer: no allocation, an `overflowed` flag when it is full. `PreparedProfile`: a profile prepared once (structure of arrays, extents, outward unit normals) for allocation-free, windowed queries that report status aggregates and witnesses: `castRay` (first hit with a tie and ambiguity report, from a bound origin or a one-off), `sweepCircle` (time of impact over vertex caps and segment faces), `contourIntrusion` (a closed contour's deepest cut into the solid), `onSurface`, and a `measureCircleClearance` overload that agrees with the unprepared one bit for bit. See [docs/TERRAIN.md](docs/TERRAIN.md) |
| Solvers | `sub0physics/solvers.hpp` | `solveBoundedPeriodicTridiagonal`: a periodic chain energy (node stiffness plus neighbour coupling) minimised under lower bounds by projected Gauss-Seidel with the wrap nodes peeled; reusable `BoundedPeriodicSolution` storage, no allocation once it has capacity; KKT-checked in `tests/` |
| Core | `sub0physics/core.hpp` | `Id<Tag, Rep>`: a strongly typed identifier with explicit construction, `<=>`, `std::hash`, an invalid sentinel and a saturating `next()` |
| Multibody | `sub0physics/multibody.hpp` | Articulated rigid-body trees in reduced coordinates: a typed description compiled by `ModelBuilder` into an immutable `Model`; `State` and a per-model `Workspace`; kinematics, site poses and Jacobians, energy; CRBA mass matrix, RNEA bias and inverse dynamics, tree-sparse LTDL factorisation; forward dynamics; force elements that bind a law to a joint coordinate or to the distance between two sites (`CoordinateElementDescription`, `TendonElementDescription`), with exact Jacobians including a tendon's geometric stiffness, site wrenches and energy; one `Stepper` shape for the default `LinearlyImplicitEuler` (stiff elements enter through their Jacobians), `SemiImplicitEuler` and the RK reference steppers (`Rk4Stepper`, `MidpointStepper`), all accounting dissipated energy and external work for the `EnergyLedger`; constraints described in the model (loop closures between sites, strict unilateral limits on a coordinate or a distance), solved at velocity level with the step's own factorised matrix and a small LCP, steps split at limit engagement, positions and velocities projected onto the closures to a requested residual (`projectConstraints`), and joint reactions and site loads from a post-solve RNEA with every force (`computeReactions`). A model whose motion stays in a coordinate plane runs in planar algebra automatically (`Formulation::Planar`, the same results to rounding); massless bodies carry extra joints; external forces at sites with optional Jacobians (`Forcing::siteForces`), an implicit-term lever (`ImplicitTerms`, including MuJoCo's `implicitfast` structure), the energy-better `LinearlyImplicitMidpoint`, and a two-phase step (`prepareLinearlyImplicitStep`, prepared site queries, `completeLinearlyImplicitStep`) for consumers that form their forcing from the state. No allocation once the workspace exists. Design: [docs/MULTIBODY.md](docs/MULTIBODY.md) |
| Terrain | `sub0physics/terrain.hpp` | `GridSpec`: paging of the plane into levels of pages on `sub0math::RegularGrid`, with a runtime cells-per-side, a placement anchor and node or cell-centre registration, so a source of any shape (a 256-cell survey page, a 40-cell fixture) is expressible exactly; ownership (`pageOf`, `cellOf`, half-open) and coverage (`cellsTouching`, `pagesTouching`, closed) are separate operations. `TileKey` and `PageFrame` (double-precision terrain frame, float page-local frame). `HeightPage`: a quantised 257 x 257 height page with an authored-node mask, a min/max pyramid and a fixed-triangulation `heightAt`. `CorridorBuilder`: oriented boxes that bound a wheel footprint over a pose pair, rasterised to grid cells. `SourceTileCache`: a bounded, deterministic tile cache (de-duplicated loads, least-recently-used residency, pins). `extractColumns` and `extractHeights`/`extractCorridor`: owned structure-of-arrays extraction of columns with spans, and of node rows and triangles with a ray cast, over the corridor's cells. `strideForTolerance` and `tightCorridor`: target-driven resolution and coverage levers. `terrain/testing.hpp`: an extruded-profile source (a 2D profile baked to pages with its vertical faces as features) and an analytic-surface catalogue (`Plane`, `Ramp`, `Step`, `RoundedBlock`, `Bridge`, `Berm`, `Rut`) with exact heights and normals for tests. Design notes: [docs/TERRAIN.md](docs/TERRAIN.md) |
| Profiling | `sub0physics/profiling.hpp` | `ProfileTask`: ITT task markers for VTune, compiled out unless `SUB0PHYSICS_ENABLE_ITT` |
| Level of detail | `sub0physics/lod.hpp` | `ModelDescriptor` and `FidelitySelector`: choose between models of one mechanism from a physical demand, with hysteresis, dwell and commit-after-transfer. See [docs/LOD.md](docs/LOD.md) |

Everything is reachable from `sub0physics/sub0physics.hpp`.

```cpp
struct Wheel                                  // a model owned by the product
{
    sub0physics::UnilateralSpringDamper tyre =
        *sub0physics::UnilateralSpringDamper::create({.stiffness = newtonsPerMetre(90000.0), .damping = newtonSecondsPerMetre(400.0)});
    std::expected<void, sub0physics::Fault> evaluate(double t, const sub0physics::StateVector<2>& x,
                                                     sub0physics::StateVector<2>& dxdt);
};

sub0physics::Rk4<2> rk4;
if (auto stepped = rk4.stepTick(wheel, tick, milliseconds(0.5), state); !stepped)
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
