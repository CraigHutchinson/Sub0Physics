# Performance

Performance here is a tested property, not a claim. Correctness gates come first; timing is measured
on named hosts; profiling is one command away. The benchmark, budget, profiler and vectorization tooling
comes from Sub0Math (`cmake/Sub0MathTooling.cmake`, `tools/perf/`, documented in its
[PERFORMANCE.md](https://github.com/CraigHutchinson/Sub0Math/blob/main/docs/PERFORMANCE.md)); this page
lists what is specific to Sub0Physics. Configure writes the scripts' location to
`<build>/sub0math_tooling_dir.txt`.

## Codified gates (every test run, every platform)

These are deterministic, so they block in CI:

| Gate | Test | Contract |
|---|---|---|
| No allocation while stepping | `Sub0Physics_AllocationTests` | RK4/midpoint steps, recording, ledgers, statistics, selection and the closed form allocate nothing; a clearance query allocates only its witness list |
| No allocation in multibody dynamics | `Sub0Physics_MultibodyAllocationTests` | Symplectic Euler steps, inverse dynamics, site kinematics and Jacobians, and energy on a 10-DOF tree allocate nothing once the `Workspace` exists |
| Evaluation count | `test_integrators.cpp` | Tick stepping costs exactly one evaluation per stage: 4 per RK4 step after the first |
| Convergence order | `test_integrators.cpp` | RK4 observes order 4, midpoint order 2; a slower method that is "faster" is not a win |
| Header self-containment | `Sub0Physics_HeaderCheck` | Every public header compiles alone |

Replacement `operator new` lives in its own test executable so it cannot affect the others.

## Timing budgets and the gate

The per-host baselines, the interleaved-round measurement, the pre-push gate and its verdicts are
Sub0Math's, and are documented once in its
[Timing budgets](https://github.com/CraigHutchinson/Sub0Math/blob/main/docs/PERFORMANCE.md#timing-budgets-per-host)
section. Sub0Physics uses them unchanged (`Sub0Physics_Bench_check` and `Sub0Physics_Bench_baseline`,
baselines in `perf/baselines/<host>.json`). What is specific here:

**Controls.** The hand-written references in `tests/bench_physics.cpp` and `tests/bench_multibody.cpp` do
not use the library, so whatever moves them is the environment. They are the controls that the gate
factors out, and each is the twin of one library call:

| Control (hand-written) | Library call it is paired with |
|---|---|
| `page height x1024, hand-written raw arrays` | `HeightPage::heightAt` |
| `clearance x256, hand-written twin` | `measureCircleClearance` (256 queries, 2048-point profile, witnesses into a fixed buffer) |
| `grid cellOf x1024, hand-written floor and range` | `GridSpec::cellOf` (ownership through `sub0math::RegularGrid`) |
| `page size x1024, 256 cells compile-time` | `page size x1024, 256 cells at run time` (the same hand-written kernel; a decision control for a runtime-sized `HeightPage`, see docs/TERRAIN.md) |
| `travel buffer x1024, hand-written scalars` | `TravelBuffer::respond` |
| `planar tree 10-DOF mass matrix, hand-written planar` | `computeMassMatrix` |
| `planar tree 10-DOF LTDL factor and solve, hand-written` | `factorMassMatrix` + `solveMassMatrix` |
| `planar tree 10-DOF forward dynamics, hand-written planar` | `forwardDynamics` |

**Pairs judged by ratio.** Each row above is a pair. The ratio of library to reference is what the gate
and the Open gaps table track, because it does not depend on the host. A reference weakened to look
better is a defect: iterate, reshape the API, or record the gap under Open gaps (Sub0Math's rule).

The remaining benchmarks (`Rk4` and `Midpoint` steps, the single-query circle clearance, surface height, the
closed form, `corridor create` and `corridor rasterise`) have no hand-written twin, so they are judged against
the host's baseline only.

## Profiling

Build with the `profile` preset: optimised code with debug info and frame pointers, and ITT task
markers (`sub0physics::ProfileTask`) that VTune shows on its timeline. Then:

```bash
cmake --preset profile -DSUB0MATH_PROFILE_FILTER="Rk4" -DSUB0MATH_PROFILE_SECONDS=5
cmake --build --preset profile --target Sub0Physics_Bench_vtune_hotspots    # where time goes
cmake --build --preset profile --target Sub0Physics_Bench_vtune_uarch       # why: front/back-end, ports, memory
cmake --build --preset profile --target Sub0Physics_Bench_advisor_roofline  # vectorization and roofline
cmake --build --preset profile --target Sub0Physics_Bench_perf_record       # Linux perf
```

Each target runs the benchmark with `--soak` so the selected hot path dominates the samples, writes
results to `build-profile/profiles/<analysis>-<filter>-<time>/`, prints the summary and tells you how to
open the GUI. Sub0Math's `tools/perf/profile.py` is the same thing without CMake. Benchmark executables
accept `--filter`, `--json`, `--soak` and `--list`.

Platform notes:

- **VTune** reads GCC/Clang DWARF and MSVC PDB symbols on Windows and Linux.
- **Advisor** on Windows gives full loop and FLOP data for MSVC or Intel compiler builds; GCC builds
  still show loops and time.
- **ITT** markers come from `intel/ittapi`, fetched at a pinned commit only when
  `SUB0PHYSICS_ENABLE_ITT=ON`; otherwise `ProfileTask` compiles to nothing.

## Vectorization and aliasing

```bash
cmake --preset vectorization
python "$(cat build-vec/sub0math_tooling_dir.txt)/vectorization_report.py" --build-dir build-vec --target Sub0Physics --target Sub0Physics_Bench
```

The report rebuilds with GCC `-fopt-info-vec-*`, Clang `-Rpass*`, or MSVC `/Qvec-report:2` and lists
remarks for this project's own files only, so you can confirm that a loop vectorized after a change
and see why another did not. The library also builds with `-fstrict-aliasing -Wstrict-aliasing=2`.

Rules that keep hot code fast and analysable:

- Fixed-size state (`std::array`) and integrator-owned scratch; no allocation, no virtual calls in steps.
- Stage loops run over contiguous arrays the integrator owns; a system writes only its output array.
  Do not alias a system's output with its input.
- Faults are values with static strings, so the error path costs no allocation.
- Interfaces carry Sub0Math quantities; models convert to plain SI doubles once, at their edge, and
  integrator state stays plain. Quantities cost nothing (Sub0Math's paired benchmarks), while the
  integrators stay generic over any state.
- Time comes from integer ticks, so the endpoint derivative is reused without approximate comparisons.

## Findings so far

| Date | Finding | Evidence | Action |
|---|---|---|---|
| 2026-10-15 | Clearance queries scanned every segment twice: 14 µs at 2048 points | nanobench | Horizontal pruning from the sorted x; 0.8 µs, brute-force agreement test |
| 2026-10-15 | Accumulated step times missed the endpoint-derivative cache on 31% of 0.5 ms steps | counting test | `stepTick`; exactly 4 evaluations per RK4 step |
| 2026-10-15 | GCC keeps `evaluate` out of line in the RK4 stages ("statement clobbers memory") | VTune hotspots, vectorization report | Forcing inlining measured no change (39.5 vs 39.8 ns/step); not applied |
| 2026-10-15 | `std::hypot` dominated dense clearance queries (11.7 µs on a 1024-segment arc in a loaded-wheel consumer) | Consumer bench | Squared distances, `sqrt` only on improvement, coordinates bounded at 1e100 m; 3.9 µs |
| 2026-10-15 | Moving the geometry interface to Sub0Math quantities regressed clearance by 52% (395 to 601 ns, best of 7) | Budget gate, then a scalar rewrite of the kernel to locate it | Fixed in Sub0Math (compile-time component expansion, no zero-initialised sums); 394 vs 402 ns, within noise |
| 2026-10-16 | Terrain page height lookup started 2.5x slower than the hand-written kernel (raw arrays, bit mask, fixed diagonal) | Paired benchmark `page height x1024` | Inline definition with reciprocal spacing (7.1 to 5.4 µs); a plain `HeightSample` return instead of `std::optional` (5.6 to 3.5 µs; `std::optional<float>` was equally slow, so the optional, not the quantity, kept the result in memory on GCC 16); two plain branches instead of a merged triangle select (3.5 to 3.0 µs). Force-inlining alone: no gain. Now 1.07x |
| 2026-10-19 | Travel buffer (`TravelBuffer::respond`, 7 outputs incl. both Jacobian entries) first measured 1.2x slower than the hand-written scalars | Paired benchmark `travel buffer x1024`; the benchmark itself fed `f32` `sub0math::Length` (float round trip), so the gap was a harness error, not the element | Fixed the benchmark to `f64` quantities: 1.00x to 1.05x over six runs (3.3 to 3.7 µs per 1024) against a constexpr-parameter reference. `std::fmax` vs compare-select: no difference; branch-free ternaries for the closing side |
| 2026-10-19 | Multibody forward dynamics on a 10-DOF planar tree: identity joint placements still paid a 3x3 product per joint, and CRBA/RNEA used general 6D products with the joint's motion subspace | Pairs `planar tree 10-DOF ...` against a hand-written planar implementation: forward dynamics 3.3x, then 2.8x | Placement-rotation product skipped when it is the identity; `S` products use that S is purely angular or purely linear |
| 2026-10-20 | Grid ownership moved onto `sub0math::RegularGrid` (`GridSpec::cellOf`, was a hand-rolled floor and range test) | Pair `grid cellOf x1024`: 0.97x to 1.07x of the hand-written floor/range kernel over three runs (3.3 to 3.9 �s per 1024) | Adopted: one set of boundary rules at no measurable cost |
| 2026-10-20 | A runtime cell count in the page lookup | Pair `page size x1024`: the same raw-array kernel with 256 as a run-time value is 1.04x to 1.28x of the constant one (about 1.1x typical) | `HeightPage` stays 256 and compile-time; other shapes get their own page type from the same `GridSpec` (docs/TERRAIN.md) |
| 2026-10-20 | `HeightPage::heightAt` after the grid change (pages built on `GridSpec`, the lookup untouched) | Pair `page height x1024`: 0.73x to 0.90x of the hand-written kernel over three runs (3.2 to 3.4 �s vs 3.7 to 4.4 �s) | None; not worse than before, and the open-gap row below stands as the conservative number |

## Open gaps

Deficiencies against the hand-written reference (rule: Sub0Math `docs/PERFORMANCE.md`).

| Since | Pair | Gap | Cause | Follow-up |
|---|---|---|---|---|
| 2026-10-16 | `page height x1024`, `HeightPage::heightAt` | about 1.07x in one process, 1.10x interleaved (pair ratio) | Per-call conversion and bounds work around an otherwise identical kernel | Terrain step G3: the batch extraction API (many samples per call, SoA out) is the real consumer; close it there or accept it with its numbers |
| 2026-10-20 | `clearance x256`, `measureCircleClearance` | 1.9x to 2.0x (about 132 to 137 vs 67 to 71 �s per 256 queries, three runs, 2048-point sinusoid) | The reference works on prepared structure-of-arrays data (per-segment extents and reciprocal squared lengths, so no division per segment), takes a window from the segment under the centre and runs four independent minima; the API divides per segment inside `sub0math::closestPoint`, reads array-of-structures points and walks with a data-dependent early-out. The witness buffer is already caller-owned (the allocation is gone), so this is the compute gap | Wave 5 (D12): `PreparedProfile` with reciprocal lengths and normals, then re-measure; the twin is the target |
| 2026-10-19 | `planar tree 10-DOF`: `computeMassMatrix`, `forwardDynamics` | 2.9x and 2.8x interleaved (276 vs 92 ns; 585 vs 216 ns); the tree-sparse LTDL factor and solve is at parity (0.92x, 69 ns) | The reference uses three-component planar algebra and a compile-time tree; the core runs the general 3D formulation, by design (docs/MULTIBODY.md section 3) | M1d, with the bike rig as the consumer: a planar joint type or a compact transform for axis-aligned revolute joints (rotation about one coordinate axis), judged on the rig's whole step. About 0.6 us per evaluation is already far inside a 0.5 ms step |
| 2026-10-20 | `travel buffer x1024`, `TravelBuffer::respond` | 1.5x to 1.7x (about 5.0 vs 3.0 µs per 1024), the same at 1f60d01, so not caused by wave 1c | Measured with GCC 16.2.0 (WinLibs), -O3, on the development host; the 2026-10-19 finding recorded 1.00x to 1.05x. The cause is not yet known: the code and the harness are unchanged, so suspect the host state or a Sub0Math change since | Re-measure at the wave 2a re-baseline (`-ffp-contract=off` and the deterministic-FP helper change both sides' code generation); if it persists, bisect Sub0Math and the fmax/compare-select choice |
