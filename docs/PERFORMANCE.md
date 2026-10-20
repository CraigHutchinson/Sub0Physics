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
| No allocation while stepping | `Sub0Physics_AllocationTests` | RK4/midpoint steps, recording, ledgers, statistics, selection and the closed form allocate nothing; clearance queries (both overloads), every `PreparedProfile` query and the bounded periodic chain solve (once the solution has capacity) allocate nothing |
| No allocation in multibody dynamics | `Sub0Physics_MultibodyAllocationTests` | Semi-implicit Euler and `Rk4Stepper` steps (with changing forcing), inverse dynamics, site kinematics and Jacobians, and energy on a 10-DOF tree allocate nothing once the `Workspace` exists; with force elements, linearly implicit Euler (dense and tree factorisation), semi-implicit Euler, RK4, `evaluateElements` and energy allocate nothing |
| Implicit stepping gate (M1b) | `test_multibody_elements.cpp` | Element Jacobians match central differences; stiff spring-damper stable at a step where RK4 diverges; energy ledger bounded; implicit Euler observed order 1 |
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
| `clearance x256, measureCircleClearance (prepared)` | The same twin, judged for the `PreparedProfile` overload |
| `chain solve 512 nodes, hand-written sweeps` | `solveBoundedPeriodicTridiagonal` (512 nodes, a pressed patch; u and gradient checked bit for bit) |
| `grid cellOf x1024, hand-written floor and range` | `GridSpec::cellOf` (ownership through `sub0math::RegularGrid`) |
| `page size x1024, 256 cells compile-time` | `page size x1024, 256 cells at run time` (the same hand-written kernel; a decision control for a runtime-sized `HeightPage`, see docs/TERRAIN.md) |
| `travel buffer x1024, hand-written scalars` | `TravelBuffer::respond` |
| `planar tree 10-DOF mass matrix, hand-written planar` | `computeMassMatrix` |
| `planar tree 10-DOF LTDL factor and solve, hand-written` | `factorMassMatrix` + `solveMassMatrix` |
| `planar tree 10-DOF forward dynamics, hand-written planar` | `forwardDynamics` |
| `planar tree 10-DOF implicit Euler step, 2 elements, hand-written planar` | `stepLinearlyImplicitEuler` (the tree with a cart spring-damper and a preloaded strut between two branches; state, dissipated energy and external work checked against the library before timing) |
| `planar tree 10-DOF elements with Jacobians, hand-written planar` | `evaluateElements` with Jacobians (world poses included on both sides) |
| `Rk4 step, contact oscillator, hand-written` | `Rk4<2>::stepTick` on the contact oscillator (raw doubles, the law inline, the endpoint derivative kept; non-finite endpoints refused; trajectory checked against the library) |
| `Midpoint step, contact oscillator, hand-written` | `Midpoint<2>::stepTick`, the same way |

**Pairs judged by ratio.** Each row above is a pair. The ratio of library to reference is what the gate
and the Open gaps table track, because it does not depend on the host. A reference weakened to look
better is a defect: iterate, reshape the API, or record the gap under Open gaps (Sub0Math's rule).

The remaining benchmarks (the single-query circle clearance, surface height, the closed form, `corridor create`
and `corridor rasterise`) have no hand-written twin, so they are judged against the host's baseline only.

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
| 2026-10-20 | `TravelBuffer::respond` at 1.5x to 1.8x: `std::fmax` compiles to a library call on MinGW GCC 16 (two `call fmax` per sample in the disassembly), while the reference's `std::max` is one `maxsd` | Disassembly of both loops; pair `travel buffer x1024` | `std::max` (same result for NaN and signed zero here): 1.02x to 1.04x (2.83 vs 2.73 µs min, 3.04 vs 2.98 µs median, 7 interleaved rounds). A precomputed reciprocal width measured no gain; not kept |
| 2026-10-20 | The M1b law response (5 common fields plus extras) pushed a station's `evaluate` past GCC's early-inlining budget ("growth 15 exceeds early-inlining-insns"): `Rk4 step, contact oscillator` 43 to 52 ns | `-fopt-info-inline-missed`, 15 interleaved rounds | `ExplicitRungeKutta::stepBetween` is flattened (GCC/Clang): 43.4 vs 43.0 ns min, 47.0 vs 46.2 median against the pre-change build; Midpoint 23.6 vs 22.5 min, 24.4 vs 24.1 median. Supersedes the 2026-10-15 "force-inlining evaluate: no gain" note, measured when GCC kept it out of line on both sides |
| 2026-10-20 | Revolute joints about z paid the general Rodrigues products | Pairs `planar tree 10-DOF ...` | Short branch for an axis of (0, 0, ±1): about 4% on mass matrix and forward dynamics; one pose pipeline (D6) does not touch the dynamics path |
| 2026-10-20 | Wave 5 (D12): `PreparedProfile` and the chain solver moved here. The recorded clearance gap (1.9x to 2.0x) did not reproduce | GCC 16, MinGW, five interleaved rounds each: unprepared `measureCircleClearance` 0.89x to 0.94x of the twin before and 0.94x to 0.96x after (unchanged code; the earlier gap predates Sub0Math's force-inlined `closestPoint`); prepared overload 0.85x to 0.90x (79 to 82 vs 90 to 93 us per 256 queries); chain solve 0.89x to 0.95x | Prepared clearance kept: window from the centre segment, four minima. A reciprocal squared length was another 10% faster (74 vs 82 us) but changes the last bit of distance and of the projection parameter, so it was not taken (docs/TERRAIN.md Exactness). The clearance gap row is removed; no gap remains |

| 2026-10-20 | M1b: the fixed-size RK step lost about 4 ns to its twins (`Rk4 step` 1.04x to 1.19x, `Midpoint step` 1.29x to 1.32x by minimum, GCC 16): the stages were members, so each step stored the endpoint and its derivative as two scalars and copied them on as 16-byte loads, which cannot be store-forwarded | Disassembly of the flattened step (`movsd` pairs then `movdqu` of the same 16 bytes); pairs over 7 interleaved rounds | The fixed size keeps its stages local to the step and commits from them; the dynamic size is unchanged. `Rk4 step` 0.97x (26.1 vs 26.7 ns min), `Midpoint step` 0.96x to 1.01x (12.8 vs 12.7 ns); 31.4 to 26.1 and 17.1 to 12.8 ns absolute. Removing the cache-state compare measured no change, so the cache contract stays |
| 2026-10-20 | M1b element evaluation (the 10-DOF tree's cart spring and branch-coupling strut, Jacobians included) first measured 2.6x its twin (281 vs 108 ns min) | Pair `planar tree 10-DOF elements with Jacobians` | The site-Hessian term as a triple product with a precomputed `u x w` per joint, and the tendon's support without prismatic joints that carry both ends (they change no entry): 190 ns, 2.0x. The rest is the 3D joint transforms (see Open gaps) |
| 2026-10-20 | Element Jacobians as dense `n x n` sums, or per-element blocks | One-off measurement, a 40-DOF tree (cart and three 13-link branches) with two branch-coupling tendons: zeroing the two 40 x 40 matrices 87 to 94 ns, `evaluateElements` with Jacobians 3.1 us, forward dynamics 4.1 to 4.3 us, the implicit step 8.8 to 9.9 us | Dense sums: the clearing is 1% of the step, an element still writes only its support, and the step needs the dense matrix anyway; blocks would add a scatter pass. Decision recorded in docs/MULTIBODY.md (M1b status) |

## Open gaps

Deficiencies against the hand-written reference (rule: Sub0Math `docs/PERFORMANCE.md`).

| Since | Pair | Gap | Cause | Follow-up |
|---|---|---|---|---|
| 2026-10-16 | `page height x1024`, `HeightPage::heightAt` | about 1.07x in one process, 1.10x interleaved (pair ratio) | Per-call conversion and bounds work around an otherwise identical kernel | Terrain step G3: the batch extraction API (many samples per call, SoA out) is the real consumer; close it there or accept it with its numbers |
| 2026-10-20 | `planar tree 10-DOF implicit Euler step, 2 elements` | 2.0x (881 vs 445 ns min, 921 vs 454 median, 1.98x to 2.05x per round, 7 rounds, GCC 16); Clang 22 `-march=x86-64-v3` 2.5x | Mostly the recorded CRBA/RNEA gap below (about 320 of the 440 ns: the step runs the same transforms, RNEA and CRBA as `forwardDynamics`), then the elements (next row); the dense `LDL^T` and assembly are written the same way on both sides | M1d with the CRBA/RNEA gap (a planar joint type or compact axis-aligned transforms), judged on the rig's whole step; 0.9 us per step is far inside a 0.5 ms step |
| 2026-10-20 | `planar tree 10-DOF elements with Jacobians` | 2.0x (190 vs 92 ns min, 7 rounds, GCC 16); Clang 3.1x | The tendon geometry runs in 3D on the general joint transforms (world rotation of each axis, 3-vector crosses) where the twin uses planar angles | With the M1d planar-joint work; it shares the cause |
| 2026-10-20 | `Rk4 step` and `Midpoint step, contact oscillator` under Clang 22 (`-march=x86-64-v3`, Windows) | 1.2x and 1.35x to 1.4x (GCC 16 is at parity, above) | Clang does not inline the flattened step's nested calls the way GCC does (`stepBetween` or `advance` stay out of line); force-inlining them measured no clear gain, so it was not kept | Re-measure on the Linux Clang CI leg; GCC is the host baseline's compiler |
| 2026-10-19 | `planar tree 10-DOF`: `computeMassMatrix`, `forwardDynamics` | 2.9x and 2.65x by median at the M1b foundations (271 vs 94 ns; 591 vs 223 ns; 3.2x and 2.9x by minimum), against 2.8x and 2.9x for the build before it, so not worsened; the tree-sparse LTDL factor and solve is at parity (1.02x median) | The reference uses three-component planar algebra and a compile-time tree; the core runs the general 3D formulation, by design (docs/MULTIBODY.md section 3) | M1d, with the bike rig as the consumer: a planar joint type or a compact transform for axis-aligned revolute joints (rotation about one coordinate axis), judged on the rig's whole step. About 0.6 us per evaluation is already far inside a 0.5 ms step |
