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
| No allocation while stepping | `Sub0Physics_AllocationTests` | RK4/midpoint steps, recording, ledgers, statistics, selection and the closed form allocate nothing; clearance queries (both overloads), every `PreparedProfile` query, the `ConvexTerrain` contact and clearance queries and the bounded periodic chain solve (once the solution has capacity) allocate nothing |
| No allocation in multibody dynamics | `Sub0Physics_MultibodyAllocationTests` | Semi-implicit Euler and `Rk4Stepper` steps (with changing forcing), inverse dynamics, site kinematics and Jacobians, and energy on a 10-DOF tree allocate nothing once the `Workspace` exists; with force elements, linearly implicit Euler (dense and tree factorisation), semi-implicit Euler, RK4, `evaluateElements` and energy allocate nothing; with constraints, implicit and semi-implicit Euler (closures, an LCP row, split steps), RK4 with closures, `projectConstraints` and `computeReactions` allocate nothing; in both formulations, steps with site forces and forcing Jacobians (implicit Euler, midpoint, semi-implicit Euler), the prepared step and its site queries, batch site kinematics with energy, `computeAppliedForce` and `computeReactionsWithForcing` allocate nothing |
| Dyadic step and swept certificate (M2d) | `test_dyadic_step.cpp`, `test_convex_terrain.cpp` | Steps aligned, never across a boundary, landing on it, events only inside floor steps, replay identical; step doubling meets its target and cost falls as it loosens; the swept bound never exceeds the true gap at 1001 instants and certifies more than the one-sided bound |
| Delassus forms (M2d) | `test_multibody_ring_contact.cpp` | The uniform form with supports is bit for bit the form without; the structured form agrees to 1e-11 relative; a kept preparation completes bit for bit as a fresh one |
| Ring substructure gate (M2a) | `test_multibody_ring.cpp` | The condensed ring agrees with its expanded oracle (bodies and tendons, dense path) at N = 8 and 16: dynamics, energy, element forces and Jacobian blocks, reactions to 1e-13, steps and runs of every stepper to 1e-10; mass properties; free flight; the ledger closes with ring damping; a modal frequency to 1e-5 of the closed form; the periodic block solver against dense solves; a run digest printed for the cross-compiler comparison |
| Planar formulation, forcing and two-phase step gate (M1d) | `test_multibody_planar.cpp` | Every query of the planar formulation equals the general one to 1e-13 and stepped trajectories to 1e-11; site forces equal `J^T F` and are applied where they act in reactions; forcing Jacobians step like the same element; `DampingOnly` equals the closed-form implicitfast step; the midpoint step keeps an oscillator's energy to 1e-12 and is order 2; prepare then complete equals the one-call step bit for bit; a massless leaf is refused |
| Constraints gate (M1c) | `test_multibody_constraints.cpp` | Slider-crank closure residual within the requested tolerance and not drifting (every stepper); a hard stop never overshot, the impact located to the requested gap and its loss in the ledger; reactions equal to independent statics; ledger channels first order |
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
| `ring contacts 64 edges, square rock, hand-written twin` | `ConvexTerrain::ringContacts` (64 casing edges, 8 mm skin, against the 40 mm square rock's pieces: the same window, height culls and closed-form capsule-slab distance on raw arrays, slab corners built from the vertex arrays; every contact checked equal before timing). `ring clearance bound 64 edges` is timed alongside, without a twin |
| `chain solve 512 nodes, hand-written sweeps` | `solveBoundedPeriodicTridiagonal` (512 nodes, a pressed patch; u and gradient checked bit for bit) |
| `grid cellOf x1024, hand-written floor and range` | `GridSpec::cellOf` (ownership through `sub0math::RegularGrid`) |
| `page size x1024, 256 cells compile-time` | `page size x1024, 256 cells at run time` (the same hand-written kernel; a decision control for a runtime-sized `HeightPage`, see docs/TERRAIN.md) |
| `page heights batch x1024, hand-written raw arrays` | `HeightPage::heightsAt` (structure-of-arrays in and out; heights and flags checked against the library before timing) |
| `height extract 128x64 cells, hand-written raw arrays` | `extractHeights` (node rows, cell states, triangles and normals of a 128 x 64-cell block; every output array checked against the library) |
| `column extract 336 columns, 4 tiles, hand-written` | `extractColumns` (the R0a shape: 24 x 14 cells over four 40-cell tiles, spans included) |
| `tile cache x1024 accesses, hand-written array LRU` | `SourceTileCache::access` and `pump` (1024 accesses over 48 tiles in a 32-tile cache; hit, miss and eviction counts checked after warm-up) |
| `travel buffer x1024, hand-written scalars` | `TravelBuffer::respond` |
| `planar tree 10-DOF mass matrix, hand-written planar` | `computeMassMatrix` |
| `planar tree 10-DOF LTDL factor and solve, hand-written` | `factorMassMatrix` + `solveMassMatrix` |
| `planar tree 10-DOF forward dynamics, hand-written planar` | `forwardDynamics` |
| `planar tree 10-DOF implicit Euler step, 2 elements, hand-written planar` | `stepLinearlyImplicitEuler` (the tree with a cart spring-damper and a preloaded strut between two branches; state, dissipated energy and external work checked against the library before timing) |
| `planar tree 10-DOF elements with Jacobians, hand-written planar` | `evaluateElements` with Jacobians (world poses included on both sides) |
| `slider-crank constrained implicit Euler step, hand-written planar` | `stepLinearlyImplicitEuler` on the slider-crank with its closure: the same algorithm in closed form (planar mass matrix, two closure rows, 2 x 2 solves, projection to the same tolerance, velocity projection, closure work and projection energy measured); state and energy channels checked against the library over 500 steps before timing |
| `two 64-node rings implicit Euler step, hand-written condensed` | `stepLinearlyImplicitEuler` on a five-coordinate tree with two 64-node rings (the same condensation, the tree and ring size hard-coded, 2 x 2 blocks unrolled; state and dissipation checked against the library and the expanded oracle over 200 steps before timing). The pair's third row, `expanded bodies (dense path)`, is the same description as ordinary bodies and tendons |
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
| 2026-10-20 | Grid ownership moved onto `sub0math::RegularGrid` (`GridSpec::cellOf`, was a hand-rolled floor and range test) | Pair `grid cellOf x1024`: 0.97x to 1.07x of the hand-written floor/range kernel over three runs (3.3 to 3.9 µs per 1024) | Adopted: one set of boundary rules at no measurable cost |
| 2026-10-20 | A runtime cell count in the page lookup | Pair `page size x1024`: the same raw-array kernel with 256 as a run-time value is 1.04x to 1.28x of the constant one (about 1.1x typical) | `HeightPage` stays 256 and compile-time; other shapes get their own page type from the same `GridSpec` (docs/TERRAIN.md) |
| 2026-10-20 | `HeightPage::heightAt` after the grid change (pages built on `GridSpec`, the lookup untouched) | Pair `page height x1024`: 0.73x to 0.90x of the hand-written kernel over three runs (3.2 to 3.4 µs vs 3.7 to 4.4 µs) | None; not worse than before, and the open-gap row below stands as the conservative number |
| 2026-10-20 | `TravelBuffer::respond` at 1.5x to 1.8x: `std::fmax` compiles to a library call on MinGW GCC 16 (two `call fmax` per sample in the disassembly), while the reference's `std::max` is one `maxsd` | Disassembly of both loops; pair `travel buffer x1024` | `std::max` (same result for NaN and signed zero here): 1.02x to 1.04x (2.83 vs 2.73 µs min, 3.04 vs 2.98 µs median, 7 interleaved rounds). A precomputed reciprocal width measured no gain; not kept |
| 2026-10-20 | The M1b law response (5 common fields plus extras) pushed a station's `evaluate` past GCC's early-inlining budget ("growth 15 exceeds early-inlining-insns"): `Rk4 step, contact oscillator` 43 to 52 ns | `-fopt-info-inline-missed`, 15 interleaved rounds | `ExplicitRungeKutta::stepBetween` is flattened (GCC/Clang): 43.4 vs 43.0 ns min, 47.0 vs 46.2 median against the pre-change build; Midpoint 23.6 vs 22.5 min, 24.4 vs 24.1 median. Supersedes the 2026-10-15 "force-inlining evaluate: no gain" note, measured when GCC kept it out of line on both sides |
| 2026-10-20 | Revolute joints about z paid the general Rodrigues products | Pairs `planar tree 10-DOF ...` | Short branch for an axis of (0, 0, ±1): about 4% on mass matrix and forward dynamics; one pose pipeline (D6) does not touch the dynamics path |
| 2026-10-20 | Wave 5 (D12): `PreparedProfile` and the chain solver moved here. The recorded clearance gap (1.9x to 2.0x) did not reproduce | GCC 16, MinGW, five interleaved rounds each: unprepared `measureCircleClearance` 0.89x to 0.94x of the twin before and 0.94x to 0.96x after (unchanged code; the earlier gap predates Sub0Math's force-inlined `closestPoint`); prepared overload 0.85x to 0.90x (79 to 82 vs 90 to 93 us per 256 queries); chain solve 0.89x to 0.95x | Prepared clearance kept: window from the centre segment, four minima. A reciprocal squared length was another 10% faster (74 vs 82 us) but changes the last bit of distance and of the projection parameter, so it was not taken (docs/TERRAIN.md Exactness). The clearance gap row is removed; no gap remains |

| 2026-10-20 | M1b: the fixed-size RK step lost about 4 ns to its twins (`Rk4 step` 1.04x to 1.19x, `Midpoint step` 1.29x to 1.32x by minimum, GCC 16): the stages were members, so each step stored the endpoint and its derivative as two scalars and copied them on as 16-byte loads, which cannot be store-forwarded | Disassembly of the flattened step (`movsd` pairs then `movdqu` of the same 16 bytes); pairs over 7 interleaved rounds | The fixed size keeps its stages local to the step and commits from them; the dynamic size is unchanged. `Rk4 step` 0.97x (26.1 vs 26.7 ns min), `Midpoint step` 0.96x to 1.01x (12.8 vs 12.7 ns); 31.4 to 26.1 and 17.1 to 12.8 ns absolute. Removing the cache-state compare measured no change, so the cache contract stays **Reverted the same day:** with local stages, the consumer stations (an expensive `evaluate` per stage) ran 1.7x slower (a loaded-wheel adaptive drop went from 116 to 202 us, A/B with only this header swapped), so the member stages are back and the oscillator pairs are an open gap again |
| 2026-10-20 | M1b element evaluation (the 10-DOF tree's cart spring and branch-coupling strut, Jacobians included) first measured 2.6x its twin (281 vs 108 ns min) | Pair `planar tree 10-DOF elements with Jacobians` | The site-Hessian term as a triple product with a precomputed `u x w` per joint, and the tendon's support without prismatic joints that carry both ends (they change no entry): 190 ns, 2.0x. The rest is the 3D joint transforms (see Open gaps) |
| 2026-10-20 | Element Jacobians as dense `n x n` sums, or per-element blocks | One-off measurement, a 40-DOF tree (cart and three 13-link branches) with two branch-coupling tendons: zeroing the two 40 x 40 matrices 87 to 94 ns, `evaluateElements` with Jacobians 3.1 us, forward dynamics 4.1 to 4.3 us, the implicit step 8.8 to 9.9 us | Dense sums: the clearing is 1% of the step, an element still writes only its support, and the step needs the dense matrix anyway; blocks would add a scatter pass. Decision recorded in docs/MULTIBODY.md (M1b status) |
| 2026-10-20 | M1c constrained step (slider-crank, closure) first measured 6x to 9x its twin; VTune: the small pivoted solves and the per-row tree solves led | Pair `slider-crank constrained ...`, VTune hotspots on a probe | Closure rows that no coordinate moves dropped at assembly (planar out-of-plane rows), one solve for equality-only sets, the projection's first iteration reusing the energy measurement's kinematics, an index-permuted pivoted LDL^T with one reciprocal per pivot. The refactor of the implicit step (shared assembly with the constrained path) left its pair unchanged: 1.92x median against 1.98x for the base build, 5 interleaved rounds on a loaded host |
| 2026-10-21 | G3 extraction. `extractHeights` first measured 2.5x its hand-written twin (97 vs 38 us for 128 x 64 cells): the page, mask and neighbour lookups were per node and per cell, and triangles were written with a capacity-checked push | Pair `height extract 128x64 cells`; a phase split (nodes alone 13 us, cells 42 us) | Nodes filled in runs that share a page with the row and bit words resolved once; cells in page runs with the page presence decided once; four `nodeKnown` bytes combined with `&`; triangles, normals and ids written through raw pointers (capacity checked once for the block). 1.0x (50 vs 50 us, three interleaved runs) |
| 2026-10-21 | `HeightPage::heightAt` per-call gap, at the batch API | Pairs `page height x1024` (single sample, 5 interleaved runs: 0.95x to 1.08x, median 1.06x) and `page heights batch x1024` (1.00x, same runs) | The lookup is one shared inline kernel (`lookup`) with the arrays hoisted by the batch; the batch API, which extraction and any many-sample caller uses, is at parity. The single-sample gap is accepted at about 1.06x |
| 2026-10-21 | `extractColumns` and the tile cache against hand-written twins | Pairs `column extract 336 columns` (1.00x: 1.49 vs 1.50 us) and `tile cache x1024 accesses` (0.57x: 4.7 vs 8.2 us) | Columns: page runs and indexed writes, as the R0a port found (PERF note in the code). Cache: the twin is R0a's array LRU (linear find, rotate); the library's hash map and intrusive list win at 32 tiles and the gap widens with the limit |
| 2026-10-21 | Ray cast against an 8192-cell patch | `castRay` 1.3 us against 45 us for testing every triangle (no hand-written twin; the brute force is the correctness reference) | Cell walk (Amanatides and Woo) over delivered cells only |
| 2026-10-21 | Resolution lever cost | `extractHeights` at stride 4: 8.3 us against 50 us at stride 1 for the same 128 x 64 block | Cost follows the cell count (s^2); see docs/TERRAIN.md Levers for the error evidence |

| 2026-10-21 | M1d: the planar pairs ran the general 3D formulation (the recorded gaps below were "the price of generality") | Pairs `planar tree 10-DOF ...`, GCC 16, three interleaved rounds against the previous build, medians | `Formulation::Planar` (three-component algebra, chosen at build): mass matrix 2.5-3.0x to 1.16-1.22x, forward dynamics 2.5-2.7x to 1.20-1.24x, implicit step 1.9-2.1x to 1.1-1.5x, elements 2.0-2.4x to 1.5-1.6x. Clang 22 `-march=x86-64-v3`: 0.87-0.90x, 1.34-1.37x, 1.50-1.56x, 1.8-2.1x. The constrained step stays at about 6.8x (its constraint machinery dominates). LTDL unchanged code |
| 2026-10-21 | M1d consumer: a bicycle rig's whole run (48 000 steps, 5 DOF, site forces, a tendon, buffers, loads every 80 steps) first measured 3.0x its hand-written twin (22.6 vs 7.4 ms) | The consumer's paired benchmark; VTune hotspots | The rig's observation and step each ran their own kinematics, and the energy a third: the two-phase step (`prepareLinearlyImplicitStep`, prepared site queries) shares one pass, as the twin does; site rotations from the half angle instead of atan2, cos and sin; a position-and-velocity site query; site forces mapped on the step's transforms; plane vectors built by a branch instead of indexed stores (store forwarding); K skipped when not implicit. Now about 2.0x (16-17 vs 8.3 ms, GCC; Clang 2.4-2.6x). The consumer's station pairs are unchanged (four interleaved rounds against the previous Sub0Physics) |
| 2026-10-27 | M2b ring contacts first measured 2.6x (GCC 16) to 3.5x (Clang 22) their twin, although only 3 (edge, piece) pairs survive the culls: the cost was per edge, not in the contact kernel | Pair `ring contacts 64 edges`; disassembly of the per-edge call: the edge was built as a `Capsule` with a 256-bit store and read back as scalars (a failed store forward per edge), and the edge routine stayed out of line | Edge coordinates passed as scalars, the per-edge routine force-inlined, input validation once per ring (exponent-bit test), and a whole-edge cull against the highest piece (the twin was given the same cull, so it stays the fastest version): GCC about 1.1x, Clang about 1.3x. The clearance bound fell from about 1.4 us to 0.5 us with the same changes |

| 2026-10-27 | M2a: one implicit Euler step of two 64-node rings on a five-coordinate tree | Pair `two 64-node rings implicit Euler step`, five interleaved rounds, loaded host: GCC 16 library 10.3-14.1 us, twin 8.0-9.1 us (about 1.3x); Clang 22 `-march=x86-64-v3` 11.2-17.6 against 8.2-12.7 us (about 1.35x). The expanded bodies (dense `LDL^T` of n = 261, plus dense element accumulation) 0.8-1.4 ms | Condensation adopted (about 100x the dense path). The consumer's rig and station pairs against the previous Sub0Physics, four interleaved rounds: rig run 15.2-17.6 ms before, 15.5-19.1 ms after (ratios to its twin 1.6-2.1x both), stations within noise: no regression |
| 2026-10-27 | M2c: the reference's step matrix on a ring is block-diagonal (`DampingOnlyMassPattern` keeps only each node's own damping block; the mass block is always diagonal), yet both condensations ran the periodic block sweep with zero neighbour blocks | The consumer's moving-wheel pair against its hand-written twin, which condenses with 2 x 2 inverses | `condenseRings` detects a block-diagonal ring block and uses the inverse pivots alone (bit for bit: the ring, contact and consumer digests are unchanged). The consumer's rig and station pairs against its main build, interleaved: within noise **A/B, not a verdict (2026-10-27):** the library's uniform per-row matrix solve and the twin's O(k^2) Schur form are both kept as candidate paths. The matrix form may map better onto SIMD or batched (many-bike, many-contact) hardware. Each gets three optimisation passes, measured on the consumer scenario and on a batched shape, before either is chosen per target (Sub0Math PERFORMANCE "Choosing between alternatives") |

## Delassus A/B (M2d, 2026-10-28)

The soft-constraint solve's `J M^-1 J^T` in two forms behind `SoftConstraintSolver::solve(..., supports, DelassusForm)`
(docs/WHEEL.md, "M2d status"), each given three optimisation passes (Sub0Math PERFORMANCE, "Choosing between
alternatives"). Measured on the library solve (one 64-node ring on the five-coordinate rig tree, k contact rows of a
patch, all pushing; `soft solve 64-node ring, k rows, <form>` in `tests/bench_delassus.cpp`), on a batched shape (both
forms as structure-of-arrays kernels over W wheels with the wheel index innermost, each wheel's patch at its own nodes;
`Delassus W wheels x 8 rows SoA, <form>`), and on the consumer's moving-wheel step (BikeDynamics docs/PERFORMANCE.md).
CrogLegion shared with other agents; medians with min-max of interleaved rounds.

**Passes.**

| Pass | Uniform (full solve per row; bit for bit the reference's arithmetic) | Structured (Schur, O(k^2)) |
|---|---|---|
| 0 | `solveMassMatrix` per row, dense products | `t = j_t - J_h^T sum B^T A^-1 j_w`, a full tree solve per row, shared-node terms |
| 1 | Products `J Y`, `J a0`, `J^T f` over each row's support (exact zeros skipped, same order) | Half a tree solve: `t^T S^-1 t = w^T D^-1 w` after one sweep of the factor |
| 2 | All rows in one coordinate-major sweep (`solveMassMatrixColumns`), rows sharing SIMD lanes | `J a0` and `J^T f` straight from the support, no index lists |
| 3 | The sweep's ring elimination skips nodes no row touches | Shared nodes by integer keys |

Library solve, GCC 16, ns at 4 / 8 / 16 rows: pass 0 uniform 1670 / 4430 / 14 730, structured 575 / 1500 / 6170; after
pass 1 about 1080 / 2700 / 8800 and 290 / 990 / 5300; pass 2 (uniform only) no gain under GCC's SSE2; after pass 3 810
/ 1930 / 7470 and 284 / 930 / 5390 (5 rounds). The pivoting (one Cholesky of the active set per pivot) dominates both
at 16 rows.

**After three passes, per target** (3 interleaved rounds; ns, 4 / 8 / 16 rows):

| Target | Uniform | Structured | Structured advantage |
|---|---|---|---|
| GCC 16 | 1300 / 2900 / 6800 | 460 / 1160 / 4600 | 2.8x / 2.5x / 1.5x |
| GCC 16 `-march=x86-64-v3` | 720 / 1530 / 6680 | 277 / 786 / 5240 | 2.6x / 1.9x / 1.3x |
| Clang 22 `-march=x86-64-v3` | 812 / 1415 / 4840 | 265 / 774 / 3520 | 3.1x / 1.8x / 1.4x |

**Batched shape** (ns for the batch; 8 rows a wheel; 5 rounds): the uniform kernel gained nothing from accumulating `D`
during its backward sweep (its pass 1), so its own passes are the library's; the order-of-magnitude asymptotic difference
(O(k N) against O(k^2) at N = 64) is the rule's exception, so the batched decision does not wait on further passes.

| Target | 1 wheel: uniform, structured | 8 wheels | 64 wheels | Per wheel at 64: uniform, structured |
|---|---|---|---|---|
| GCC 16 | 2450, 161 | 14 700, 1170 | 148 000, 13 300 | 2.3 us, 0.21 us (11x) |
| GCC 16 `-march=x86-64-v3` | 2450, 164 | 11 300, 1150 | 99 900, 13 100 | 1.6 us, 0.21 us (7.6x) |
| Clang 22 `-march=x86-64-v3` | 2130, 194 | 18 500, 1550 | 99 000, 14 700 | 1.5 us, 0.23 us (6.7x) |

AVX2 helps the uniform kernel most (2.3 to 1.6 us a wheel at 64 wheels): it is fixed-shape and streams. The structured
kernel gathers each wheel's two nodes per lane, and got slower when the wheels' patches stopped sharing nodes (9 to 13
us at 64 wheels), yet stays 7 to 11 times ahead.

**Decision per target:** the structured form on every target measured, scalar and batched. The uniform form is not
deleted: it is the reference path's arithmetic (its digest depends on it), it handles rows without ring supports, and it
is the candidate if a target with wide fixed-shape throughput (GPU) or much smaller rings changes the balance. On the
consumer's whole step the difference is small (BikeDynamics: 1.1x to 1.2x on flat ground, within noise on the rock),
because the forward pass's other work dominates.

## Open gaps

Deficiencies against the hand-written reference (rule: Sub0Math `docs/PERFORMANCE.md`).

| Since | Pair | Gap | Cause | Follow-up |
|---|---|---|---|---|
| 2026-10-27 | `ring contacts 64 edges, square rock` | GCC 16 0.93x to 1.19x, median 1.09x (119-143 vs 130-165 ns); Clang 22 1.19x to 1.36x, median 1.31x (129-138 vs 154-188 ns); `-march=x86-64-v3`, 5 interleaved rounds, loaded host | The library writes a fuller record (axis parameter, capsule point, crossing flag, quantities; 112 bytes against the twin's 72), keeps found/measured counts and the smallest gap, and validates its inputs | Re-measure when M2c binds contacts to elements (the record is what the element needs); accept or trim then |
| 2026-10-27 | The consumer's moving-wheel reference step (one 64-node ring, soft-constraint contact rows) against its hand-written twin | about 2.9x to 3.2x after M2d (GCC 16 5.45 against 1.90 us a step on flat ground; GCC x86-64-v3 5.3 against 1.6; Clang 22 7.0 against 2.4); was 2.6x to 4.5x | M2d: the Delassus A/B (above) took the per-row solve's cost down on both forms; the reference keeps the uniform form (its bits). What remains: the forward pass's energy and acceleration observation every state, the general tree and element code | A forward pass without observation for consumers that do not record every state; the reference may switch to the structured form only with a re-pinned digest and its L2 bounds (an owner decision) |
| 2026-10-27 | `two 64-node rings implicit Euler step` | about 1.3x (GCC 16), 1.35x (Clang 22) | The library evaluates the laws through `respond` (potential and the full response per coordinate), the ring's right-hand side in a separate pass, generic periodic-solve loops over `std::array` blocks and the general tree code; the twin fuses these | Profile when contact (M2c) is in the step; the step is far inside the frame budget |
| 2026-10-21 | `page height x1024`, `HeightPage::heightAt` (single sample) | about 1.06x median (0.95x to 1.08x over 5 interleaved runs); the batch `heightsAt` pair is at 1.00x | Per-call conversion and bounds work around an otherwise identical kernel | Accepted with these numbers: callers with many samples use `heightsAt` (G3), which extraction does |
| 2026-10-20 | `planar tree 10-DOF implicit Euler step, 2 elements` | 1.1x to 1.5x since M1d (GCC 16, three rounds; was 2.0x); Clang 22 `-march=x86-64-v3` 1.5x | The planar formulation closed the transform cost; the rest is the generic element loop and the forcing/assembly checks per step | The consumer's whole rig run is the judge (about 2.0x, below); next: the elements (next row) |
| 2026-10-20 | `slider-crank constrained implicit Euler step` | about 6.8x after M1d (unchanged: the planar formulation does not reach the constraint machinery), earlier about 6.5x by median (5.0x to 6.9x per round, 5 interleaved rounds, GCC 16, loaded host: 680 to 770 vs 100 to 155 ns); Clang 22 `-march=x86-64-v3` 6.0x (520 vs 86 ns) | The unconstrained step is already about 117 ns of it (the CRBA/RNEA gap below on a 3-DOF tree); the rest is the general machinery for a few rows: a full tree solve per row for `J A^-1 J^T`, a pivoted LDL^T that detects dependent rows, 3D site Jacobians, and four kinematics passes against the twin's closed forms | With M1d, on the rig's whole step (closures on a 10-20 DOF tree, where the fixed costs amortise); first candidates in docs/MULTIBODY.md (M1c post-parity ideas: Delassus from one triangular sweep). 0.7 us per step is far inside a 0.5 ms step |
| 2026-10-20 | `planar tree 10-DOF elements with Jacobians` | 1.5-1.6x since M1d (GCC 16; was 2.0x); Clang 1.8-2.1x | The tendon geometry now runs in planar algebra, but through the general element loop: dense `n x n` clearing of K and D, 3-vector storage of planar columns, the law dispatch | A planar column type and support-only clearing; measure on a consumer with more tendons before changing the dense-sum decision (M1b) |
| 2026-10-20 | `Rk4 step` and `Midpoint step, contact oscillator` (GCC 16) | 1.12x and 1.27x (7 interleaved rounds) | Member stages: endpoint stores that cannot be store-forwarded (Findings, 2026-10-20). Local stages close it on the oscillator, but make stations whose `evaluate` is large 1.7x slower; the cause of that is not yet found | Find a form that is at parity on both the twin and a station-sized system (add a station-sized twin), then re-apply; the stations are the consumer that matters |
| 2026-10-20 | `Rk4 step` and `Midpoint step, contact oscillator` under Clang 22 (`-march=x86-64-v3`, Windows) | 1.2x and 1.35x to 1.4x (GCC 16 is at parity, above) | Clang does not inline the flattened step's nested calls the way GCC does (`stepBetween` or `advance` stay out of line); force-inlining them measured no clear gain, so it was not kept | Re-measure on the Linux Clang CI leg; GCC is the host baseline's compiler |
| 2026-10-19 | `planar tree 10-DOF`: `computeMassMatrix`, `forwardDynamics` | 1.16-1.22x and 1.20-1.24x since M1d (GCC 16; were 2.9x and 2.65x); Clang 0.9x and 1.35x | The planar formulation runs the same algorithms as the twin; what remains is the run-time tree and joint-type branches against the twin's compile-time ones | Accepted unless a consumer shows it; the rig run (2.0x) is dominated by other work |
| 2026-10-21 | The consumer's bicycle rig whole run against its hand-written twin | about 2.0x (GCC 16), 2.4-2.6x (Clang 22) | The twin computes the mass matrix and bias once per state and reuses them for the energy, the step and the loads, with closed-form kinematics and tyre Jacobian rows; the library still recomputes the forward dynamics for the loads on a second workspace, maps site forces through the tree, and pays the generic element loop and the tree factorisation's index chasing (VTune: factor and solve about 15%, site queries and checks about 10%) | Reuse of a prepared state for reactions; a dense factor below a size threshold |

## Wheel-model study (2026-10-28)

The consumer's explicit nodal ring (its FIDELITY_LEVERS 4.11, PERFORMANCE "Findings (wheel-model study)") evaluates a
64-node casing in about 2 to 3 us with no factorisation, against about 20 us for the implicit step; at the Interactive
setting a moving-wheel run is 20 to 90 times cheaper than the improved model. The functions added here
(`casing_models.hpp`) are scalar closed forms called per contact or per switch; they have no measurable cost of their own.

