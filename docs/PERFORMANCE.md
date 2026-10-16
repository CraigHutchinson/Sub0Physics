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
| Evaluation count | `test_integrators.cpp` | Tick stepping costs exactly one evaluation per stage: 4 per RK4 step after the first |
| Convergence order | `test_integrators.cpp` | RK4 observes order 4, midpoint order 2; a slower method that is "faster" is not a win |
| Header self-containment | `Sub0Physics_HeaderCheck` | Every public header compiles alone |

Replacement `operator new` lives in its own test executable so it cannot affect the others.

## Timing budgets (per host)

```bash
cmake --preset profile
cmake --build --preset profile --target Sub0Physics_Bench_check      # fails if >10% slower than baseline
cmake --build --preset profile --target Sub0Physics_Bench_baseline   # record this host's baseline
```

Baselines live in `perf/baselines/<host>.json` (median ns/op from nanobench). Timings from one machine
say nothing about another, so each host has its own file; a host without one passes and is told how
to adopt the gate. Regenerate a baseline deliberately and say why in the commit. CI records
`github-ubuntu` with a 25% tolerance and reports drift in the job summary without blocking, because
shared runners are noisy.

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
| 2026-10-15 | `std::hypot` dominated dense clearance queries (11.7 µs on a 1024-segment arc in BikeDynamics) | BikeDynamics bench | Squared distances, `sqrt` only on improvement, coordinates bounded at 1e100 m; 3.9 µs |
| 2026-10-15 | Moving the geometry interface to Sub0Math quantities regressed clearance by 52% (395 to 601 ns, best of 7) | Budget gate, then a scalar rewrite of the kernel to locate it | Fixed in Sub0Math (compile-time component expansion, no zero-initialised sums); 394 vs 402 ns, within noise |
| 2026-10-16 | Terrain page height lookup started 2.5x slower than the hand-written kernel (raw arrays, bit mask, fixed diagonal) | Paired benchmark `page height x1024` | Inline definition with reciprocal spacing (7.1 to 5.4 µs); a plain `HeightSample` return instead of `std::optional` (5.6 to 3.5 µs; `std::optional<float>` was equally slow, so the optional, not the quantity, kept the result in memory on GCC 16); two plain branches instead of a merged triangle select (3.5 to 3.0 µs). Force-inlining alone: no gain. Now 1.07x |

## Open gaps

Deficiencies against the hand-written reference (rule: Sub0Math `docs/PERFORMANCE.md`).

| Since | Pair | Gap | Cause | Follow-up |
|---|---|---|---|---|
| 2026-10-16 | `page height x1024`, `HeightPage::heightAt` | about 1.07x in one process, 1.10x interleaved (pair ratio) | Per-call conversion and bounds work around an otherwise identical kernel | Terrain step G3: the batch extraction API (many samples per call, SoA out) is the real consumer; close it there or accept it with its numbers |
