# Multibody core (design, step M1 of ADR 0003)

- **Status:** Design for implementation (2026-10-19)
- **Decision:** [ADR 0003](https://github.com/CraigHutchinson/.github/blob/main/docs/decisions/0003-dynamics-engine.md).
  Sub0Physics owns the dynamics core; MuJoCo 3.16.0 is a test-only reference.

This document fixes the contracts that the implementation steps (M1a–M1d below) build to, so they can
be delegated and reviewed against one specification.

## 1. Scope

The core models articulated bodies with reduced coordinates:

- a tree of rigid bodies;
- force elements that may be stiff and supply their own Jacobians;
- a few constraints: loop closures, strict unilateral limits and attachments.

The first consumer is the bicycle rig: a planar hardtail or full-suspension frame, fork, swing arm and
guided rider. Then come the two spatial wheels and the connected trail bike (M2–M3), and later 3D (M4).

It is sized for 10–40 degrees of freedom. Dense matrices are the right choice at that size: they are
simpler, deterministic and fast. Nothing here is a general contact-rich robotics engine. Tyre and soil
contact are **force elements** computed by our own terrain code. They are not constraints solved by
the core.

## 2. Model description (ours, not MJCF)

A `ModelBuilder` takes typed descriptions (Sub0Math quantities and roles) and compiles them into an
immutable `Model` of flat arrays in topological order. The description is the asset that must not be
locked in: a bicycle-rig consumer builds its models with it, and the test-only MuJoCo adapter compiles the
same description to `mjSpec`.

| Item | Contents |
|---|---|
| `Body` | Name, mass, centre of mass, inertia about the COM (positive definite), parent joint |
| `Joint` | Type (`Prismatic`, `Revolute`; later `Planar`, `Spherical`, `Free`), parent body, frame in the parent, axis, optional spring reference; generalized coordinates are ordered by the tree |
| `Site` | A named frame on a body, for attachments, element endpoints and load reporting |
| `ForceElement` | One of the element kinds in §4, attached to coordinates or sites |
| `Constraint` | Loop closure between two sites (point or frame equality), or a unilateral limit on a coordinate or distance |
| `Gravity` | A vector, typed `Acceleration` |

The builder validates all of this and returns `std::expected<Model, Fault>`: tree connectivity,
positive mass and inertia, finite values, unique names, and every reference resolved. There is no XML
and no runtime parsing.

## 3. State and dynamics

- **State** holds `q` and `qd` (double), the tick, per-element state (histories, buffer state) and
  constraint state (active set, warm start). The whole state is a value that can be snapshotted, which
  is what replay and restart need.
- **Equations:** `M(q) qdd + c(q, qd) = tau(q, qd, t) + J(q)^T lambda`.
  - `M` comes from CRBA (dense, n ≤ ~40), and `c` (Coriolis, centrifugal, gravity) from RNEA.
  - `M` is factorised by Featherstone's LTDL (`M = LᵀDL`, eliminating from the leaves): the ordering of
    `LDLᵀ` that creates no fill-in on a tree, so the cost follows the tree's depth. Landed in M1a.
  - ABA only if measurement ever asks for it.
- **Spatial algebra:** spatial vectors, spatial inertias and Plücker transforms, in double, through a
  6D vocabulary added to Sub0Math when M1a first needs it, following ADR 0002's "grow by use".
  Planar models use the same 3D formulation with planar joints; there is no separate 2D engine.

## 4. Force elements: the extension point

Each element kind is a type with this contract (a C++ concept, not a virtual hierarchy, so hot loops
inline):

```cpp
struct ElementContribution            // written into preallocated workspace
{
    // generalized force on the coordinates the element touches (sparse: indices + values)
    // optional Jacobians d tau / d q and d tau / d qd on the same index set (dense small block)
    // potential energy, dissipated power, external power
    // wrench at each attached site (for reactions and loads)
};

template<typename E>
concept ForceElement = requires(const E e, const ElementInputs& in, ElementContribution& out) {
    { e.evaluate(in, out) } noexcept;          // never allocates
    { e.hasJacobians() } -> std::same_as<bool>;
};
```

Inputs are the state of the touched coordinates or the site kinematics: positions, velocities,
Jacobians.

**First element kinds:**

| Kind | Notes |
|---|---|
| Coordinate spring-damper | On one joint coordinate; linear, with a reference |
| `TravelBuffer` | Progressive buffer on a coordinate; ported from BikeLab with exact analytic Jacobians. **Landed** as a stand-alone law (`element/travel_buffer.hpp`: force, potential, loss power and both Jacobian entries, built by the validating `TravelBuffer::create`); binding it to a coordinate as a `ForceElement` is M1b |
| Spatial tendon spring-damper | Two-site length, like BikeLab's shock: force along the line, Jacobian through the site Jacobians |
| External wrench | Prescribed or scripted, at a site |
| Tyre and contact (M2) | Over terrain extracts, with Jacobians from the contact law |

**Why Jacobians matter.** Stiff elements (tyre casing, soil, buffers) enter the implicit integrator
through their own Jacobians, so they no longer force tiny explicit steps. That is the main capability
MuJoCo could not give us: its custom forces miss its implicit derivatives.

## 5. Integration

- **Default: linearly implicit Euler on velocities**, including the element Jacobians:
  `(M − h·D − h²·K) Δqd = h·(tau − c + Jᵀλ)`, where `D = Σ dtau/dqd` and `K = Σ dtau/dq` over the
  elements that supply them. Then `qd += Δqd` and `q` is advanced on the manifold for each joint
  type. This is the structure of MuJoCo's `implicitfast`, extended with our elements' derivatives.
- **Higher-order variants** (implicit midpoint, SDIRK2) are added only when convergence measurement
  asks for them. The existing explicit RK family remains for non-stiff reference runs.
- **Steps run on integer ticks** (`stepTick`, as now). Events are located within a step where required
  (first contact, limit engagement), never stepped over. The step is split at the event.
- **Fidelity selection** (`FidelitySelector`) chooses between element models and integrators per
  mechanism, as today.

## 6. Constraints

Constraints here are few by design:

- **Loop closures:** solved at velocity level with the factorised `M`, using a small dense
  `J M⁻¹ Jᵀ` solve (a handful of rows), followed by position projection onto the constraint manifold
  each accepted step. No Baumgarte tuning. Residuals are recorded and gated.
- **Strict unilateral limits:** travel ends, ultimate stops and rider stops. They are solved as a
  small LCP by direct pivoting (exact for a few rows), with event location at engagement. They are
  never overshot silently: an accepted state inside limits is an invariant, and a violation is a fault
  or a rejected step. Progressive buffers stay force elements; the hard limit behind them is a
  constraint.
- **Attachments with release** (rider contacts, D5): constraints that can deactivate on a force
  threshold, with the transition logged.

## 7. Observability by construction

- **Joint reaction wrenches** come from a post-solve RNEA that includes every force: elements,
  constraint impulses and gravity. There are no omitted categories, unlike MuJoCo's post-constraint
  RNE, which skips flex contacts.
- **Each element reports the wrench at its sites.** Loads at any named site come from the same
  numbers that moved the bodies.
- **The energy ledger** (`EnergyLedger`) holds kinetic, gravitational and element potential energy,
  dissipated energy and external work, each as a channel, closed over every accepted step within a
  declared budget.

## 8. Determinism, precision, memory and performance

- **Determinism:**
  - fixed evaluation order;
  - no parallel reductions inside a step;
  - floating-point contraction disabled for the core's translation units
    (`-ffp-contract=off`, `/fp:precise`);
  - one set of transcendental functions for state evolution: our own, or audited `std::` ones,
    decided at M3, where cross-platform hashes are the gate.
- **Precision:** state and dynamics in double. Terrain data stays float and is converted at the
  contact boundary (terrain design §3).
- **Memory:** the `Model` is immutable and holds flat arrays; a `Workspace` is preallocated per model.
  A step never allocates, and an allocation test gates this.
- **Performance:** paired benchmarks against hand-written kernels for a fixed small tree (CRBA, RNEA,
  the `LDLᵀ` solve, a whole step). The rules in Sub0Math `docs/PERFORMANCE.md` apply: the reference
  is never weakened, gaps are closed or recorded, and iterated code carries `// PERF:` notes.

## 9. Delivery

| Step | Delivers | Gate |
|---|---|---|
| M1a | Spatial algebra (Sub0Math); `ModelBuilder`, `Model`, `State`, prismatic and revolute joints; CRBA, RNEA, LTDL factorisation; explicit stepping | Pendulum, double pendulum and slider analytic checks; energy conservation with no dissipation; allocation gate; benchmarks |
| M1b | `ForceElement` concept; coordinate spring-damper, tendon spring-damper, `TravelBuffer`; linearly implicit integrator with element Jacobians | Stiff spring-damper stable at large `h`, with energy and dissipation accounted; Jacobians checked against finite differences; convergence order measured |
| M1c | Loop closures, strict limits (LCP), event location, post-solve reactions, energy ledger | Slider-crank closure residuals; a limit never overshot; reactions match independent statics; the energy ledger closes |
| M1d | The consumer's bike-rig description (hardtail and full suspension, buffers, rider guide) on the core | Parity against the M0 MuJoCo reference and BikeLab goldens within the declared tolerances; timestep advantage measured |

Each step ships through the pre-push gate, with goldens and provenance as usual.


### M1a status (2026-10-19)

**Landed.** Sub0Math `linear/`: `SpatialMotion`, `SpatialForce` (`dot`, `cross` as crm/crf),
`SpatialInertia`, `PluckerTransform`. Sub0Physics `multibody/`: the description (`BodyDescription`,
`JointDescription`, `SiteDescription`, `RotationalInertia`, `BodyId`, `SiteId`), `ModelBuilder` →
`std::expected<Model, Fault>`, `Model`, `State`, `Workspace`; kinematics (`updateKinematics`, site pose,
velocity and translational Jacobian, `computeEnergy`); `computeMassMatrix` (CRBA),
`computeBiasForces` and `inverseDynamics` (RNEA), `factorMassMatrix` / `solveMassMatrix`,
`forwardDynamics`; `stepSymplecticEuler` (a semi-implicit Euler step). The analytic gates, allocation gate and paired benchmarks
are in `tests/test_multibody*.cpp` and `tests/bench_multibody.cpp`.

**Deviations and decisions:**

- The factorisation is Featherstone's LTDL, `M = Lᵀ D L`, eliminating from the leaves: the ordering of
  `LDLᵀ` that creates no fill-in on a tree. It is used for every size; no dense fallback was needed
  (parity with a hand-written version, 69 ns at 10 DOF).
- A body's frame is its joint's frame; centre of mass, inertia and sites are given in it. One
  coordinate per joint, so `dofCount() == bodyCount()` until multi-DOF joints arrive.
- Spatial vectors hold plain SI doubles (like `Mat`); quantities and roles sit at the description.
- The explicit stepping is semi-implicit Euler. With a configuration-dependent `M` it is not
  symplectic, so energy is bounded to first order rather than conserved: about 1e-3 of `m g l` at
  `h` = 0.1 ms on the double pendulum, halving with `h`. Conservation to round-off is shown with the
  existing `Rk4` through an `OdeSystem` adapter in the tests (9e-11 over 20 s at 1 ms).
- `Workspace` exposes its arrays as fields, sized once (the Pinocchio `Data` pattern); the site and
  energy queries recompute kinematics per call (O(n)). A workspace is bound to the model it was built for
  (`Model::id`): every public entry point returns `InvalidInput` for another model's workspace or a state
  that is not sized for the model, and the queries that used to return plain values (`computeMassMatrix`,
  `computeBiasForces`, `updateKinematics`, `computeSiteKinematics`, `computeEnergy`, `solveMassMatrix`) return
  `std::expected`. The mass matrix is zero outside ancestor pairs from construction and only ancestor entries
  are ever written.
- `ModelBuilder` also rejects inertias whose principal moments break the triangle inequality, as MuJoCo does.
- The multibody sources build with `-ffp-contract=off` (`/fp:precise` on MSVC), per section 8.
- Gaps against hand-written planar code (CRBA 2.9x, forward dynamics 2.8x) are recorded in
  `docs/PERFORMANCE.md` with an M1d follow-up.

## 10. Open questions, to be answered by measurement

- **Higher-order implicit scheme:** whether one is needed, and which, for the tyre ring (M2).
- **LCP method:** direct pivoting against projected Gauss-Seidel once contacts or limits outnumber a
  handful (expected never for one bike).
- **Transcendentals:** whether `std::sin`/`std::cos` are good enough for cross-platform hashes, or
  whether we need our own (M3 gate).
