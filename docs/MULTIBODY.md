# Multibody core (design, step M1 of ADR 0003)

- **Status:** M1a to M1d landed (M1d 2026-10-21); M2 designed in [WHEEL.md](WHEEL.md); M2a (ring substructure,
  condensation, expanded oracle) landed 2026-10-27, see [WHEEL.md, M2a status](WHEEL.md#m2a-status-2026-10-27); M2c (contact
  elements, the reference soft-constraint law, `ImplicitTerms::DampingOnlyMassPattern`) landed 2026-10-27, see
  [M2c status](#m2c-status-2026-10-27)
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
| `Joint` | Type (`Prismatic`, `Revolute`; later `Spherical`, `Free`), parent body, frame in the parent, axis, optional spring reference; generalized coordinates are ordered by the tree. A body with several joints is a chain with massless outer bodies (M1d). A planar mechanism needs no planar joint type: the model is evaluated in planar algebra when it qualifies (M1d, `Formulation`) |
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
  Planar models use the same model description; when a model qualifies, the same algorithms run in three-component
  planar algebra (M1d, `Formulation::Planar`). There is no separate 2D engine or description.

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
| `TravelBuffer` | Progressive buffer on a coordinate; ported from BikeLab with exact analytic Jacobians. **Landed** as a stand-alone law (`element/travel_buffer.hpp`, a `LawResponse` with per-buffer forces, built by the validating `TravelBuffer::create`); binding it to a coordinate as a `ForceElement` is M1b |
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
  - one set of transcendental functions for state evolution, decided (2026-10-27): our own, `sub0math::det`
    (fdlibm-derived sin, cos, exp, log; `pow` by squaring or exp-log), never the C library. The ring run digest
    matched on every Windows build and differed on Linux (`abd8d2ab7e62fbc7`) and macOS (`f67fc6a4aa623f6d`)
    while libm `sin`/`cos` drove the joints; the cause was the C libraries, not the compiler. Call sites: the
    joint transforms (revolute and planar), the corridor builder, the soft-constraint impedance power and the
    step-size controller's power law. `std::sqrt` stays (IEEE-exact). Not on the state path and still libm: the
    reference damped oscillator, the principal-moment validation (`acos`), `tan` in the berm test surface and the
    convergence-order `log`. A new state-path call to a libm transcendental is a bug; add the function to
    Sub0Math instead.
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
`forwardDynamics`; a semi-implicit Euler step (now `stepSemiImplicitEuler`, see M1b foundations). The analytic gates, allocation gate and paired benchmarks
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
  existing `Rk4` through an `OdeSystem` adapter (9e-11 over 20 s at 1 ms); since M1b foundations that
  adapter is the library's `ModelOde` and `Rk4Stepper`.
- `Workspace` exposes its arrays as fields, sized once (the Pinocchio `Data` pattern); the site and
  energy queries recompute kinematics per call (O(n)). A workspace is bound to the model it was built for
  (`Model::id`): every public entry point returns `InvalidInput` for another model's workspace or a state
  that is not sized for the model, and the queries that used to return plain values (`computeMassMatrix`,
  `computeBiasForces`, `updateKinematics`, `computeSiteKinematics`, `computeEnergy`, `solveMassMatrix`) return
  `std::expected`. The mass matrix is zero outside ancestor pairs from construction and only ancestor entries
  are ever written.
- `ModelBuilder` also rejects inertias whose principal moments break the triangle inequality, as MuJoCo does.
- The whole library builds with `-ffp-contract=off` (`/fp:precise` on MSVC), per section 8, and the option is PUBLIC so consumers compile the inline headers the same way (Sub0Math `sub0math_target_deterministic_fp`, 2026-10-20). A CI leg builds with FMA hardware enabled (`-march=x86-64-v3`) to prove it; without the option, bit-exact tests failed there and on Apple Silicon.
- Gaps against hand-written planar code (CRBA 2.9x, forward dynamics 2.8x) are recorded in
  `docs/PERFORMANCE.md` with an M1d follow-up.

### M1b contracts (2026-10-20)

M1b builds on these contracts, which a cross-repository review fixed before more elements land.

**Laws and elements.**

- A **law** (`element/`) is a scalar constitutive function of one coordinate and its rate.
  - It is a value type built by a validating `static create(const Spec&) -> std::expected<Law, Fault>`.
  - It answers `respond(Coordinate q, Rate qd) -> LawResponse<Coordinate>`:

    ```cpp
    template<typename Coordinate>          // f64::Length (linear) or f64::Angle (rotational)
    struct LawResponse
    {
        GeneralizedForce<Coordinate> force; // acts along +q: Force for Length, Torque for Angle
        decltype(force / q) dForceDPosition; // exact partial derivative
        decltype(force / qd) dForceDVelocity;
        f64::Energy potential;               // stored energy, relative to the law's reference
        f64::Power dissipationPower;         // removed power, never negative
    };
    ```
  - **Sign:** `force` is the generalized force on the law's own coordinate. For conservative and dissipative
    laws, `force · qd = −dV/dt − dissipationPower`, and a test checks this identity for every law. A
    spring on a compression coordinate therefore reports a *negative* force when compressed.
    `TravelBuffer` already follows this rule. `AxialSpringDamper` and `UnilateralSpringDamper`, which report
    a positive "spreading" force, migrate to it, gain their exact Jacobians and rename `storedEnergy` and
    `dissipationRate` to `potential` and `dissipationPower`. Extra diagnostic fields (`topForce`,
    `springForce`) may stay alongside the common ones.
- An **element** (`multibody/`, M1b) binds a law to the model.
  - **Coordinate element:** the law on one joint coordinate, with an optional reference offset.
  - **Tendon element:** the law on the distance `l(q)` between two sites. `tau = Jₗᵀ F`, and the stiffness is
    `Jₗᵀ (dF/dl) Jₗ` plus the geometric term `F · d²l/dq²`. The geometric term is included, because
    otherwise the implicit step loses its stability margin for preloaded tendons.
  - Elements write `tau`, the `dtau/dq` and `dtau/dqd` blocks, energy, and the site wrenches into the
    workspace (§4). They never allocate. The law exists once; elements never copy its formula.

**Energy.**

- One `EnergyAccount {kinetic, potential, dissipated, externalWork}` in `f64::Energy`.
- `computeEnergy` fills `kinetic` and gravitational `potential`. Elements add their `potential`, and the
  stepper integrates `dissipationPower` and external power into `dissipated` and `externalWork`.
- `EnergyLedger` records an `EnergyAccount`, so there are no raw-double argument lists. Its scale policy is
  unchanged. Stations that use it today migrate in the same change.

**Stepping.**

- Every public stepper has the shape
  `step(const Model&, const Forcing&, std::uint64_t tick, f64::Time h, State&, Workspace&) -> std::expected<void, Fault>`.
  - It is **transactional**: on a fault, `State` is unchanged.
  - Time is `tick · h`, derived inside the stepper, never passed separately.
- `State` stores `[q; qd]` in one contiguous buffer and exposes the `q()` and `qd()` spans. Integrators work
  on spans of it.
- The explicit Runge–Kutta family keeps its fixed-size form for small station systems, which are measured
  and perf-tuned. It gains `f64::Time h` at its public boundary and a span-based form that can drive a
  `Model` (the RK reference run moves from the tests into the library).
- The semi-implicit Euler step, the linearly implicit Euler step (§5) and RK all satisfy one `Stepper`
  concept, so `FidelitySelector` can switch between them. The planned fused small-system RK is a
  `Stepper`, not a separate family.

**One pose pipeline.**

- Each joint's motion is evaluated once per call: one sin/cos pair of q (Rodrigues) builds
  `parentToBody`.
- World `bodyPose` is composed from those transforms through Sub0Math `toRigidTransform(PluckerTransform)`.
- The quaternion half-angle path is removed, so kinematics, Jacobians, energy and dynamics all come from
  one representation.

### M1b foundations status (2026-10-20)

**Landed** (the four contracts above; consolidation D2, D4, D5, D6):

- **Laws:** `LawResponse<Coordinate>`, `CoordinateRate`, `GeneralizedForce` and the `Law` concept in
  `element/law_response.hpp`. `AxialSpringDamper`, `UnilateralSpringDamper` and `TravelBuffer` are built by
  validating `create(Spec)` and return responses derived from `LawResponse<f64::Length>`, keeping their
  diagnostic extras (`springForce`/`damperForce`, `compression`/`normalForce`, `topForce`/`bottomForce`).
  A typed test over every law checks `force · qd = −dV/dt − dissipationPower` along trajectories through
  every branch, and both Jacobian entries against one-sided finite differences (at a kink the entry must
  match the side its branch selects).
- **Energy:** `EnergyAccount` (`verification/energy_account.hpp`); `computeEnergy` returns it with
  `kinetic` and gravitational `potential`; `EnergyLedger` is built from an initial account and records
  accounts. `MechanicalEnergy` is gone.
- **Stepping:** `State` is `values = [q; qd]` with `q()`/`qd()` spans (`State::fromCoordinates` for set-up);
  `Forcing{tau}`, the `Stepper` concept, `stepSemiImplicitEuler` and `SemiImplicitEuler` (`stepSymplecticEuler`
  is removed, no alias was needed); `ExplicitRungeKutta` takes `f64::Time` and has a `std::dynamic_extent`
  span form; `ModelOde` and `Rk4Stepper`/`MidpointStepper` (`multibody/runge_kutta_stepper.hpp`) are the RK
  reference run, satisfying `Stepper`.
- **One pose pipeline:** the quaternion half-angle path is gone. `updateKinematics` composes
  `worldToBody` from `parentToBody` and sets `bodyPose = toRigidTransform(worldToBody)`; site velocities,
  Jacobians and energy read the world transforms. `Model` stores each placement once (`treeTransforms()`;
  `jointPlacements()` is removed).
- **Ids (D8):** `BodyId` and `SiteId` are `Id<Tag, std::uint32_t>` (`core/id.hpp`); `WorldBody` is the one
  sentinel, `BodyId::invalid()`. A default-constructed id is now invalid rather than body 0.

**Deviations and decisions:**

- `State` keeps its `tick`; a stepper refuses (`InvalidInput`) a `tick` argument that differs from
  `state.tick`, and sets `state.tick = tick + 1` on success. The argument stays because the contract fixes the
  signature; the check stops the two from disagreeing.
- `Forcing` holds only `tau`, constant over the step. Elements join it with the element binding (rest of M1b).
- `RungeKuttaStepper` is bound to the model it was built for, like a `Workspace`, and keeps a copy of the
  forcing: when `tau` changes between ticks the cached first stage is discarded, so a change of forcing at a
  tick boundary keeps the method's order (tested).
- `EnergyAccount::potential` holds gravitational plus element potential, so no separate elastic channel.
- The fixed-size RK step flattens its call tree (`SUB0PHYSICS_FLATTEN`, GCC/Clang): the richer law
  response pushed a station's `evaluate` past GCC's early-inlining budget (RK4 contact step 43 to 52 ns);
  flattened, it is back at or below the old time.
- Revolute joints about the z axis take a short branch of the joint transform (about 4% on CRBA and forward
  dynamics); the CRBA/forward-dynamics gap is otherwise unchanged (docs/PERFORMANCE.md).

### M1b status (2026-10-20)

**Landed** (the M1b gate, section 9):

- **Elements** (`multibody/force_element.hpp`): `CoordinateElementDescription<Law>` (a law on one joint coordinate,
  optional reference, defaulting to the joint's `springReference`, and a `LawSense`) and
  `TendonElementDescription<Law>` (a law on the distance between two sites, with a reference length and a sense:
  extension or compression). `ModelBuilder::addElement` collects them and `build` compiles and validates them
  (`std::expected`): resolved bodies and sites, a law on a length for a prismatic joint or a tendon and a law on an
  angle for a revolute joint, finite references, unique names. `evaluateElements` writes into the `Workspace`
  without allocating: `elementForce`, `elementStiffness` (dtau/dq) and `elementDamping` (dtau/dqd) as dense
  `n x n` sums, `elementPotential`, `elementDissipationPower` and two `ElementSiteWrench` per tendon. The
  `ForceElement` concept is satisfied by `CoordinateElement<L>` and `TendonElement<L>`; the laws are a closed
  compile-time list (`CoordinateLaws`, `TendonLaws`), so evaluation is a fixed-order loop over typed vectors with no
  virtual dispatch, and the law's formula exists only in the law.
- **Tendon Jacobians, exact:** `dtau/dq = J_l^T (dF/dl) J_l + F d2l/dq2 + J_l^T (dF/dl') (d2l/dq2 qd)^T`, the
  material, geometric and rate-transport terms; `d2l/dq2` comes from the site Hessians (`dJ_k/dq_j = w_j x J_k`
  for a revolute ancestor j). Checked against central differences on a 6-DOF spatial tree with skew axes, rotated
  placements, a world anchor and every law (largest error 1.3e-11 of the largest entry for dtau/dq, 8e-10 for
  dtau/dqd), against the closed form of a pendulum anchored to the world (1e-12), against `-dV/dq` at rest rates, and
  for a pure preload, where the whole stiffness is geometric (an implementation without the term would report
  zero). Site wrenches reproduce `tau_e` through the site Jacobians to 1e-12.
- **Energy:** `State` carries `dissipated` and `externalWork`; every stepper adds the elements' dissipated energy and
  the applied forces' work `tau . (q1 - q0)`; `computeEnergy` adds the elements' potential and reports the flows, so
  an `EnergyLedger` fed its accounts closes over a run. Over 1 s of the 10-DOF tree with both elements and applied
  forces the relative residual is 1.2e-7 for RK4 at 1 ms, and 1.9e-3, 1.0e-3, 5.4e-4 for implicit Euler at 2, 1 and
  0.5 ms (order 0.93): the residual is implicit Euler's own numerical damping.
- **`stepLinearlyImplicitEuler` / `LinearlyImplicitEuler`** (`multibody/implicit_euler.hpp`):
  `(M - h D - h^2 K~) dqd = h (tau + tau_e - c + h K~ qd)`, then `qd += dqd`, `q += h qd`. Gate evidence: a 1 kg
  slider on 1e6 N/m at h = 10 ms (h omega = 10) is stable and its balance never rises, while `Rk4Stepper` reaches
  1e12 J within 100 steps; the 10-DOF tree with a 1e6 N/m strut between branches and a 1e6 N/m cart spring at
  h = 5 ms is stable while RK4 diverges; the step matches the exact backward-Euler map of a linear oscillator to
  1e-12 and an independently assembled and solved system on both factorisation paths to 1e-11; measured order
  1.16, 1.10, 1.05 for h = 2 to 0.25 ms against an RK4 reference (10-DOF tree with elements, 0.5 s).
- **Allocation gate:** implicit Euler (dense and tree paths), semi-implicit Euler and RK4 with elements, element
  evaluation and energy allocate nothing (`Sub0Physics_MultibodyAllocationTests`).

**Deviations and decisions:**

- **`Forcing` stays the external input** (`tau`). The elements belong to the `Model`, so every stepper,
  `forwardDynamics` (`qdd = M^-1 (tau + tau_e - c)`), `inverseDynamics` (`tau = M qdd + c - tau_e`, its inverse) and
  `computeEnergy` include them without being told; that is how the steppers "include element forces". Prescribed
  external wrenches at sites (section 4) are not in M1b.
- **`LinearSpringDamper<Coordinate>`** is a new law: the coordinate spring-damper of section 4 on a length or an angle
  (the other laws are on lengths, and a revolute joint needs a law on an angle). The joint's `springReference`
  is its default reference.
- **Sites may be fixed to the world** (`SiteDescription::body = WorldBody`), for tendon anchors; `computeSiteKinematics`
  and `computeSiteJacobian` report them as fixed.
- **Dense `n x n` accumulation**, not per-element blocks: the implicit step needs the dense matrix anyway, and
  zeroing two 40 x 40 matrices costs 90 ns against a 9.4 us implicit step on a 40-DOF tree with two tendons (one-off
  measurement, docs/PERFORMANCE.md). An element adds only on its support.
- **The step uses the symmetric part of K.** K is exact in the workspace (the rate-transport term of a damped tendon
  is not symmetric); the step is a W-method, which keeps first order with any Jacobian, so the symmetric part costs
  no accuracy and allows a symmetric factorisation.
- **Factorisation:** while every tendon's ends lie on one chain (`Model::elementsFollowTree`) the matrix keeps the
  tree sparsity and Featherstone's LTDL factorises it with no fill-in; a tendon between branches couples them and a
  dense `LDL^T` is used (n <= 40). No pivoting: the matrix is symmetric positive definite whenever it is
  factorised.
- **Indefinite K:** every law here is passive (`dForceDPosition <= 0`, `dForceDVelocity <= 0`), so
  `M - h D - h^2 K_material` is positive definite. The geometric term (a compressed strut softens transverse motion)
  can make the matrix indefinite when the instability it describes is faster than `1/h`; the step then falls back to
  the material stiffness only (`Workspace::geometricStiffnessDropped`), keeping stiff elements implicit, and faults
  only if that fails. Dropping all of K (the first design) made a stiff strut explicit and diverged.
- **Dissipation is integrated with the power at the start of the step** (left rectangle, first order like the
  step). Starting from rest it telescopes to at most the true discrete dissipation of backward Euler, so the balance
  never rises above its start on a linear system (tested); RK integrates the dissipated energy as one more state, to
  its own order. `ModelOde` accepts `[q; qd]` or `[q; qd; E]`.
- **Common prismatic ancestors are dropped from a tendon's support:** they move both ends equally and change no
  entry (1.4x on the element evaluation).
- **M is not differentiated** (as in MuJoCo's `implicitfast`): the step is implicit in the elements only.

### M1c status (2026-10-20)

**Landed** (the M1c gate, section 9):

- **Description** (`multibody/constraint.hpp`): `LoopClosureDescription` (point or frame equality between two sites),
  `CoordinateLimitDescription` (lower and/or upper stop on a joint coordinate) and `DistanceLimitDescription` (minimum
  and/or maximum distance between two sites), each with a restitution coefficient for limits; `ConstraintSettings`
  holds the solver's targets. `ModelBuilder::addConstraint` and `setConstraintSettings`; `build` validates (resolved
  sites on two different bodies, ordered finite bounds, restitution in [0, 1], positive tolerances, unique names) and
  compiles a fixed **constraint layout**: every closure's rows (3 or 6), then one row per one-sided limit. Per-row
  results (`Workspace::constraintImpulse`, `constraintForce`) follow it.
- **Velocity-level solve with the step's own matrix.** The free step is the stepper's (`A = M` for semi-implicit
  Euler, `A = M - h D - h^2 K~` for linearly implicit Euler, already factorised). The rows `J` of the closures and of
  the limits engaged at the start give `G = J A^-1 J^T` (a few rows, dense), and a mixed LCP gives the impulses `p`:
  closure rows `J qd+ = 0`; limit rows `J qd+ >= -e min(J qd-, 0)`, `p >= 0`, complementary. `qd+ = qd_free + A^-1 J^T p`.
  The LCP is solved by principal pivoting (Murty's least-index rule) over an LDL^T with symmetric diagonal pivoting
  that drops dependent rows, so a planar mechanism's out-of-plane closure rows need no mask (rows that no coordinate
  moves are dropped at assembly, the rest by the pivot tolerance). Exact for a few rows; equality-only sets take one
  solve.
- **Projection, every accepted step:** Gauss-Newton on the closures (and on engaged limits that are violated),
  `dq = -A^-1 J^T (J A^-1 J^T)^-1 phi`, until `max |phi| <= closureTolerance`; coordinate limits still violated are set
  to their bound; then the velocities onto the closures' tangent space. No Baumgarte terms. `projectConstraints` does
  the same in the mass metric, for assembling a mechanism from an approximate configuration and for the RK stepper.
- **Strict limits and event location.** A limit is engaged when its gap is within `eventTolerance`. A trial step that
  leaves a limit that was not engaged with a negative gap is split: the sub-step size is found by Illinois regula
  falsi on the smallest free gap, re-running the stepper's own trial (its matrix depends on the sub-step size), until
  the gap lands in `[0, eventTolerance]`; the rest of the tick then runs with the limit engaged, where the LCP gives
  the impact impulse. An accepted state is inside its limits (coordinate limits exactly, distance limits within
  `closureTolerance`); a state given outside them is refused (InvalidInput), and failing a target within its guard
  rejects the step (OutsideDomain), state unchanged.
- **Ticks.** A split step still advances exactly one tick, from `tick h` to `(tick + 1) h`. Its sub-steps divide the
  integration only: their sizes sum to h (the last is what remains), nothing between them is exposed, and time is
  still derived from the integer tick. `ConstraintReport::substeps` and `engagementGap` report what happened.
- **Steppers.** `stepLinearlyImplicitEuler` (default) and `stepSemiImplicitEuler` honour closures and limits.
  `Rk4Stepper`/`MidpointStepper` honour closures (acceleration level in every stage through `forwardDynamics`, then
  `projectConstraints`; the cached endpoint derivative is dropped, so RK4 with closures costs 5 evaluations per step)
  and refuse limits (InvalidInput: an impact needs event location and an impulse). `forwardDynamics` adds the
  acceleration-level forces of the closures and of engaged, non-separating limits (`J qdd + Jdot qd = 0`, with
  `Jdot qd` from the zero-acceleration RNEA's body accelerations) and writes them to `constraintForce`.
  `inverseDynamics` stays the tree's (constraint forces are an input there).
- **Reactions** (`computeReactions`): the constrained accelerations from `forwardDynamics`, then RNEA with every
  force: gravity, applied forces, elements (coordinate elements in their joints, tendon site wrenches), constraint
  forces (closure and distance-limit site wrenches, coordinate limits in their joints). It writes `jointReactions`
  (the wrench each joint transmits to its body, world axes, about the joint origin) and `constraintSiteWrenches`;
  a joint reaction's component on its motion axis equals the generalized force on that coordinate (tested).
- **Energy.** Limit impulses' work (never positive) is dissipation and goes to `State::dissipated`. The closures'
  impulse work (`State::closureWork`, zero for an exact ideal constraint) and the projection's energy change
  (`State::projectionEnergy`, measured with two energy evaluations when `measureProjectionEnergy`) are reported
  apart from the balance, so a ledger's residual can be read against them.
- **Allocation:** constrained implicit and semi-implicit steps with closures, an LCP row and split steps, RK4 with
  closures, `projectConstraints` and `computeReactions` allocate nothing (`Sub0Physics_MultibodyAllocationTests`).

**Gate evidence** (`tests/test_multibody_constraints.cpp`, GCC 16 and Clang 22 `-march=x86-64-v3`):

- *Slider-crank closure residual:* 10 s at 1 ms with implicit and semi-implicit Euler (3 s with RK4): worst residual
  9.998e-11 m against the requested 1e-10, velocity residual 1.8e-15 m/s, the second half no worse than the first.
  With projection off (`maxProjectionIterations = 0`) the velocity-level closure alone drifts to 7.1e-3 m in 5 s.
- *A limit never overshot:* a 1 m pendulum released from 1.2 rad onto a stop at -0.5 rad never goes below it (exact
  comparison every step; implicit Euler at 1 ms, semi-implicit at 0.5 ms). Engagement located to a gap in
  [0, 1e-9] in 3 or 4 trials. Impact rate -3.903 rad/s against the exact -3.894 (the step's own error); rebound 0 for
  restitution 0 and -0.499 of the impact rate for 0.5; impact dissipation 2.529 J against the exact 2.527 J
  (restitution 0) and 1.896 J against 1.895 J (`(1 - e^2)` of it); ledger relative residual 1.8e-3 and 1.3e-3. A rope
  (maximum-distance limit) on a falling point mass never exceeds 0.5 m. Split steps are deterministic (two runs
  identical) and keep integer ticks.
- *Reactions match independent statics:* a slider-crank held against its end stop by a crank torque and a slider
  force: the closure force, the stop force and all three joint reactions match planar free-body statics to 1e-9
  relative (closure 19.2687, -1.2401 N; stop -40.7313 N), and the accelerations are below 1e-9.
- *The energy ledger closes:* RK4 with the closure, 2 s at 1 ms with applied forces: relative residual 2.1e-8, the
  projection's energy -3.4e-10 J. Linearly implicit Euler, 2 s without forcing at 2, 1 and 0.5 ms: balance residual
  -0.250, -0.135, -0.073 J; closure work -0.139, -0.069, -0.033 J; projection -0.289, -0.153, -0.081 J (initial kinetic
  energy 0.16 J, peak about 0.4 J). Every channel is first order (each halves with h), so it is the scheme's own
  error, measured and attributed, not a leak.

**Target-driven levers** (docs/LOD.md, "Continuous before discrete"). Every lever is deterministic (state,
`ConstraintSettings` and tick only) and is set at model build; a model built with other settings is another model, so
switching in-run means a `decide`/`commit` with a state copy, which needs no transfer (the state shape is the same).

| Lever | Target | Achieved against requested | Cost against target |
|---|---|---|---|
| `closureTolerance` | max closure residual after projection | 1e-5 and 1e-8 requested: 1.7e-9 achieved; 1e-11 and 1e-14: 5.6e-17 (assembling from 2 cm off) | 3, 3, 4, 4 Gauss-Newton iterations (quadratic convergence: cost grows with log log of the target); a step typically needs 1 |
| `eventTolerance` | gap at engagement (and the engagement distance) | 1e-3 and 1e-6: 1.8e-9; 1e-9 and 1e-12: 0 | 3, 3, 4, 4 trial sub-steps per engagement |
| `solverTolerance` | relative pivot and sign tolerance of the small solves | dependent rows dropped (planar closures); LCP sign tests | none measurable at a few rows |
| `maxProjectionIterations`, `maxEventIterations`, `maxEventsPerStep` | guards, not targets | a target missed within its guard rejects the step | - |
| `measureProjectionEnergy` | observability | - | about 70 ns of a 550 ns slider-crank step |

**Deviations and decisions:**

- **The velocity-level solve uses the step matrix `A`, not `M`.** For semi-implicit Euler they are the same; for
  linearly implicit Euler `J A^-1 J^T` is what keeps the impulses consistent with the implicit elements (a stiff
  element and a closure on the same coordinates do not fight). The projection uses the same factor as its metric.
- **Reactions are instantaneous, at a state** (acceleration level), not the step's impulses: a step's impulses include
  impacts and are reported per row in `constraintImpulse`. "The same numbers that moved the bodies" holds for the
  forces at that state; a step-consistent reaction would need the implicit elements' effective forces too.
- **Active set and warm start are not stored in `State`:** engagement is recomputed from the gap each step (direct
  pivoting needs no warm start), so a snapshot needs nothing more.
- **Missed grazes:** a limit entered and left within one (sub-)step is not detected; event location looks only at the
  end of each trial.
- **First-order energy loss with closures.** The velocity-level solve plus projection is a first-order constrained
  Euler method: on the slider-crank about 0.07 J/s at 1 ms (against a peak kinetic energy of about 0.4 J) leaves
  through the closures' impulse work and the projection. Evaluating the rows at the predicted end of the step instead moved the
  loss from the projection into the closure work without reducing it (measured, not kept). See the post-parity ideas.
- **Frame closures** use the orientation error `vee(skew(R1 R2^T))`, whose rate is `w1 - w2` on the manifold; the
  projection converges quadratically near it. A straight chain welded at both ends is a singular configuration (the
  along-chain row vanishes); the pivoted solve drops it and the chain settles by about 1e-5 rad.

**Post-parity improvement ideas** (docs/LOD.md; none is needed for the M1c gate):

- A constrained symplectic scheme (SHAKE/RATTLE-type: the impulse chosen to satisfy the position constraint at the end
  of the step, then a velocity projection), as another stepper behind the same interface, to remove the first-order
  energy loss above. Judge it on the rig's ledger and its cost.
- The Delassus matrix from the factor's first triangular sweep (`G = Y D^-1 Y^T`, `Y = L^-T J^T`) and one full solve
  to apply `J^T p`, instead of a full solve per row: fewer tree sweeps once closures have more rows.
- Projection with a reused `J A^-1 J^T` across the Gauss-Newton iterations (a quasi-Newton lever: cheaper per
  iteration, linear convergence), selected by the same residual target.
- An `eventTolerance` and `closureTolerance` scaled from the step's own error estimate, so a coarse step does not pay
  for a fine target.

### M1d status (2026-10-21)

**Landed** (the M1d gate, section 9). The consumer is a bicycle rig (a frame on heave and pitch, a fork, a swingarm with
a shock between two eyes, a guided rider, progressive end buffers, two vertical tyre springs over moving platens), built
by its consumer with `ModelBuilder` from typed descriptions; nothing here names it.

- **Planar formulation** (`Formulation`, `PlanarModel`, `multibody/planar.hpp`). `build` selects it automatically when
  every revolute axis lies along one coordinate axis k (either sense), every prismatic axis and placement rotation is
  in its plane, and gravity has no k component; the tests are exact, so a model built from in-plane numbers qualifies.
  `ModelBuilder::setFormulation` forces `General` (or requires `Planar`). Transforms, RNEA, CRBA, the tendons with their
  exact Jacobians, energy, site kinematics and site forces then run in three-component planar algebra (RBDA 2.15) and
  write `Workspace::planar`; out-of-plane offsets are kept (a tendon's length includes the constant k separation of its
  sites). The constraint solver and `computeReactions` still read the general fields, which the planar path fills from
  its own results without trigonometry (`expandPlanarKinematics`). One set of algorithms per formulation, chosen per
  call; there is still no separate 2D engine or model description.
- **Massless bodies.** A body with mass and inertia exactly zero only carries a joint: a body with several joints in
  MJCF is a chain of them (the rig's frame: heave, then pitch). `build` checks that the mass matrix at the zero
  configuration still factorises.
- **Site forces and forcing Jacobians** (`Forcing::siteForces`, `Forcing::stiffness`, `Forcing::damping`). External
  forces at sites, held over the step: mapped through `J^T` at the start of the step in every stepper
  (`computeAppliedForce`), and applied where they act in `computeReactionsWithForcing`, so a tyre load reaches the
  joints through the bodies. A state-dependent input (a contact computed by the consumer) may supply its generalized
  Jacobians, which the linearly implicit steps treat like an element's. `forwardDynamicsWithForcing`.
- **Implicit-term lever** (`ImplicitTerms`: `Full`, `MaterialStiffness`, `DampingOnly`). `DampingOnly` is MuJoCo's
  `implicitfast` structure (only the dampers implicit); with the forces held at the start of the step it reproduces
  MuJoCo's scheme, which is how the consumer's engine-parity reference runs.
- **Linearly implicit midpoint** (`stepLinearlyImplicitMidpoint`, `LinearlyImplicitMidpoint`): the same linearisation
  with weight 1/2 (the trapezoid rule on velocities), sharing the implicit Euler code path, which stays bit for bit.
  A-stable, no numerical damping on a linear oscillator (energy kept to 1e-12 where backward Euler loses it), order 2
  measured for constant M. Constraints are refused (InvalidInput) until a RATTLE-type form exists.
- **Two-phase step** (`prepareLinearlyImplicitStep`, `computePreparedSiteKinematics`, `computePreparedSiteMotion`,
  `completeLinearlyImplicitStep`): a consumer that forms its forcing from the state reads sites and energy from the
  step's own kinematics pass instead of a second and third one. prepare then complete equals the one-call step bit for
  bit (tested, both formulations, both schemes); any other entry point withdraws a preparation, and complete refuses a
  different state. `computeSiteKinematics` also has a batch form that can return the energy from the same pass.

**Gate evidence.**

- *Parity with the M0 MuJoCo reference and BikeLab's goldens* (the consumer's parity tests and its MuJoCo
  reference documentation): with `DampingOnly` and forces held, all six golden scenarios
  (hardtail, full suspension, compliant rider with buffers, rider stops, and the two travel-limit terminations)
  agree on every one of 60 channels, the settled state and the run summary to 1.7e-12 of the channel's scale (worst:
  a near-zero fixture load against its 1 N floor; typically 1e-16 to 1e-13), with identical termination instants, under
  M0's 1e-10 bounds. The native MuJoCo reference run with the same MJCF at 0.5, 2 and 10 ms agrees with the core at the
  same steps to 2.8e-12.
- *Timestep advantage, measured* (the consumer's fidelity study): against a converged run, MuJoCo's scheme is
  already stable to 10 ms on these scenarios (5 ms with the weak shock buffer; checked on native MuJoCo), so making the
  springs implicit buys no stability here: the stiffest mode, the wheel on its tyre, has `w h` = 1.7 at 10 ms. It buys
  accuracy only with the midpoint step: BikeLab's own accuracy at 4 to 8 times its step on the smooth scenarios
  (whole runs 1.9 to 3.7 times faster), none where buffers engage (their law kinks; event location is the lever, as for
  the stations). Linearly implicit Euler is never as accurate as BikeLab's scheme at the same step: it damps the
  tyre mode at about `w^2 h / 2` per second.
- *Planar formulation agrees with the general one*: every query (mass matrix, bias, forward and inverse dynamics,
  element force and Jacobians, energy, site kinematics, reactions) to 1e-13 relative, stepped trajectories to 1e-11
  (implicit, semi-implicit, constrained), on the 10-DOF tree, the tree with elements, the slider-crank and a rig-like
  tree in the xz plane with a massless carrier, rotated placements, a -y axis and out-of-plane offsets
  (`tests/test_multibody_planar.cpp`).
- *Determinism*: the consumer's whole-run digest over every channel of the six scenarios and all three steppers is
  identical from GCC 16, GCC 16 `-march=x86-64-v3` and Clang 22 `-march=x86-64-v3`.
- *Allocation*: site forces, forcing Jacobians, the midpoint and prepared steps, the site-motion query and reactions
  with forcing allocate nothing, in both formulations (`Sub0Physics_MultibodyAllocationTests`).

**Performance** (docs/PERFORMANCE.md). On the 10-DOF tree, GCC 16, three interleaved rounds against the previous
build: mass matrix 2.5-3.0x to 1.16-1.22x of the hand-written planar twin, forward dynamics 2.5-2.7x to 1.20-1.24x,
implicit step 1.9-2.1x to 1.1-1.5x, elements 2.0-2.4x to 1.5-1.6x; Clang 22: mass matrix 0.9x, forward dynamics 1.35x,
implicit step 1.5x, elements 1.8-2.1x. The constrained step is unchanged (about 6.8x): its general constraint
machinery dominates, not the transforms. On the consumer's whole rig run against a hand-written twin of the same
run: first 3.0x, now about 2.0x (GCC; Clang 2.4-2.6x).

**Deviations and decisions:**

- **The rig needs no loop closure.** BikeLab's linkage is a single pivot with the shock as a spatial tendon between two
  eyes, so the consumer uses a tendon element, whose length is checked against its linkage oracle to 1e-15. The
  closures of M1c are not exercised by the rig; their energy question (below) stays open for four-bar linkages.
- **No fixed joint.** A welded body (the hardtail's rear, `fusestatic="false"` in MJCF) is merged into its parent by
  the consumer, which reports the weld's load from the parent's acceleration. A fixed joint would break one coordinate
  per body; not added for one consumer.
- **The general kinematic fields are not written by a planar model's dynamics** (documented on `Workspace`); the
  functions that read them (constraints, reactions) fill them. Bit-identity between the formulations is not a goal:
  each is deterministic, and they agree to rounding.
- **Forcing Jacobians take the dense factorisation**: they may couple any coordinates. At the rig's 4 or 5 DOF this
  costs nothing measurable; a sparsity hint can come with a larger consumer.
- **Site forces are mapped at the start of the step** (like MuJoCo's `xfrc_applied`), so a time-dependent input is the
  consumer's to evaluate at the time its scheme needs (the rig's improved steppers take the platen at the end or the
  middle of the step; held at the start, the midpoint step was first order).
- **Midpoint keeps `M` and `c` at the start of the step**, like the Euler form: second order only while M is constant;
  on the rig the configuration-dependent terms are below the converged truth's own agreement (3e-4 mm).

**Energy (the open question of M1c).** On the rig, implicit Euler's numerical damping is material: the run's ledger
residual over its energy scale is 6.6e-5 to 5.3e-3 at 125 us, against 3.3e-6 to 6.5e-4 for MuJoCo's scheme and
5.6e-7 to 2.6e-4 for the midpoint step (the consumer's fidelity study data). The midpoint step is the
energy-better model for unconstrained mechanisms; with closures, a RATTLE-type midpoint (impulse chosen to satisfy the
position constraint at the end of the step) is the next prototype, to be judged on a four-bar linkage. It does not
block M1d.

**Post-parity levers this unit exposes** (docs/LOD.md):

| Lever | Kind | Status |
|---|---|---|
| Formulation (planar or general) | Not a fidelity lever | Same results to rounding; chosen at build |
| Implicit terms (`Full`, `MaterialStiffness`, `DampingOnly`) | Discrete | Measured on the rig: DampingOnly (MuJoCo's) is more accurate than Full at equal step for lightly damped modes; Full matters only past the explicit stability limit |
| Implicit Euler or midpoint | Discrete | Midpoint dominates on smooth phases (4-8x step); equal at buffer kinks |
| Step size with error control | Continuous | Not yet for multibody: the midpoint and Euler steps form an embedded pair (their difference estimates Euler's error at no extra solve), the natural controller; with event location at the buffers' and tyres' switches |
| Jacobian and factorisation reuse across steps | Continuous | Not yet: the W-method keeps first order with a stale Jacobian, so `prepare` could be skipped for k steps in Interactive and Ambient profiles |
| Geometric stiffness on or off | Discrete | `MaterialStiffness`; no measurable effect on the rig (its one tendon is a linear shock) |

### M2–M3 plan (proposed, 2026-10-21)

**Design only; nothing implemented.** The deformable wheel and its assembly under a moving tree are designed in
[WHEEL.md](WHEEL.md). In brief:

- **Casing:** a ring of material nodes with two displacement coordinates each, attached to a host body as a
  structured *substructure*. Its coordinates are generalized coordinates; the step eliminates them first (periodic
  block-tridiagonal, a Schur complement onto the host's motion), then factorises the tree as now. Two rings of 64
  nodes make n ≈ 264, beyond §1's dense sizing: a throwaway prototype measured 4.9 µs for the condensed solve against
  1.4 ms for a dense `LDLᵀ`. Generic bodies and tendons (the same description expanded) remain the test oracle at small
  N. Cheaper casing models (rigid ring with an enveloping footprint, modal, quasi-static) are further models of the
  same mechanism.
- **Contact:** the reference engine's soft-constraint law as the parity model (ADR 0003, decision 4); a compliant
  contact element with exact Jacobians and exact Coulomb friction as bounded rows as the improved model; event-located
  onset with a certified swept bound.
- **Drive and brakes:** a generalized torque, and a bounded row that gives an exact wheel lock.

| Step | Delivers | Gate |
|---|---|---|
| M2a | Ring substructure, condensation, expanded oracle | Oracle agreement to rounding; mass and inertia; free flight; Jacobians; allocation; paired benchmark |
| M2b | Capsule-edge contact geometry and witnesses | Contact lists against the reference engine on its own states |
| M2c | Soft-constraint law; the moving wheel | Open-loop force parity, then whole runs within bounds set by the reference's own sensitivity |
| M2d | Improved contact, swept certification, error-controlled midpoint | Converged-run comparison, ledgers, no tunnelling |
| M3a–M3d | Connected vehicle frictionless, traction, drive and brakes, profiles | As M2c and M2d on the connected scenarios; cross-compiler digest |

**Proposed contract changes, for decision before M2 builds on them** (WHEEL.md §10): substructures; bounded
velocity rows (brakes, friction against rigid surfaces), which amends §1's "contact is force elements only" for
friction while normal contact and soil stay elements; a contact element (a law on a signed gap); block storage for
element Jacobians when substructures are present; velocity-only rows admitted in the midpoint step.

## 10. Open questions, to be answered by measurement

- **Higher-order implicit scheme:** whether one is needed, and which, for the tyre ring (M2).
- **LCP method:** direct pivoting against projected Gauss-Seidel once contacts or limits outnumber a
  handful (expected never for one bike). M1c uses principal pivoting (Murty) over a pivoted LDL^T; its cost at a
  few rows is not measurable against the rest of the step.
- **Constrained energy behaviour:** the rig has no closures (M1d status); implicit Euler's own numerical damping is
  material on it, and the linearly implicit midpoint step removes it for unconstrained mechanisms. Whether a
  RATTLE-type midpoint is needed for closures is to be measured on a four-bar linkage when one is modelled.
- **Error control for multibody steps:** done in M2d (step doubling on a dyadic lattice, events bisected; M2d status).
- **Transcendentals:** decided, see section 8: our own (`sub0math::det`). Remaining to verify: the ring and ring-contact
  digests match on Linux and macOS CI.

### M2a status (2026-10-27)

**Landed:** substructures (WHEEL.md section 10, item 1) and block storage (item 4). A model may attach deformable rings
to bodies; their coordinates follow the tree's (`Model::dofCount` = `treeDofCount` + 2 per node), the Workspace's dense
matrices are `treeDofCount` square, and every step, forward and inverse dynamics, energy and reactions condense the
rings onto their hosts before the tree factorisation. Models without rings are unchanged bit for bit (the consumer's
parity tests pass unchanged). Details, gate evidence and deviations: [WHEEL.md, M2a status](WHEEL.md#m2a-status-2026-10-27).

### M2c status (2026-10-27)

**Landed:** contact elements (WHEEL.md section 10, item 3) and the reference engine's soft-constraint law, with these
contract changes (details, derivation and evidence in [WHEEL.md, M2c status](WHEEL.md#m2c-status-2026-10-27)):

- `Forcing::ringLoads`: contact element loads on rings, held over the step like `tau`; their force is mapped through
  the host's chain, and the linearly implicit steps take their exact Jacobians through the ring's block storage, the
  condensation's border and the tree's matrix (the chain's pairs only, so LTDL keeps its sparsity). Dense forcing
  Jacobians stay refused with rings.
- `ImplicitTerms::DampingOnlyMassPattern`: MuJoCo's `implicitfast` exactly. The engine gathers its velocity derivative
  into the mass matrix's sparsity, so dampers between branches (a ring's neighbouring nodes, a tendon across branches)
  stay explicit; `DampingOnly` is the same for models whose dampers follow the tree (the rig: its parity is unchanged
  bit for bit).
- `multibody/soft_constraint.hpp`: the soft-constraint law written once for every row kind, and one exact solve of a
  step's rows with the mass matrix (contacts now; soft equality closures in M4's reference path). Nothing of it lives in
  `State`: the solve is direct, so there is no warm start to keep.
- `multibody/ring_contact.hpp`: a contact point's row on a ring, the node inverse weights the reference scales its
  regularisation by, and the contact element's gap kinematics with exact first and second derivatives.
- The ring's condensation skips the periodic sweep when its block is block-diagonal (always for the mass matrix), bit
  for bit.

### M2d status (2026-10-28)

**Landed** (details and evidence in [WHEEL.md, M2d status](WHEEL.md#m2d-status-2026-10-28)): the error-controlled
multibody step section 10 asked for, on a dyadic lattice (`integrator/dyadic_step.hpp`: `DyadicClock`,
`DyadicStepController`, `advanceDyadically`), with step doubling of the linearly implicit midpoint step as the estimate
rather than the Euler-midpoint pair (the two do not share a matrix: `M - h D - h^2 K` against `M - h/2 D - h^2/4 K`);
`PreparedStepUse::Keep` lets one preparation serve the coarse step and the first half step. Contract changes:

- `completeLinearlyImplicitStep(..., PreparedStepUse use = Consume)`: Keep leaves the preparation valid for the same
  starting values (bit for bit a fresh preparation).
- `SoftConstraintSolver::solve(..., supports, DelassusForm)`: row supports make the uniform form's products sparse (bit
  for bit) and enable the structured form.
- `UniformRingSpec::foundationViscosity` (default zero).

Section 10's open question "error control for multibody steps" is answered for the ring; "higher-order implicit scheme"
stays open: the midpoint step is not L-stable, and the improved wheel shows where that matters (an unresolved stiff
contact mode is not damped, WHEEL.md M2d finding 1).

### M3a consumer status (2026-10-28)

**No core change.** The consumer's connected trail bike (BikeDynamics docs/CONNECTED_BIKE.md, "M3a status") runs on the
M2a to M2d pieces as they are: one tree of up to eight coordinates (massless carriers for the frame's travel and heave,
the frame on pitch with a centre-of-mass offset, the fork, the swingarm with a shock tendon, the guided rider, two hubs)
and two 64-node rings on hubs deep in the tree (264 coordinates). Its reference stepping (both rings' soft rows solved
together with `M`, `DampingOnlyMassPattern`) reproduces MuJoCo's step to 7e-13 of the step's change on every engine state
of ten runs when given the engine's constraint force, and its whole-run digest is identical on the three builds.

**Contract note found by the consumer** (documentation, not a defect): `computeReactionsWithForcing` treats `Forcing::tau`
as forces inside the joints, so a contact on a ring passed as its generalized force `J^T f` (the reference path's form)
reaches the reactions as joint forces, and the ring nodes' Newton-Euler sees no external force: the reactions omit ring
contact loads, as MuJoCo's `mj_rnePostConstraint` omits flex contacts. Consumers remove the ring's ground wrench where
they need a subtree's interface load, as BikeLab does. "No force category is left out" holds for site forces, elements,
closures and limits; a ring contact load given as a site-like external force on the ring would remove the correction
(a candidate contract change, not needed for parity).
