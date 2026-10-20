# Multibody core (design, step M1 of ADR 0003)

- **Status:** M1a and M1b landed (2026-10-20); M1c next
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

## 10. Open questions, to be answered by measurement

- **Higher-order implicit scheme:** whether one is needed, and which, for the tyre ring (M2).
- **LCP method:** direct pivoting against projected Gauss-Seidel once contacts or limits outnumber a
  handful (expected never for one bike).
- **Transcendentals:** whether `std::sin`/`std::cos` are good enough for cross-platform hashes, or
  whether we need our own (M3 gate).
