# Level of detail

How Sub0Physics lets a product run cheaper models where they are good enough, and detailed ones
where they are needed, without losing correctness or replayability. The design follows BikeLab's
[physics fidelity architecture](https://github.com/CraigHutchinson/BikeLab/blob/main/docs/PHYSICS_FIDELITY.md).

## Principles

- **Fidelity belongs to a mechanism, not the world.** Each mechanism (a contact, a damper, a tyre)
  can have several models. A `ModelDescriptor` names one by id and revision, with a relative cost.
  There is no global "LOD 0..3" integer that claims to be more accurate everywhere.
- **Select on physical need, not camera distance.** The product defines a scalar *demand* for each
  mechanism (contact activity, strain rate, predicted error, closing speed). A `FidelityBand` says
  over what range of demand a model is valid.
- **Hysteresis and dwell stop chatter.** A model may be entered only at or below `enterBelow` and stays
  valid up to `leaveAbove`. Refinement is immediate when the active model leaves its envelope.
  Coarsening waits for `minimumDwellTicks`.
- **Decide, transfer, commit.** `decide` proposes a change; the product maps state between the models;
  `commit` applies it only if that transfer succeeded. A failed transfer leaves the old model running.
  A stale decision cannot be committed.
- **Nothing valid is a reported outcome.** `NoValidModel` keeps the active model and tells the caller,
  who must report loss of validity or stop. It is never silently ignored.
- **Deterministic.** Decisions depend only on demand and the integer tick, so a replay reproduces them.

## Tiers provided today

| Mechanism | Detailed | Cheaper | Exactness of the cheaper tier |
|---|---|---|---|
| Time integration | `Rk4` (order 4, 4 evaluations per step) | `Midpoint` (order 2, 2 evaluations) | Exact while acceleration is constant, e.g. ballistic flight between contacts |

The first consumer is a loaded wheel: midpoint while airborne, RK4 from just before
contact. That halves the cost of flight phases with no loss of accuracy there.

## Verification of a switch

A product that switches models must show, per mechanism:

1. The reduced model agrees with the detailed one within declared per-channel limits on held-out
   scenarios (`ErrorStatistics`), and its energy balance holds (`EnergyLedger`).
2. Both are in their asymptotic range (`observedOrder` near the nominal order).
3. Repeated switches near the threshold do not inject energy or duplicate reactions.

## Reference models, levers and quality profiles

- **A port is a reference, not a ceiling.**
  - A mechanism ported for parity with an upstream implementation becomes that mechanism's *reference model*. It is
    frozen, and changes only when its upstream changes.
  - Once it is stable, a post-parity review asks whether the port's choices limit us: upstream arithmetic
    semantics, fixed resolutions, fixed steps, cold-started iterations, explicit methods where an implicit
    one would allow larger steps.
  - Improvements become further models behind the same interface. The reference is never edited for speed.
- **Levers.** A lever is a deterministic parameter of one model that trades error for cost: step size, node count,
  iteration tolerance, warm start, decimation, arithmetic mode, Jacobian reuse. A lever parameterises one
  implementation; it never copies a formula.
- **Evidence for a lever.**
  - An error-against-cost frontier against the reference on held-out scenarios (`ErrorStatistics` per
    channel, plus a benchmark).
  - An energy balance at every setting (`EnergyLedger`).
  - Determinism across compilers.
  - Whether the lever can switch in-run (`decide`/`commit` with a state transfer) or only at setup.
- **Continuous before discrete: pay for what you need.**
  - Discrete LOD switches a model (a tier). Continuous LOD varies a lever smoothly: step size, solver
    tolerance or iteration count, resolution, decimation.
  - Prefer a continuous lever driven by a *target* over a fixed setting:
    - an **error estimator** (embedded Runge–Kutta pair, step doubling, solver residual, a certified clearance
      bound);
    - a deterministic **controller** that sets the lever each step to meet a requested tolerance, or to
      spend a compute budget where the marginal error per unit cost is highest.
  - Compute then concentrates where the physics demands it (contact onset, impacts, stiff transients) and
    falls away where it does not (flight, steady rolling).
  - Discrete models remain for genuinely different formulations, and a `FidelitySelector` switches them on
    demand.
  - **Evidence for a controlled lever:** the achieved error meets the requested target across the scenario set
    (measured against the reference); cost falls smoothly as the target loosens; and decisions depend only on
    state, configuration and tick, so a replay reproduces them.
- **Quality profiles.**
  - A product names profiles (for example offline reference, engineering analysis, interactive, ambient). Each maps
    every mechanism to a model and, preferably, an error target or compute budget for its controlled levers,
    rather than fixed lever values. The combination is validated together as compatible.
  - A profile is not a global LOD integer: within a profile, a `FidelitySelector` can still raise one mechanism on
    physical demand.

## Levers provided

Each is deterministic (a function of configuration, state and tick) and its default reproduces the existing result
bit for bit. *Discrete* levers select between models; *continuous* ones are driven toward a target, so compute is
spent only where it is needed.

| Lever | Kind | Where | Target and controller |
|---|---|---|---|
| Norm arithmetic (`NormArithmetic::Compensated` or `Fast`) | Discrete | `geometry/norm.hpp`, `PreparedProfile` | None: `Fast` is a guarded `sqrt(a^2 + b^2)`, safe for every finite input |
| Chain solver (`solveBoundedPeriodicTridiagonal` or `...ActiveSet`) | Discrete | `solver/` | The active-set solve is exact to rounding: no tolerance to tune |
| Chain residual tolerance | Continuous | `BoundedPeriodicOptions::tolerance` | A residual target; with node stiffness `a` it bounds every node's error by about `tolerance / a` |
| Warm start | Discrete (on or off) | `BoundedPeriodicOptions::warmStart`, `BoundedPeriodicActiveSetOptions::warmStart` | Start from the caller's guess (the previous solution of a slowly changing problem) |
| Integrator (RK4, midpoint, embedded 3(2) or 5(4), linearly implicit Euler) | Discrete | `integrator/` | Selected per mechanism; linearly implicit Euler needs the system's Jacobian (`OdeSystemWithJacobian`) |
| Step size | Continuous | `advanceAdaptively` | An error tolerance per component (`ErrorTolerance`); a deterministic controller (`StepSizeController`) sets each sub-step inside the caller's fixed tick |
| Event location | Discrete (on or off) | `advanceAdaptively` event function | Steps end on a sign change of the event function (a contact or branch switch) within a time tolerance, so no step crosses a kink in the force law |
| Extraction resolution (`HeightRequest::stride`) | Continuous | `terrain/height_extract.hpp` | A height-error tolerance; `strideForTolerance` picks the coarsest power-of-two stride whose estimated error `M s^2 / 4` meets it (docs/TERRAIN.md Levers) |
| Corridor tightness (`tightCorridor`) | Continuous | `terrain/corridor_builder.hpp` | A relative cell-saving target; doubles the sub-interval count until a doubling gains less. Containment holds at every count |
| Quadrature slots | Discrete (per slot) | `LinearlyImplicitEuler::setQuadrature` | Integrals nothing depends on use the trapezoid rule on the endpoint derivatives, so energy ledgers measure the method, not the bookkeeping |
| Multibody implicit terms (`ImplicitTerms::Full`, `MaterialStiffness`, `DampingOnly`, `DampingOnlyMassPattern`) | Discrete | `multibody/implicit_euler.hpp` | None. `DampingOnly` is MuJoCo's `implicitfast` structure and the more accurate at equal step for lightly damped modes; `Full` is needed only past the explicit stability limit (docs/MULTIBODY.md, M1d status). `DampingOnlyMassPattern` is `implicitfast` exactly (dampers between branches explicit); a parity setting, not a fidelity lever (M2c) |
| Multibody step (linearly implicit Euler or midpoint) | Discrete | `stepLinearlyImplicitEuler`, `stepLinearlyImplicitMidpoint` | None yet. The midpoint step has no numerical damping and order 2 for constant M; the two share their solve, so they are an embedded pair for a future step-size controller |
| Multibody formulation (planar or general) | Not a fidelity lever | `ModelBuilder::setFormulation` | Same results to rounding; the planar one is cheaper and chosen automatically |

## Levers of the deformable ring (recorded, not yet studied)

The ring substructure (docs/WHEEL.md, M2a status) exposes two levers of the nodal ring and reserves the cheaper models:

| Lever or model | Kind | Status |
|---|---|---|
| Node count N | Discrete | Any N >= 3 through the description; the canonical tables cover divisors of 128. No error estimate yet: judge against a converged ring (the 32/64/128 sensitivity) |
| Non-uniform spacing | Continuous (spacing profile) | Per-node and per-edge coefficients make a fixed non-uniform ring a description, not a code change; between runs only, after the mass, first-moment and spin-inertia checks (WHEEL.md section 2.5) |
| Quasi-static, modal and rigid-ring casing models | Discrete (models of the same mechanism) | Future cheaper models behind the substructure (WHEEL.md section 2.4); they need state transfers |
| Contact law (reference soft constraint, contact element) | Discrete (models of the same contact) | M2c: the reference rows (parity) and the contact element with exact Jacobians (the improved model's carrier); the element's law, its parameters and event-located onset are M2d's, judged against a converged run (continuous, target-driven: step and event tolerance first) |
| The soft solve's Delassus matrix (one mass solve per row, or the condensed Schur complement) | Cost only (same result) | M2c: one solve per row; the O(k^2) Schur form is the recorded lever (docs/PERFORMANCE.md, open gap) |

## Not provided yet

Shadow references running beside the live model, asynchronous comparison windows, and state
adapters between models with different state shapes. They arrive with their first real consumer.
