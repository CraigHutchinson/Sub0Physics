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

## Not provided yet

Shadow references running beside the live model, asynchronous comparison windows, and state
adapters between models with different state shapes. They arrive with their first real consumer.
