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

## Not provided yet

Shadow references running beside the live model, asynchronous comparison windows, and state
adapters between models with different state shapes. They arrive with their first real consumer.
