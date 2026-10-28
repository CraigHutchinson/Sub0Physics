# Deformable wheel (design, steps M2–M3 of ADR 0003)

- **Status:** design, 2026-10-21. M2d (the improved model's dyadic error-controlled step, swept certificate,
  foundation viscosity and the Delassus A/B) landed 2026-10-28, see [M2d status](#m2d-status-2026-10-28). M2c (the contact element, the reference soft-constraint law and their place in the
  substructure's step) landed 2026-10-27, see [M2c status](#m2c-status-2026-10-27). M2b (contact geometry, §3.1, and the certified skip of §3.4) is implemented
  (2026-10-27): casing edges against canonical convex pieces with exact witnesses and the rotation-aware clearance
  bound, in `geometry/convex_terrain.hpp` on Sub0Math's closed-form capsule contacts; see
  [TERRAIN.md](TERRAIN.md#casing-contact-convexterrain). M2a (the ring substructure, its condensation and its oracle) landed
  2026-10-27, see [M2a status](#m2a-status-2026-10-27). The rest fixes the core-side contracts for the
  spatial deformable wheel (M2) and for assembling several of them under a moving tree (M3), so the steps can
  be delegated and reviewed against one specification.
- **Builds on:** [MULTIBODY.md](MULTIBODY.md) (tree, elements, constraints, steppers), [LOD.md](LOD.md)
  (models, levers, profiles), [TERRAIN.md](TERRAIN.md) (prepared profiles, witnesses, corridors, the periodic
  chain solver).
- **Provenance:** BikeLab's moving-wheel and connected-vehicle experiments (iterations 15–19: `MOVING_WHEEL.md`,
  `TRAIL_BIKE.md`, `ADAPTIVE_WHEEL.md`, and the MJCF built by `dist/moving-wheel.js` and `dist/trail-bike.js`) at
  BikeLab `e73f6cc`. They run on MuJoCo 3.16.0 flex bodies. Here MuJoCo stays a test-only reference.

## 1. What is received

BikeLab's wheel is a ring of N material nodes (64 by default; 32 and 128 supported) on a hub that spins freely
on a carrier:

| Part | Reference mechanism (BikeLab MJCF) |
|---|---|
| Hub | A body on a hinge about the axle, relative to its carrier. Mass and inertia of the rigid hub |
| Casing nodes | N point masses, children of the hub, each on two slide joints along the hub's x and z axes: displacements `u`, `v` in the hub frame from canonical rest positions on the centreline radius `R = outer radius − skin`. Node mass `tyreMass / N`, inertia `m skin² / 2` |
| Foundation | A coordinate spring on each slide, stiffness `k ds` with `k` in N/m² and `ds = 2πR/N`: energy `Σ ½ k ds (u² + v²)` |
| Edges | A closed chain of N capsules of radius `skin` between neighbouring nodes. Extension energy `Σ ½ (EA/L0)(L − L0)²`, `L0 = 2R sin(π/N)`, and a dashpot `viscosity / L0`. Tension and compression; no bending, inflation or pressure law |
| Contact | Native capsule-edge contact against convex terrain prisms (one per polyline segment after collinear merging). Soft, regularised (`solref`, `solimp`); `condim=1` frictionless, or elliptic `condim=3` with one homogeneous coefficient |
| Drive and brake | A generalized torque on the hub hinge; a brake as the hinge's regularised friction loss whose capacity is the command |
| Rest geometry | Canonical binary64 circle coordinates and chord ratios for the supported N, not runtime trigonometry (an ulp of `sin` differed between runtimes and contact amplified it) |

The step is MuJoCo's `implicitfast` at 31.25 µs. That step was chosen for collision-normal accuracy, not
stability: the ring's own admission bound (`h sqrt((k ds + 4 EA/L0) / m) ≤ 0.55`, m the node mass) allows about 250 µs.

What BikeLab records as open, and this design must not hide: sharp-corner casing-mesh convergence (64 to 128
nodes changes a square-edge impulse by about 4.6%), complete swept collision of the deforming casing, exact
static friction and exact wheel lock (regularised constraints creep), and reactions inside the wheel (MuJoCo's
post-constraint RNE omits flex contacts, so only complete-wheel cuts are corrected).

## 2. Casing representation

### 2.1 The structure of the step

Put the nodes' displacements into the generalized coordinates, after the tree's, as BikeLab does. The
linearly implicit step solves with `A = M − h D − h² K` (MULTIBODY §5). For a ring on a host body (the hub):

- **Mass.** A node's world position is `x_i = p + R(θ)(r_i + w_i)`, with `w_i = (u_i, v_i)`. Its kinetic energy couples
  `ẇ_i` to the host body's motion, and through it to every ancestor coordinate, but never to another node: each
  node is a leaf. The mass matrix is an *arrowhead*: the tree block, a 2×2 block per node, and a border from each
  node to the ancestors.
- **Foundation.** Coordinate springs: diagonal 2×2 blocks.
- **Edges.** The length `|r_i + w_i − r_j − w_j|` does not depend on the host's pose (the rotation cancels). Edges
  couple only neighbouring nodes: a *periodic block-tridiagonal* ring, with no border term.
- **Contact.** A contact force depends on world positions, so its Jacobian couples the touched nodes (one, or two
  for a point inside an edge) and the host's motion.

Every coupling between a ring and the rest of the system passes through the **host body's motion**: three
components in the planar formulation, six in general. The step matrix is therefore a tree, plus per ring a
periodic block-tridiagonal block, bordered by a matrix of rank at most 3 (6) onto the host.

Element Jacobians break the tree sparsity that the LTDL factorisation relies on (MULTIBODY M1b: a tendon between
branches falls back to a dense `LDLᵀ`). Dense is the wrong tool here: n ≈ 6 + 2(1 + 2·64) = 264 for a two-wheeled
vehicle, against the 10–40 that §8 of MULTIBODY sized dense matrices for.

### 2.2 Options

| Option | What it is | Step cost (two rings, N = 64) | Parity with the flex reference | Generality | Verdict |
|---|---|---|---|---|---|
| (a) Nodes as model bodies | 2N bodies (or N bodies on a two-joint chain) per ring, edges as tendon elements, a general sparse factorisation with a fill-reducing ordering | Asymptotically as (b2) with the right ordering, but every node pays spatial-body kinematics, a tendon's site Hessians and the dense `n × n` element accumulation of M1b (two 264² matrices: about 1.1 MB to clear per step). Dense `LDLᵀ` of the whole matrix: **1.4 ms** (prototype, below) | Coordinates map one to one | Highest: any body can be "deformable" | Rejected as the runtime path. Kept as the **test oracle** at small N (8, 16): the same description expanded into generic bodies and tendons, factorised densely |
| (b1) Ring as a co-simulated subsystem | Its own state and integrator, coupled to the hub by interface forces each step (or subcycled) | Cheap per step | Coordinates map; dynamics coupled at first order | Medium | **Rejected.** It reintroduces the held-force, first-order coupling that ADR 0003 (N4) left MuJoCo to avoid, and needs its own stability argument for subcycling |
| (b2) Ring as a structured substructure in one monolithic step | The ring's coordinates are in `State`; the step eliminates them first (static condensation, a Schur complement onto the host's motion), then factorises the tree as now | **4.9 µs** for both rings in the prototype (forward elimination and tree solve; back-substitution adds about as much again). Linear in N | Coordinates map one to one, so trajectories compare channel by channel | A substructure kind attached to one body; the same pattern serves a 3D casing (block-banded), a modal ring and a quasi-static ring | **Recommended** for every ring that has dynamic casing state |
| (c) Modal ring | Fourier modes of the ring (radial and tangential, k ≤ K) as coordinates; exact eigenmodes of the uniform unloaded ring, so mass and stiffness are diagonal | O(K) plus contact | None with the nodal reference | Uniform rings only; adaptivity is in K | A **reduced model**, not a reference: a localised contact (an edge, a rock corner) excites modes up to N/2, so K ≪ N/2 smooths the footprint. Useful in flight and for the global casing modes. Its natural partner is (d) for the footprint (the "rigid ring plus enveloping footprint" family in the tyre literature) |
| (d) Quasi-static ring | No casing coordinates. Each step the footprint solver (`solveBoundedPeriodicTridiagonal` or its active-set form) finds the static deflection under the hub's pose; the result is a force element on the hub with Jacobians from the active-set sensitivity | One footprint solve per wheel per step (tens of µs at 64 nodes, measured families in the consumer), but steps of milliseconds: no casing inertia, so the stiffest mode is the hub on its tyre | None (no casing dynamics) | Any profile | The **cheapest deformable tier**: correct static enveloping of steps and corners, no casing waves, no inertia transfer at impact |

**Prototype** (throwaway, not in the repository): two rings of N periodic 2×2 blocks with a dense six-wide
border onto a six-coordinate tree, solved by static condensation against a dense `LDLᵀ` of the assembled matrix,
GCC 16 `-O2 -ffp-contract=off`, CrogLegion, unloaded:

| N | n | Dense `LDLᵀ` + solve | Condensation (forward + tree) | Tree-coordinate agreement |
|---|---|---|---|---|
| 32 | 134 | 193 µs | 3.4 µs | 3e-14 relative |
| 64 | 262 | 1405 µs | 4.9 µs | 1.2e-15 |
| 128 | 518 | 12 600 µs | 11.1 µs | 2.0e-15 |

The dense kernel is naive; a blocked one would be several times faster and still two orders of magnitude behind.
The border here is six wide; through the host's planar motion it is three wide, which is cheaper.

### 2.3 Recommendation: one substructure, several casing models

- **Description.** A `DeformableRingDescription` attaches to a host body: node count, canonical rest coordinates
  (a table, never computed with runtime trigonometry), per-node mass and integration weight, the foundation law,
  the edge law (extension and viscosity coefficients per edge), the skin radius and the contact law selection.
  It validates like every description (`std::expected<…, Fault>`). The test-only MuJoCo adapter compiles it to
  node bodies, slide joints and a dim-1 flex, exactly as BikeLab's MJCF does, so the description drives both
  engines (ADR 0003, decision 2).
- **State.** The ring's `w` and `ẇ` are generalized coordinates, after the tree's, in material-node order. Their
  ordering never depends on where the ring touches the ground.
- **Laws exist once.** The foundation is `LinearSpringDamper` on each slide; an edge is a law on the distance
  between two nodes, with the tendon element's exact Jacobian (material, geometric and rate-transport terms,
  MULTIBODY M1b). The substructure's kernel evaluates the same laws over structure-of-arrays storage; it does not
  copy their formulas.
- **Solve.** Per ring, factorise the periodic block-tridiagonal block (wrap handled by one fill column, as the
  prototype), eliminate it onto the host's motion, add the resulting symmetric 3×3 (6×6) matrix to the host's
  composite inertia in the tree factorisation, then back-substitute. The scalar bounded periodic solver in
  `solver/` and this block solver share one periodic elimination over a block type (scalar or 2×2, later 3×3), so
  there is one implementation of the periodic sweep.
- **Element accumulation.** The dense `n × n` element matrices of M1b stay for trees without substructures. A
  substructure accumulates into its own block storage (diagonal and neighbour blocks, border rows). The decision in
  M1b ("dense sums, an element writes only its support") is re-measured at n = 264 before it is kept anywhere else.
- **Oracle.** The same description expands into generic bodies and tendon elements (option a). At N = 8 and 16
  the structured and the expanded models must agree to rounding: kinematics, energy, forces, Jacobians, one step and
  a short run.

### 2.4 The model ladder

Fidelity belongs to the mechanism (LOD.md). Each casing model is a `ModelDescriptor`; the ones with the same state
shape switch freely, the others need a state transfer.

| Model | Casing state | Step set by | Use | Transfer to the nodal ring |
|---|---|---|---|---|
| **Nodal ring, reference** | 2 per node | The reference contact law's time constant and collision accuracy (31.25 µs in the reference) | Parity with the flex reference | — |
| **Nodal ring, improved** | 2 per node, any N, non-uniform spacing allowed (§2.5) | The lowest casing modes that matter and contact onset (event-located), with an error-controlled step | Analysis, detailed interactive | Same shape: none |
| **Rigid ring with enveloping footprint** | 3 (the belt's rigid motion relative to the hub) | The belt-on-sidewall mode | Interactive, ambient on rough ground | Project nodal displacement onto modes k = 0, 1 (momentum and energy kept; the remainder goes to a transfer-loss ledger, never into damping) |
| **Modal ring** | 2 per mode, K modes | The K-th mode | Flight; global casing modes for comfort or noise channels | Fourier projection; exact for the retained modes |
| **Quasi-static ring** | None | The hub on its tyre (milliseconds) | Ambient; static analysis | Reconstruct the static deflection; casing momentum is lost and logged |
| **Rigid disc on a radial spring** | None | The same | Distant vehicles | As above |

The enveloping footprint of the rigid-ring and quasi-static models is the same solver, so it exists once.

### 2.5 Material coordinates and non-uniform rings

BikeLab's adaptive-wheel review is adopted as the data contract:

- **Material identity.** A node has a stable id, a rest arc coordinate, canonical rest position, integration weight,
  mass and intrinsic inertia. An edge has its endpoint ids, rest arc interval, chord rest length and coefficients.
  Contacts and any material history (friction state, coating, wear) refer to these ids and the ring's topology
  revision, never to storage slots.
- **Measure-aware coefficients.** With weights `w_i = (Δs_prev + Δs_next)/2`: foundation `k w_i`; edge stiffness
  `EA / L0_e`, viscosity `c / L0_e`, per edge. The substructure stores coefficients per edge, so a non-uniform ring
  needs no workaround (the flex reference stores them per flex, which is what blocked this in BikeLab).
- **Mass consistency.** Arc-weighted lumped masses keep the total but can move the rest centre of mass of an
  asymmetrically sampled circle. A non-uniform ring must pass independent mass, first-moment and spin-inertia checks;
  a rotating imbalance must never appear as tyre behaviour.
- **Order of delivery.** Uniform rings first; then fixed non-uniform rings chosen *between runs* (a lever, judged
  against uniform controls); in-run refinement only after conservation of mass, centre of mass, momentum and angular
  momentum across a topology change is shown, with transfer energy in its own ledger (LOD.md "state adapters").

## 3. Contact

### 3.1 Geometry

- **Pairs.** A casing edge is a capsule (segment plus skin). Terrain is a prepared profile whose collinear runs are
  merged (author subdivision must not add support), split into convex pieces. One contact per (edge, convex piece)
  pair that overlaps or lies within the activation margin: this is the flex reference's multiplicity, and it is
  independent of the terrain's subdivision.
- **Closest features, not averages.** The distance between a segment and a convex polygon is computed in closed
  form. Its witness (TERRAIN.md) is a vertex with the true approach direction, or an edge with its face normal. At an
  exposed convex corner the admissible normals span the neighbouring faces' cone; the witness picks the one the
  geometry gives and never averages. A point inside a casing edge distributes to its two nodes by barycentric weights
  (for a dim-1 flex the reference's normalised inverse-distance weights are the same numbers).
- **Candidates.** A periodic interval hierarchy over material arc (ADAPTIVE_WHEEL) bounds groups of edges by their
  actual deformed positions; only groups whose bound meets the terrain window are tested. Measure the plain scan
  first: at N = 64 it may be cheaper than the hierarchy.
- **Certified skip.** A group whose clearance at the start of the step exceeds its largest possible movement in the
  step (node displacement plus the host's rotation, bounded as `CorridorBuilder` bounds rotation) cannot touch: no
  contact work, no detail needed.

### 3.2 Normal law

| Option | Description | Verdict |
|---|---|---|
| Soft-constraint emulation | The reference engine's regularised constraint rows: an impedance from `solimp`, stiffness and damping from `solref`, a regulariser from a diagonal approximation of `J M⁻¹ Jᵀ`, solved as its convex problem. The problem is strictly convex, so its solution is unique, and any solver converged to a tight tolerance reproduces it | **Reference model** (ADR 0003, decision 4: "an equivalent regularised contact law as one selectable law"). Lives here as a generic law, selected per contact |
| Compliant contact element | A scalar law on the signed gap along the witness normal (`UnilateralSpringDamper`, or a Hunt–Crossley variant), bound to the contact by a contact element: the analogue of the tendon element, with the gap's Jacobian, `Jₙᵀ (dF/dg) Jₙ` and the geometric term from the normal's rotation at a vertex witness | **Improved model.** Implicit through its Jacobians, no solver iterations, and the compliance is a declared tread property instead of a numerical regulariser. Keeps MULTIBODY §1's rule that tyre and soil contact are force elements, which the soil working set needs (its deformation changes the surface within the step) |
| Rigid unilateral rows | Exact non-penetration of the nodes as LCP rows | **Rejected for casing nodes.** A light node meeting rigid ground loses its normal velocity at every touchdown, so dissipation scales with the mesh, not with a material. Kept where a rigid stop is the physics (a rim strike, a travel end) |

### 3.3 Friction

| Option | Description | Verdict |
|---|---|---|
| Regularised Coulomb (soft constraint, elliptic cone) | Part of the soft-constraint emulation. Creeps under sustained load; its initial slip stiffness is set by the regulariser | **Reference model** for parity |
| Regularised Coulomb as a force element | A smooth slip-velocity law with Jacobians | Cheapest; the same creep. A lever for ambient models |
| **Exact Coulomb as bounded rows** | One tangential row per contact at velocity level, `|f_t| ≤ μ f_n` with `f_n` from the normal element, solved with the step's own factorised matrix by the M1c machinery (bounded rows in a mixed LCP) | **Improved model.** In the planar formulation the friction cone is an interval, so the LCP is exact, with no linearisation. A node in sticking contact has exactly zero slip, so the slip stiffness at small slip comes from the casing's tangential compliance (a brush-like adhesion and sliding zone emerges from the ring), not from a regulariser |

Notes:

- `f_n` for the bound is taken at the start of the step (first order), or iterated with the step to a residual
  target (a continuous lever). The friction work is the rows' impulse times velocity: never positive, counted once,
  as dissipation.
- In 3D the cone is a disc: a linearised pyramid with a facet-count lever, or an exact second-order-cone projection.
  That choice is M4's.
- Soil interfaces are shear-displacement laws with history (force elements). Per contact, the terrain material
  selects the tangential model, so rock and hardpack can use exact Coulomb while soil uses its own law.
- Rows per wheel are few: a 64-node casing (3.5 cm edges) on flat ground is expected to have a handful of edges in
  contact (to be counted on the reference runs). The Delassus matrix comes from
  one triangular sweep per row through the condensed factor (MULTIBODY M1c, post-parity ideas). Direct pivoting is
  the default; projected Gauss–Seidel with a residual target is the lever for many rows (3D, wide patches).

### 3.4 Swept contact and events

- **Onset is an event.** A contact law kinks where the gap reaches zero; a stepper that steps over the kink loses an
  order (as the station and rig studies showed). The step is split at the first contact of any edge group, found by
  regula falsi on the group's gap (as M1c locates limits), to an `eventTolerance`.
- **No tunnelling by construction.** The certified-skip bound of §3.1 also proves no contact was missed: the gap is
  1-Lipschitz in the movement of the casing, and the movement bound includes rotation. Where the bound cannot certify,
  the step is split. This covers the deforming casing within the step's motion, which the reference's discrete
  checks do not.
- **Rim strike** stays a rigid swept query of the rim circle (`PreparedProfile::sweepCircle`) along the hub's path,
  with the same rotation-aware bound for curved paths.

### 3.5 Sharp corners and faceted ground

- Forces act along per-contact witness normals. A vehicle's support frame (the load-weighted mean normal) is the
  consumer's, for terrain extraction and prediction; it never replaces the contact normals.
- Corner contacts are tested for admissibility (the witness normal lies in the corner's cone) and for witness
  position on the authored terrain, at declared tolerances.
- Mesh sensitivity at a corner is a casing-model property. Refining terrain or contact search does not repair it;
  the levers are N and non-uniform spacing near the patch (§2.5), each judged against a converged nodal reference.

## 4. Drive, brakes and command inputs

- **Drive** is a generalized torque on the hub's hinge, relative to its carrier. The equal and opposite reaction on
  the carrier is automatic in reduced coordinates. Its work is `τ · Δθ_relative`, already the stepper's external work.
- **Brake, reference:** the reference engine's regularised friction loss on the hinge, capacity from the command,
  emulated with the soft-constraint law. It creeps and can show small positive instantaneous power.
- **Brake, improved:** a bounded velocity row on the hinge, `|τ_b| ≤ C`. When the bound is not reached the relative
  spin rate after the step is exactly zero (to the solve's rounding): an exact wheel lock. Over capacity it slips with
  `τ_b = −C sign(ω)`. Its work is dissipation, counted once. Brake and tyre-friction rows are solved together, so a
  locked wheel sliding on the ground is one consistent complementarity problem.
- **Bounded rows are velocity-level only.** They need no position projection, so they can be admitted in the linearly
  implicit midpoint step, which refuses closures until a RATTLE-type form exists (MULTIBODY M1d).
- **Inputs.** Row bounds (brake capacities) and generalized torques belong to `Forcing`, held over the step like
  `tau`. A consumer's command schedule maps the step's tick to a `Forcing`; the core keeps no command state, so a
  rolled-back step needs nothing restored beyond `State`.

## 5. Integration, determinism and state

- **Steppers.** The reference uses the soft-constraint emulation with MuJoCo's scheme (`ImplicitTerms::DampingOnly`,
  forces held, constraint forces solved in the forward pass and held through the velocity update: to be confirmed
  against the pinned source, see the consumer's plan). Improved models use the linearly implicit midpoint step on
  smooth phases, event splitting at contact and stick–slip switches, and an error-controlled step (the Euler–midpoint
  embedded pair).
- **Unresolved casing modes.** The highest nodal modes (a node on its contact compliance, axial edge modes) are mesh
  artefacts. An L-stable step may damp them if the observed channels (vehicle motion, wheel loads, slip, work) stay
  within the profile's budget against a converged run. The global casing modes must be resolved.
- **What `State` holds.** `q`, `qd`, the tick, the energy flows, and anything that changes the next step's bits. The
  reference path is warm-start free (direct solves to tolerance), so nothing else. An iterative lever that warm-starts
  keeps its start in `State`, so a snapshot replays bit for bit. Persistent material history (friction or coating
  state) is `State`, keyed by material id.
- **Transactional.** A step that faults or leaves a declared support domain leaves `State` unchanged (MULTIBODY
  M1b). Support limits (casing strain, node displacement, overlap, rim clearance) are the consumer's guards over the
  accepted state; the core reports what it measured.
- **Arithmetic.** Rest geometry from canonical tables; the hub's rotation needs one sine and cosine per step. The
  decision is made (MULTIBODY §8): `sub0math::det::sincos`, bit-identical on every platform.
- **Bounded memory.** Contact rows and candidate lists have fixed capacities set at build. Overflow is a fault, never
  truncation.

## 6. Observability

- **Reactions include everything.** `computeReactions` adds the substructure's forces at the host and the contact
  forces at their points, so the bearing load, the carrier-on-wheel wrench and the load at any ancestor joint come
  from the same numbers that moved the bodies. Loads inside the casing (per node, per edge) are available too; the
  reference engine supports only complete-wheel cuts.
- **Per contact:** material ids (edge, nodes, weights), terrain piece and witness kind, point, normal, gap, normal and
  tangential force, slip velocity, and the row's state (stick, slip, separating).
- **Energy channels:** foundation and edge potential; edge and contact damping; friction and brake dissipation; drive
  work; the soft-constraint law's signed constraint work (it stores and releases energy, so it is not all
  dissipation); transfer losses of a model switch. Each enters the ledger once.

## 7. Levers

| Lever | Kind | Target and controller | Switchable in-run |
|---|---|---|---|
| Casing model (§2.4) | Discrete | On physical demand: contact activity, closing speed, predicted error | With a state transfer (decide, transfer, commit) |
| Node count N | Discrete | None until an error estimate exists (the footprint study's lesson) | Between runs; in-run only after the §2.5 gates |
| Non-uniform spacing | Continuous (spacing profile) | A local error indicator near the patch | Between runs first |
| Step size | Continuous | Error tolerance per channel; Euler–midpoint pair | Every step |
| Event tolerance (contact onset, stick–slip) | Continuous | Gap at engagement | Every step |
| Friction model | Discrete | None | Setup |
| Friction bound iteration | Continuous | Residual target on `f_n` | Every step |
| Row solver (pivoting, projected Gauss–Seidel) | Discrete, then continuous | Complementarity residual | Every step |
| Jacobian and factor reuse | Continuous (reuse count) | The W-method keeps first order with a stale Jacobian | Every step |
| Certified skip | Discrete (on or off) | The clearance bound | Every step |

Quality profiles map these to targets, not fixed values (LOD.md).

## 8. Cost estimates (to be measured)

Per step, two planar rings, N = 64, scalar double, one core. The condensation figure is the prototype's; the others
are estimates from the measured rig (about 0.35 µs per step for a five-coordinate tree) and the footprint and
clearance studies.

| Part | Estimate |
|---|---|
| Tree (six coordinates, planar): kinematics, CRBA, RNEA, factor | 0.5–1 µs |
| Node kinematics (world positions of 128 nodes) | < 1 µs |
| Ring elements (256 foundation springs, 128 edges with Jacobians) | 2–4 µs |
| Contact search and geometry, with the certified skip | 1–2 µs |
| Condensation, tree solve and back-substitution | 5–8 µs |
| Bounded rows (10–20 friction and brake rows) | 3–10 µs |
| **Whole step** | **about 12–25 µs** |

Per simulated second: at the reference's 31.25 µs (32 000 steps) about 0.4–0.8 s, offline. At a 125–250 µs average
step (improved models with event splitting) about 50–200 ms. Ladder models below the nodal ring are expected at
1–20 ms per simulated second. A 3D casing several nodes wide (blocks of 3W) costs about `N (3W)³`: roughly ten
times the planar figure at W = 3.

## 9. Verification

| Check | Applies to |
|---|---|
| Structured model against its expanded oracle at N = 8, 16: kinematics, energy, forces, Jacobians, one step, a short run, to rounding | Every substructure change |
| Mass, first moment, spin inertia (`hub + m (R² + skin²/2)` for the uniform ring) | Every ring, uniform or not |
| Free flight: linear and angular momentum, spin; energy without dissipation | Every stepper |
| Element Jacobians against finite differences, and `F·u = −dV/dt − dissipation` for every law | Edges, foundation, contact elements |
| Contact geometry against the reference engine's contact list on its own states (points, normals, gaps, multiplicity) | Reference contact |
| Constraint forces against the reference engine's on its own states (open loop, no chaos) | Soft-constraint emulation |
| Corner witness admissibility and position; no averaged normals | All contact |
| No tunnelling across speeds and step sizes; onset located within `eventTolerance` | Swept contact |
| Exact stick: a loaded wheel on an incline below the friction limit holds with zero creep; a braked wheel under capacity holds zero relative spin | Exact rows |
| Ledger closes at every step size; friction and brake work never positive | All |
| Allocation gate; paired benchmark against a hand-written ring step | All |

## 10. Proposed changes to the multibody contracts

**Accepted by the architect, 2026-10-21.** Each item lands in the delivery step that first needs it: substructures
and block storage in M2a; contact elements in M2c; bounded velocity rows in M3b and M3c; velocity-only rows in the
midpoint step when a constrained midpoint run is first needed. Owner decisions: planar first with 3D-ready storage,
and an Interactive budget of 2 ms per 16.7 ms frame.

1. **Substructures.** A model may attach structured substructures (first kind: the deformable ring) to a body. Their
   coordinates are generalized coordinates; the step eliminates them before the tree factorisation.
2. **Bounded velocity rows** join closures and limits in the constraint layer: per-row bounds from `Forcing` (brake
   capacity) or from an element's force (`μ f_n`). §1's rule stays for normal contact and soil (force elements);
   friction against rigid surfaces and brakes may be rows.
3. **Contact elements**: a scalar law on a signed gap along a witness normal, with exact Jacobians including the
   normal's rotation, as the tendon element is a law on a distance.
4. **Element accumulation** into the substructure's block storage when present; dense `n × n` only for small trees.
5. **Velocity-only rows in the midpoint step.**

## 11. Open questions

- Whether planar-first (exact Coulomb as an LCP, parity with the planar references) or 3D-first (finite width,
  lateral modes) is the first casing. Recommended: planar first, with storage and laws written for three components.
- The improved normal law's form and parameters (linear, Hunt–Crossley, a tread layer in series with the node), and
  whether a massless contact layer removes the node-on-contact mode that sets the step. M2d: revision 1 is linear and
  keeps the reference row's transient; switches of the witness feature dominate its cost, so a law distributed along
  the edge (no force jump between nodes) is the next candidate (M2d status).
- Whether `f_n` must be iterated for the friction bound, or the start-of-step value is within budget.
- Whether the hierarchy pays for itself at N = 64.
- The cost and accuracy of the rigid-ring model's state transfer at contact onset (the most demanding switch).

## M2a status (2026-10-27)

**Landed:** the ring substructure (section 10, items 1 and 4), planar, with its oracle.

- **Description** (`multibody/deformable_ring.hpp`): `DeformableRingDescription` on a host body: rest positions (3D,
  host frame, one out-of-plane offset), two orthonormal in-plane displacement axes, per-node mass and inertia, a
  foundation law per node (on each coordinate), an edge law and rest length per edge, skin radius, material ids per
  node and edge and a topology revision. `ModelBuilder::addRing` validates it (`std::expected`). The laws are
  `LinearSpringDamper<Length>` and are evaluated through `respond`, so the formula exists once. `uniformRing` builds
  the reference's ring from the canonical tables (`canonicalRingDirection`, `canonicalChordRatio`: the quarter-sine
  table, reflected, and the 128-node chord; N a divisor of 128, so the oracle sizes 8 and 16 come from the same table).
- **Conventions transcribed from the reference MJCF** (BikeLab `e73f6cc`, `wheelParts` in `dist/moving-wheel.js`, shared
  by `dist/trail-bike.js`): node i rests at `R (sin t_i, -cos t_i)` along (x, z) of the hub, `t_i = 2 pi i / N`
  (node 0 at the bottom); slides along the hub's x then z; foundation = joint stiffness `k ds` (no joint damping),
  relaxed at zero displacement; edges = flex edge `stiffness EA/L` and `damping c/L` with `L = R chordRatio(N)`, a
  linear spring on `l - l0` plus a dashpot on the edge rate (MuJoCo's flex edge passive force); node mass `M/N`, inertia
  `(M/N) skin^2/2` (it rotates with the hub). The law's reference `l0` defaults to the distance between the canonical
  rest nodes; MuJoCo computes `flexedge_length0` from world positions at `qpos0`, which may differ by an ulp, so the
  description takes it per edge and M2c can supply MuJoCo's.
- **Storage.** The ring's coordinates follow the tree's in `State` (`Model::dofCount` = `treeDofCount` + 2 per node).
  The Workspace's dense matrices are `treeDofCount` square; ring Jacobians live in `RingWorkspace` blocks (2 x 2
  diagonal and neighbour blocks of dtau/dw and dtau/dwd, exact, including the edges' geometric and rate-transport
  terms), the mass coupling to the host's planar motion in a per-node 2 x 3 border. Blocks are templated on the
  component count (`Block<Dim>`), so the spatial ring adds Dim = 3.
- **Dynamics.** Hooks in the planar algorithms: node kinematics after the tree's transforms, Newton-Euler over the
  nodes between the tree's passes (the nodes' force reaches the host), the ring's rigid inertia added to the host's
  composite and its border, the host's motion Jacobian over its chain. **Condensation:** per ring the step matrix's
  block `A_ww = m I - h D - h^2 K~` is factorised by the periodic block LDL^T (`solver/periodic_block_tridiagonal.hpp`,
  one fill column for the wrap), `Z = A_ww^-1 B` and `B^T Z` (3 x 3) are formed, `J^T (B^T Z) J` is removed from the tree's
  matrix over the host's ancestors (so the tree factorisation keeps its sparsity), the right-hand side is condensed the
  same way, the tree is solved, and `w = y - Z J dqd_tree`. The mass matrix's condensation serves forward dynamics,
  `solveMassMatrix`, semi-implicit Euler and RK; the step matrix's serves both linearly implicit steps, all three
  `ImplicitTerms`, the geometric fallback and the prepared step. Reactions add the nodes' inertia and weight to their
  host as an external wrench.
- **Not yet with rings** (refused, `InvalidInput`): constraints (build), the general formulation (build), forcing
  Jacobians (step). Contact elements and their border terms are M2c.

**Gate evidence** (`tests/test_multibody_ring.cpp`, GCC 16, GCC 16 and Clang 22 `-march=x86-64-v3`):

- *Oracle agreement at N = 8 and 16* (the same description as massless carriers, node bodies, coordinate elements and
  tendons, factorised densely): inverse dynamics (M and c) and forward dynamics to 1e-15 relative, `solveMassMatrix`,
  kinetic and potential energy (identical bits at these states), element forces, potential and dissipation to 1e-15,
  stiffness and damping Jacobians (block against dense) to 1e-16, joint reactions at every tree joint to 1e-13
  relative; one step to 1e-13 and 300 steps of implicit Euler (Full, DampingOnly), the implicit midpoint step and
  semi-implicit Euler to 1e-10, RK4 100 steps to 1e-10; prepared and one-call steps bit for bit.
- *Mass properties:* spin inertia `hub + m (R^2 + skin^2/2)` and the slide's mass to 1e-14; zero rest first moment.
- *Free flight:* (no gravity, no viscosity, midpoint at 20 us, 40 ms): linear momentum drift 1.6e-5 N s (the step
  holds M at the start, first order), energy 3.3e-5 relative.
- *Energy ledger with ring damping:* RK4 (dissipation integrated as a state) at 10 us over 20 ms: 4.16 J dissipated,
  residual 6.5e-10 J; the midpoint step's residual 1.6e-3, 4.4e-4, 1.3e-4 J at 40, 20, 10 us (order 1.7; the
  dissipation is integrated at the start of each step).
- *Modal frequency against the closed form:* uniform 16-node ring on a fixed hub, mode j = 3 of the edge branch
  `w^2 = (kf + ke (s^2 |1+z|^2 + c^2 |1-z|^2)) / m`: 859.542 rad/s closed form and measured (zero crossings, midpoint
  step at `w h` = 0.002) to 1e-5.
- *Allocation:* every step, forward and inverse dynamics, `solveMassMatrix`, elements, energy and reactions with two
  64-node rings allocate nothing (`Sub0Physics_MultibodyAllocationTests`).
- *Determinism:* the digest of three 2000-step runs (two 64-node rings, each implicit scheme) is `761d96d625139bc2`
  on all three builds.
- *Benchmark* (docs/PERFORMANCE.md): two 64-node rings on a five-coordinate tree, one implicit Euler step: library
  10.3-14.1 us against 8.0-9.1 us for the hand-written twin (GCC 16, five interleaved rounds, loaded host: about 1.3x),
  Clang 22 `-march=x86-64-v3` 11.2-17.6 against 8.2-12.7 us (about 1.35x); the expanded bodies take the dense path at
  0.8-1.4 ms, about 100 times the condensed step. The whole step is inside the section 8 estimate (12-25 us with
  contact and rows still to come).

**Deviations and decisions:**

- **No shared scalar sweep.** Section 2.3 asked for one periodic elimination over a block type shared with the scalar
  solver in `solver/`; the existing chain solvers (projected Gauss-Seidel, and an active-set method whose free runs are
  open chains) contain no periodic direct elimination, so the block solver is the one implementation and its scalar
  instantiation (Dim = 1) is tested for when a periodic scalar solve is needed.
- **Laws on the nodes are `LinearSpringDamper` only.** Other length laws can join through the same `respond` when a
  model needs them.
- **The step uses the symmetric part of the ring's K**, as for the tree (a W-method); the exact blocks stay in the
  workspace.
- **The border is the mass coupling only.** The ring's own elements do not depend on the host's pose; contact (M2c)
  will add its host-motion terms to the same border.
- **Levers recorded, not studied:** node count and a non-uniform spacing are levers of the nodal ring (per-node and
  per-edge coefficients already allow a non-uniform ring; the section 2.5 mass, first-moment and spin-inertia checks
  apply); the quasi-static, modal and rigid-ring models of section 2.4 are future cheaper models behind the same
  substructure (docs/LOD.md).
- **Not done here:** the MuJoCo adapter's ring (consumer side, test-only) and Jacobian checks by finite differences of
  the ring alone (the blocks are checked against the expanded tendons, whose Jacobians are checked against central
  differences in M1b).

## M2c status (2026-10-27)

**Landed:** the contact element (section 10, item 3) and the reference soft-constraint law (section 3.2), with the
ring-node border terms in the substructure's blocks. The consumer's moving wheel runs on them at parity with the
reference engine (its docs/CONNECTED_BIKE.md and docs/MUJOCO_REFERENCE.md).

**Gap kinematics, one implementation** (`multibody/ring_contact.hpp`). A contact acts at a point of a casing edge
along a unit normal. `ringContactJacobian` gives its row `J = n^T dp/dq` (the host's chain through each joint's motion
at the point, the two nodes' coordinates weighted): with the normal held, the gap's gradient (the witnesses' own motion
drops out by the envelope theorem). The reference attaches a flex contact to the node bodies with normalised
inverse-distance weights from the contact point; the element uses barycentric weights at its axis witness. Both are the
same row with different weights.

**The reference law, written once for every row kind** (`multibody/soft_constraint.hpp`). `softConstraintRow` turns a
row's position (a contact's gap, an equality's violation), its rate and the diagonal approximation of `J M^-1 J^T` into
the reference acceleration and the regularisation; `SoftConstraintSolver` solves all rows of a step together, contacts
unilateral (`SoftRowKind::Unilateral`) and equalities always active (`SoftRowKind::Equality`, for the soft closures M4
will need). Transcribed from MuJoCo 3.16.0 at commit
[a892750](https://github.com/google-deepmind/mujoco/tree/a892750516336d8358949f4518fcc4fbbd85beae), file
`src/engine/engine_core_constraint.c`:

| Step | MuJoCo | Here |
|---|---|---|
| Time constant raised to `2 h` (non-discrete integrators) | getsolparam, [:1780](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1780) | `softConstraintRow` |
| Impedance `d(g)` from solimp (flat when `d0 = dwidth`) | getimpedance, [:1851-1900](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1851) | same |
| `R = max(mjMINVAL, (1 - d) diagApprox / d)`, `K = 1 / (dwidth^2 tc^2 dr^2)`, `B = 2 / (dwidth tc)` | mj_makeImpedance, [:1930](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1930), [:1947](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1947), [:1956](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1956) | same |
| `aref = -B v - K d g` | mj_referenceConstraint, [:3195](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L3195) | same |
| `diagApprox`: the contact's body weights times `body_invweight0` (the world side weighs zero) | mj_diagApprox, [:1679-1697](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L1679); `body_invweight0` at qpos0, `engine_setconst.c` [:1250-1257](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_setconst.c#L1250) | `ringNodeInverseWeights` (the mean of the translational `J M^-1 J^T` diagonal over three axes, once at the reference configuration) |
| Row cost `1/2 D jar^2` for `jar < 0`, zero when satisfied; force `-D jar` | mj_constraintUpdate_impl, [:3239](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L3239), [:3288-3300](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_core_constraint.c#L3288) | The mixed LCP `f >= 0`, `w = (J M^-1 J^T + R) f + J a0 - aref >= 0`, `f w = 0`: its unique solution is the cost's minimiser. Murty's principal pivoting solves it exactly where the engine's Newton stops at its tolerance (1e-10) |

**How it combines with `implicitfast`** (the open risk of the consumer's CONNECTED_BIKE section 9, settled from the
source and by the open-loop check). The constraint forces are solved in the forward pass with `M` (the free
acceleration `a0 = M^-1 qfrc_smooth`), then held: the velocity update solves `(M - h D) dqd = h (qfrc_smooth +
qfrc_constraint)` ([engine_forward.c:1819-1845](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_forward.c#L1819)).
So the parity path is the equivalent regularised row solved at the start of the step and passed as a held generalized
force, not a force element with Jacobians (an implicit element would differ from the engine at first order). One detail
mattered: the engine **gathers its velocity derivative into the mass matrix's sparsity** (`mapD2M`,
[engine_forward.c:1842](https://github.com/google-deepmind/mujoco/blob/a892750516336d8358949f4518fcc4fbbd85beae/src/engine/engine_forward.c#L1842)),
so a flex edge's damping between neighbouring nodes (sibling bodies, no mass-matrix entry) is dropped from its implicit
matrix while each node's own block stays. `ImplicitTerms::DampingOnlyMassPattern` is that rule (a damper between
branches stays explicit). With plain `DampingOnly`, one step from the engine's state differed by 0.3% of the step's
velocity change in free flight; with the mass pattern, by 3e-14 m/s (rounding).

**The contact element** (`RingContactGeometry`, `measureRingContact`, `ringContactLoad`). A scalar law on the signed
gap along the witness normal (the law's coordinate is the penetration), applied at the casing's axis witness, with its
exact Jacobians in local coordinates (the host's planar motion in its frame and the edge's four node coordinates):
`K = J^T (dF/dg) J + F H + (dF/dgdot) J^T (H v)^T`, `D = J^T (dF/dgdot) J`. The gap's Hessian `H` holds the host's
rotation (the point's second derivative), and at a vertex witness the rotation of the normal and the witness's sliding
along the edge (`H = P^T (I - n n^T) P / rho + n^T P'' - phi_qt phi_tq / phi_tt`, the envelope over the axis
parameter). The host's chain adds its own curvature (`F C(grad g)` and its rate transport), mapped exactly by
`ringLoadsToGeneralized` (a dense oracle path). The steps take loads through `Forcing::ringLoads` (their work is
external work, as for every forcing): the chain's block into the tree's matrix (ancestor pairs, so the tree's sparsity
is kept), the nodes' blocks into per-node and per-edge contact blocks (`RingWorkspace::contactDiagonal`,
`contactUpper`), and the node-to-host coupling into the condensation's border (`RingCondensation::border`: the mass
coupling plus the loads' terms, the place section 10 reserved for them).

**Gate evidence** (`tests/test_multibody_ring_contact.cpp`, GCC 16, GCC 16 and Clang 22 `-march=x86-64-v3`):

- *The law:* `softConstraintRow` against the source's formulas term by term (flat and graded impedance, the `2 h` bound).
- *The solve:* 40 random sets of 1 to 6 rows on a 16-node ring (some with an equality row): forces non-negative on
  unilateral rows, and no perturbation of the forces lowers the primal cost.
- *The row:* `J qd` against a central difference of the attached point (off the axis, rotating with the host) to 1e-8.
- *The element's Jacobians:* central differences of the generalized force (geometry re-derived at each state) on a ring
  whose host sits on a chain with two revolute joints, for a face witness and two vertex witnesses with interior axis
  points: about 1e-9 of the stiffness's scale (bound 1e-8); dropping the chain's curvature term gives 4e-3.
- *The structured step with loads:* against a dense solve of the same step (dense `M` by inverse dynamics, the loads'
  dense Jacobians) to 1e-10, for `Full` and `DampingOnly`; prepared and one-call steps bit for bit; a later step
  without loads is unaffected by the earlier contact blocks.
- *The mass pattern:* the condensed ring and its expanded oracle agree over 300 steps of `DampingOnlyMassPattern`
  (1e-10; the oracle drops the cross-branch tendon damping by the same rule), and the mode differs from `DampingOnly`
  exactly when the edges are viscous.
- *A ring dropped on contact elements* (64 nodes, elastic, implicit midpoint at 200 us): energy plus contact potential
  within 6.6e-4 of the start (the step crosses each onset's kink; event location is M2d); run digest
  `11241068545413185626` on all three builds.
- *Allocation:* rows, the soft solve, loads and every stepper with loads allocate nothing.
- *Consumer parity* (the reference stepper on these pieces): open loop on the engine's states, one step to rounding
  (median 5e-13 of the step's change), geometry and forces at the engine's GJK/EPA tolerance; whole runs within the
  engine's own collision-tolerance envelope.

**Deviations and decisions:**

- **Parity law as rows, not an element.** Held rows reproduce `implicitfast` exactly; an implicit force element would
  not. The element is the improved model's.
- **Exact LCP instead of Newton.** The same minimiser; the engine stops at 1e-10, which the open-loop check sees as up
  to 1e-5 relative in the generalized constraint force (median).
- **One mass-matrix solve per row** for `J M^-1 J^T`; a hand-written twin forms it from the condensed Schur complement
  in O(k^2) (docs/PERFORMANCE.md, open gap). A block-diagonal ring block (the mass matrix always, the reference's step
  matrix) now skips the periodic sweep, bit for bit.
- **The structured step drops one term of the element's exact K:** the rate transport through the host chain's
  curvature on the nodes' rows (it does not factor through the host's motion). It is zero for hosts on prismatic chains
  and for a law without damping; the step is a W-method, so only stability could depend on it. The dense map keeps it.
- **Deep overlap at a vertex** (the casing's axis past the obstacle) omits the normal's rotation terms; outside the
  element's support.
- **Found, not caused here:** the M2a ring digest (`ring run digest`) is `761d96d625139bc2` on GCC 16 with and without
  `-march=x86-64-v3` but `b96b9f5285379502` on Clang 22 `-march=x86-64-v3` with Sub0Math f7316fe, already at this
  branch's base (b7c004b). The new contact digest and the consumer's moving-wheel digest agree on all three builds.
- **Not here:** the improved model's law and parameters, event-located onset and swept certification (M2d), friction
  rows (M3b).

## M2d status (2026-10-28)

**Landed:** the pieces of the improved moving wheel (section 3.2's contact element model, section 3.4's swept contact and
events, section 5's error-controlled step). The consumer's improved moving wheel (its docs/CONNECTED_BIKE.md, "M2d
status") runs on them beside the frozen reference; the reference path is unchanged bit for bit.

- **Error-controlled steps on a dyadic lattice** (`integrator/dyadic_step.hpp`; the consumer's time-base decision D4).
  Time is an integer count of fine units (`base / 2^S`); a step of level k is `base * 2^k`, starts aligned to its own
  length, never crosses a boundary the caller names (a command tick, the end) and lands on it. `DyadicStepController`
  changes the level by whole factors of two from fixed thresholds (`ratio * 2^(q+1) <= safety^(q+1)` grows), so it
  evaluates no power function. `advanceDyadically` drives any system that proposes a step with an error ratio and an
  event flag: an error above target is retried at the predicted level, an event (a switch inside the step) is bisected,
  both down to `minLevel`, where the step is accepted (the switch is then located to the floor step). Faults retry
  smaller. Decisions depend on the reports, the configuration and the position only.
- **Step doubling from one preparation** (`PreparedStepUse::Keep`): `completeLinearlyImplicitStep` may leave the
  preparation in place, so the coarse step and the first half step share it. Each completion is bit for bit what a
  fresh preparation gives (tested with ring contact loads, both schemes).
- **The swept two-sided certificate** (`ConvexTerrain::certifySweptRingClearance`, `SweptStepMotion`). Along a step's
  straight interpolant in the coordinates (hub centre and angle and hub-frame node displacements each linear), every
  casing point is within `tau R` of its start and `(1 - tau) R` of its end, so each edge's gap is at least
  `(g0 + g1 - R) / 2` at every instant, `R = hub travel + node travel + lever |dtheta|` (rotation uncapped so it stays
  proportional). Conservative in a 1001-instant property test on four fixtures, and it certifies twice as many random
  steps as the one-sided `certifyRingClearance` (26 against 12 of 160). Per-edge bounds are returned.
- **Hub-relative foundation viscosity** (`UniformRingSpec::foundationViscosity`, N s/m^2 per unit centreline length, a
  dashpot `c_f ds` beside each slide). Zero is the reference ring bit for bit; a non-zero value is an improved-model
  description value (the consumer's review decision D5), never a numerical fix.
- **The Delassus A/B** (`DelassusForm`, `SoftRowSupport`): the soft-constraint solve forms `J M^-1 J^T` either by a
  full mass-matrix solve per row (uniform, the reference's bits) or from the condensed Schur complement in O(k^2) per
  pair (structured: `t = j_t - J_h^T sum B^T A^-1 j_w`, `D_rc = t_r^T S^-1 t_c + shared-node terms`, with `t^T S^-1 t`
  from one half sweep of the tree factor). Both behind one call; three optimisation passes each; numbers per target in
  docs/PERFORMANCE.md ("Delassus A/B").

**The improved contact element** (the consumer's choice, recorded here because it answers section 11's open question
on the law). Each (edge, piece) pair in overlap carries a unilateral spring-damper on its gap at the casing's axis
witness (`measureRingContact`, `ringContactLoad` with the law's exact derivatives), so the contact is implicit through
its Jacobians in the linearly implicit midpoint step (`ImplicitTerms::Full`). Three findings from building it:

1. **The reference row is not a spring.** The soft row drives the gap's acceleration to `aref` (a critically damped
   response of time constant `tc`), but its force balances slow loads through `J a = d aref + (1 - d) J a0`, which is
   `1 / (1 - d)` (100 at `d = 0.99`) times stiffer than the transient: about 2.3e6 N/m and 4.6e3 N s/m per contact on
   this ring, a node-on-contact mode near 1.4 kHz, heavily overdamped. A force element matched to that is a stiff,
   unresolved mode that the midpoint step (A-stable, not L-stable) does not damp: 19 000 steps on flat ground and 0.25% of
   the energy created at a 1e-5 m target. Revision 1 of the improved law keeps the row's *transient* instead (`k = m_c
   / (d tc^2)`, `c = 2 m_c / tc` on the mid-edge effective mass `m_c = 2 m`), a declared tread compliance whose mode the
   step resolves. Its static deflection is larger than the reference's (a model-form difference, reported at L3, not a
   numerical one).
2. **The gap's minimising feature is a switch.** The gap of a capsule against a convex piece is a minimum over
   features: on a face the force jumps from one end node of the edge to the other as the edge rolls past level, at a
   vertex it turns onto the next one. Not locating those cost first order: the energy created fell only as the 0.3th
   power of the tolerance. Keying the switch on the witness feature and the axis end (as well as the pair and its
   pushing state) gave a tenfold drop; what remains scales with the event floor (first order in the floor step).
3. **Events dominate the cost.** On rolling contact every edge switches several times; each switch costs a bisection
   to the floor and a regrowth. A smooth (distributed) tread law along the edge, whose force does not jump between
   nodes, is the recorded next lever (section 11).

**Global casing modes.** The ring's translation and rotation against the hub change no edge length, so no edge viscosity
damps them; the reference ring has no foundation damping. The midpoint step adds no numerical damping: in the
consumer's controlled case (free flight under a constant hub force) the mode energy's late peak is 1.01 of its early
peak (the reference scheme's 0.97). Foundation viscosity damps them when the model asks for it.

**Not done here:** an L-stable step for unresolved modes (TR-BDF2 or a massless contact layer), the distributed tread
law, friction rows (M3b).
