#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/multibody/constraint.hpp"
#include "sub0physics/multibody/deformable_ring.hpp"
#include "sub0physics/multibody/force_element.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/planar.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/solver/periodic_block_tridiagonal.hpp"

namespace sub0physics::multibody
{

/**
 * Preallocated scratch of the constraint solver, sized from the model's constraint layout. Internal; read the
 * results from Workspace::constraintReport, constraintImpulse and constraintForce.
 */
struct ConstraintScratch
{
    /**
     * Sizes the scratch.
     *
     * @param dof Degrees of freedom.
     * @param rows Rows of the constraint layout (closure rows plus one per one-sided limit).
     * @param limits One-sided limits.
     * @param model The model, for the scratch states.
     */
    ConstraintScratch(std::size_t dof, std::size_t rows, std::size_t limits, const Model& model)
        : jacobian(rows * dof, 0.0), solved(rows * dof, 0.0), delassus(rows * rows, 0.0), factor(rows * rows, 0.0),
          rhs(rows, 0.0), target(rows, 0.0), impulse(rows, 0.0), slack(rows, 0.0), residual(rows, 0.0), bias(rows, 0.0),
          vector(rows, 0.0), stepJacobian(rows * dof, 0.0), stepLayoutRow(rows, 0), layoutRow(rows, 0), subset(rows, 0),
          permutation(rows, 0), active(rows, 0), gap(limits, 0.0),
          engaged(limits, 0), elementForce(dof, 0.0), velocity(dof, 0.0), bodyExternal(dof), start(State::atRest(model)),
          trial(State::atRest(model))
    {
    }

    std::vector<double> jacobian;     ///< rows x dof, row-major, the assembled rows.
    std::vector<double> solved;       ///< rows x dof: `A^-1 J^T`, stored by row.
    std::vector<double> delassus;     ///< rows x rows: `J A^-1 J^T`.
    std::vector<double> factor;       ///< rows x rows: the pivoted factor.
    std::vector<double> rhs;          ///< rows.
    std::vector<double> target;       ///< rows.
    std::vector<double> impulse;      ///< rows.
    std::vector<double> slack;        ///< rows.
    std::vector<double> residual;     ///< rows.
    std::vector<double> bias;         ///< rows: `Jdot qd`.
    std::vector<double> vector;       ///< rows.
    std::vector<double> stepJacobian; ///< rows x dof: the rows at the sub-step's start.
    std::vector<std::size_t> stepLayoutRow; ///< rows: the layout row of each row at the sub-step's start.
    std::vector<std::size_t> layoutRow;  ///< rows: the layout row of each assembled row.
    std::vector<std::size_t> subset;     ///< rows: the rows of the LCP's active set.
    std::vector<std::size_t> permutation; ///< rows.
    std::vector<std::uint8_t> active;     ///< rows: assembled rows in the LCP's active set.
    std::vector<double> gap;          ///< limits: gap at the current state.
    std::vector<std::uint8_t> engaged; ///< limits: engaged at the sub-step's start.
    std::vector<double> elementForce; ///< dof: the elements' force at the sub-step's start.
    std::vector<double> velocity;     ///< dof.
    std::vector<sub0math::SpatialForced> bodyExternal; ///< dof: external wrenches on each body, body frame.
    State start;     ///< The sub-step's start.
    State trial;     ///< A trial end of the sub-step.
    double startDissipationPower = 0.0; ///< Elements' dissipation power at the sub-step's start [W].
    std::size_t stepRows = 0;           ///< Rows at the sub-step's start.
    std::size_t stepBilateral = 0;      ///< Of which closure rows (the first ones).
};


/**
 * The planar formulation's kinematics (Formulation::Planar), in the plane's axes: written by the dynamics, element,
 * energy and site queries of a planar model instead of the general fields. Empty for a general model.
 */
struct PlanarKinematics
{
    /**
     * Sizes the arrays for a model.
     *
     * @param bodies Bodies of a planar model, or zero for a general one.
     */
    explicit PlanarKinematics(std::size_t bodies)
        : parentToBody(bodies), worldToBody(bodies), velocity(bodies), acceleration(bodies), force(bodies),
          composite(bodies)
    {
    }

    std::vector<PlanarTransform> parentToBody; ///< Parent frame to body frame, at q.
    std::vector<PlanarTransform> worldToBody;  ///< World frame to body frame.
    std::vector<PlanarMotion> velocity;        ///< Each body's velocity, in its frame.
    std::vector<PlanarMotion> acceleration;    ///< Each body's acceleration (gravity included as -g), in its frame.
    std::vector<PlanarForce> force;            ///< Net force each joint transmits, in its body's frame.
    std::vector<PlanarInertia> composite;      ///< Inertia of each body's subtree, in its frame.
};

/**
 * One condensation of a ring's block of a step matrix onto its host's motion (docs/WHEEL.md section 2.3): the factor
 * of the ring's periodic block-tridiagonal block `A_ww`, `Z = A_ww^-1 B` for the border B (the ring's coupling to the
 * host's planar motion) and the host-motion matrix `B^T Z` that is removed from the tree's matrix. Internal.
 */
struct RingCondensation
{
    /**
     * Sizes the arrays.
     *
     * @param nodes Ring nodes.
     */
    explicit RingCondensation(std::size_t nodes)
        : inversePivot(nodes), next(nodes), wrap(nodes), coupling(nodes), border(nodes)
    {
    }

    /** @return The factor's views. */
    [[nodiscard]] PeriodicBlockFactor<PlanarRingComponents> factor() noexcept
    {
        return {inversePivot, next, wrap};
    }

    std::vector<Block<PlanarRingComponents>> inversePivot; ///< The factor's inverse pivots.
    std::vector<Block<PlanarRingComponents>> next;         ///< Multipliers onto the next node.
    std::vector<Block<PlanarRingComponents>> wrap;         ///< Multipliers onto the last node.
    std::vector<BlockColumns<PlanarRingComponents, 3>> coupling; ///< Z = A_ww^-1 B, per node (2 x 3: w, x, y of the host's motion).
    std::array<double, 9> hostMatrix{};                    ///< B^T Z, row-major 3 x 3 in the host frame (w, x, y).
    std::vector<BlockColumns<PlanarRingComponents, 3>> border; ///< The border B this condensation used: the mass coupling, plus the contact loads' terms for a step.
    bool blockDiagonal = false; ///< Every neighbour block was zero: `inversePivot` alone is the factor (no periodic sweep).
};

/**
 * Per-ring scratch and results (Workspace::rings), sized once from the ring. Ring quantities are in the host frame's
 * plane axes (positions, node motion) or in the ring's displacement axes (generalized forces and the element blocks).
 */
struct RingWorkspace
{
    /**
     * Sizes the arrays for a ring.
     *
     * @param ring The compiled ring.
     */
    explicit RingWorkspace(const DeformableRing& ring)
        : positionX(ring.nodeCount, 0.0), positionY(ring.nodeCount, 0.0), velocity(ring.nodeCount),
          stiffnessDiagonal(ring.nodeCount), stiffnessUpper(ring.nodeCount), dampingDiagonal(ring.nodeCount),
          dampingUpper(ring.nodeCount), border(ring.nodeCount), hostJacobian(ring.hostChain.size()),
          stepDiagonal(ring.nodeCount), stepUpper(ring.nodeCount), solution(ring.nodeCount), mass(ring.nodeCount),
          step(ring.nodeCount), contactDiagonal(ring.nodeCount), contactUpper(ring.nodeCount),
          contactBorder(ring.nodeCount)
    {
    }

    std::vector<double> positionX;       ///< Each node's position along the plane's i axis, host frame [m].
    std::vector<double> positionY;       ///< Along j [m].
    std::vector<PlanarMotion> velocity;  ///< Each node's velocity in its frame (the host's axes at the node).
    PlanarForce hostForce{};             ///< From the last Newton-Euler pass: the nodes' net force on the host, host frame (their inertia and gravity, with qdd when given).

    // Element Jacobians in block storage (docs/WHEEL.md section 10, item 4), in the displacement axes: exact, like
    // Workspace::elementStiffness. Block (i, i + 1) and block (i + 1, i) of dtau/dw are both `stiffnessUpper[i]`
    // (an edge's Jacobian blocks are G and -G; not transposed, as a damped edge's G is not symmetric).
    std::vector<Block<PlanarRingComponents>> stiffnessDiagonal; ///< dtau_i/dw_i.
    std::vector<Block<PlanarRingComponents>> stiffnessUpper;    ///< dtau_i/dw_{i+1} = dtau_{i+1}/dw_i.
    std::vector<Block<PlanarRingComponents>> dampingDiagonal;   ///< dtau_i/dwd_i (symmetric).
    std::vector<Block<PlanarRingComponents>> dampingUpper;      ///< dtau_i/dwd_{i+1} = dtau_{i+1}/dwd_i (symmetric).

    std::vector<BlockColumns<PlanarRingComponents, 3>> border; ///< The mass matrix's coupling of each node's coordinates to the host's planar motion (w, x, y; host frame).
    std::vector<PlanarMotion> hostJacobian; ///< Per DeformableRing::hostChain entry: the host's motion (host frame) per unit rate of that body's joint.

    std::vector<Block<PlanarRingComponents>> stepDiagonal; ///< Scratch: the step matrix's ring blocks.
    std::vector<Block<PlanarRingComponents>> stepUpper;    ///< Scratch.
    mutable std::vector<BlockColumns<PlanarRingComponents, 1>> solution; ///< Scratch: one right-hand side (a solve with a const workspace uses it).
    RingCondensation mass; ///< The mass matrix's condensation (forward dynamics, solveMassMatrix).
    RingCondensation step; ///< The step matrix's condensation (linearly implicit steps).

    // The contact loads' terms of the last step matrix (Forcing::ringLoads), already scaled: `-h D - h^2 K~`.
    std::vector<Block<PlanarRingComponents>> contactDiagonal; ///< Per node.
    std::vector<Block<PlanarRingComponents>> contactUpper;    ///< Per edge: the block (e, e + 1); its transpose is (e + 1, e).
    std::vector<BlockColumns<PlanarRingComponents, 3>> contactBorder; ///< Per node: added to the border.
    bool contactBlocks = false; ///< Whether the contact terms are set (they are zero otherwise).
};

/**
 * Preallocated scratch and results for one model's kinematics and dynamics.
 *
 * Sized once from the model; the algorithms write into it and never resize it, so a step does not
 * allocate. For a model with Formulation::Planar the dynamics, element, energy and site queries write the compact
 * `planar` kinematics instead of the general kinematic fields (parentToBody, bodyVelocity, bodyAcceleration,
 * bodyForce, compositeInertia, worldToBody, bodyPose); updateKinematics, computeReactions and the constrained steps
 * write the general ones. Fields hold the results of the last algorithm that wrote them (documented on each
 * function) and are meant to be read, not resized. Matrices are row-major `treeDofCount x treeDofCount`
 * (the tree's coordinates; equal to dofCount for a model without rings); a ring's coordinates live in its
 * RingWorkspace's block storage.
 * Use one workspace per thread.
 *
 * A workspace is bound to the model it was built for (not merely to its size): every algorithm
 * refuses a different model with an InvalidInput fault. The mass matrix is zero outside the
 * ancestor pairs of that model's tree, and the algorithms only ever write the ancestor entries.
 */
struct Workspace
{
    /**
     * Sizes a workspace for a model.
     *
     * @param model The model it will serve.
     */
    explicit Workspace(const Model& model)
        : modelId(model.id()),
          dof(model.dofCount()),
          treeDof(model.treeDofCount()),
          parentToBody(model.bodyCount()),
          bodyVelocity(model.bodyCount()),
          bodyAcceleration(model.bodyCount()),
          bodyForce(model.bodyCount()),
          compositeInertia(model.bodyCount()),
          worldToBody(model.bodyCount()),
          bodyPose(model.bodyCount()),
          massMatrix(treeDof * treeDof, 0.0),
          massFactor(treeDof * treeDof, 0.0),
          bias(dof),
          qdd(dof),
          elementForce(dof, 0.0),
          elementStiffness(treeDof * treeDof, 0.0),
          elementDamping(treeDof * treeDof, 0.0),
          elementSiteWrenches(2 * model.tendonCount()),
          implicitMatrix(treeDof * treeDof, 0.0),
          implicitRhs(dof, 0.0),
          elementVectors(3 * treeDof),
          elementScalars(treeDof * (treeDof + 2), 0.0),
          constraintImpulse(model.constraintRowCount(), 0.0),
          constraintForce(model.constraintRowCount(), 0.0),
          constraintSiteWrenches(model.constraintWrenchCount()),
          jointReactions(model.bodyCount()),
          appliedForce(dof, 0.0),
          preparedValues(2 * dof, 0.0),
          planar(model.formulation() == Formulation::Planar ? model.bodyCount() : 0U),
          constraintScratch(dof, model.constraintRowCount(), model.limits().size(), model)
    {
        rings.reserve(model.rings().size());
        for (const DeformableRing& ring : model.rings())
        {
            rings.emplace_back(ring);
        }
    }

    /**
     * Checks that this workspace was built for a model.
     *
     * The model's id fixes the topology and therefore every array size, so two comparisons suffice as
     * long as the arrays are not resized, which their documentation forbids.
     *
     * @param model The model about to be used with it.
     * @return True when the model is the one (or a copy of the one) the workspace was built for.
     */
    [[nodiscard]] bool isBuiltFor(const Model& model) const noexcept
    {
        return modelId == model.id() && dof == model.dofCount();
    }

    std::uint64_t modelId; ///< Id of the model this workspace was built for.
    std::size_t dof;       ///< Number of generalized coordinates of the model.
    std::size_t treeDof;   ///< Number of the tree's coordinates: the size of the dense matrices below.

    std::vector<sub0math::PluckerTransformd> parentToBody; ///< Transform from parent frame to body frame, at q.
    std::vector<sub0math::SpatialMotiond> bodyVelocity;     ///< Spatial velocity of each body, in its frame.
    std::vector<sub0math::SpatialMotiond> bodyAcceleration; ///< Spatial acceleration (gravity included as -g).
    std::vector<sub0math::SpatialForced> bodyForce;         ///< Net spatial force transmitted by each joint.
    std::vector<sub0math::SpatialInertiad> compositeInertia; ///< Inertia of each body's subtree, in its frame.
    std::vector<sub0math::PluckerTransformd> worldToBody; ///< Transform from the world frame to each body frame.
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> bodyPose; ///< Each body frame in the world: `toRigidTransform(worldToBody)`.

    std::vector<double> massMatrix; ///< The tree's block of M(q) from CRBA (a ring's inertia included in its host's); symmetric; zero except ancestor pairs (zeroed at construction). A ring's border is in RingWorkspace::border and its own block is diagonal (the node masses).
    std::vector<double> massFactor; ///< Tree-sparse factor of M (lower triangle: L below the diagonal, D on it).
    std::vector<double> bias;       ///< c(q, qd): Coriolis, centrifugal and gravity forces from RNEA.
    std::vector<double> qdd;        ///< Generalized accelerations from forward dynamics.

    // Force elements (evaluateElements): sums over every element of the model at the last evaluated state.
    std::vector<double> elementForce;     ///< Generalized force tau of the elements [N or N m].
    std::vector<double> elementStiffness; ///< dtau/dq, row-major; exact (not symmetric for a damped tendon); written only when Jacobians are asked for; the material part only after an implicit step's fallback (geometricStiffnessDropped).
    std::vector<double> elementDamping;   ///< dtau/dqd, row-major, symmetric; written only when Jacobians are asked for.
    double elementPotential = 0.0;        ///< Stored energy of the elements [J].
    double elementDissipationPower = 0.0; ///< Power the elements remove [W]; never negative.
    std::vector<ElementSiteWrench> elementSiteWrenches; ///< Two per tendon (first end, then second), in the model's order.

    // The linearly implicit Euler step (stepLinearlyImplicitEuler).
    std::vector<double> implicitMatrix; ///< `M - h D - h^2 K`, then its dense LDL^T factor (lower triangle).
    std::vector<double> implicitRhs;    ///< Right-hand side, then the velocity change.
    bool geometricStiffnessDropped = false; ///< Whether the last implicit step fell back to the material stiffness because `M - h D - h^2 K` was not positive definite.

    // Scratch for the tendon geometry.
    std::vector<sub0math::Vec3d> elementVectors; ///< 3 dofCount vectors.
    std::vector<double> elementScalars;          ///< dofCount (dofCount + 2) values.

    // Constraints (docs/MULTIBODY.md section 6), in the model's constraint layout (Model::constraintRowCount):
    // every closure's rows, then one row per one-sided limit.
    std::vector<double> constraintImpulse; ///< Impulse of each row over the last step [N s or N m s], summed over its sub-steps; closure rows in world axes, limit rows along the gap (never negative).
    std::vector<double> constraintForce;   ///< Force of each row from the last forwardDynamics or computeReactions [N or N m]: the acceleration-level solve.
    std::vector<ConstraintSiteWrench> constraintSiteWrenches; ///< From computeReactions: two per closure and per distance limit (first site, then second).
    std::vector<JointReaction> jointReactions; ///< From computeReactions: the wrench each joint transmits to its body.
    ConstraintReport constraintReport{};    ///< What the last step or projection achieved against the requested targets.

    std::vector<double> appliedForce; ///< The last step's applied generalized force: Forcing::tau plus the site forces mapped through their Jacobians [N or N m].
    std::vector<double> preparedValues; ///< The state a prepared step was prepared at (prepareLinearlyImplicitStep); internal.
    bool stepPrepared = false;          ///< Whether a prepared step is waiting for completeLinearlyImplicitStep; internal.
    ImplicitTerms preparedTerms = ImplicitTerms::Full; ///< The prepared step's implicit terms; internal.
    PlanarKinematics planar;          ///< The planar formulation's kinematics (Formulation::Planar); empty otherwise.
    ConstraintScratch constraintScratch;    ///< Solver scratch; internal.
    std::vector<RingWorkspace> rings;       ///< One per ring of the model, in its order.
};

} // namespace sub0physics::multibody
