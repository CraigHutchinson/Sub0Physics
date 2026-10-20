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
#include "sub0physics/multibody/force_element.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"

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
 * Preallocated scratch and results for one model's kinematics and dynamics.
 *
 * Sized once from the model; the algorithms write into it and never resize it, so a step does not
 * allocate. Fields hold the results of the last algorithm that wrote them (documented on each
 * function) and are meant to be read, not resized. Matrices are row-major `dofCount x dofCount`.
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
          parentToBody(model.bodyCount()),
          bodyVelocity(model.bodyCount()),
          bodyAcceleration(model.bodyCount()),
          bodyForce(model.bodyCount()),
          compositeInertia(model.bodyCount()),
          worldToBody(model.bodyCount()),
          bodyPose(model.bodyCount()),
          massMatrix(dof * dof, 0.0),
          massFactor(dof * dof, 0.0),
          bias(dof),
          qdd(dof),
          elementForce(dof, 0.0),
          elementStiffness(dof * dof, 0.0),
          elementDamping(dof * dof, 0.0),
          elementSiteWrenches(2 * model.tendonCount()),
          implicitMatrix(dof * dof, 0.0),
          implicitRhs(dof, 0.0),
          elementVectors(3 * dof),
          elementScalars(dof * (dof + 2), 0.0),
          constraintImpulse(model.constraintRowCount(), 0.0),
          constraintForce(model.constraintRowCount(), 0.0),
          constraintSiteWrenches(model.constraintWrenchCount()),
          jointReactions(model.bodyCount()),
          constraintScratch(dof, model.constraintRowCount(), model.limits().size(), model)
    {
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

    std::vector<sub0math::PluckerTransformd> parentToBody; ///< Transform from parent frame to body frame, at q.
    std::vector<sub0math::SpatialMotiond> bodyVelocity;     ///< Spatial velocity of each body, in its frame.
    std::vector<sub0math::SpatialMotiond> bodyAcceleration; ///< Spatial acceleration (gravity included as -g).
    std::vector<sub0math::SpatialForced> bodyForce;         ///< Net spatial force transmitted by each joint.
    std::vector<sub0math::SpatialInertiad> compositeInertia; ///< Inertia of each body's subtree, in its frame.
    std::vector<sub0math::PluckerTransformd> worldToBody; ///< Transform from the world frame to each body frame.
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> bodyPose; ///< Each body frame in the world: `toRigidTransform(worldToBody)`.

    std::vector<double> massMatrix; ///< M(q) from CRBA; symmetric; zero except ancestor pairs (zeroed at construction).
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
    ConstraintScratch constraintScratch;    ///< Solver scratch; internal.
};

} // namespace sub0physics::multibody
