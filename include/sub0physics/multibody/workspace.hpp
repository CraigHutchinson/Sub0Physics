#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/multibody/model.hpp"

namespace sub0physics::multibody
{

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
          bodyPose(model.bodyCount()),
          massMatrix(dof * dof, 0.0),
          massFactor(dof * dof, 0.0),
          bias(dof),
          qdd(dof)
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
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> bodyPose; ///< Each body frame in the world.

    std::vector<double> massMatrix; ///< M(q) from CRBA; symmetric; zero except ancestor pairs (zeroed at construction).
    std::vector<double> massFactor; ///< Tree-sparse factor of M (lower triangle: L below the diagonal, D on it).
    std::vector<double> bias;       ///< c(q, qd): Coriolis, centrifugal and gravity forces from RNEA.
    std::vector<double> qdd;        ///< Generalized accelerations from forward dynamics.
};

} // namespace sub0physics::multibody
