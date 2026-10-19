#pragma once

#include <cstddef>
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
 */
struct Workspace
{
    /**
     * Sizes a workspace for a model.
     *
     * @param model The model it will serve.
     */
    explicit Workspace(const Model& model)
        : dof(model.dofCount()),
          parentToBody(model.bodyCount()),
          bodyVelocity(model.bodyCount()),
          bodyAcceleration(model.bodyCount()),
          bodyForce(model.bodyCount()),
          compositeInertia(model.bodyCount()),
          bodyPose(model.bodyCount()),
          massMatrix(dof * dof),
          massFactor(dof * dof),
          bias(dof),
          qdd(dof),
          scratch(dof)
    {
    }

    std::size_t dof; ///< Number of generalized coordinates of the model.

    std::vector<sub0math::PluckerTransformd> parentToBody; ///< Transform from parent frame to body frame, at q.
    std::vector<sub0math::SpatialMotiond> bodyVelocity;     ///< Spatial velocity of each body, in its frame.
    std::vector<sub0math::SpatialMotiond> bodyAcceleration; ///< Spatial acceleration (gravity included as -g).
    std::vector<sub0math::SpatialForced> bodyForce;         ///< Net spatial force transmitted by each joint.
    std::vector<sub0math::SpatialInertiad> compositeInertia; ///< Inertia of each body's subtree, in its frame.
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> bodyPose; ///< Each body frame in the world.

    std::vector<double> massMatrix; ///< M(q) from CRBA; symmetric.
    std::vector<double> massFactor; ///< Tree-sparse factor of M (lower triangle: L below the diagonal, D on it).
    std::vector<double> bias;       ///< c(q, qd): Coriolis, centrifugal and gravity forces from RNEA.
    std::vector<double> qdd;        ///< Generalized accelerations from forward dynamics.
    std::vector<double> scratch;    ///< Internal scratch, one entry per coordinate.
};

} // namespace sub0physics::multibody
