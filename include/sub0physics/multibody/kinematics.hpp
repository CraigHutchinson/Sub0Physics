#pragma once

#include <expected>
#include <span>

#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/description.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

namespace sub0physics::multibody
{

/** Pose and velocity of a site in the world frame. */
struct SiteKinematics
{
    sub0math::RigidTransform<double, sub0math::f64::Length> pose;   ///< Site frame in the world.
    sub0math::Vec3<sub0math::f64::Velocity> velocity;               ///< Velocity of the site's origin.
    sub0math::Vec3<sub0math::f64::AngularVelocity> angularVelocity; ///< Angular velocity of its body.
};

/** Kinetic and gravitational potential energy of a model in a state. */
struct MechanicalEnergy
{
    sub0math::f64::Energy kinetic;   ///< Sum of 1/2 v . I v over the bodies.
    sub0math::f64::Energy potential; ///< Gravitational: -sum m g . c, zero with every centre of mass at the world origin.

    /** @return kinetic + potential. */
    [[nodiscard]] sub0math::f64::Energy total() const noexcept
    {
        return kinetic + potential;
    }
};

/**
 * Computes each body's transform from its parent, spatial velocity and world pose.
 *
 * Writes `parentToBody`, `bodyVelocity` and `bodyPose`.
 *
 * @param model The model.
 * @param state Coordinates and rates, sized for the model.
 * @param workspace Workspace built for the model.
 * @return Nothing, or InvalidInput if the workspace was built for another model or the state is not
 *         sized for this one.
 */
[[nodiscard]] std::expected<void, Fault> updateKinematics(const Model& model, const State& state,
                                                          Workspace& workspace) noexcept;

/**
 * Computes a site's pose and velocity.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param site A site of the model.
 * @param workspace Workspace built for the model; its kinematics are updated.
 * @return The site's world pose and velocity, or InvalidInput for a site that is not in the model, a
 *         workspace built for another model or a state that is not sized for this one.
 */
[[nodiscard]] std::expected<SiteKinematics, Fault> computeSiteKinematics(const Model& model, const State& state,
                                                                         SiteId site, Workspace& workspace) noexcept;

/**
 * Computes the translational Jacobian of a site's origin: its world velocity is `J qd`.
 *
 * @param model The model.
 * @param state Coordinates (rates are not used).
 * @param site A site of the model.
 * @param workspace Workspace built for the model; its kinematics are updated.
 * @param jacobian Output, row-major 3 x dofCount, in m per unit coordinate.
 * @return Nothing, or InvalidInput if `jacobian` is not 3 x dofCount, the site is not in the model, the
 *         workspace was built for another model or the state is not sized for this one.
 */
[[nodiscard]] std::expected<void, Fault> computeSiteJacobian(const Model& model, const State& state, SiteId site,
                                                             Workspace& workspace, std::span<double> jacobian) noexcept;

/**
 * Computes the kinetic and gravitational potential energy.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param workspace Workspace built for the model; its kinematics are updated.
 * @return The energies, or InvalidInput for a workspace built for another model or a state that is not
 *         sized for this one.
 */
[[nodiscard]] std::expected<MechanicalEnergy, Fault> computeEnergy(const Model& model, const State& state,
                                                                  Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
