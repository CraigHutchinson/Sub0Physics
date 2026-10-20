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
#include "sub0physics/verification/energy_account.hpp"

namespace sub0physics::multibody
{

/** Pose and velocity of a site in the world frame. */
struct SiteKinematics
{
    sub0math::RigidTransform<double, sub0math::f64::Length> pose;   ///< Site frame in the world.
    sub0math::Vec3<sub0math::f64::Velocity> velocity;               ///< Velocity of the site's origin.
    sub0math::Vec3<sub0math::f64::AngularVelocity> angularVelocity; ///< Angular velocity of its body.
};

/**
 * Computes each body's transform from its parent, spatial velocity and world pose.
 *
 * Writes `parentToBody`, `bodyVelocity`, `worldToBody` and `bodyPose`. There is one pose pipeline: each
 * joint's transform is evaluated once (one sin/cos pair of its coordinate), world transforms are composed
 * from them, and each world pose is `toRigidTransform(worldToBody)`.
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
 * Kinetic energy is the sum of `v . I v / 2` over the bodies; gravitational potential is `-sum m g . c`,
 * zero with every centre of mass at the world origin. Elements add their own potential, and a run
 * accumulates `dissipated` and `externalWork`; this call leaves those two zero.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param workspace Workspace built for the model; its kinematics are updated.
 * @return The account (kinetic and gravitational potential), or InvalidInput for a workspace built for
 *         another model or a state that is not sized for this one.
 */
[[nodiscard]] std::expected<EnergyAccount, Fault> computeEnergy(const Model& model, const State& state,
                                                                Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
