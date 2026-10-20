#include <cstddef>
#include <expected>
#include <span>

#include <sub0math/geometry/point.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "multibody_joint.hpp"
#include "sub0physics/multibody/kinematics.hpp"

namespace sub0physics::multibody
{

namespace
{

using Pose = sub0math::RigidTransform<double, sub0math::f64::Length>;

sub0math::Vec3d plain(const sub0math::Vec3<sub0math::f64::Length>& v) noexcept
{
    return {v.x.si(), v.y.si(), v.z.si()};
}

// Kinematics without the precondition checks, for the entry points that have made them. One pose
// pipeline (docs/MULTIBODY.md, M1b): each joint's transform is evaluated once, world transforms are
// composed from them, and the world pose is recovered from the world transform.
void updateUnchecked(const Model& model, const State& state, Workspace& workspace) noexcept
{
    detail::updateTransforms<true>(model, state, workspace);
    const auto parents = model.parents();
    for (std::size_t i = 0; i < parents.size(); ++i)
    {
        const std::size_t p = parents[i];
        workspace.worldToBody[i] =
            p == NoParent ? workspace.parentToBody[i] : workspace.parentToBody[i] * workspace.worldToBody[p];
        workspace.bodyPose[i] = sub0math::toRigidTransform<sub0math::f64::Length>(workspace.worldToBody[i]);
    }
}

} // namespace

std::expected<void, Fault> updateKinematics(const Model& model, const State& state, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    updateUnchecked(model, state, workspace);
    return {};
}

std::expected<SiteKinematics, Fault> computeSiteKinematics(const Model& model, const State& state, SiteId site,
                                                           Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return std::unexpected(checked.error());
    }
    if (site.value() >= model.siteCount())
    {
        return invalidInput("site is not a site of the model");
    }
    updateUnchecked(model, state, workspace);
    const std::size_t body = model.siteBody(site).value();
    const Pose& frame = model.siteFrame(site);
    const sub0math::SpatialMotiond& v = workspace.bodyVelocity[body];
    const sub0math::Mat3d& toBody = workspace.worldToBody[body].rotation;
    // Velocity of the body point at the site, in body axes, then rotated into the world.
    const sub0math::Vec3d local = v.linear + cross(v.angular, plain(frame.translation));
    const sub0math::Vec3d world = detail::bodyToWorld(toBody, local);
    const sub0math::Vec3d omega = detail::bodyToWorld(toBody, v.angular);
    using sub0math::f64::AngularVelocity;
    using sub0math::f64::Velocity;
    return SiteKinematics{workspace.bodyPose[body] * frame,
            {Velocity::fromSi(world.x), Velocity::fromSi(world.y), Velocity::fromSi(world.z)},
            {AngularVelocity::fromSi(omega.x), AngularVelocity::fromSi(omega.y), AngularVelocity::fromSi(omega.z)}};
}

std::expected<void, Fault> computeSiteJacobian(const Model& model, const State& state, SiteId site,
                                               Workspace& workspace, std::span<double> jacobian) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    const std::size_t n = model.dofCount();
    if (jacobian.size() != 3 * n)
    {
        return invalidInput("site Jacobian must be 3 x dofCount");
    }
    if (site.value() >= model.siteCount())
    {
        return invalidInput("site is not a site of the model");
    }
    updateUnchecked(model, state, workspace);
    for (double& e : jacobian)
    {
        e = 0.0;
    }
    const std::size_t body = model.siteBody(site).value();
    const sub0math::PluckerTransformd& bodyTransform = workspace.worldToBody[body];
    const sub0math::Vec3d point =
        bodyTransform.translation + detail::bodyToWorld(bodyTransform.rotation, plain(model.siteFrame(site).translation));
    const auto parents = model.parents();
    const auto subspace = model.motionSubspaces();
    // Each ancestor joint moves the site by its world motion axis, taken at the site.
    for (std::size_t j = body; j != NoParent; j = parents[j])
    {
        const sub0math::PluckerTransformd& x = workspace.worldToBody[j];
        const sub0math::Vec3d omega = detail::bodyToWorld(x.rotation, subspace[j].angular);
        const sub0math::Vec3d linear = detail::bodyToWorld(x.rotation, subspace[j].linear);
        const sub0math::Vec3d column = linear + cross(omega, point - x.translation);
        const std::size_t c = model.coordinates()[j];
        jacobian[c] = column.x;
        jacobian[n + c] = column.y;
        jacobian[2 * n + c] = column.z;
    }
    return {};
}

std::expected<EnergyAccount, Fault> computeEnergy(const Model& model, const State& state,
                                                  Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return std::unexpected(checked.error());
    }
    updateUnchecked(model, state, workspace);
    const auto inertias = model.inertias();
    const sub0math::Vec3d g = model.baseAcceleration().linear * -1.0;
    double kinetic = 0.0;
    double potential = 0.0;
    for (std::size_t i = 0; i < inertias.size(); ++i)
    {
        const sub0math::SpatialMotiond& v = workspace.bodyVelocity[i];
        kinetic += 0.5 * dot(v, inertias[i] * v);
        const sub0math::PluckerTransformd& x = workspace.worldToBody[i];
        const sub0math::Vec3d centre = detail::bodyToWorld(x.rotation, inertias[i].centreOfMass()) + x.translation;
        potential -= inertias[i].mass * dot(g, centre);
    }
    return EnergyAccount{.kinetic = sub0math::f64::Energy::fromSi(kinetic), .potential = sub0math::f64::Energy::fromSi(potential)};
}

} // namespace sub0physics::multibody
