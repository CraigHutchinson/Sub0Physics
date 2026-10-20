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
#include "multibody_planar.hpp"
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
    detail::composeWorldTransforms(model, workspace);
    for (std::size_t i = 0; i < model.bodyCount(); ++i)
    {
        workspace.bodyPose[i] = sub0math::toRigidTransform<sub0math::f64::Length>(workspace.worldToBody[i]);
    }
}

// Kinetic and potential energy of a planar model whose planar transforms, velocities and world transforms are in
// the workspace; the elements' energy terms are evaluated (planar).
std::expected<detail::EnergyTerms, Fault> planarEnergyTerms(const Model& model, const State& state,
                                                            Workspace& workspace) noexcept
{
    const PlanarModel& pm = model.planar();
    const PlanarKinematics& k = workspace.planar;
    double kinetic = 0.0;
    double potential = 0.0;
    for (std::size_t i = 0; i < pm.bodies.size(); ++i)
    {
        const PlanarInertia& in = pm.bodies[i].inertia;
        const PlanarMotion& v = k.velocity[i];
        kinetic += 0.5 * dot(v, in * v);
        // m times the world centre of mass: m o + R h.
        const PlanarTransform& x = k.worldToBody[i];
        const double cx = in.m * x.x + (x.c * in.hx - x.s * in.hy);
        const double cy = in.m * x.y + (x.s * in.hx + x.c * in.hy);
        potential -= pm.gravityX * cx + pm.gravityY * cy;
    }
    if (model.elementCount() > 0)
    {
        if (auto added = detail::evaluateElementsUnchecked(model, state, workspace, false, true, true); !added)
        {
            return std::unexpected(added.error());
        }
        potential += workspace.elementPotential;
    }
    return detail::EnergyTerms{kinetic, potential};
}

// A site's kinematics from the planar kinematics in the workspace.
SiteKinematics planarSiteKinematics(const Model& model, const Workspace& workspace, SiteId site) noexcept
{
    using sub0math::f64::AngularVelocity;
    using sub0math::f64::Length;
    using sub0math::f64::Velocity;
    const PlanarModel& pm = model.planar();
    const std::size_t s = site.value();
    const Pose& frame = model.siteFrame(site);
    const std::size_t body = model.siteBodies()[s];
    if (body == NoParent)
    {
        return SiteKinematics{frame, {}, {}};
    }
    const PlanarTransform& x = workspace.planar.worldToBody[body];
    const PlanarMotion& v = workspace.planar.velocity[body];
    const detail::PlanarPoint p = detail::planarToWorld(x, pm.siteX[s], pm.siteY[s]);
    // Velocity of the body point at the site, in body axes (v + w k x r), then rotated into the world.
    const double lx = v.x - v.w * pm.siteY[s];
    const double ly = v.y + v.w * pm.siteX[s];
    const sub0math::Vec3d position = detail::planarVectorToWorld(pm, p.x, p.y, pm.siteNormal[s]);
    const sub0math::Vec3d velocity = detail::planarVectorToWorld(pm, x.c * lx - x.s * ly, x.s * lx + x.c * ly, 0.0);
    const sub0math::Vec3d omega = detail::planarVectorToWorld(pm, 0.0, 0.0, v.w);
    // The site's orientation: its body's (a rotation about k by the world angle), then its own frame.
    const double half = 0.5 * std::atan2(x.s, x.c);
    sub0math::Quatd rotation{std::cos(half), 0.0, 0.0, 0.0};
    const double sh = std::sin(half);
    rotation.x = pm.normal == 0 ? sh : 0.0;
    rotation.y = pm.normal == 1 ? sh : 0.0;
    rotation.z = pm.normal == 2 ? sh : 0.0;
    return SiteKinematics{{rotation * frame.rotation,
                           {Length::fromSi(position.x), Length::fromSi(position.y), Length::fromSi(position.z)}},
                          {Velocity::fromSi(velocity.x), Velocity::fromSi(velocity.y), Velocity::fromSi(velocity.z)},
                          {AngularVelocity::fromSi(omega.x), AngularVelocity::fromSi(omega.y),
                           AngularVelocity::fromSi(omega.z)}};
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
    const Pose& frame = model.siteFrame(site);
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return SiteKinematics{frame, {}, {}}; // fixed in the world
    }
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

std::expected<void, Fault> computeSiteKinematics(const Model& model, const State& state,
                                                std::span<const SiteId> sites, std::span<SiteKinematics> out,
                                                Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (out.size() != sites.size())
    {
        return invalidInput("site kinematics need one output per site");
    }
    for (const SiteId site : sites)
    {
        if (site.value() >= model.siteCount())
        {
            return invalidInput("site is not a site of the model");
        }
    }
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, workspace);
        detail::planarWorldTransforms(model, workspace);
        for (std::size_t s = 0; s < sites.size(); ++s)
        {
            out[s] = planarSiteKinematics(model, workspace, sites[s]);
        }
        return {};
    }
    updateUnchecked(model, state, workspace);
    for (std::size_t s = 0; s < sites.size(); ++s)
    {
        const SiteId site = sites[s];
        const Pose& frame = model.siteFrame(site);
        const std::size_t body = model.siteBodies()[site.value()];
        if (body == NoParent)
        {
            out[s] = SiteKinematics{frame, {}, {}};
            continue;
        }
        const sub0math::SpatialMotiond& v = workspace.bodyVelocity[body];
        const sub0math::Mat3d& toBody = workspace.worldToBody[body].rotation;
        const sub0math::Vec3d local = v.linear + cross(v.angular, plain(frame.translation));
        const sub0math::Vec3d world = detail::bodyToWorld(toBody, local);
        const sub0math::Vec3d omega = detail::bodyToWorld(toBody, v.angular);
        using sub0math::f64::AngularVelocity;
        using sub0math::f64::Velocity;
        out[s] = SiteKinematics{workspace.bodyPose[body] * frame,
                                {Velocity::fromSi(world.x), Velocity::fromSi(world.y), Velocity::fromSi(world.z)},
                                {AngularVelocity::fromSi(omega.x), AngularVelocity::fromSi(omega.y),
                                 AngularVelocity::fromSi(omega.z)}};
    }
    return {};
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
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return {}; // fixed in the world: no coordinate moves it
    }
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
    std::expected<detail::EnergyTerms, Fault> terms;
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, workspace);
        detail::planarWorldTransforms(model, workspace);
        terms = planarEnergyTerms(model, state, workspace);
    }
    else
    {
        updateUnchecked(model, state, workspace);
        terms = detail::energyTerms(model, state, workspace);
    }
    if (!terms)
    {
        return std::unexpected(terms.error());
    }
    return EnergyAccount{.kinetic = sub0math::f64::Energy::fromSi(terms->kinetic),
                         .potential = sub0math::f64::Energy::fromSi(terms->potential),
                         .dissipated = state.dissipated,
                         .externalWork = state.externalWork};
}

namespace detail
{

std::expected<EnergyTerms, Fault> energyTerms(const Model& model, const State& state, Workspace& workspace) noexcept
{
    const auto inertias = model.inertias();
    const sub0math::Vec3d g = model.baseAcceleration().linear * -1.0;
    double kinetic = 0.0;
    double potential = 0.0;
    for (std::size_t i = 0; i < inertias.size(); ++i)
    {
        const sub0math::SpatialMotiond& v = workspace.bodyVelocity[i];
        kinetic += 0.5 * dot(v, inertias[i] * v);
        const sub0math::PluckerTransformd& x = workspace.worldToBody[i];
        // m c in world axes as m o + R h: no division by the mass, so a massless body adds exactly nothing.
        const sub0math::Vec3d moment = detail::bodyToWorld(x.rotation, inertias[i].firstMoment) +
                                       x.translation * inertias[i].mass;
        potential -= dot(g, moment);
    }
    if (model.elementCount() > 0)
    {
        if (auto added = detail::evaluateElementsUnchecked(model, state, workspace, false); !added)
        {
            return std::unexpected(added.error());
        }
        potential += workspace.elementPotential;
    }
    return EnergyTerms{kinetic, potential};
}

} // namespace detail

} // namespace sub0physics::multibody
