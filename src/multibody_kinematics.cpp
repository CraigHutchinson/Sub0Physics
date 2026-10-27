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
#include "multibody_ring.hpp"
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

// Kinetic and gravitational energy of a planar model from the planar velocities and world transforms in the workspace.
detail::EnergyTerms planarMotionEnergy(const Model& model, const Workspace& workspace) noexcept
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
    if (!model.rings().empty())
    {
        detail::ringMotionEnergy(model, workspace, kinetic, potential);
    }
    return {kinetic, potential};
}

// The same in the general formulation, from bodyVelocity and worldToBody.
detail::EnergyTerms generalMotionEnergy(const Model& model, const Workspace& workspace) noexcept
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
    return {kinetic, potential};
}

// A site's kinematics from the general kinematics in the workspace (bodyVelocity, worldToBody, bodyPose).
SiteKinematics generalSiteKinematics(const Model& model, const Workspace& workspace, SiteId site) noexcept
{
    const Pose& frame = model.siteFrame(site);
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return SiteKinematics{frame, {}, {}};
    }
    const sub0math::SpatialMotiond& v = workspace.bodyVelocity[body];
    const sub0math::Mat3d& toBody = workspace.worldToBody[body].rotation;
    const sub0math::Vec3d local = v.linear + cross(v.angular, plain(frame.translation));
    const sub0math::Vec3d world = detail::bodyToWorld(toBody, local);
    const sub0math::Vec3d omega = detail::bodyToWorld(toBody, v.angular);
    using sub0math::f64::AngularVelocity;
    using sub0math::f64::Velocity;
    return SiteKinematics{workspace.bodyPose[body] * frame,
                          {Velocity::fromSi(world.x), Velocity::fromSi(world.y), Velocity::fromSi(world.z)},
                          {AngularVelocity::fromSi(omega.x), AngularVelocity::fromSi(omega.y),
                           AngularVelocity::fromSi(omega.z)}};
}

// Kinetic and potential energy of a planar model whose planar transforms, velocities and world transforms are in
// the workspace; the elements' energy terms are evaluated (planar).
std::expected<detail::EnergyTerms, Fault> planarEnergyTerms(const Model& model, const State& state,
                                                            Workspace& workspace) noexcept
{
    detail::EnergyTerms terms = planarMotionEnergy(model, workspace);
    double& potential = terms.potential;
    const double kinetic = terms.kinetic;
    if (model.hasForceElements())
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
    // The site's orientation: its body's (a rotation about k by the world angle), then its own frame. The half-angle
    // quaternion from the cosine and sine, choosing the well-conditioned root.
    // PERF: half-angle from (c, s) beats atan2, cos and sin per site.
    double cw = 0.0;
    double sh = 0.0;
    if (x.c > -0.5)
    {
        cw = std::sqrt(0.5 * (1.0 + x.c));
        sh = x.s / (2.0 * cw);
    }
    else
    {
        sh = std::copysign(std::sqrt(0.5 * (1.0 - x.c)), x.s);
        cw = x.s / (2.0 * sh);
    }
    sub0math::Quatd rotation{cw, 0.0, 0.0, 0.0};
    rotation.x = pm.normal == 0 ? sh : 0.0;
    rotation.y = pm.normal == 1 ? sh : 0.0;
    rotation.z = pm.normal == 2 ? sh : 0.0;
    // Most sites are not rotated in their body: skip the product then.
    const sub0math::Quatd& own = frame.rotation;
    const bool unrotated = own.w == 1.0 && own.x == 0.0 && own.y == 0.0 && own.z == 0.0;
    return SiteKinematics{{unrotated ? rotation : rotation * own,
                           {Length::fromSi(position.x), Length::fromSi(position.y), Length::fromSi(position.z)}},
                          {Velocity::fromSi(velocity.x), Velocity::fromSi(velocity.y), Velocity::fromSi(velocity.z)},
                          {AngularVelocity::fromSi(omega.x), AngularVelocity::fromSi(omega.y),
                           AngularVelocity::fromSi(omega.z)}};
}

EnergyAccount accountOf(const detail::EnergyTerms& terms, const State& state) noexcept
{
    return EnergyAccount{.kinetic = sub0math::f64::Energy::fromSi(terms.kinetic),
                         .potential = sub0math::f64::Energy::fromSi(terms.potential),
                         .dissipated = state.dissipated,
                         .externalWork = state.externalWork};
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
                                                Workspace& workspace, EnergyAccount* energy) noexcept
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
        if (energy != nullptr)
        {
            const auto terms = planarEnergyTerms(model, state, workspace);
            if (!terms)
            {
                return std::unexpected(terms.error());
            }
            *energy = accountOf(*terms, state);
        }
        return {};
    }
    updateUnchecked(model, state, workspace);
    if (energy != nullptr)
    {
        const auto terms = detail::energyTerms(model, state, workspace);
        if (!terms)
        {
            return std::unexpected(terms.error());
        }
        *energy = accountOf(*terms, state);
    }
    for (std::size_t s = 0; s < sites.size(); ++s)
    {
        out[s] = generalSiteKinematics(model, workspace, sites[s]);
    }
    return {};
}

namespace
{

// The checks both prepared site queries make.
std::expected<void, Fault> checkPrepared(const Model& model, const State& state, std::span<const SiteId> sites,
                                         std::size_t outputs, const Workspace& workspace) noexcept
{
    if (!workspace.isBuiltFor(model) || !state.isSizedFor(model) || !workspace.stepPrepared ||
        !std::ranges::equal(state.values, workspace.preparedValues))
    {
        return invalidInput("the workspace holds no step prepared at this state");
    }
    if (outputs != sites.size())
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
    return {};
}

void preparedEnergy(const Model& model, const State& state, const Workspace& workspace, EnergyAccount* energy) noexcept
{
    if (energy != nullptr)
    {
        // The prepared step evaluated the elements at this state: their potential is in the workspace.
        const detail::EnergyTerms terms =
            detail::isPlanar(model) ? planarMotionEnergy(model, workspace) : generalMotionEnergy(model, workspace);
        *energy = accountOf({terms.kinetic, terms.potential + workspace.elementPotential}, state);
    }
}

} // namespace

std::expected<void, Fault> computePreparedSiteKinematics(const Model& model, const State& state,
                                                        std::span<const SiteId> sites, std::span<SiteKinematics> out,
                                                        const Workspace& workspace, EnergyAccount* energy) noexcept
{
    if (auto checked = checkPrepared(model, state, sites, out.size(), workspace); !checked)
    {
        return checked;
    }
    const bool planar = detail::isPlanar(model);
    for (std::size_t s = 0; s < sites.size(); ++s)
    {
        out[s] = planar ? planarSiteKinematics(model, workspace, sites[s]) : generalSiteKinematics(model, workspace, sites[s]);
    }
    preparedEnergy(model, state, workspace, energy);
    return {};
}

std::expected<void, Fault> computePreparedSiteMotion(const Model& model, const State& state,
                                                    std::span<const SiteId> sites, std::span<SiteMotion> out,
                                                    const Workspace& workspace, EnergyAccount* energy) noexcept
{
    if (auto checked = checkPrepared(model, state, sites, out.size(), workspace); !checked)
    {
        return checked;
    }
    using sub0math::f64::Length;
    using sub0math::f64::Velocity;
    if (detail::isPlanar(model))
    {
        const PlanarModel& pm = model.planar();
        for (std::size_t s = 0; s < sites.size(); ++s)
        {
            const std::size_t i = sites[s].value();
            const std::size_t body = model.siteBodies()[i];
            if (body == NoParent)
            {
                const auto& t = model.siteFrame(sites[s]).translation;
                out[s] = {t, {}};
                continue;
            }
            const PlanarTransform& x = workspace.planar.worldToBody[body];
            const PlanarMotion& v = workspace.planar.velocity[body];
            const detail::PlanarPoint p = detail::planarToWorld(x, pm.siteX[i], pm.siteY[i]);
            const double lx = v.x - v.w * pm.siteY[i];
            const double ly = v.y + v.w * pm.siteX[i];
            const sub0math::Vec3d position = detail::planarVectorToWorld(pm, p.x, p.y, pm.siteNormal[i]);
            const sub0math::Vec3d velocity =
                detail::planarVectorToWorld(pm, x.c * lx - x.s * ly, x.s * lx + x.c * ly, 0.0);
            out[s] = {{Length::fromSi(position.x), Length::fromSi(position.y), Length::fromSi(position.z)},
                      {Velocity::fromSi(velocity.x), Velocity::fromSi(velocity.y), Velocity::fromSi(velocity.z)}};
        }
    }
    else
    {
        for (std::size_t s = 0; s < sites.size(); ++s)
        {
            const SiteKinematics k = generalSiteKinematics(model, workspace, sites[s]);
            out[s] = {k.pose.translation, k.velocity};
        }
    }
    preparedEnergy(model, state, workspace, energy);
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
    EnergyTerms terms = generalMotionEnergy(model, workspace);
    if (model.hasForceElements())
    {
        if (auto added = detail::evaluateElementsUnchecked(model, state, workspace, false); !added)
        {
            return std::unexpected(added.error());
        }
        terms.potential += workspace.elementPotential;
    }
    return terms;
}

} // namespace detail

} // namespace sub0physics::multibody
