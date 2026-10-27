#pragma once

#include <cmath>
#include <cstddef>

#include <sub0math/config.hpp>
#include <sub0math/linear/mat.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/linear/vec.hpp>

#include "multibody_ring.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/planar.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// The planar formulation (docs/MULTIBODY.md, M1d status): the same algorithms as the general path (transforms,
// RNEA, CRBA, energy, site kinematics) in three-component planar algebra, for a model with
// Formulation::Planar. They write Workspace::planar and never the general kinematic fields; expandPlanarKinematics
// fills those for the code that reads them (constraints, reactions).

namespace sub0physics::multibody::detail
{

// Whether the model runs the planar formulation.
[[nodiscard]] SUB0MATH_FORCE_INLINE bool isPlanar(const Model& model) noexcept
{
    return model.formulation() == Formulation::Planar;
}

// Transform from the parent's frame to the body's frame at coordinate q.
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarTransform planarJoint(const PlanarBody& b, double q) noexcept
{
    const PlanarTransform& t = b.placement;
    if (b.revolute)
    {
        const double c = std::cos(q);
        const double s = b.axisSign * std::sin(q);
        if (b.unrotated)
        {
            return {c, s, t.x, t.y};
        }
        return {c * t.c - s * t.s, s * t.c + c * t.s, t.x, t.y};
    }
    // The axis is fixed in the joint frame, which the placement rotates: origin r_T + R(theta_T) a q.
    const double ax = b.axisX * q;
    const double ay = b.axisY * q;
    if (b.unrotated)
    {
        return {t.c, t.s, t.x + ax, t.y + ay};
    }
    return {t.c, t.s, t.x + t.c * ax - t.s * ay, t.y + t.s * ax + t.c * ay};
}

// S qd: the joint's motion in the body frame.
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarMotion planarSubspace(const PlanarBody& b, double rate) noexcept
{
    return b.revolute ? PlanarMotion{b.axisSign * rate, 0.0, 0.0} : PlanarMotion{0.0, b.axisX * rate, b.axisY * rate};
}

// S . f
[[nodiscard]] SUB0MATH_FORCE_INLINE double planarSubspaceDot(const PlanarBody& b, const PlanarForce& f) noexcept
{
    return b.revolute ? b.axisSign * f.n : b.axisX * f.x + b.axisY * f.y;
}

// Fills planar.parentToBody and, when requested, planar.velocity.
template<bool WithVelocity>
SUB0MATH_FORCE_INLINE void planarTransforms(const Model& model, const State& state, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const PlanarBody* bodies = model.planar().bodies.data();
    PlanarKinematics& k = ws.planar;
    const std::size_t n = parents.size();
    const double* q = state.values.data();
    const double* qd = q + model.dofCount();
    for (std::size_t i = 0; i < n; ++i)
    {
        k.parentToBody[i] = planarJoint(bodies[i], q[i]);
        if constexpr (WithVelocity)
        {
            const PlanarMotion joint = planarSubspace(bodies[i], qd[i]);
            const std::size_t p = parents[i];
            k.velocity[i] = p == NoParent ? joint : apply(k.parentToBody[i], k.velocity[p]) + joint;
        }
    }
    if (!model.rings().empty())
    {
        ringTransforms(model, state, ws, WithVelocity);
    }
}

// Composes each body's world transform (world to body) from planar.parentToBody.
SUB0MATH_FORCE_INLINE void planarWorldTransforms(const Model& model, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    PlanarKinematics& k = ws.planar;
    for (std::size_t i = 0; i < parents.size(); ++i)
    {
        const std::size_t p = parents[i];
        k.worldToBody[i] = p == NoParent ? k.parentToBody[i] : compose(k.parentToBody[i], k.worldToBody[p]);
    }
}

// Recursive Newton-Euler over transforms and velocities already in planar: tau = M qdd + c, or c without qdd.
// planar.acceleration ends as each body's acceleration (gravity included as -g), planar.force as its joint's
// transmitted force.
template<bool WithAcceleration>
SUB0MATH_FORCE_INLINE void planarNewtonEuler(const Model& model, const State& state, const double* qdd, double* tau,
                                             Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const PlanarModel& pm = model.planar();
    const PlanarBody* bodies = pm.bodies.data();
    PlanarKinematics& k = ws.planar;
    const std::size_t n = parents.size();
    const double* qd = state.values.data() + model.dofCount();
    const PlanarMotion base{0.0, -pm.gravityX, -pm.gravityY};
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::size_t p = parents[i];
        const PlanarMotion& v = k.velocity[i];
        const PlanarMotion joint = planarSubspace(bodies[i], qd[i]);
        PlanarMotion a = apply(k.parentToBody[i], p == NoParent ? base : k.acceleration[p]) + crossMotion(v, joint);
        if constexpr (WithAcceleration)
        {
            a = a + planarSubspace(bodies[i], qdd[i]);
        }
        k.acceleration[i] = a;
        const PlanarInertia& in = bodies[i].inertia;
        k.force[i] = in * a + crossForce(v, in * v);
    }
    if (!model.rings().empty())
    {
        ringNewtonEuler(model, state, WithAcceleration ? qdd : nullptr, tau, ws);
    }
    for (std::size_t i = n; i-- > 0;)
    {
        tau[i] = planarSubspaceDot(bodies[i], k.force[i]);
        if (const std::size_t p = parents[i]; p != NoParent)
        {
            k.force[p] = k.force[p] + applyTranspose(k.parentToBody[i], k.force[i]);
        }
    }
}

// Composite-rigid-body algorithm over transforms already in planar; writes massMatrix's ancestor entries.
SUB0MATH_FORCE_INLINE void planarCompositeRigidBody(const Model& model, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const PlanarBody* bodies = model.planar().bodies.data();
    PlanarKinematics& k = ws.planar;
    const std::size_t n = parents.size();
    double* h = ws.massMatrix.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        k.composite[i] = bodies[i].inertia;
    }
    if (!model.rings().empty())
    {
        ringCompositeInertia(model, ws);
    }
    for (std::size_t i = n; i-- > 0;)
    {
        if (const std::size_t p = parents[i]; p != NoParent)
        {
            k.composite[p] = k.composite[p] + applyTranspose(k.parentToBody[i], k.composite[i]);
        }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        const PlanarBody& b = bodies[i];
        const PlanarInertia& c = k.composite[i];
        // I S: {sigma I, -sigma hy, sigma hx} for a revolute joint, {hx ay - hy ax, m ax, m ay} for a prismatic one.
        PlanarForce f = b.revolute ? PlanarForce{b.axisSign * c.i, -b.axisSign * c.hy, b.axisSign * c.hx}
                                   : PlanarForce{c.hx * b.axisY - c.hy * b.axisX, c.m * b.axisX, c.m * b.axisY};
        h[i * n + i] = planarSubspaceDot(b, f);
        for (std::size_t j = i; parents[j] != NoParent;)
        {
            f = applyTranspose(k.parentToBody[j], f);
            j = parents[j];
            const double value = planarSubspaceDot(bodies[j], f);
            h[i * n + j] = value;
            h[j * n + i] = value;
        }
    }
    if (!model.rings().empty())
    {
        ringHostJacobians(model, ws);
    }
}

// A body point in world plane coordinates.
struct PlanarPoint
{
    double x;
    double y;
};

// World plane position of a point with body-frame plane coordinates (x, y) on a body with world transform t.
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarPoint planarToWorld(const PlanarTransform& t, double x, double y) noexcept
{
    return {t.x + t.c * x - t.s * y, t.y + t.s * x + t.c * y};
}

// A site's world plane position; needs planar.worldToBody.
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarPoint planarSitePoint(const Model& model, const Workspace& ws,
                                                               std::size_t site) noexcept
{
    const PlanarModel& pm = model.planar();
    const std::size_t body = model.siteBodies()[site];
    if (body == NoParent)
    {
        return {pm.siteX[site], pm.siteY[site]};
    }
    return planarToWorld(ws.planar.worldToBody[body], pm.siteX[site], pm.siteY[site]);
}

// A plane vector (along i and j) and a k component as a world vector.
[[nodiscard]] SUB0MATH_FORCE_INLINE sub0math::Vec3d planarVectorToWorld(const PlanarModel& pm, double x, double y,
                                                                       double normal) noexcept
{
    // PERF: a branch per plane beats indexed stores (store forwarding).
    switch (pm.normal)
    {
    case 0:
        return {normal, x, y}; // (i, j) = (y, z)
    case 1:
        return {y, normal, x}; // (i, j) = (z, x)
    default:
        return {x, y, normal}; // (i, j) = (x, y)
    }
}

// Fills the general kinematic fields (parentToBody, worldToBody, bodyVelocity and bodyAcceleration) from planar,
// without trigonometry, for the code that reads them (constraints, reactions). The k offsets come from the model.
inline void expandPlanarKinematics(const Model& model, Workspace& ws, bool withAcceleration) noexcept
{
    const PlanarModel& pm = model.planar();
    const std::size_t i = pm.axisI;
    const std::size_t j = pm.axisJ;
    const std::size_t kk = pm.normal;
    const auto tree = model.treeTransforms();
    const auto parents = model.parents();
    const PlanarKinematics& k = ws.planar;
    for (std::size_t b = 0; b < parents.size(); ++b)
    {
        const PlanarTransform& t = k.parentToBody[b];
        sub0math::Mat3d e{};
        e(kk, kk) = 1.0;
        e(i, i) = t.c;
        e(j, j) = t.c;
        e(i, j) = t.s;
        e(j, i) = -t.s;
        sub0math::Vec3d r{};
        r[i] = t.x;
        r[j] = t.y;
        r[kk] = tree[b].translation[kk];
        ws.parentToBody[b] = {e, r};
        const std::size_t p = parents[b];
        ws.worldToBody[b] = p == NoParent ? ws.parentToBody[b] : ws.parentToBody[b] * ws.worldToBody[p];
        const PlanarMotion& v = k.velocity[b];
        sub0math::SpatialMotiond v3{};
        v3.angular[kk] = v.w;
        v3.linear[i] = v.x;
        v3.linear[j] = v.y;
        ws.bodyVelocity[b] = v3;
        if (withAcceleration)
        {
            const PlanarMotion& a = k.acceleration[b];
            sub0math::SpatialMotiond a3{};
            a3.angular[kk] = a.w;
            a3.linear[i] = a.x;
            a3.linear[j] = a.y;
            // Gravity along k is zero for a planar model, so the k components are exact.
            ws.bodyAcceleration[b] = a3;
        }
    }
}

} // namespace sub0physics::multibody::detail
