#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <tuple>

#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "multibody_joint.hpp"
#include "multibody_planar.hpp"
#include "multibody_ring.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/force_element.hpp"

namespace sub0physics::multibody
{

namespace
{

sub0math::Vec3d sitePosition(const Model& model, const Workspace& ws, SiteId site) noexcept
{
    const auto& offset = model.siteFrame(site).translation;
    const sub0math::Vec3d local{offset.x.si(), offset.y.si(), offset.z.si()};
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return local;
    }
    const sub0math::PluckerTransformd& x = ws.worldToBody[body];
    return x.translation + detail::bodyToWorld(x.rotation, local);
}

template<typename Lists, typename F>
void forEachElement(const Lists& lists, F&& f)
{
    std::apply([&](const auto&... list) { (std::ranges::for_each(list, f), ...); }, lists);
}

} // namespace

template<ElementLaw L>
bool CoordinateElement<L>::evaluate(const ElementInputs& in, ElementOutputs& out) const noexcept
{
    using Coordinate = LawCoordinate<L>;
    const std::size_t n = in.model.treeDofCount(); // the Jacobians' stride
    const double q = in.state.values[coordinate];
    const double qd = in.state.values[in.model.dofCount() + coordinate];
    const auto r = law.respond(Coordinate::fromSi(sign * (q - reference)), CoordinateRate<Coordinate>::fromSi(sign * qd));
    out.force[coordinate] += sign * r.force.si();
    out.potential += r.potential.si();
    out.dissipationPower += r.dissipationPower.si();
    if (!out.stiffness.empty())
    {
        if (out.wantStiffness)
        {
            out.stiffness[coordinate * n + coordinate] += r.dForceDPosition.si();
        }
        out.damping[coordinate * n + coordinate] += r.dForceDVelocity.si();
    }
    return true;
}

template<TendonLaw L>
bool TendonElement<L>::evaluate(const ElementInputs& in, ElementOutputs& out) const noexcept
{
    const Model& model = in.model;
    const Workspace& ws = in.workspace;
    const std::size_t n = model.treeDofCount(); // the Jacobians' stride and the scratch's length
    const std::size_t m = support.size();

    // Over the support: the first site's Jacobian column J1, the difference Jd = J1 - J2, u x w for a
    // revolute joint's world axis w (zero for a prismatic one), and the length Jacobian Jl = u . Jd.
    // The planar formulation stores these as (i, j, 0) in the plane's axes: every product below is a dot product
    // between them, which a change of axes does not change.
    sub0math::Vec3d* j1 = out.vectors.data();
    sub0math::Vec3d* jd = j1 + n;
    sub0math::Vec3d* uw = jd + n;
    double* jl = out.scalars.data();
    double* g = jl + n;
    double* hessian = g + n;
    const double* qd = in.state.values.data() + model.dofCount();
    double rate = 0.0;
    double l = 0.0;
    sub0math::Vec3d u{};      // in the axes the columns use
    sub0math::Vec3d uWorld{}; // in world axes, for the site wrenches
    if (in.planar)
    {
        const PlanarModel& pm = model.planar();
        const detail::PlanarPoint p1 = detail::planarSitePoint(model, ws, first.value());
        const detail::PlanarPoint p2 = detail::planarSitePoint(model, ws, second.value());
        const double dx = p1.x - p2.x;
        const double dy = p1.y - p2.y;
        const double dk = pm.siteNormal[first.value()] - pm.siteNormal[second.value()];
        l = std::sqrt(dx * dx + dy * dy + dk * dk);
        if (!(l > 0.0) || !std::isfinite(l))
        {
            return false;
        }
        const double inverse = 1.0 / l;
        u = {dx * inverse, dy * inverse, 0.0};
        uWorld = detail::planarVectorToWorld(pm, dx * inverse, dy * inverse, dk * inverse);
        for (std::size_t a = 0; a < m; ++a)
        {
            const std::size_t body = support[a];
            const PlanarBody& b = pm.bodies[body];
            const PlanarTransform& x = ws.planar.worldToBody[body];
            sub0math::Vec3d c1{};
            sub0math::Vec3d c2{};
            if (b.revolute)
            {
                // w = sigma k: w x (p - o) = sigma (-(p - o).y, (p - o).x); u x w = sigma (u.y, -u.x).
                if (onFirst[a] != 0U)
                {
                    c1 = {-b.axisSign * (p1.y - x.y), b.axisSign * (p1.x - x.x), 0.0};
                }
                if (onSecond[a] != 0U)
                {
                    c2 = {-b.axisSign * (p2.y - x.y), b.axisSign * (p2.x - x.x), 0.0};
                }
                uw[a] = {b.axisSign * u.y, -b.axisSign * u.x, 0.0};
            }
            else
            {
                const sub0math::Vec3d w{x.c * b.axisX - x.s * b.axisY, x.s * b.axisX + x.c * b.axisY, 0.0};
                c1 = onFirst[a] != 0U ? w : sub0math::Vec3d{};
                c2 = onSecond[a] != 0U ? w : sub0math::Vec3d{};
                uw[a] = {};
            }
            j1[a] = c1;
            jd[a] = c1 - c2;
            jl[a] = u.x * jd[a].x + u.y * jd[a].y;
            rate += jl[a] * qd[supportCoordinates[a]];
        }
    }
    else
    {
        const sub0math::Vec3d p1 = sitePosition(model, ws, first);
        const sub0math::Vec3d p2 = sitePosition(model, ws, second);
        const sub0math::Vec3d d = p1 - p2;
        l = std::sqrt(dot(d, d));
        if (!(l > 0.0) || !std::isfinite(l))
        {
            return false;
        }
        u = d * (1.0 / l);
        uWorld = u;
        const auto types = model.jointTypes();
        const auto axes = model.jointAxes();
        for (std::size_t a = 0; a < m; ++a)
        {
            const std::size_t body = support[a];
            const sub0math::PluckerTransformd& x = ws.worldToBody[body];
            const sub0math::Vec3d w = detail::bodyToWorld(x.rotation, axes[body]);
            const bool revolute = types[body] == JointType::Revolute;
            const sub0math::Vec3d c1 = onFirst[a] != 0U ? (revolute ? cross(w, p1 - x.translation) : w) : sub0math::Vec3d{};
            const sub0math::Vec3d c2 = onSecond[a] != 0U ? (revolute ? cross(w, p2 - x.translation) : w) : sub0math::Vec3d{};
            j1[a] = c1;
            jd[a] = c1 - c2;
            uw[a] = revolute ? cross(u, w) : sub0math::Vec3d{};
            jl[a] = dot(u, jd[a]);
            rate += jl[a] * qd[supportCoordinates[a]];
        }
    }
    const double inverseLength = 1.0 / l;

    // The law sees s = sign (l - reference); F = sign * force is the force conjugate to l.
    const auto r = law.respond(sub0math::f64::Length::fromSi(sign * (l - reference)),
                               sub0math::f64::Velocity::fromSi(sign * rate));
    const double force = sign * r.force.si();
    for (std::size_t a = 0; a < m; ++a)
    {
        out.force[supportCoordinates[a]] += jl[a] * force;
    }
    out.potential += r.potential.si();
    out.dissipationPower += r.dissipationPower.si();
    out.siteWrenches[wrenchSlot] = {id, first, uWorld * force, {}};
    out.siteWrenches[wrenchSlot + 1] = {id, second, uWorld * -force, {}};
    if (out.stiffness.empty())
    {
        return true;
    }
    const double dForceDLength = r.dForceDPosition.si();
    const double dForceDRate = r.dForceDVelocity.si();
    if (!out.geometricStiffness)
    {
        // The material terms only: J_l^T (dF/dl) J_l, negative semidefinite for a law whose force never grows
        // with its coordinate.
        for (std::size_t a = 0; a < m; ++a)
        {
            const std::size_t row = supportCoordinates[a] * n;
            for (std::size_t b = 0; b < m; ++b)
            {
                const std::size_t column = supportCoordinates[b];
                if (out.wantStiffness)
                {
                    out.stiffness[row + column] += dForceDLength * jl[a] * jl[b];
                }
                out.damping[row + column] += dForceDRate * jl[a] * jl[b];
            }
        }
        return true;
    }

    // d2l/dq_a dq_b = (Jd_a . Jd_b - Jl_a Jl_b) / l + u . (dJ1_b/dq_a - dJ2_b/dq_a), with a the ancestor (the
    // support is in topological order, so on one chain the earlier entry is the ancestor) and
    // dJ_b/dq_a = w_a x J_b for a revolute a: u . (w_a x J_b) = J_b . (u x w_a). J2 = J1 - Jd.
    // PERF: triple product via u x w; common prisms pruned (1.4x).
    for (std::size_t a = 0; a < m; ++a)
    {
        for (std::size_t b = a; b < m; ++b)
        {
            double h = (dot(jd[a], jd[b]) - jl[a] * jl[b]) * inverseLength;
            if (onFirst[a] != 0U && onFirst[b] != 0U)
            {
                h += dot(uw[a], j1[b]);
            }
            if (onSecond[a] != 0U && onSecond[b] != 0U)
            {
                h -= dot(uw[a], j1[b] - jd[b]);
            }
            hessian[a * m + b] = h;
            hessian[b * m + a] = h;
        }
    }
    // g = (d2l/dq2) qd: how the rate l' = Jl qd moves with q.
    for (std::size_t a = 0; a < m; ++a)
    {
        double sum = 0.0;
        for (std::size_t b = 0; b < m; ++b)
        {
            sum += hessian[a * m + b] * qd[supportCoordinates[b]];
        }
        g[a] = sum;
    }
    for (std::size_t a = 0; a < m; ++a)
    {
        const std::size_t row = supportCoordinates[a] * n;
        const double materialA = dForceDLength * jl[a];
        const double dampingA = dForceDRate * jl[a];
        for (std::size_t b = 0; b < m; ++b)
        {
            const std::size_t column = supportCoordinates[b];
            out.stiffness[row + column] += materialA * jl[b] + force * hessian[a * m + b] + dampingA * g[b];
            out.damping[row + column] += dampingA * jl[b];
        }
    }
    return true;
}

// The closed list of element types (force_element.hpp, CoordinateLaws and TendonLaws).
template struct CoordinateElement<LinearSpringDamper<sub0math::f64::Length>>;
template struct CoordinateElement<LinearSpringDamper<sub0math::f64::Angle>>;
template struct CoordinateElement<AxialSpringDamper>;
template struct CoordinateElement<UnilateralSpringDamper>;
template struct CoordinateElement<TravelBuffer>;
template struct TendonElement<LinearSpringDamper<sub0math::f64::Length>>;
template struct TendonElement<AxialSpringDamper>;
template struct TendonElement<UnilateralSpringDamper>;
template struct TendonElement<TravelBuffer>;

namespace detail
{

std::expected<void, Fault> evaluateElementsUnchecked(const Model& model, const State& state, Workspace& ws,
                                                     bool jacobians, bool geometricStiffness, bool planar,
                                                     bool stiffness) noexcept
{
    std::ranges::fill(ws.elementForce, 0.0);
    if (jacobians)
    {
        if (stiffness)
        {
            // PERF: skip stiffness fill and assembly when only damping is implicit.
            std::ranges::fill(ws.elementStiffness, 0.0);
        }
        std::ranges::fill(ws.elementDamping, 0.0);
    }
    ElementOutputs out{.force = ws.elementForce,
                       .stiffness = jacobians ? std::span<double>(ws.elementStiffness) : std::span<double>{},
                       .damping = jacobians ? std::span<double>(ws.elementDamping) : std::span<double>{},
                       .wantStiffness = stiffness,
                       .geometricStiffness = geometricStiffness,
                       .siteWrenches = ws.elementSiteWrenches,
                       .vectors = ws.elementVectors,
                       .scalars = ws.elementScalars};
    const ElementInputs in{model, state, ws, planar};
    bool inDomain = true;
    const auto evaluate = [&](const auto& element) { inDomain = element.evaluate(in, out) && inDomain; };
    forEachElement(model.coordinateElements(), evaluate);
    forEachElement(model.tendonElements(), evaluate);
    if (!model.rings().empty())
    {
        inDomain = detail::ringElements(model, state, ws, jacobians, geometricStiffness, out.potential,
                                        out.dissipationPower) && inDomain;
    }
    ws.elementPotential = out.potential;
    ws.elementDissipationPower = out.dissipationPower;
    if (!inDomain)
    {
        return std::unexpected(Fault{FaultKind::OutsideDomain, "tendon sites or ring nodes coincide or are not finite"});
    }
    return {};
}

} // namespace detail

std::expected<void, Fault> evaluateElements(const Model& model, const State& state, Workspace& workspace,
                                            ElementJacobians jacobians) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    const bool planar = detail::isPlanar(model);
    if (model.tendonCount() > 0)
    {
        if (planar)
        {
            detail::planarTransforms<false>(model, state, workspace);
            detail::planarWorldTransforms(model, workspace);
        }
        else
        {
            detail::updateTransforms<false>(model, state, workspace);
            detail::composeWorldTransforms(model, workspace);
        }
    }
    return detail::evaluateElementsUnchecked(model, state, workspace, jacobians == ElementJacobians::Compute, true,
                                             planar);
}

} // namespace sub0physics::multibody
