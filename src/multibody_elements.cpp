#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <tuple>

#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "multibody_joint.hpp"
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
    const std::size_t n = in.model.dofCount();
    const double q = in.state.values[coordinate];
    const double qd = in.state.values[n + coordinate];
    const auto r = law.respond(Coordinate::fromSi(sign * (q - reference)), CoordinateRate<Coordinate>::fromSi(sign * qd));
    out.force[coordinate] += sign * r.force.si();
    out.potential += r.potential.si();
    out.dissipationPower += r.dissipationPower.si();
    if (!out.stiffness.empty())
    {
        out.stiffness[coordinate * n + coordinate] += r.dForceDPosition.si();
        out.damping[coordinate * n + coordinate] += r.dForceDVelocity.si();
    }
    return true;
}

template<TendonLaw L>
bool TendonElement<L>::evaluate(const ElementInputs& in, ElementOutputs& out) const noexcept
{
    const Model& model = in.model;
    const Workspace& ws = in.workspace;
    const std::size_t n = model.dofCount();
    const std::size_t m = support.size();
    const sub0math::Vec3d p1 = sitePosition(model, ws, first);
    const sub0math::Vec3d p2 = sitePosition(model, ws, second);
    const sub0math::Vec3d d = p1 - p2;
    const double l = std::sqrt(dot(d, d));
    if (!(l > 0.0) || !std::isfinite(l))
    {
        return false;
    }
    const sub0math::Vec3d u = d * (1.0 / l);

    // Columns of the two site Jacobians, the joints' world axes and the length Jacobian, over the support.
    sub0math::Vec3d* j1 = out.vectors.data();
    sub0math::Vec3d* j2 = j1 + n;
    sub0math::Vec3d* axis = j2 + n;
    double* jl = out.scalars.data();
    double* g = jl + n;
    double* hessian = g + n;
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const auto coordinates = model.coordinates();
    const double* qd = in.state.values.data() + n;
    double rate = 0.0;
    for (std::size_t a = 0; a < m; ++a)
    {
        const std::size_t body = support[a];
        const sub0math::PluckerTransformd& x = ws.worldToBody[body];
        const sub0math::Vec3d w = detail::bodyToWorld(x.rotation, axes[body]);
        const bool revolute = types[body] == JointType::Revolute;
        axis[a] = w;
        j1[a] = onFirst[a] != 0U ? (revolute ? cross(w, p1 - x.translation) : w) : sub0math::Vec3d{};
        j2[a] = onSecond[a] != 0U ? (revolute ? cross(w, p2 - x.translation) : w) : sub0math::Vec3d{};
        jl[a] = dot(u, j1[a] - j2[a]);
        rate += jl[a] * qd[coordinates[body]];
    }

    // The law sees s = sign (l - reference); F = sign * force is the force conjugate to l.
    const auto r = law.respond(sub0math::f64::Length::fromSi(sign * (l - reference)),
                               sub0math::f64::Velocity::fromSi(sign * rate));
    const double force = sign * r.force.si();
    for (std::size_t a = 0; a < m; ++a)
    {
        out.force[coordinates[support[a]]] += jl[a] * force;
    }
    out.potential += r.potential.si();
    out.dissipationPower += r.dissipationPower.si();
    out.siteWrenches[wrenchSlot] = {id, first, u * force, {}};
    out.siteWrenches[wrenchSlot + 1] = {id, second, u * -force, {}};
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
            const std::size_t row = coordinates[support[a]] * n;
            for (std::size_t b = 0; b < m; ++b)
            {
                const std::size_t column = coordinates[support[b]];
                out.stiffness[row + column] += dForceDLength * jl[a] * jl[b];
                out.damping[row + column] += dForceDRate * jl[a] * jl[b];
            }
        }
        return true;
    }

    // d2l/dq_a dq_b = (Jd_a . Jd_b - Jl_a Jl_b) / l + u . (dJ1_b/dq_a - dJ2_b/dq_a), with a the ancestor
    // (support is in topological order, so on one chain the earlier entry is the ancestor).
    for (std::size_t a = 0; a < m; ++a)
    {
        const sub0math::Vec3d da = j1[a] - j2[a];
        const bool revolute = types[support[a]] == JointType::Revolute;
        for (std::size_t b = a; b < m; ++b)
        {
            double h = (dot(da, j1[b] - j2[b]) - jl[a] * jl[b]) / l;
            if (revolute)
            {
                if (onFirst[a] != 0U && onFirst[b] != 0U)
                {
                    h += dot(u, cross(axis[a], j1[b]));
                }
                if (onSecond[a] != 0U && onSecond[b] != 0U)
                {
                    h -= dot(u, cross(axis[a], j2[b]));
                }
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
            sum += hessian[a * m + b] * qd[coordinates[support[b]]];
        }
        g[a] = sum;
    }
    for (std::size_t a = 0; a < m; ++a)
    {
        const std::size_t row = coordinates[support[a]] * n;
        for (std::size_t b = 0; b < m; ++b)
        {
            const std::size_t column = coordinates[support[b]];
            out.stiffness[row + column] +=
                dForceDLength * jl[a] * jl[b] + force * hessian[a * m + b] + dForceDRate * jl[a] * g[b];
            out.damping[row + column] += dForceDRate * jl[a] * jl[b];
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
                                                     bool jacobians, bool geometricStiffness) noexcept
{
    std::ranges::fill(ws.elementForce, 0.0);
    if (jacobians)
    {
        std::ranges::fill(ws.elementStiffness, 0.0);
        std::ranges::fill(ws.elementDamping, 0.0);
    }
    ElementOutputs out{.force = ws.elementForce,
                       .stiffness = jacobians ? std::span<double>(ws.elementStiffness) : std::span<double>{},
                       .damping = jacobians ? std::span<double>(ws.elementDamping) : std::span<double>{},
                       .geometricStiffness = geometricStiffness,
                       .siteWrenches = ws.elementSiteWrenches,
                       .vectors = ws.elementVectors,
                       .scalars = ws.elementScalars};
    const ElementInputs in{model, state, ws};
    bool inDomain = true;
    const auto evaluate = [&](const auto& element) { inDomain = element.evaluate(in, out) && inDomain; };
    forEachElement(model.coordinateElements(), evaluate);
    forEachElement(model.tendonElements(), evaluate);
    ws.elementPotential = out.potential;
    ws.elementDissipationPower = out.dissipationPower;
    if (!inDomain)
    {
        return std::unexpected(Fault{FaultKind::OutsideDomain, "tendon sites coincide or are not finite"});
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
    if (model.tendonCount() > 0)
    {
        detail::updateTransforms<false>(model, state, workspace);
        detail::composeWorldTransforms(model, workspace);
    }
    return detail::evaluateElementsUnchecked(model, state, workspace, jacobians == ElementJacobians::Compute);
}

} // namespace sub0physics::multibody
