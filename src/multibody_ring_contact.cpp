#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include "multibody_planar.hpp"
#include "multibody_ring.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/ring_contact.hpp"

// Contact on a ring's casing (docs/WHEEL.md, "M2c status"): the row of a contact point, the reference engine's soft
// constraint law and its exact solve, and the contact element's gap kinematics with exact first and second
// derivatives. Planar: every vector here is in the model's plane coordinates (i, j).

namespace sub0physics::multibody
{

namespace
{

constexpr std::size_t C = PlanarRingComponents;
constexpr std::size_t Local = RingContactLoad::Coordinates;

// A vector in the plane (i, j).
struct V2
{
    double i;
    double j;
};

constexpr V2 operator+(V2 a, V2 b) noexcept { return {a.i + b.i, a.j + b.j}; }
constexpr V2 operator-(V2 a, V2 b) noexcept { return {a.i - b.i, a.j - b.j}; }
constexpr V2 operator*(double s, V2 a) noexcept { return {s * a.i, s * a.j}; }
constexpr double dot(V2 a, V2 b) noexcept { return a.i * b.i + a.j * b.j; }
// k x a: a quarter turn about the plane's normal.
constexpr V2 perp(V2 a) noexcept { return {-a.j, a.i}; }

V2 inPlane(const PlanarModel& pm, const sub0math::Vec3d& v) noexcept
{
    return {v[pm.axisI], v[pm.axisJ]};
}

// A host-frame vector in world plane axes.
constexpr V2 rotate(const PlanarTransform& t, V2 a) noexcept
{
    return {t.c * a.i - t.s * a.j, t.s * a.i + t.c * a.j};
}

bool validRing(const Model& model, std::size_t ring, std::size_t edge) noexcept
{
    return detail::isPlanar(model) && ring < model.rings().size() && edge < model.rings()[ring].nodeCount;
}

// The world-plane velocity of point p per unit rate of body b's joint (needs planar.worldToBody).
V2 jointPointVelocity(const Model& model, const Workspace& ws, std::size_t body, V2 p) noexcept
{
    const PlanarBody& b = model.planar().bodies[body];
    const PlanarTransform& x = ws.planar.worldToBody[body];
    if (b.revolute)
    {
        return b.axisSign * perp(p - V2{x.x, x.y});
    }
    return rotate(x, V2{b.axisX, b.axisY});
}

// The local kinematics of a point at axis parameter t on edge (a, b): host-frame positions, world vectors.
struct EdgeKinematics
{
    std::size_t a;     // node e
    std::size_t b;     // node e + 1
    V2 localA;         // host frame
    V2 localB;
    V2 axis0;          // R a_0 in world
    V2 axis1;          // R a_1 in world
    PlanarTransform host;
};

EdgeKinematics edgeKinematics(const DeformableRing& ring, const RingWorkspace& rw, const Workspace& ws,
                              std::size_t edge) noexcept
{
    const std::size_t b = edge + 1U == ring.nodeCount ? 0U : edge + 1U;
    const PlanarTransform& host = ws.planar.worldToBody[ring.host];
    return {edge,
            b,
            {rw.positionX[edge], rw.positionY[edge]},
            {rw.positionX[b], rw.positionY[b]},
            rotate(host, V2{ring.axis0[0], ring.axis0[1]}),
            rotate(host, V2{ring.axis1[0], ring.axis1[1]}),
            host};
}

// d(world point)/d(local coordinate), the local coordinates (theta, v_i, v_j, w_a0, w_a1, w_b0, w_b1) of the chart
// p = o0 + R0 (v + Rot(theta) l): see ringContactLoad.
struct PointDerivatives
{
    std::array<V2, Local> first{};
    // Non-zero second derivatives: (theta, theta) and (theta, node coordinate).
    V2 thetaTheta{};
    std::array<V2, 4> thetaNode{};
};

PointDerivatives pointDerivatives(const EdgeKinematics& k, double t) noexcept
{
    const V2 l = (1.0 - t) * k.localA + t * k.localB;
    PointDerivatives d;
    d.first[0] = rotate(k.host, perp(l));
    d.first[1] = rotate(k.host, V2{1.0, 0.0});
    d.first[2] = rotate(k.host, V2{0.0, 1.0});
    d.first[3] = (1.0 - t) * k.axis0;
    d.first[4] = (1.0 - t) * k.axis1;
    d.first[5] = t * k.axis0;
    d.first[6] = t * k.axis1;
    d.thetaTheta = -1.0 * rotate(k.host, l);
    d.thetaNode = {(1.0 - t) * perp(k.axis0), (1.0 - t) * perp(k.axis1), t * perp(k.axis0), t * perp(k.axis1)};
    return d;
}

// The gap, its local gradient and Hessian (row-major 7 x 7).
struct GapDerivatives
{
    double gap = 0.0;
    std::array<double, Local> gradient{};
    std::array<double, Local * Local> hessian{};
};

GapDerivatives gapDerivatives(const EdgeKinematics& k, const RingContactGeometry& g, const PlanarModel& pm) noexcept
{
    const double t = g.axisT;
    const PointDerivatives d = pointDerivatives(k, t);
    const detail::PlanarPoint wa = detail::planarToWorld(k.host, k.localA.i, k.localA.j);
    const detail::PlanarPoint wb = detail::planarToWorld(k.host, k.localB.i, k.localB.j);
    const V2 p = (1.0 - t) * V2{wa.x, wa.y} + t * V2{wb.x, wb.y};
    const V2 w = inPlane(pm, g.witnessPoint);
    GapDerivatives out;
    V2 n = inPlane(pm, g.normal);
    double rho = 0.0;
    if (g.witness == ContactWitness::Vertex)
    {
        const V2 diff = p - w;
        rho = std::sqrt(dot(diff, diff));
        const double sign = dot(diff, n) >= 0.0 ? 1.0 : -1.0;
        if (rho > 0.0)
        {
            n = (sign / rho) * diff;
        }
        out.gap = sign * rho - g.skin;
        if (sign < 0.0)
        {
            rho = 0.0; // the axis crossed: no curvature terms (docs: deep overlap is outside the element's support)
        }
    }
    else
    {
        out.gap = dot(n, p - w) - g.skin;
    }
    for (std::size_t a = 0; a < Local; ++a)
    {
        out.gradient[a] = dot(n, d.first[a]);
    }
    // n . p'' (the host's rotation).
    out.hessian[0] = dot(n, d.thetaTheta);
    for (std::size_t c = 0; c < 4; ++c)
    {
        const double v = dot(n, d.thetaNode[c]);
        out.hessian[3 + c] = v;
        out.hessian[(3 + c) * Local] = v;
    }
    if (g.witness == ContactWitness::Vertex && rho > 0.0)
    {
        // The normal turns: P^T (I - n n^T) P / rho, with (I - n n^T) x = (n_perp . x) n_perp.
        const V2 m = perp(n);
        std::array<double, Local> tangential{};
        for (std::size_t a = 0; a < Local; ++a)
        {
            tangential[a] = dot(m, d.first[a]);
        }
        for (std::size_t a = 0; a < Local; ++a)
        {
            for (std::size_t c = 0; c < Local; ++c)
            {
                out.hessian[a * Local + c] += tangential[a] * tangential[c] / rho;
            }
        }
        // The casing's witness slides along the edge when it is interior: minus phi_qt phi_tq / phi_tt.
        if (t > 0.0 && t < 1.0)
        {
            const detail::PlanarPoint ea = wa;
            const detail::PlanarPoint eb = wb;
            const V2 e{eb.x - ea.x, eb.y - ea.y};
            const double et = dot(m, e);
            const double phiTT = et * et / rho;
            if (phiTT > 0.0)
            {
                // d e / d local: theta: R k x (l_b - l_a); node a: -R a; node b: +R a.
                std::array<V2, Local> de{};
                de[0] = rotate(k.host, perp(k.localB - k.localA));
                de[3] = -1.0 * k.axis0;
                de[4] = -1.0 * k.axis1;
                de[5] = k.axis0;
                de[6] = k.axis1;
                std::array<double, Local> phiTQ{};
                for (std::size_t a = 0; a < Local; ++a)
                {
                    phiTQ[a] = et * tangential[a] / rho + dot(n, de[a]);
                }
                for (std::size_t a = 0; a < Local; ++a)
                {
                    for (std::size_t c = 0; c < Local; ++c)
                    {
                        out.hessian[a * Local + c] -= phiTQ[a] * phiTQ[c] / phiTT;
                    }
                }
            }
        }
    }
    return out;
}

// The chart velocity: the host's motion in its frame and the two nodes' rates.
std::array<double, Local> localVelocity(const Model& model, const State& state, const Workspace& ws,
                                        const DeformableRing& ring, std::size_t a, std::size_t b) noexcept
{
    const PlanarMotion& v = ws.planar.velocity[ring.host];
    const double* wd = state.values.data() + model.dofCount() + ring.firstCoordinate;
    return {v.w, v.x, v.y, wd[C * a], wd[C * a + 1U], wd[C * b], wd[C * b + 1U]};
}

// The host origin's world-plane velocity per unit rate of each chain joint, and the chain's curvature:
// d2 o / dq_j dq_k = [j revolute] w_j k x (d o / dq_k) for j an ancestor of k or k itself (the chain is host first,
// so a larger chain index is an ancestor). Host-frame components (R0^T).
void chainCurvature(const Model& model, const Workspace& ws, const DeformableRing& ring,
                    std::span<double> curvatureI, std::span<double> curvatureJ) noexcept
{
    const std::span<const std::size_t> chain = ring.hostChain;
    const std::size_t m = chain.size();
    const PlanarTransform& host = ws.planar.worldToBody[ring.host];
    const V2 o{host.x, host.y};
    for (std::size_t a = 0; a < m; ++a)
    {
        for (std::size_t c = 0; c < m; ++c)
        {
            // The ancestor of the pair is the one with the larger chain index.
            const std::size_t ancestor = std::max(a, c);
            const std::size_t other = std::min(a, c);
            const PlanarBody& body = model.planar().bodies[chain[ancestor]];
            V2 second{0.0, 0.0};
            if (body.revolute)
            {
                second = body.axisSign * perp(jointPointVelocity(model, ws, chain[other], o));
            }
            // In host-frame components.
            curvatureI[a * m + c] = host.c * second.i + host.s * second.j;
            curvatureJ[a * m + c] = -host.s * second.i + host.c * second.j;
        }
    }
}

} // namespace

std::expected<void, Fault> ringContactJacobian(const Model& model, const Workspace& workspace,
                                               const RingContactPoint& contact, std::span<double> row) noexcept
{
    if (!workspace.isBuiltFor(model) || !validRing(model, contact.ring, contact.edge) || row.size() != model.dofCount())
    {
        return invalidInput("ring contact Jacobian needs a planar model, a ring and edge in it, and dofCount entries");
    }
    std::ranges::fill(row, 0.0);
    const PlanarModel& pm = model.planar();
    const DeformableRing& ring = model.rings()[contact.ring];
    const V2 p = inPlane(pm, contact.point);
    const V2 n = inPlane(pm, contact.normal);
    for (const std::size_t body : ring.hostChain)
    {
        row[model.coordinates()[body]] = dot(n, jointPointVelocity(model, workspace, body, p));
    }
    const PlanarTransform& host = workspace.planar.worldToBody[ring.host];
    const double along0 = dot(n, rotate(host, V2{ring.axis0[0], ring.axis0[1]}));
    const double along1 = dot(n, rotate(host, V2{ring.axis1[0], ring.axis1[1]}));
    const std::size_t a = contact.edge;
    const std::size_t b = a + 1U == ring.nodeCount ? 0U : a + 1U;
    row[ring.firstCoordinate + C * a] += contact.weight0 * along0;
    row[ring.firstCoordinate + C * a + 1U] += contact.weight0 * along1;
    row[ring.firstCoordinate + C * b] += contact.weight1 * along0;
    row[ring.firstCoordinate + C * b + 1U] += contact.weight1 * along1;
    return {};
}

std::expected<void, Fault> ringNodeInverseWeights(const Model& model, const Workspace& workspace, std::size_t ring,
                                                  std::span<double> out)
{
    if (!workspace.isBuiltFor(model) || !validRing(model, ring, 0) || out.size() != model.rings()[ring].nodeCount)
    {
        return invalidInput("ring node inverse weights need a planar model, one of its rings and one entry per node");
    }
    const DeformableRing& r = model.rings()[ring];
    const RingWorkspace& rw = workspace.rings[ring];
    const PlanarTransform& host = workspace.planar.worldToBody[r.host];
    const V2 a0 = rotate(host, V2{r.axis0[0], r.axis0[1]});
    const V2 a1 = rotate(host, V2{r.axis1[0], r.axis1[1]});
    std::vector<double> row(model.dofCount(), 0.0);
    std::vector<double> x(model.dofCount(), 0.0);
    for (std::size_t node = 0; node < r.nodeCount; ++node)
    {
        const detail::PlanarPoint w = detail::planarToWorld(host, rw.positionX[node], rw.positionY[node]);
        const V2 p{w.x, w.y};
        // The mean of the diagonal over three world axes; the axis normal to the plane has no motion.
        double sum = 0.0;
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
            std::ranges::fill(row, 0.0);
            for (const std::size_t body : r.hostChain)
            {
                const V2 v = jointPointVelocity(model, workspace, body, p);
                row[model.coordinates()[body]] = axis == 0 ? v.i : v.j;
            }
            row[r.firstCoordinate + C * node] = axis == 0 ? a0.i : a0.j;
            row[r.firstCoordinate + C * node + 1U] = axis == 0 ? a1.i : a1.j;
            x = row;
            if (auto solved = solveMassMatrix(model, workspace, x); !solved)
            {
                return std::unexpected(solved.error());
            }
            for (std::size_t k = 0; k < row.size(); ++k)
            {
                sum += row[k] * x[k];
            }
        }
        out[node] = sum / 3.0;
    }
    return {};
}

std::expected<void, Fault> measureRingContact(const Model& model, const State& state, const Workspace& workspace,
                                              const RingContactGeometry& geometry, double& gap,
                                              double& gapRate) noexcept
{
    if (!workspace.isBuiltFor(model) || !state.isSizedFor(model) || !validRing(model, geometry.ring, geometry.edge))
    {
        return invalidInput("ring contact needs a planar model, a sized state, and a ring and edge of the model");
    }
    const DeformableRing& ring = model.rings()[geometry.ring];
    const EdgeKinematics k = edgeKinematics(ring, workspace.rings[geometry.ring], workspace, geometry.edge);
    const GapDerivatives d = gapDerivatives(k, geometry, model.planar());
    const auto v = localVelocity(model, state, workspace, ring, k.a, k.b);
    gap = d.gap;
    gapRate = 0.0;
    for (std::size_t a = 0; a < Local; ++a)
    {
        gapRate += d.gradient[a] * v[a];
    }
    return {};
}

std::expected<void, Fault> ringContactLoad(const Model& model, const State& state, const Workspace& workspace,
                                           const RingContactGeometry& geometry, const ContactLawValue& law,
                                           RingContactLoad& load) noexcept
{
    if (!workspace.isBuiltFor(model) || !state.isSizedFor(model) || !validRing(model, geometry.ring, geometry.edge))
    {
        return invalidInput("ring contact needs a planar model, a sized state, and a ring and edge of the model");
    }
    const DeformableRing& ring = model.rings()[geometry.ring];
    const EdgeKinematics k = edgeKinematics(ring, workspace.rings[geometry.ring], workspace, geometry.edge);
    const GapDerivatives d = gapDerivatives(k, geometry, model.planar());
    const auto v = localVelocity(model, state, workspace, ring, k.a, k.b);
    // H v: the rate transport of the gap rate.
    std::array<double, Local> hv{};
    for (std::size_t a = 0; a < Local; ++a)
    {
        double sum = 0.0;
        for (std::size_t c = 0; c < Local; ++c)
        {
            sum += d.hessian[a * Local + c] * v[c];
        }
        hv[a] = sum;
    }
    load.ring = geometry.ring;
    load.edge = geometry.edge;
    load.gradient = d.gradient;
    load.normalForce = law.force;
    load.dForceDGapRate = law.dForceDGapRate;
    for (std::size_t a = 0; a < Local; ++a)
    {
        load.force[a] = law.force * d.gradient[a];
        for (std::size_t c = 0; c < Local; ++c)
        {
            load.stiffness[a * Local + c] = law.dForceDGap * d.gradient[a] * d.gradient[c] +
                                            law.force * d.hessian[a * Local + c] +
                                            law.dForceDGapRate * d.gradient[a] * hv[c];
            load.damping[a * Local + c] = law.dForceDGapRate * d.gradient[a] * d.gradient[c];
        }
    }
    return {};
}

namespace detail
{

void ringLoadBlocks(const Model& model, const State& state, const Workspace& ws, const RingContactLoad& load,
                    bool structured, RingLoadBlocks& out) noexcept
{
    constexpr std::size_t Cap = RingLoadBlocks::Capacity;
    const std::size_t n = model.dofCount();
    const double* qd = state.values.data() + n;
    const DeformableRing& ring = model.rings()[load.ring];
    const RingWorkspace& rw = ws.rings[load.ring];
    const std::span<const std::size_t> chain = ring.hostChain;
    const std::size_t m = chain.size();
    const std::size_t a = load.edge;
    const std::size_t b = a + 1U == ring.nodeCount ? 0U : a + 1U;
    out.chain = m;
    out.count = m + 4U;
    const std::size_t count = out.count;
    for (std::size_t i = 0; i < m; ++i)
    {
        out.index[i] = model.coordinates()[chain[i]];
    }
    out.index[m] = ring.firstCoordinate + C * a;
    out.index[m + 1U] = ring.firstCoordinate + C * a + 1U;
    out.index[m + 2U] = ring.firstCoordinate + C * b;
    out.index[m + 3U] = ring.firstCoordinate + C * b + 1U;
    // L: local (7) -> touched generalized (count): chain columns through the host's motion, nodes one to one.
    std::array<double, Local * Cap> map{};
    for (std::size_t i = 0; i < m; ++i)
    {
        const PlanarMotion& sm = rw.hostJacobian[i];
        map[0 * Cap + i] = sm.w;
        map[1 * Cap + i] = sm.x;
        map[2 * Cap + i] = sm.y;
    }
    for (std::size_t r = 0; r < 4; ++r)
    {
        map[(3 + r) * Cap + m + r] = 1.0;
    }
    const auto lift = [&](const std::array<double, Local * Local>& local, std::array<double, RingLoadBlocks::Size>& g) {
        g.fill(0.0);
        std::array<double, Local * Cap> right{}; // local * L
        for (std::size_t r = 0; r < Local; ++r)
        {
            for (std::size_t j = 0; j < count; ++j)
            {
                double sum = 0.0;
                for (std::size_t s = 0; s < Local; ++s)
                {
                    sum += local[r * Local + s] * map[s * Cap + j];
                }
                right[r * Cap + j] = sum;
            }
        }
        for (std::size_t i = 0; i < count; ++i)
        {
            for (std::size_t j = 0; j < count; ++j)
            {
                double sum = 0.0;
                for (std::size_t r = 0; r < Local; ++r)
                {
                    sum += map[r * Cap + i] * right[r * Cap + j];
                }
                g[i * Cap + j] = sum;
            }
        }
    };
    for (std::size_t i = 0; i < count; ++i)
    {
        double sum = 0.0;
        for (std::size_t r = 0; r < Local; ++r)
        {
            sum += map[r * Cap + i] * load.force[r];
        }
        out.force[i] = sum;
    }
    lift(load.damping, out.damping);
    lift(load.stiffness, out.stiffness);
    // The host chain's curvature: F C(grad_v g) on the chain block, and the rate transport through it,
    // f_gdot (L^T grad)(C(grad_v) qd)^T (in the structured form only on the chain's rows, where it factors).
    std::array<double, RingLoadBlocks::MaxChain * RingLoadBlocks::MaxChain> ci{};
    std::array<double, RingLoadBlocks::MaxChain * RingLoadBlocks::MaxChain> cj{};
    chainCurvature(model, ws, ring, std::span<double>(ci.data(), m * m), std::span<double>(cj.data(), m * m));
    std::array<double, RingLoadBlocks::MaxChain> cqd{};
    for (std::size_t i = 0; i < m; ++i)
    {
        double sum = 0.0;
        for (std::size_t j = 0; j < m; ++j)
        {
            const double curvature = load.gradient[1] * ci[i * m + j] + load.gradient[2] * cj[i * m + j];
            out.stiffness[i * Cap + j] += load.normalForce * curvature;
            sum += curvature * qd[out.index[j]];
        }
        cqd[i] = sum;
    }
    const std::size_t rows = structured ? m : count;
    for (std::size_t i = 0; i < rows; ++i)
    {
        double grad = 0.0;
        for (std::size_t r = 0; r < Local; ++r)
        {
            grad += map[r * Cap + i] * load.gradient[r];
        }
        for (std::size_t j = 0; j < m; ++j)
        {
            out.stiffness[i * Cap + j] += load.dForceDGapRate * grad * cqd[j];
        }
    }
}

void ringLoadStepTerms(const Model& model, const State& state, Workspace& ws, std::span<const RingContactLoad> loads,
                       double hd, double h2, double dt, double thetaStep, double* tree, double* b) noexcept
{
    constexpr std::size_t Cap = RingLoadBlocks::Capacity;
    const std::size_t nt = model.treeDofCount();
    const double* qd = state.values.data() + model.dofCount();
    if (!loads.empty())
    {
        planarWorldTransforms(model, ws); // the chain's curvature reads them; a one-call step composes them only for tendons
    }
    for (RingWorkspace& rw : ws.rings)
    {
        if (rw.contactBlocks)
        {
            std::ranges::fill(rw.contactDiagonal, Block<C>{});
            std::ranges::fill(rw.contactUpper, Block<C>{});
            std::ranges::fill(rw.contactBorder, BlockColumns<C, 3>{});
            rw.contactBlocks = false;
        }
    }
    for (const RingContactLoad& load : loads)
    {
        RingLoadBlocks g;
        ringLoadBlocks(model, state, ws, load, true, g);
        const std::size_t m = g.chain;
        const std::size_t count = g.count;
        // A = -hd D - h2 K~ over the touched coordinates.
        std::array<double, RingLoadBlocks::Size> step{};
        for (std::size_t i = 0; i < count; ++i)
        {
            for (std::size_t j = 0; j < count; ++j)
            {
                const double k = 0.5 * (g.stiffness[i * Cap + j] + g.stiffness[j * Cap + i]);
                const double d = 0.5 * (g.damping[i * Cap + j] + g.damping[j * Cap + i]);
                step[i * Cap + j] = -hd * d - h2 * k;
            }
        }
        // The host chain's block (ancestor pairs, so the tree's sparsity is kept), lower triangle.
        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < m; ++j)
            {
                if (g.index[i] >= g.index[j])
                {
                    tree[g.index[i] * nt + g.index[j]] += step[i * Cap + j];
                }
            }
        }
        if (h2 != 0.0)
        {
            for (std::size_t i = 0; i < count; ++i)
            {
                double sum = 0.0;
                for (std::size_t j = 0; j < count; ++j)
                {
                    sum += 0.5 * (g.stiffness[i * Cap + j] + g.stiffness[j * Cap + i]) * qd[g.index[j]];
                }
                b[g.index[i]] += dt * thetaStep * sum;
            }
        }
        // The nodes' blocks and the border (in the host's motion, from the local matrices).
        const DeformableRing& ring = model.rings()[load.ring];
        RingWorkspace& rw = ws.rings[load.ring];
        rw.contactBlocks = true;
        const std::size_t e = load.edge;
        const std::size_t f = e + 1U == ring.nodeCount ? 0U : e + 1U;
        for (std::size_t r = 0; r < C; ++r)
        {
            for (std::size_t c = 0; c < C; ++c)
            {
                rw.contactDiagonal[e][r * C + c] += step[(m + r) * Cap + m + c];
                rw.contactDiagonal[f][r * C + c] += step[(m + 2 + r) * Cap + m + 2 + c];
                rw.contactUpper[e][r * C + c] += step[(m + r) * Cap + m + 2 + c];
            }
        }
        for (std::size_t node = 0; node < 2; ++node)
        {
            BlockColumns<C, 3>& border = node == 0 ? rw.contactBorder[e] : rw.contactBorder[f];
            for (std::size_t r = 0; r < C; ++r)
            {
                const std::size_t w = 3 + 2 * node + r;
                for (std::size_t s = 0; s < 3; ++s)
                {
                    const double k = 0.5 * (load.stiffness[w * Local + s] + load.stiffness[s * Local + w]);
                    const double d = 0.5 * (load.damping[w * Local + s] + load.damping[s * Local + w]);
                    border[r * 3 + s] += -hd * d - h2 * k;
                }
            }
        }
    }
}

} // namespace detail

std::expected<void, Fault> ringLoadsToGeneralized(const Model& model, const State& state, const Workspace& workspace,
                                                  std::span<const RingContactLoad> loads, std::span<double> tau,
                                                  std::span<double> stiffness, std::span<double> damping) noexcept
{
    const std::size_t n = model.dofCount();
    if (!workspace.isBuiltFor(model) || !state.isSizedFor(model) || tau.size() != n ||
        (!stiffness.empty() && stiffness.size() != n * n) || (!damping.empty() && damping.size() != n * n))
    {
        return invalidInput("ring loads map needs a sized state, tau of dofCount and dense dofCount^2 Jacobians");
    }
    for (const RingContactLoad& load : loads)
    {
        if (!validRing(model, load.ring, load.edge) ||
            model.rings()[load.ring].hostChain.size() > detail::RingLoadBlocks::MaxChain)
        {
            return invalidInput("ring load names a ring or edge that is not in the model, or a host chain too long");
        }
        detail::RingLoadBlocks blocks;
        detail::ringLoadBlocks(model, state, workspace, load, false, blocks);
        for (std::size_t i = 0; i < blocks.count; ++i)
        {
            tau[blocks.index[i]] += blocks.force[i];
            for (std::size_t j = 0; j < blocks.count; ++j)
            {
                const std::size_t at = blocks.index[i] * n + blocks.index[j];
                if (!stiffness.empty())
                {
                    stiffness[at] += blocks.stiffness[i * detail::RingLoadBlocks::Capacity + j];
                }
                if (!damping.empty())
                {
                    damping[at] += blocks.damping[i * detail::RingLoadBlocks::Capacity + j];
                }
            }
        }
    }
    return {};
}

} // namespace sub0physics::multibody
