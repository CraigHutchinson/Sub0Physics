#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <numbers>
#include <span>

#include <sub0math/units/quantities.hpp>

#include "multibody_planar.hpp"
#include "multibody_ring.hpp"
#include "sub0physics/multibody/deformable_ring.hpp"
#include "sub0physics/solver/periodic_block_tridiagonal.hpp"

namespace sub0physics::multibody
{

namespace
{

constexpr std::size_t C = PlanarRingComponents;
using Block2 = Block<C>;

// sin(i pi / 64), i = 0..32: canonical binary64 values, reflected for the rest of the circle (docs/WHEEL.md section 1,
// "Rest geometry"). Ported from the reference's moving-wheel model (upstream/bikelab.json).
constexpr std::array<double, 33> QuarterSine{
    0.0,
    0.049067674327418015,
    0.0980171403295606,
    0.14673047445536175,
    0.19509032201612825,
    0.24298017990326387,
    0.29028467725446233,
    0.33688985339222005,
    0.3826834323650898,
    0.4275550934302821,
    0.47139673682599764,
    0.5141027441932217,
    0.5555702330196022,
    0.5956993044924334,
    0.6343932841636455,
    0.6715589548470183,
    0.7071067811865475,
    0.7409511253549591,
    0.7730104533627369,
    0.8032075314806448,
    0.8314696123025452,
    0.8577286100002721,
    0.8819212643483549,
    0.9039892931234433,
    0.9238795325112867,
    0.9415440651830208,
    0.9569403357322089,
    0.970031253194544,
    0.9807852804032304,
    0.989176509964781,
    0.9951847266721968,
    0.9987954562051724,
    1.0};

// 2 sin(pi / 128): the one supported chord ratio that is not twice a table entry.
constexpr double ChordRatio128 = 0.049082457045824576;

// sin(i pi / 64) for any integer i, by reflection of the quarter table.
constexpr double ringSine(std::size_t i) noexcept
{
    i %= 128U;
    if (i <= 32U)
    {
        return QuarterSine[i];
    }
    if (i <= 64U)
    {
        return QuarterSine[64U - i];
    }
    if (i <= 96U)
    {
        return -QuarterSine[i - 64U];
    }
    return -QuarterSine[128U - i];
}

constexpr bool supportedNodeCount(std::size_t nodes) noexcept
{
    return nodes >= 4U && nodes <= 128U && 128U % nodes == 0U;
}

// The symmetric part of a block.
constexpr Block2 symmetric(const Block2& b) noexcept
{
    const double off = 0.5 * (b[1] + b[2]);
    return {b[0], off, off, b[3]};
}

// The plane coordinates (i, j) of a node's displacement w along the ring's axes.
struct PlaneVector
{
    double x;
    double y;
};

PlaneVector alongAxes(const DeformableRing& ring, double w0, double w1) noexcept
{
    return {ring.axis0[0] * w0 + ring.axis1[0] * w1, ring.axis0[1] * w0 + ring.axis1[1] * w1};
}

} // namespace

std::expected<CanonicalRingDirection, Fault> canonicalRingDirection(std::size_t nodes, std::size_t index) noexcept
{
    if (!supportedNodeCount(nodes))
    {
        return invalidInput("canonical ring geometry needs a node count that divides 128, at least 4");
    }
    const std::size_t a = (index % nodes) * (128U / nodes);
    return CanonicalRingDirection{ringSine(a), ringSine(a + 32U)};
}

std::expected<double, Fault> canonicalChordRatio(std::size_t nodes) noexcept
{
    if (!supportedNodeCount(nodes))
    {
        return invalidInput("canonical ring geometry needs a node count that divides 128, at least 4");
    }
    // 2 sin(pi / N) = 2 sin((64 / N) pi / 64); doubling is exact.
    return nodes == 128U ? ChordRatio128 : 2.0 * QuarterSine[64U / nodes];
}

std::expected<DeformableRingDescription, Fault> uniformRing(const UniformRingSpec& spec)
{
    using sub0math::f64::Length;
    const std::size_t n = spec.nodes;
    const auto chord = canonicalChordRatio(n);
    if (!chord)
    {
        return std::unexpected(chord.error());
    }
    const double r = spec.centrelineRadius.si();
    const double skin = spec.skinRadius.si();
    const double mass = spec.ringMass.si();
    const double k = spec.foundationStiffness.si();
    const double ea = spec.extensionStiffness.si();
    const double c = spec.edgeViscosity.si();
    const bool finite = std::isfinite(r) && std::isfinite(skin) && std::isfinite(mass) && std::isfinite(k) &&
                        std::isfinite(ea) && std::isfinite(c);
    if (!finite || !(r > 0.0) || !(skin >= 0.0) || !(mass > 0.0) || !(k >= 0.0) || !(ea >= 0.0) || !(c >= 0.0))
    {
        return invalidInput("uniform ring needs finite values: positive radius and mass, the rest not negative");
    }
    // The reference's arithmetic, in its order: ds = 2 pi R / N, L = R chord, m = M / N, I = m skin^2 / 2.
    const double ds = 2.0 * std::numbers::pi * r / static_cast<double>(n);
    const double chordLength = r * *chord;
    const double nodeMass = mass / static_cast<double>(n);
    const double nodeInertia = nodeMass * (skin * skin) / 2.0;
    const auto foundation = LinearSpringDamper<Length>::create({.stiffness = sub0math::f64::Stiffness::fromSi(k * ds),
                                                                .damping = sub0math::f64::Damping::fromSi(0.0)});
    const auto edge = LinearSpringDamper<Length>::create({.stiffness = sub0math::f64::Stiffness::fromSi(ea / chordLength),
                                                          .damping = sub0math::f64::Damping::fromSi(c / chordLength)});
    if (!foundation || !edge)
    {
        return invalidInput("uniform ring coefficients are not finite");
    }
    DeformableRingDescription d;
    d.name = spec.name;
    d.host = spec.host;
    d.displacementAxes = spec.displacementAxes;
    d.restPositions.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const CanonicalRingDirection t = *canonicalRingDirection(n, i);
        const double along0 = r * t.sine;
        const double along1 = -r * t.cosine;
        const auto& a0 = spec.displacementAxes[0];
        const auto& a1 = spec.displacementAxes[1];
        d.restPositions.push_back({Length::fromSi(a0.x * along0 + a1.x * along1),
                                   Length::fromSi(a0.y * along0 + a1.y * along1),
                                   Length::fromSi(a0.z * along0 + a1.z * along1)});
    }
    d.nodeMass.assign(n, sub0math::f64::Mass::fromSi(nodeMass));
    d.nodeInertia.assign(n, sub0math::f64::MomentOfInertia::fromSi(nodeInertia));
    d.foundation.assign(n, *foundation);
    d.edgeLaws.assign(n, *edge);
    d.edgeRestLength.reserve(n);
    for (std::size_t e = 0; e < n; ++e)
    {
        const auto& p = d.restPositions[e];
        const auto& q = d.restPositions[(e + 1U) % n];
        const double dx = q.x.si() - p.x.si();
        const double dy = q.y.si() - p.y.si();
        const double dz = q.z.si() - p.z.si();
        d.edgeRestLength.push_back(Length::fromSi(std::sqrt(dx * dx + dy * dy + dz * dz)));
    }
    d.skinRadius = spec.skinRadius;
    return d;
}

namespace detail
{

void ringTransforms(const Model& model, const State& state, Workspace& ws, bool withVelocity) noexcept
{
    const auto rings = model.rings();
    const std::size_t n = model.dofCount();
    const double* q = state.values.data();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        const double* w = q + ring.firstCoordinate;
        const double* wd = w + n;
        const PlanarMotion host = withVelocity ? ws.planar.velocity[ring.host] : PlanarMotion{};
        for (std::size_t i = 0; i < ring.nodeCount; ++i)
        {
            const PlaneVector d = alongAxes(ring, w[C * i], w[C * i + 1U]);
            const double x = ring.restX[i] + d.x;
            const double y = ring.restY[i] + d.y;
            rw.positionX[i] = x;
            rw.positionY[i] = y;
            if (withVelocity)
            {
                const PlaneVector v = alongAxes(ring, wd[C * i], wd[C * i + 1U]);
                rw.velocity[i] = apply(PlanarTransform{1.0, 0.0, x, y}, host) + PlanarMotion{0.0, v.x, v.y};
            }
        }
    }
}

void ringNewtonEuler(const Model& model, const State& state, const double* qdd, double* tau, Workspace& ws) noexcept
{
    const auto rings = model.rings();
    const std::size_t n = model.dofCount();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        const double* wd = state.values.data() + n + ring.firstCoordinate;
        const PlanarMotion& hostAcceleration = ws.planar.acceleration[ring.host];
        PlanarForce total{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < ring.nodeCount; ++i)
        {
            const double x = rw.positionX[i];
            const double y = rw.positionY[i];
            const PlanarTransform toNode{1.0, 0.0, x, y};
            const PlaneVector rate = alongAxes(ring, wd[C * i], wd[C * i + 1U]);
            const PlanarMotion joint{0.0, rate.x, rate.y};
            const PlanarMotion& v = rw.velocity[i];
            PlanarMotion a = apply(toNode, hostAcceleration) + crossMotion(v, joint);
            if (qdd != nullptr)
            {
                const double* wdd = qdd + ring.firstCoordinate;
                const PlaneVector accel = alongAxes(ring, wdd[C * i], wdd[C * i + 1U]);
                a = a + PlanarMotion{0.0, accel.x, accel.y};
            }
            const PlanarInertia in{ring.mass[i], 0.0, 0.0, ring.inertia[i]};
            const PlanarForce f = in * a + crossForce(v, in * v);
            double* t = tau + ring.firstCoordinate + C * i;
            t[0] = ring.axis0[0] * f.x + ring.axis0[1] * f.y;
            t[1] = ring.axis1[0] * f.x + ring.axis1[1] * f.y;
            total = total + applyTranspose(toNode, f);
        }
        rw.hostForce = total;
        ws.planar.force[ring.host] = ws.planar.force[ring.host] + total;
    }
}

void ringCompositeInertia(const Model& model, Workspace& ws) noexcept
{
    const auto rings = model.rings();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        PlanarInertia total{0.0, 0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < ring.nodeCount; ++i)
        {
            const double m = ring.mass[i];
            const double x = rw.positionX[i];
            const double y = rw.positionY[i];
            total = total + applyTranspose(PlanarTransform{1.0, 0.0, x, y}, PlanarInertia{m, 0.0, 0.0, ring.inertia[i]});
            // The node's momentum per unit host motion, projected on its axes: m a_k . (-y w + v_x, x w + v_y).
            BlockColumns<C, 3>& b = rw.border[i];
            b = {m * (ring.axis0[1] * x - ring.axis0[0] * y), m * ring.axis0[0], m * ring.axis0[1],
                 m * (ring.axis1[1] * x - ring.axis1[0] * y), m * ring.axis1[0], m * ring.axis1[1]};
        }
        ws.planar.composite[ring.host] = ws.planar.composite[ring.host] + total;
    }
}

void ringHostJacobians(const Model& model, Workspace& ws) noexcept
{
    const auto rings = model.rings();
    const PlanarBody* bodies = model.planar().bodies.data();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        const std::span<const std::size_t> chain = ring.hostChain;
        for (std::size_t c = 0; c < chain.size(); ++c)
        {
            // The joint's motion in its own frame, carried down the chain to the host.
            PlanarMotion m = planarSubspace(bodies[chain[c]], 1.0);
            for (std::size_t d = c; d-- > 0;)
            {
                m = apply(ws.planar.parentToBody[chain[d]], m);
            }
            rw.hostJacobian[c] = m;
        }
    }
}

bool ringElements(const Model& model, const State& state, Workspace& ws, bool jacobians, bool geometricStiffness,
                  double& potential, double& dissipationPower) noexcept
{
    using sub0math::f64::Length;
    using sub0math::f64::Velocity;
    const auto rings = model.rings();
    const std::size_t n = model.dofCount();
    bool inDomain = true;
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        const std::size_t nodes = ring.nodeCount;
        const double* w = state.values.data() + ring.firstCoordinate;
        const double* wd = w + n;
        double* force = ws.elementForce.data() + ring.firstCoordinate;
        // Foundation: the node's law on each of its coordinates.
        for (std::size_t i = 0; i < nodes; ++i)
        {
            double stiffness = 0.0;
            double damping = 0.0;
            for (std::size_t k = 0; k < C; ++k)
            {
                const std::size_t j = C * i + k;
                const auto res = ring.foundation[i].respond(Length::fromSi(w[j]), Velocity::fromSi(wd[j]));
                force[j] = res.force.si();
                potential += res.potential.si();
                dissipationPower += res.dissipationPower.si();
                stiffness = res.dForceDPosition.si();
                damping = res.dForceDVelocity.si();
            }
            if (jacobians)
            {
                rw.stiffnessDiagonal[i] = {stiffness, 0.0, 0.0, stiffness};
                rw.dampingDiagonal[i] = {damping, 0.0, 0.0, damping};
            }
        }
        // Edges: the law on the distance between nodes e and e + 1, in the displacement axes.
        for (std::size_t e = 0; e < nodes; ++e)
        {
            const std::size_t f = e + 1U == nodes ? 0U : e + 1U;
            const double dx = (ring.rest0[f] + w[C * f]) - (ring.rest0[e] + w[C * e]);
            const double dy = (ring.rest1[f] + w[C * f + 1U]) - (ring.rest1[e] + w[C * e + 1U]);
            const double vx = wd[C * f] - wd[C * e];
            const double vy = wd[C * f + 1U] - wd[C * e + 1U];
            const double l = std::sqrt(dx * dx + dy * dy);
            if (!(l > 0.0) || !std::isfinite(l))
            {
                inDomain = false;
                continue;
            }
            const double inverse = 1.0 / l;
            const double ux = dx * inverse;
            const double uy = dy * inverse;
            const double rate = ux * vx + uy * vy;
            const auto res = ring.edgeLaws[e].respond(Length::fromSi(l - ring.edgeRestLength[e]), Velocity::fromSi(rate));
            const double tension = res.force.si(); // conjugate to l: on node f along u, on node e against it
            force[C * f] += tension * ux;
            force[C * f + 1U] += tension * uy;
            force[C * e] -= tension * ux;
            force[C * e + 1U] -= tension * uy;
            potential += res.potential.si();
            dissipationPower += res.dissipationPower.si();
            if (!jacobians)
            {
                continue;
            }
            const double kl = res.dForceDPosition.si();
            const double cl = res.dForceDVelocity.si();
            // G = dtau_f/dw_f: material k u u^T, geometric F (I - u u^T) / l, rate transport c u ((I - u u^T) v / l)^T.
            Block2 g{kl * ux * ux, kl * ux * uy, kl * uy * ux, kl * uy * uy};
            if (geometricStiffness)
            {
                const double p = tension * inverse;
                const double tx = (vx - rate * ux) * inverse; // (I - u u^T) v / l
                const double ty = (vy - rate * uy) * inverse;
                g[0] += p * (1.0 - ux * ux) + cl * ux * tx;
                g[1] += -p * ux * uy + cl * ux * ty;
                g[2] += -p * uy * ux + cl * uy * tx;
                g[3] += p * (1.0 - uy * uy) + cl * uy * ty;
            }
            const Block2 d{cl * ux * ux, cl * ux * uy, cl * uy * ux, cl * uy * uy};
            for (std::size_t k = 0; k < C * C; ++k)
            {
                rw.stiffnessDiagonal[e][k] += g[k];
                rw.stiffnessDiagonal[f][k] += g[k];
                rw.stiffnessUpper[e][k] = -g[k];
                rw.dampingDiagonal[e][k] += d[k];
                rw.dampingDiagonal[f][k] += d[k];
                rw.dampingUpper[e][k] = -d[k];
            }
        }
    }
    return inDomain;
}

void ringMotionEnergy(const Model& model, const Workspace& ws, double& kinetic, double& potential) noexcept
{
    const auto rings = model.rings();
    const PlanarModel& pm = model.planar();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        const RingWorkspace& rw = ws.rings[r];
        const PlanarTransform& host = ws.planar.worldToBody[ring.host];
        for (std::size_t i = 0; i < ring.nodeCount; ++i)
        {
            const PlanarInertia in{ring.mass[i], 0.0, 0.0, ring.inertia[i]};
            const PlanarMotion& v = rw.velocity[i];
            kinetic += 0.5 * dot(v, in * v);
            const PlanarPoint p = planarToWorld(host, rw.positionX[i], rw.positionY[i]);
            potential -= pm.gravityX * (ring.mass[i] * p.x) + pm.gravityY * (ring.mass[i] * p.y);
        }
    }
}

bool condenseRings(const Model& model, Workspace& ws, RingSystem system, double hd, double h2, double* tree) noexcept
{
    const auto rings = model.rings();
    const std::size_t nt = model.treeDofCount();
    const bool withElements = hd != 0.0 || h2 != 0.0;
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        RingWorkspace& rw = ws.rings[r];
        RingCondensation& cond = system == RingSystem::Mass ? rw.mass : rw.step;
        const std::size_t nodes = ring.nodeCount;
        const bool contact = system == RingSystem::Step && rw.contactBlocks;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            const double m = ring.mass[i];
            Block2 diagonal{m, 0.0, 0.0, m};
            Block2 upper{};
            if (withElements)
            {
                const Block2 kd = symmetric(rw.stiffnessDiagonal[i]);
                const Block2 ku = symmetric(rw.stiffnessUpper[i]);
                for (std::size_t k = 0; k < C * C; ++k)
                {
                    diagonal[k] -= hd * rw.dampingDiagonal[i][k] + h2 * kd[k];
                    upper[k] = -(hd * rw.dampingUpper[i][k] + h2 * ku[k]);
                }
            }
            if (contact)
            {
                for (std::size_t k = 0; k < C * C; ++k)
                {
                    diagonal[k] += rw.contactDiagonal[i][k];
                    upper[k] += rw.contactUpper[i][k];
                }
            }
            rw.stepDiagonal[i] = diagonal;
            rw.stepUpper[i] = upper;
            BlockColumns<C, 3>& border = cond.border[i];
            border = rw.border[i];
            if (contact)
            {
                for (std::size_t k = 0; k < 3 * C; ++k)
                {
                    border[k] += rw.contactBorder[i][k];
                }
            }
        }
        if (!factorPeriodicBlockTridiagonal<C>(rw.stepDiagonal, rw.stepUpper, cond.factor()))
        {
            return false;
        }
        std::span<BlockColumns<C, 3>> z = cond.coupling;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            z[i] = cond.border[i];
        }
        solvePeriodicBlockTridiagonal<C, 3>(cond.inversePivot, cond.next, cond.wrap, z);
        // B^T Z, symmetric.
        std::array<double, 9> h{};
        for (std::size_t i = 0; i < nodes; ++i)
        {
            const BlockColumns<C, 3>& b = cond.border[i];
            for (std::size_t a = 0; a < 3; ++a)
            {
                for (std::size_t c = a; c < 3; ++c)
                {
                    h[a * 3 + c] += b[a] * z[i][c] + b[3 + a] * z[i][3 + c];
                }
            }
        }
        for (std::size_t a = 0; a < 3; ++a)
        {
            for (std::size_t c = 0; c < a; ++c)
            {
                h[a * 3 + c] = h[c * 3 + a];
            }
        }
        cond.hostMatrix = h;
        // tree -= J^T H J over the host's chain (ancestor pairs: the tree's sparsity is kept), lower triangle.
        const std::span<const std::size_t> chain = ring.hostChain;
        for (std::size_t a = 0; a < chain.size(); ++a)
        {
            const PlanarMotion& ja = rw.hostJacobian[a];
            const double ha0 = h[0] * ja.w + h[1] * ja.x + h[2] * ja.y;
            const double ha1 = h[3] * ja.w + h[4] * ja.x + h[5] * ja.y;
            const double ha2 = h[6] * ja.w + h[7] * ja.x + h[8] * ja.y;
            for (std::size_t b = a; b < chain.size(); ++b)
            {
                const PlanarMotion& jb = rw.hostJacobian[b];
                // chain[b] is an ancestor of chain[a] (or the same body), so its index is not larger.
                tree[model.coordinates()[chain[a]] * nt + model.coordinates()[chain[b]]] -=
                    jb.w * ha0 + jb.x * ha1 + jb.y * ha2;
            }
        }
    }
    return true;
}

void ringStepRhs(const Model& model, const State& state, const Workspace& ws, const double* tau,
                 const double* elementForce, double dt, double thetaStep, double* b) noexcept
{
    const auto rings = model.rings();
    const std::size_t n = model.dofCount();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        const RingWorkspace& rw = ws.rings[r];
        const std::size_t nodes = ring.nodeCount;
        const std::size_t first = ring.firstCoordinate;
        const double* wd = state.values.data() + n + first;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            std::array<double, C> rate{};
            if (thetaStep != 0.0)
            {
                const std::size_t previous = i == 0U ? nodes - 1U : i - 1U;
                const std::size_t next = i + 1U == nodes ? 0U : i + 1U;
                const Block2 kd = symmetric(rw.stiffnessDiagonal[i]);
                const Block2 kn = symmetric(rw.stiffnessUpper[i]);
                const Block2 kp = symmetric(rw.stiffnessUpper[previous]);
                for (std::size_t a = 0; a < C; ++a)
                {
                    double sum = 0.0;
                    for (std::size_t c = 0; c < C; ++c)
                    {
                        sum += kd[a * C + c] * wd[C * i + c] + kn[a * C + c] * wd[C * next + c] +
                               kp[a * C + c] * wd[C * previous + c];
                    }
                    rate[a] = sum;
                }
            }
            for (std::size_t a = 0; a < C; ++a)
            {
                const std::size_t j = first + C * i + a;
                const double force = tau == nullptr ? 0.0 : tau[j];
                const double element = elementForce == nullptr ? 0.0 : elementForce[j];
                b[j] = dt * (force + element - ws.bias[j] + thetaStep * rate[a]);
            }
        }
    }
}

void ringForward(const Model& model, const Workspace& ws, RingSystem system, double* x) noexcept
{
    const auto rings = model.rings();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        const RingWorkspace& rw = ws.rings[r];
        const RingCondensation& cond = system == RingSystem::Mass ? rw.mass : rw.step;
        const std::size_t nodes = ring.nodeCount;
        double* xw = x + ring.firstCoordinate;
        const std::span<BlockColumns<C, 1>> y = rw.solution;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            y[i] = {xw[C * i], xw[C * i + 1U]};
        }
        solvePeriodicBlockTridiagonal<C, 1>(cond.inversePivot, cond.next, cond.wrap, y);
        // g = B^T y, then the tree rows lose J^T g.
        double g0 = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            const BlockColumns<C, 3>& b = cond.border[i];
            xw[C * i] = y[i][0];
            xw[C * i + 1U] = y[i][1];
            g0 += b[0] * y[i][0] + b[3] * y[i][1];
            g1 += b[1] * y[i][0] + b[4] * y[i][1];
            g2 += b[2] * y[i][0] + b[5] * y[i][1];
        }
        const std::span<const std::size_t> chain = ring.hostChain;
        for (std::size_t c = 0; c < chain.size(); ++c)
        {
            const PlanarMotion& j = rw.hostJacobian[c];
            x[model.coordinates()[chain[c]]] -= j.w * g0 + j.x * g1 + j.y * g2;
        }
    }
}

void ringBackward(const Model& model, const Workspace& ws, RingSystem system, double* x) noexcept
{
    const auto rings = model.rings();
    for (std::size_t r = 0; r < rings.size(); ++r)
    {
        const DeformableRing& ring = rings[r];
        const RingWorkspace& rw = ws.rings[r];
        const RingCondensation& cond = system == RingSystem::Mass ? rw.mass : rw.step;
        // The host's motion from the tree's solution, then w = y - Z dV.
        double v0 = 0.0;
        double v1 = 0.0;
        double v2 = 0.0;
        const std::span<const std::size_t> chain = ring.hostChain;
        for (std::size_t c = 0; c < chain.size(); ++c)
        {
            const PlanarMotion& j = rw.hostJacobian[c];
            const double rate = x[model.coordinates()[chain[c]]];
            v0 += j.w * rate;
            v1 += j.x * rate;
            v2 += j.y * rate;
        }
        double* xw = x + ring.firstCoordinate;
        for (std::size_t i = 0; i < ring.nodeCount; ++i)
        {
            const BlockColumns<C, 3>& z = cond.coupling[i];
            xw[C * i] -= z[0] * v0 + z[1] * v1 + z[2] * v2;
            xw[C * i + 1U] -= z[3] * v0 + z[4] * v1 + z[5] * v2;
        }
    }
}

} // namespace detail

} // namespace sub0physics::multibody
