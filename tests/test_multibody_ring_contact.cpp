// Contact on a ring's casing (docs/WHEEL.md, "M2c status"): the reference soft-constraint law against its source's
// formulas, the soft-constraint solve against brute-force enumeration of active sets, the contact row against a
// finite difference of the attached point, the contact element's exact Jacobians against central differences of its
// generalized force (face and vertex witnesses, the host on a chain with revolute joints), and the structured step
// with contact loads against a dense solve of the same step.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include <doctest.h>

#include "sub0physics/multibody.hpp"
#include "test_ring_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using ringmodels::RingParameters;
using ringmodels::TreeDof;

namespace
{

double maxAbs(std::span<const double> v)
{
    double m = 0.0;
    for (const double x : v)
    {
        m = std::max(m, std::abs(x));
    }
    return m;
}

double maxDifference(std::span<const double> a, std::span<const double> b)
{
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        m = std::max(m, std::abs(a[i] - b[i]));
    }
    return m;
}

// World plane (i, j) position of node `k` of ring `r` at the workspace's kinematics.
std::array<double, 2> nodeWorld(const Model& model, const Workspace& ws, std::size_t r, std::size_t k)
{
    const PlanarTransform& t = ws.planar.worldToBody[model.rings()[r].host];
    const double lx = ws.rings[r].positionX[k];
    const double ly = ws.rings[r].positionY[k];
    return {t.x + t.c * lx - t.s * ly, t.y + t.s * lx + t.c * ly};
}

sub0math::Vec3d toWorld(const Model& model, double i, double j)
{
    sub0math::Vec3d v{};
    v[model.planar().axisI] = i;
    v[model.planar().axisJ] = j;
    return v;
}

void prepare(const Model& model, const State& s, Workspace& ws)
{
    REQUIRE(prepareLinearlyImplicitStep(model, s, ws).has_value());
}

// The contact element's geometry at a state: a face of fixed normal under the casing's node `edge` (t = 0), or a
// vertex whose nearest axis point is interior.
RingContactGeometry geometryAt(const Model& model, const Workspace& ws, std::size_t edge, ContactWitness witness,
                               std::array<double, 2> fixedPoint, std::array<double, 2> faceNormal, double skin)
{
    RingContactGeometry g;
    g.ring = 0;
    g.edge = edge;
    g.witness = witness;
    g.skin = skin;
    const std::size_t b = (edge + 1) % model.rings()[0].nodeCount;
    const auto pa = nodeWorld(model, ws, 0, edge);
    const auto pb = nodeWorld(model, ws, 0, b);
    if (witness == ContactWitness::Face)
    {
        g.axisT = 0.0;
        g.witnessPoint = toWorld(model, fixedPoint[0], fixedPoint[1]);
        g.normal = toWorld(model, faceNormal[0], faceNormal[1]);
        return g;
    }
    const double ex = pb[0] - pa[0];
    const double ey = pb[1] - pa[1];
    const double t = ((fixedPoint[0] - pa[0]) * ex + (fixedPoint[1] - pa[1]) * ey) / (ex * ex + ey * ey);
    g.axisT = std::clamp(t, 0.0, 1.0);
    const double px = pa[0] + g.axisT * ex;
    const double py = pa[1] + g.axisT * ey;
    const double rho = std::hypot(px - fixedPoint[0], py - fixedPoint[1]);
    g.witnessPoint = toWorld(model, fixedPoint[0], fixedPoint[1]);
    g.normal = toWorld(model, (px - fixedPoint[0]) / rho, (py - fixedPoint[1]) / rho);
    return g;
}

// A linear law in the penetration, unclipped (so central differences see one smooth law).
ContactLawValue linearLaw(double gap, double rate, double k, double c)
{
    return {.force = -k * gap - c * rate, .dForceDGap = -k, .dForceDGapRate = -c};
}

// The element's generalized force at a state, its geometry re-derived there.
std::vector<double> elementForce(const Model& model, const State& s, Workspace& ws, std::size_t edge,
                                 ContactWitness witness, std::array<double, 2> point, std::array<double, 2> normal,
                                 double k, double c, RingContactLoad* loadOut = nullptr)
{
    prepare(model, s, ws);
    const RingContactGeometry g = geometryAt(model, ws, edge, witness, point, normal, 0.008);
    double gap = 0.0;
    double rate = 0.0;
    REQUIRE(measureRingContact(model, s, ws, g, gap, rate).has_value());
    RingContactLoad load;
    REQUIRE(ringContactLoad(model, s, ws, g, linearLaw(gap, rate, k, c), load).has_value());
    std::vector<double> tau(model.dofCount(), 0.0);
    REQUIRE(ringLoadsToGeneralized(model, s, ws, std::span<const RingContactLoad>(&load, 1), tau, {}, {}).has_value());
    if (loadOut != nullptr)
    {
        *loadOut = load;
    }
    return tau;
}

} // namespace

TEST_CASE("soft-constraint row follows the reference's formulas")
{
    // The moving wheel's contact: solref (0.001, 1), solimp (0.99, 0.99, 0.001): a flat impedance of 0.99.
    const SoftConstraintParameters flat{.timeConstant = 0.001, .dampingRatio = 1.0, .impedance = {0.99, 0.99, 0.001, 0.5, 2.0}};
    const SoftConstraintRow r = softConstraintRow(flat, -2e-4, -0.3, 1.7, 3.125e-5);
    CHECK(r.impedance == 0.5 * (0.99 + 0.99));
    CHECK(r.stiffness == 1 / (0.99 * 0.99 * 0.001 * 0.001 * 1.0 * 1.0));
    CHECK(r.damping == 2 / (0.99 * 0.001));
    CHECK(r.regularisation == (1 - r.impedance) * 1.7 / r.impedance);
    CHECK(r.referenceAcceleration == -r.damping * -0.3 - r.stiffness * r.impedance * -2e-4);
    // The time constant is raised to 2 h.
    const SoftConstraintRow slow = softConstraintRow(flat, -2e-4, -0.3, 1.7, 1e-3);
    CHECK(slow.damping == 2 / (0.99 * 0.002));
    // A graded impedance (the default solimp), below and above the midpoint, and saturated.
    const SoftConstraintParameters graded{};
    const double low = softConstraintRow(graded, -0.0002, 0.0, 1.0, 1e-3).impedance;
    const double high = softConstraintRow(graded, -0.0008, 0.0, 1.0, 1e-3).impedance;
    const double full = softConstraintRow(graded, -0.002, 0.0, 1.0, 1e-3).impedance;
    CHECK(low == doctest::Approx(0.9 + 2.0 * 0.2 * 0.2 * 0.05).epsilon(1e-14));
    CHECK(high == doctest::Approx(0.9 + (1 - 2.0 * 0.2 * 0.2) * 0.05).epsilon(1e-14));
    CHECK(full == 0.95);
}

TEST_CASE("soft-constraint solve is the unique minimiser: brute force over active sets")
{
    const Model model = ringmodels::ringRig({.nodes = 8}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    const State s = ringmodels::randomState(model, 3);
    prepare(model, s, ws);
    REQUIRE(factorMassMatrix(model, ws).has_value());
    SoftConstraintSolver solver(model, 8);
    for (int trial = 0; trial < 40; ++trial)
    {
        const std::size_t k = 1 + static_cast<std::size_t>(trial % 6);
        std::vector<double> jac(k * n, 0.0);
        std::vector<double> aref(k);
        std::vector<double> reg(k);
        std::vector<SoftRowKind> kinds(k, SoftRowKind::Unilateral);
        for (std::size_t r = 0; r < k; ++r)
        {
            RingContactPoint p;
            p.edge = static_cast<std::size_t>(trial + 3 * static_cast<int>(r)) % 8;
            p.weight0 = 0.5 + 0.4 * u(rng);
            p.weight1 = 1.0 - p.weight0;
            const auto a = nodeWorld(model, ws, 0, p.edge);
            p.point = toWorld(model, a[0] + 0.01 * u(rng), a[1] + 0.01 * u(rng));
            const double angle = 3.0 * u(rng);
            p.normal = toWorld(model, std::cos(angle), std::sin(angle));
            REQUIRE(ringContactJacobian(model, ws, p, std::span<double>(jac.data() + r * n, n)).has_value());
            aref[r] = 50.0 * u(rng);
            reg[r] = 0.01 * (1.5 + u(rng));
            if (trial % 5 == 4 && r == 0)
            {
                kinds[r] = SoftRowKind::Equality;
            }
        }
        std::vector<double> a0(n);
        for (double& v : a0)
        {
            v = 10.0 * u(rng);
        }
        std::vector<double> force(k);
        std::vector<double> generalized(n);
        const auto solved = solver.solve(model, ws, jac, aref, reg, a0, kinds, force, generalized);
        REQUIRE(solved.has_value());
        // The primal cost at a = a0 + M^-1 J^T f, against every active set's own stationary point.
        const auto cost = [&](std::span<const double> f) {
            std::vector<double> x(n, 0.0);
            for (std::size_t r = 0; r < k; ++r)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    x[i] += jac[r * n + i] * f[r];
                }
            }
            std::vector<double> y = x;
            REQUIRE(solveMassMatrix(model, ws, y).has_value());
            double value = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                value += 0.5 * x[i] * y[i]; // 1/2 (a - a0)^T M (a - a0)
            }
            for (std::size_t r = 0; r < k; ++r)
            {
                double jar = -aref[r];
                for (std::size_t i = 0; i < n; ++i)
                {
                    jar += jac[r * n + i] * (a0[i] + y[i]);
                }
                if (kinds[r] == SoftRowKind::Equality || jar < 0.0)
                {
                    value += 0.5 * jar * jar / reg[r];
                }
            }
            return value;
        };
        const double best = cost(force);
        // Complementarity: f_i = -(J_i a - aref_i) / R_i on active rows.
        for (std::size_t r = 0; r < k; ++r)
        {
            if (kinds[r] == SoftRowKind::Unilateral)
            {
                CHECK(force[r] >= 0.0);
            }
        }
        // Perturbations of the forces never lower the cost (strict convexity in f-space of the primal).
        for (int probe = 0; probe < 20; ++probe)
        {
            std::vector<double> f = force;
            for (double& v : f)
            {
                v += 1e-3 * (1.0 + std::abs(v)) * u(rng);
            }
            CHECK(cost(f) >= best - 1e-12 * std::max(1.0, std::abs(best)));
        }
    }
}

TEST_CASE("contact row is the attached point's rate along the normal")
{
    const Model model = ringmodels::ringRig({.nodes = 16}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    const State s = ringmodels::randomState(model, 5);
    prepare(model, s, ws);
    RingContactPoint p;
    p.edge = 15; // the wrap edge (15, 0)
    p.weight0 = 0.3;
    p.weight1 = 0.7;
    const auto a = nodeWorld(model, ws, 0, 15);
    const auto b = nodeWorld(model, ws, 0, 0);
    const PlanarTransform h0 = ws.planar.worldToBody[model.rings()[0].host];
    // A point off the axis, attached to the host's rotation (the reference's attachment).
    const double ox = 0.3 * a[0] + 0.7 * b[0] + 0.004;
    const double oy = 0.3 * a[1] + 0.7 * b[1] - 0.003;
    p.point = toWorld(model, ox, oy);
    p.normal = toWorld(model, 0.6, 0.8);
    std::vector<double> row(n);
    REQUIRE(ringContactJacobian(model, ws, p, row).has_value());
    // Offset in the host frame, carried with the host: p(q) = sum w x_c(q) + R(q) d.
    const double dx = ox - (0.3 * a[0] + 0.7 * b[0]);
    const double dy = oy - (0.3 * a[1] + 0.7 * b[1]);
    const double lx = h0.c * dx + h0.s * dy;
    const double ly = -h0.s * dx + h0.c * dy;
    const auto pointAt = [&](double eps) {
        State t = s;
        for (std::size_t i = 0; i < n; ++i)
        {
            t.q()[i] += eps * s.qd()[i];
        }
        prepare(model, t, ws);
        const auto pa = nodeWorld(model, ws, 0, 15);
        const auto pb = nodeWorld(model, ws, 0, 0);
        const PlanarTransform& h = ws.planar.worldToBody[model.rings()[0].host];
        return std::array<double, 2>{0.3 * pa[0] + 0.7 * pb[0] + h.c * lx - h.s * ly,
                                     0.3 * pa[1] + 0.7 * pb[1] + h.s * lx + h.c * ly};
    };
    const double eps = 1e-6;
    const auto plus = pointAt(eps);
    const auto minus = pointAt(-eps);
    const double rate = (0.6 * (plus[0] - minus[0]) + 0.8 * (plus[1] - minus[1])) / (2 * eps);
    double jqd = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        jqd += row[i] * s.qd()[i];
    }
    CHECK(std::abs(jqd - rate) < 1e-8 * std::max(1.0, std::abs(rate)));
}

TEST_CASE("contact element Jacobians are exact: central differences of its generalized force")
{
    // The rear hub sits on a chain with two revolute joints (frame pitch, hub spin) and two slides.
    const Model model = ringmodels::ringRig({.nodes = 16}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    struct Case
    {
        ContactWitness witness;
        std::size_t edge;
    };
    for (const Case cs : {Case{ContactWitness::Face, 1}, Case{ContactWitness::Vertex, 0}, Case{ContactWitness::Vertex, 14}})
    {
        CAPTURE(static_cast<int>(cs.witness));
        CAPTURE(cs.edge);
        const State s = ringmodels::randomState(model, 11);
        prepare(model, s, ws);
        // An obstacle near the edge: a face through a point below node `edge`, or a vertex below the edge's middle.
        const std::size_t b = (cs.edge + 1) % 16;
        const auto pa = nodeWorld(model, ws, 0, cs.edge);
        const auto pb = nodeWorld(model, ws, 0, b);
        const std::array<double, 2> mid{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])};
        // The outward direction of the casing at the edge: from the hub.
        const PlanarTransform& h = ws.planar.worldToBody[model.rings()[0].host];
        double ox = mid[0] - h.x;
        double oy = mid[1] - h.y;
        const double ol = std::hypot(ox, oy);
        ox /= ol;
        oy /= ol;
        const std::array<double, 2> point = cs.witness == ContactWitness::Face
                                                ? std::array<double, 2>{pa[0] + 0.003 * ox, pa[1] + 0.003 * oy}
                                                : std::array<double, 2>{mid[0] + 0.004 * ox + 0.002 * oy, mid[1] + 0.004 * oy - 0.002 * ox};
        const std::array<double, 2> normal{-ox, -oy};
        const double k = 2e5;
        const double c = 150.0;
        RingContactLoad load;
        (void)elementForce(model, s, ws, cs.edge, cs.witness, point, normal, k, c, &load);
        std::vector<double> tau(n, 0.0);
        std::vector<double> stiffness(n * n, 0.0);
        std::vector<double> damping(n * n, 0.0);
        prepare(model, s, ws);
        REQUIRE(ringLoadsToGeneralized(model, s, ws, std::span<const RingContactLoad>(&load, 1), tau, stiffness, damping)
                    .has_value());
        const double scale = maxAbs(stiffness);
        REQUIRE(scale > 0.0);
        double worstK = 0.0;
        double worstD = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            const double eps = 1e-7;
            State plus = s;
            State minus = s;
            plus.q()[j] += eps;
            minus.q()[j] -= eps;
            const auto fp = elementForce(model, plus, ws, cs.edge, cs.witness, point, normal, k, c);
            const auto fm = elementForce(model, minus, ws, cs.edge, cs.witness, point, normal, k, c);
            State vp = s;
            State vm = s;
            vp.qd()[j] += eps;
            vm.qd()[j] -= eps;
            const auto gp = elementForce(model, vp, ws, cs.edge, cs.witness, point, normal, k, c);
            const auto gm = elementForce(model, vm, ws, cs.edge, cs.witness, point, normal, k, c);
            for (std::size_t i = 0; i < n; ++i)
            {
                worstK = std::max(worstK, std::abs((fp[i] - fm[i]) / (2 * eps) - stiffness[i * n + j]));
                worstD = std::max(worstD, std::abs((gp[i] - gm[i]) / (2 * eps) - damping[i * n + j]));
            }
        }
        CHECK(worstK < 1e-8 * scale); // central differences: about 1e-9 measured; dropping a curvature term gives 4e-3
        CHECK(worstD < 1e-6 * std::max(1.0, maxAbs(damping)));
        // The generalized force is the load's force mapped: J^T F with J the contact row at the casing witness.
        prepare(model, s, ws);
        const RingContactGeometry g = geometryAt(model, ws, cs.edge, cs.witness, point, normal, 0.008);
        RingContactPoint p{.ring = 0, .edge = cs.edge, .weight0 = 1.0 - g.axisT, .weight1 = g.axisT};
        const auto qa = nodeWorld(model, ws, 0, cs.edge);
        const auto qb = nodeWorld(model, ws, 0, b);
        p.point = toWorld(model, (1 - g.axisT) * qa[0] + g.axisT * qb[0], (1 - g.axisT) * qa[1] + g.axisT * qb[1]);
        p.normal = g.normal;
        std::vector<double> row(n);
        REQUIRE(ringContactJacobian(model, ws, p, row).has_value());
        double worstF = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            worstF = std::max(worstF, std::abs(row[i] * load.normalForce - tau[i]));
        }
        CHECK(worstF < 1e-12 * std::max(1.0, maxAbs(tau)));
    }
}

TEST_CASE("structured step with contact loads equals the dense step")
{
    // No ring elements (free nodes), so the step matrix is M plus the loads' terms and a dense solve is the oracle.
    RingParameters p{.nodes = 16};
    p.foundation = 0.0;
    p.extension = 0.0;
    p.viscosity = 0.0;
    const Model model = ringmodels::ringRig(p, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    const double h = 1e-3;
    for (const ImplicitTerms terms : {ImplicitTerms::Full, ImplicitTerms::DampingOnly})
    {
        const State s = ringmodels::randomState(model, 21);
        // Two loads on neighbouring edges of the ring, the second on the wrap edge; no gap-rate dependence in K's
        // chain curvature term (c = 0 for one load), which the structured step keeps only on the chain's rows.
        std::array<RingContactLoad, 2> loads{};
        prepare(model, s, ws);
        {
            const auto pa = nodeWorld(model, ws, 0, 15);
            const auto pb = nodeWorld(model, ws, 0, 0);
            const RingContactGeometry g = geometryAt(model, ws, 15, ContactWitness::Vertex,
                                                     {0.5 * (pa[0] + pb[0]) + 0.003, 0.5 * (pa[1] + pb[1]) - 0.002},
                                                     {0.0, 0.0}, 0.008);
            double gap = 0.0;
            double rate = 0.0;
            REQUIRE(measureRingContact(model, s, ws, g, gap, rate).has_value());
            REQUIRE(ringContactLoad(model, s, ws, g, linearLaw(gap, rate, 3e5, 0.0), loads[0]).has_value());
            const auto qa = nodeWorld(model, ws, 0, 3);
            const RingContactGeometry f = geometryAt(model, ws, 3, ContactWitness::Face, {qa[0], qa[1] - 0.005},
                                                     {0.0, 1.0}, 0.008);
            REQUIRE(measureRingContact(model, s, ws, f, gap, rate).has_value());
            REQUIRE(ringContactLoad(model, s, ws, f, linearLaw(gap, rate, 2e5, 0.0), loads[1]).has_value());
        }
        std::vector<double> tauApplied(n, 0.0);
        std::vector<double> stiffness(n * n, 0.0);
        std::vector<double> damping(n * n, 0.0);
        REQUIRE(ringLoadsToGeneralized(model, s, ws, loads, tauApplied, stiffness, damping).has_value());
        // Dense M and c by inverse dynamics.
        std::vector<double> mass(n * n);
        std::vector<double> bias(n);
        std::vector<double> zero(n, 0.0);
        REQUIRE(inverseDynamics(model, s, zero, bias, ws).has_value());
        for (std::size_t j = 0; j < n; ++j)
        {
            std::vector<double> e(n, 0.0);
            e[j] = 1.0;
            std::vector<double> col(n);
            REQUIRE(inverseDynamics(model, s, e, col, ws).has_value());
            for (std::size_t i = 0; i < n; ++i)
            {
                mass[i * n + j] = col[i] - bias[i];
            }
        }
        const bool withK = terms == ImplicitTerms::Full;
        std::vector<double> a(n * n);
        std::vector<double> rhs(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            double kqd = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                const double ks = withK ? 0.5 * (stiffness[i * n + j] + stiffness[j * n + i]) : 0.0;
                const double ds = 0.5 * (damping[i * n + j] + damping[j * n + i]);
                a[i * n + j] = mass[i * n + j] - h * ds - h * h * ks;
                kqd += ks * s.qd()[j];
            }
            rhs[i] = h * (tauApplied[i] - bias[i] + h * kqd);
        }
        // Gaussian elimination with partial pivoting.
        for (std::size_t col = 0; col < n; ++col)
        {
            std::size_t pivot = col;
            for (std::size_t r = col + 1; r < n; ++r)
            {
                if (std::abs(a[r * n + col]) > std::abs(a[pivot * n + col]))
                {
                    pivot = r;
                }
            }
            for (std::size_t c = 0; c < n; ++c)
            {
                std::swap(a[col * n + c], a[pivot * n + c]);
            }
            std::swap(rhs[col], rhs[pivot]);
            for (std::size_t r = col + 1; r < n; ++r)
            {
                const double f = a[r * n + col] / a[col * n + col];
                for (std::size_t c = col; c < n; ++c)
                {
                    a[r * n + c] -= f * a[col * n + c];
                }
                rhs[r] -= f * rhs[col];
            }
        }
        for (std::size_t r = n; r-- > 0;)
        {
            double sum = rhs[r];
            for (std::size_t c = r + 1; c < n; ++c)
            {
                sum -= a[r * n + c] * rhs[c];
            }
            rhs[r] = sum / a[r * n + r];
        }
        State stepped = s;
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{.tau = {}, .ringLoads = loads}, s.tick, sub0math::f64::Time::fromSi(h),
                                          stepped, ws, terms)
                    .has_value());
        std::vector<double> change(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            change[i] = stepped.qd()[i] - s.qd()[i];
        }
        CHECK(maxDifference(change, rhs) < 1e-10 * maxAbs(rhs));
        // The prepared step and the one-call step agree bit for bit with loads too.
        State prepared = s;
        REQUIRE(prepareLinearlyImplicitStep(model, prepared, ws, terms).has_value());
        REQUIRE(completeLinearlyImplicitStep(model, Forcing{.tau = {}, .ringLoads = loads}, prepared.tick,
                                             sub0math::f64::Time::fromSi(h), prepared, ws)
                    .has_value());
        CHECK(prepared == stepped);
        // Without loads the step is unchanged by the earlier contact terms (they are cleared).
        State plain = s;
        State plainAgain = s;
        Workspace fresh(model);
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{}, s.tick, sub0math::f64::Time::fromSi(h), plain, ws, terms).has_value());
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{}, s.tick, sub0math::f64::Time::fromSi(h), plainAgain, fresh, terms)
                    .has_value());
        CHECK(plain == plainAgain);
    }
}
