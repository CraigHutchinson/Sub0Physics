// The deformable ring substructure (docs/WHEEL.md section 2.3, step M2a): the periodic block solver, the canonical
// tables, and the condensed ring against its expanded oracle (the same description as ordinary bodies, coordinate
// elements and tendons, factorised densely): dynamics, energy, element forces and Jacobians, steps, runs and
// reactions to rounding; mass properties; free flight; the energy ledger; modal frequencies against the closed form.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <random>
#include <span>
#include <vector>

#include <doctest.h>

#include "sub0physics/multibody.hpp"
#include "sub0physics/solver/periodic_block_tridiagonal.hpp"
#include "test_ring_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using ringmodels::RingParameters;
using ringmodels::TreeDof;

namespace
{

// A random symmetric positive definite periodic block-tridiagonal system, solved densely and by the block solver.
template<std::size_t Dim>
double periodicSolveError(std::size_t n, unsigned seed)
{
    std::mt19937_64 rng(seed);
    const auto u = [](std::mt19937_64& g) { return ringmodels::uniformSigned(g); };
    std::vector<Block<Dim>> diagonal(n);
    std::vector<Block<Dim>> upper(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        for (std::size_t i = 0; i < Dim; ++i)
        {
            for (std::size_t j = 0; j <= i; ++j)
            {
                const double v = i == j ? 20.0 + u(rng) : u(rng);
                diagonal[k][i * Dim + j] = v;
                diagonal[k][j * Dim + i] = v;
            }
        }
        for (double& v : upper[k])
        {
            v = 4.0 * u(rng);
        }
    }
    const std::size_t m = n * Dim;
    std::vector<double> a(m * m, 0.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t r = k * Dim;
        const std::size_t c = ((k + 1) % n) * Dim;
        for (std::size_t i = 0; i < Dim; ++i)
        {
            for (std::size_t j = 0; j < Dim; ++j)
            {
                a[(r + i) * m + r + j] += diagonal[k][i * Dim + j];
                a[(r + i) * m + c + j] += upper[k][i * Dim + j];
                a[(c + j) * m + r + i] += upper[k][i * Dim + j];
            }
        }
    }
    constexpr std::size_t Columns = 2;
    std::vector<BlockColumns<Dim, Columns>> x(n);
    std::vector<double> b(m * Columns);
    for (std::size_t k = 0; k < n; ++k)
    {
        for (double& v : x[k])
        {
            v = u(rng);
        }
        for (std::size_t i = 0; i < Dim; ++i)
        {
            for (std::size_t c = 0; c < Columns; ++c)
            {
                b[(k * Dim + i) * Columns + c] = x[k][i * Columns + c];
            }
        }
    }
    std::vector<Block<Dim>> pivots(n);
    std::vector<Block<Dim>> next(n);
    std::vector<Block<Dim>> wrap(n);
    REQUIRE(factorPeriodicBlockTridiagonal<Dim>(diagonal, upper, {pivots, next, wrap}));
    solvePeriodicBlockTridiagonal<Dim, Columns>(pivots, next, wrap, x);
    // Residual A x - b, relative.
    double worst = 0.0;
    for (std::size_t row = 0; row < m; ++row)
    {
        for (std::size_t c = 0; c < Columns; ++c)
        {
            double sum = 0.0;
            for (std::size_t col = 0; col < m; ++col)
            {
                sum += a[row * m + col] * x[col / Dim][(col % Dim) * Columns + c];
            }
            worst = std::max(worst, std::abs(sum - b[row * Columns + c]));
        }
    }
    return worst;
}

double relativeDifference(std::span<const double> a, std::span<const double> b)
{
    double scale = 0.0;
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        scale = std::max(scale, std::abs(b[i]));
        worst = std::max(worst, std::abs(a[i] - b[i]));
    }
    return worst / std::max(scale, 1e-300);
}

struct Pair
{
    Model condensed;
    Model expanded;
};

Pair rigPair(const RingParameters& p, bool secondRing = true)
{
    return {ringmodels::ringRig(p, false, secondRing), ringmodels::ringRig(p, true, secondRing)};
}

// The ring's dense dtau/dw (or dtau/dwd) from the block storage, in the model's ring coordinates.
std::vector<double> denseRingJacobian(const Model& model, const Workspace& ws, bool damping)
{
    const std::size_t n = model.dofCount();
    const std::size_t t = model.treeDofCount();
    const std::size_t m = n - t;
    std::vector<double> dense(m * m, 0.0);
    for (std::size_t r = 0; r < model.rings().size(); ++r)
    {
        const DeformableRing& ring = model.rings()[r];
        const RingWorkspace& rw = ws.rings[r];
        const std::size_t base = ring.firstCoordinate - t;
        const std::size_t nodes = ring.nodeCount;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            const std::size_t j = (i + 1) % nodes;
            const auto& d = damping ? rw.dampingDiagonal[i] : rw.stiffnessDiagonal[i];
            const auto& up = damping ? rw.dampingUpper[i] : rw.stiffnessUpper[i];
            for (std::size_t a = 0; a < 2; ++a)
            {
                for (std::size_t b = 0; b < 2; ++b)
                {
                    dense[(base + 2 * i + a) * m + base + 2 * i + b] += d[a * 2 + b];
                    dense[(base + 2 * i + a) * m + base + 2 * j + b] += up[a * 2 + b];
                    dense[(base + 2 * j + a) * m + base + 2 * i + b] += up[a * 2 + b];
                }
            }
        }
    }
    return dense;
}

// The ring rows and columns of a dense dofCount x dofCount matrix.
std::vector<double> ringPart(const std::vector<double>& full, std::size_t n, std::size_t t)
{
    const std::size_t m = n - t;
    std::vector<double> out(m * m);
    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < m; ++j)
        {
            out[i * m + j] = full[(t + i) * n + t + j];
        }
    }
    return out;
}

double totalEnergy(const EnergyAccount& e)
{
    return e.kinetic.si() + e.potential.si();
}

} // namespace

TEST_CASE("periodic block-tridiagonal solve agrees with a dense solve")
{
    for (const std::size_t n : {std::size_t{3}, std::size_t{4}, std::size_t{7}, std::size_t{64}})
    {
        CHECK(periodicSolveError<1>(n, 1) < 1e-13);
        CHECK(periodicSolveError<2>(n, 2) < 1e-13);
        CHECK(periodicSolveError<3>(n, 3) < 1e-13);
    }
    // An indefinite block is refused.
    std::vector<Block<2>> diagonal(4, Block<2>{1.0, 0.0, 0.0, -1.0});
    std::vector<Block<2>> upper(4, Block<2>{});
    std::vector<Block<2>> f(4);
    std::vector<Block<2>> g(4);
    std::vector<Block<2>> h(4);
    CHECK_FALSE(factorPeriodicBlockTridiagonal<2>(diagonal, upper, {f, g, h}));
}

TEST_CASE("canonical ring tables are exact and symmetric")
{
    // The reference's chord ratios for 32, 64 and 128 nodes, bit for bit.
    CHECK(*canonicalChordRatio(32) == 0.1960342806591212);
    CHECK(*canonicalChordRatio(64) == 0.09813534865483603);
    CHECK(*canonicalChordRatio(128) == 0.049082457045824576);
    CHECK_FALSE(canonicalChordRatio(48).has_value());
    CHECK_FALSE(canonicalRingDirection(12, 0).has_value());
    for (const std::size_t n : {std::size_t{8}, std::size_t{16}, std::size_t{64}, std::size_t{128}})
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            const auto d = *canonicalRingDirection(n, i);
            const auto mirror = *canonicalRingDirection(n, n - i);
            CHECK(d.sine == -mirror.sine); // exact reflection symmetry
            CHECK(d.cosine == mirror.cosine);
            CHECK(std::abs(d.sine * d.sine + d.cosine * d.cosine - 1.0) < 4e-16);
            const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
            CHECK(std::abs(d.sine - std::sin(angle)) < 1e-15);
        }
    }
}

TEST_CASE("ring descriptions are validated")
{
    ModelBuilder builder;
    const auto hubs = ringmodels::addTree(builder, 9.81);
    auto ring = *uniformRing(ringmodels::ringSpec("ring", hubs[0], {}));
    SUBCASE("mismatched sizes")
    {
        ring.nodeMass.pop_back();
        builder.addRing(ring);
        CHECK_FALSE(builder.build().has_value());
    }
    SUBCASE("axis out of the plane")
    {
        ring.displacementAxes[1] = {0.0, 1.0, 0.0};
        builder.addRing(ring);
        CHECK_FALSE(builder.build().has_value());
    }
    SUBCASE("duplicate material ids")
    {
        ring.nodeIds.assign(ring.restPositions.size(), 7U);
        builder.addRing(ring);
        CHECK_FALSE(builder.build().has_value());
    }
    SUBCASE("general formulation")
    {
        builder.addRing(ring);
        builder.setFormulation(Formulation::General);
        CHECK_FALSE(builder.build().has_value());
    }
    SUBCASE("valid")
    {
        builder.addRing(ring);
        const auto model = builder.build();
        REQUIRE(model.has_value());
        CHECK(model->dofCount() == TreeDof + 2 * 16);
        CHECK(model->treeDofCount() == TreeDof);
        CHECK(model->findRing("ring").has_value());
        CHECK(model->formulation() == Formulation::Planar);
    }
    CHECK_FALSE(uniformRing(ringmodels::ringSpec("bad", hubs[0], {.nodes = 12})).has_value());
}

TEST_CASE("ring mass, first moment and spin inertia")
{
    // Undeformed and uniform: the hub's composite gains the ring's mass at its centre and m (R^2 + skin^2 / 2) of spin
    // inertia (docs/WHEEL.md section 9), read from the tree's mass matrix (the spin coordinate's diagonal entry).
    for (const std::size_t nodes : {std::size_t{8}, std::size_t{64}})
    {
        RingParameters p{.nodes = nodes};
        const Model model = ringmodels::ringRig(p, false, false);
        Workspace ws(model);
        State s = State::atRest(model);
        REQUIRE(computeMassMatrix(model, s, ws).has_value());
        const double spin = ws.massMatrix[3 * TreeDof + 3];
        const double expected = 0.04 + p.mass * (p.radius * p.radius + p.skin * p.skin / 2.0);
        CHECK(std::abs(spin - expected) < 1e-14 * expected);
        // The travel coordinate's entry: frame, both hubs, ring.
        CHECK(std::abs(ws.massMatrix[0] - (20.0 + 3.0 + p.mass)) < 1e-13);
        // Rest first moment about the hub: zero for the symmetric tables.
        double mx = 0.0;
        double my = 0.0;
        const DeformableRing& ring = model.rings()[0];
        for (std::size_t i = 0; i < nodes; ++i)
        {
            mx += ring.mass[i] * ring.restX[i];
            my += ring.mass[i] * ring.restY[i];
        }
        CHECK(std::abs(mx) < 1e-16);
        CHECK(std::abs(my) < 1e-16);
    }
}

TEST_CASE("condensed ring agrees with its expanded oracle: dynamics, energy, elements and reactions")
{
    for (const std::size_t nodes : {std::size_t{8}, std::size_t{16}})
    {
        CAPTURE(nodes);
        const Pair pair = rigPair({.nodes = nodes});
        const Model& c = pair.condensed;
        const Model& e = pair.expanded;
        REQUIRE(c.dofCount() == e.dofCount());
        CHECK(c.formulation() == Formulation::Planar);
        CHECK(e.formulation() == Formulation::Planar);
        Workspace wc(c);
        Workspace we(e);
        const std::size_t n = c.dofCount();
        for (unsigned seed = 1; seed <= 3; ++seed)
        {
            const State s = ringmodels::randomState(c, seed);

            // Inverse dynamics with random accelerations: M qdd + c - tau_e, every coordinate.
            std::mt19937_64 rng(seed + 100);
            const auto u = [](std::mt19937_64& g) { return ringmodels::uniformSigned(g); };
            std::vector<double> qdd(n);
            for (double& v : qdd)
            {
                v = 50.0 * u(rng);
            }
            std::vector<double> tc(n);
            std::vector<double> te(n);
            REQUIRE(inverseDynamics(c, s, qdd, tc, wc).has_value());
            REQUIRE(inverseDynamics(e, s, qdd, te, we).has_value());
            CHECK(relativeDifference(tc, te) < 1e-12);

            // Forward dynamics (the mass condensation) with applied forces on every coordinate.
            std::vector<double> tau(n);
            for (double& v : tau)
            {
                v = 10.0 * u(rng);
            }
            REQUIRE(forwardDynamics(c, s, tau, wc).has_value());
            REQUIRE(forwardDynamics(e, s, tau, we).has_value());
            CHECK(relativeDifference(wc.qdd, we.qdd) < 1e-11);
            // ... and solveMassMatrix after it.
            std::vector<double> xc = tau;
            std::vector<double> xe = tau;
            REQUIRE(solveMassMatrix(c, wc, xc).has_value());
            REQUIRE(solveMassMatrix(e, we, xe).has_value());
            CHECK(relativeDifference(xc, xe) < 1e-11);

            // Energy.
            const auto ec = computeEnergy(c, s, wc);
            const auto ee = computeEnergy(e, s, we);
            REQUIRE(ec.has_value());
            REQUIRE(ee.has_value());
            CHECK(std::abs(ec->kinetic.si() - ee->kinetic.si()) < 1e-12 * ee->kinetic.si());
            CHECK(std::abs(ec->potential.si() - ee->potential.si()) < 1e-12 * std::abs(ee->potential.si()) + 1e-12);

            // Element forces, power and Jacobians (block storage against the dense oracle).
            REQUIRE(evaluateElements(c, s, wc).has_value());
            REQUIRE(evaluateElements(e, s, we).has_value());
            CHECK(relativeDifference(wc.elementForce, we.elementForce) < 1e-12);
            CHECK(std::abs(wc.elementPotential - we.elementPotential) < 1e-12 * we.elementPotential);
            CHECK(std::abs(wc.elementDissipationPower - we.elementDissipationPower) <
                  1e-12 * we.elementDissipationPower);
            const auto kc = denseRingJacobian(c, wc, false);
            const auto ke = ringPart(we.elementStiffness, n, TreeDof);
            CHECK(relativeDifference(kc, ke) < 1e-12);
            const auto dc = denseRingJacobian(c, wc, true);
            const auto de = ringPart(we.elementDamping, n, TreeDof);
            CHECK(relativeDifference(dc, de) < 1e-12);

            // Reactions at every tree joint (the ring's inertia and weight reach the hubs).
            REQUIRE(computeReactions(c, s, tau, wc).has_value());
            REQUIRE(computeReactions(e, s, tau, we).has_value());
            for (std::size_t b = 0; b < TreeDof; ++b)
            {
                const auto& rc = wc.jointReactions[b];
                const auto& re = we.jointReactions[b];
                const double scale = std::max(1.0, sub0math::length(re.force) + sub0math::length(re.moment));
                CHECK(sub0math::length(rc.force - re.force) < 1e-11 * scale);
                CHECK(sub0math::length(rc.moment - re.moment) < 1e-11 * scale);
            }
        }
    }
}

TEST_CASE("condensed ring steps agree with the expanded oracle over a run")
{
    for (const std::size_t nodes : {std::size_t{8}, std::size_t{16}})
    {
        CAPTURE(nodes);
        const Pair pair = rigPair({.nodes = nodes});
        const Model& c = pair.condensed;
        const Model& e = pair.expanded;
        Workspace wc(c);
        Workspace we(e);
        const auto h = sub0math::f64::Time::fromSi(1e-4);
        struct Case
        {
            const char* name;
            int scheme; // 0 implicit Euler, 1 midpoint, 2 semi-implicit Euler
            ImplicitTerms terms;
            double tolerance;
        };
        for (const Case k : {Case{"implicit Euler, full", 0, ImplicitTerms::Full, 1e-10},
                             Case{"implicit Euler, damping only", 0, ImplicitTerms::DampingOnly, 1e-10},
                             Case{"implicit Euler, damping only in the mass pattern", 0,
                                  ImplicitTerms::DampingOnlyMassPattern, 1e-10},
                             Case{"implicit midpoint", 1, ImplicitTerms::Full, 1e-10},
                             Case{"semi-implicit Euler", 2, ImplicitTerms::Full, 1e-10}})
        {
            CAPTURE(k.name);
            State sc = ringmodels::randomState(c, 7);
            State se = sc;
            const double h0 = k.scheme == 2 ? 2e-5 : h.si(); // the explicit step needs a smaller step for the ring
            const auto step = sub0math::f64::Time::fromSi(h0);
            for (std::uint64_t t = 0; t < 300; ++t)
            {
                const auto run = [&](const Model& m, State& s, Workspace& w) {
                    switch (k.scheme)
                    {
                    case 0:
                        return stepLinearlyImplicitEuler(m, {}, t, step, s, w, k.terms);
                    case 1:
                        return stepLinearlyImplicitMidpoint(m, {}, t, step, s, w, k.terms);
                    default:
                        return stepSemiImplicitEuler(m, {}, t, step, s, w);
                    }
                };
                REQUIRE(run(c, sc, wc).has_value());
                REQUIRE(run(e, se, we).has_value());
                if (t == 0)
                {
                    CHECK(relativeDifference(sc.values, se.values) < 1e-13);
                }
            }
            CHECK(relativeDifference(sc.values, se.values) < k.tolerance);
            CHECK(std::abs(sc.dissipated.si() - se.dissipated.si()) < 1e-10 * se.dissipated.si());
            const auto ec = computeEnergy(c, sc, wc);
            const auto ee = computeEnergy(e, se, we);
            CHECK(std::abs(totalEnergy(*ec) - totalEnergy(*ee)) < 1e-9 * std::abs(totalEnergy(*ee)));
        }
        // RK4 drives the ring through forwardDynamics.
        State sc = ringmodels::randomState(c, 9);
        State se = sc;
        Rk4Stepper rc(c, wc);
        Rk4Stepper re(e, we);
        const auto step = sub0math::f64::Time::fromSi(1e-5);
        for (std::uint64_t t = 0; t < 100; ++t)
        {
            REQUIRE(rc.step(c, {}, t, step, sc, wc).has_value());
            REQUIRE(re.step(e, {}, t, step, se, we).has_value());
        }
        CHECK(relativeDifference(sc.values, se.values) < 1e-10);
    }
}

TEST_CASE("damping only in the mass pattern drops exactly the ring's neighbour damping")
{
    // The reference engine's implicitfast (ImplicitTerms::DampingOnlyMassPattern): the edges' dampers stay implicit on
    // each node's own block and explicit between neighbours. It differs from DampingOnly by O(h c / m) on the step,
    // and equals it without edge viscosity.
    const auto h = sub0math::f64::Time::fromSi(1e-4);
    for (const double viscosity : {0.5, 0.0})
    {
        CAPTURE(viscosity);
        RingParameters p{.nodes = 16};
        p.viscosity = viscosity;
        const Model c = ringmodels::ringRig(p, false);
        Workspace w(c);
        State a = ringmodels::randomState(c, 3);
        State b = a;
        REQUIRE(stepLinearlyImplicitEuler(c, {}, 0, h, a, w, ImplicitTerms::DampingOnly).has_value());
        REQUIRE(stepLinearlyImplicitEuler(c, {}, 0, h, b, w, ImplicitTerms::DampingOnlyMassPattern).has_value());
        if (viscosity == 0.0)
        {
            CHECK(a == b);
        }
        else
        {
            CHECK(relativeDifference(a.values, b.values) > 1e-9);
        }
    }
}

TEST_CASE("prepared ring step equals the one-call step bit for bit")
{
    const Model model = ringmodels::ringRig({.nodes = 16}, false);
    Workspace a(model);
    Workspace b(model);
    State sa = ringmodels::randomState(model, 3);
    State sb = sa;
    const auto h = sub0math::f64::Time::fromSi(1e-4);
    for (std::uint64_t t = 0; t < 20; ++t)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {}, t, h, sa, a).has_value());
        REQUIRE(prepareLinearlyImplicitStep(model, sb, b).has_value());
        REQUIRE(completeLinearlyImplicitStep(model, {}, t, h, sb, b).has_value());
    }
    CHECK(sa == sb);
}

TEST_CASE("free flight of a spinning deformed ring keeps momentum, and energy without damping")
{
    // No gravity, no edge viscosity: the midpoint step keeps a quadratic invariant of a linear system only, so energy
    // is checked to the scheme's order; linear momentum along travel and the generalized momentum of the uncoupled
    // directions are exact for any implicit step (no external force and the tree has no springs).
    RingParameters p{.nodes = 16, .viscosity = 0.0, .gravity = 0.0};
    const Model model = ringmodels::ringRig(p, false);
    Workspace ws(model);
    State s = ringmodels::randomState(model, 11, 0.003);
    const auto momentum = [&](const State& st) {
        // Linear momentum along x and z: the first two rows of M qd (the slides carry everything).
        REQUIRE(computeMassMatrix(model, st, ws).has_value());
        std::vector<double> px(2, 0.0);
        const std::size_t t = model.treeDofCount();
        for (std::size_t r = 0; r < 2; ++r)
        {
            for (std::size_t j = 0; j < t; ++j)
            {
                px[r] += ws.massMatrix[r * t + j] * st.qd()[j];
            }
            for (std::size_t k = 0; k < model.rings().size(); ++k)
            {
                const DeformableRing& ring = model.rings()[k];
                const RingWorkspace& rw = ws.rings[k];
                // Border rows map the ring rates onto the host's motion; the slides see the host's linear motion
                // through the host Jacobian's columns 0 and 1.
                double fx = 0.0;
                double fy = 0.0;
                double fw = 0.0;
                for (std::size_t i = 0; i < ring.nodeCount; ++i)
                {
                    const auto& b = rw.border[i];
                    const double w0 = st.qd()[ring.firstCoordinate + 2 * i];
                    const double w1 = st.qd()[ring.firstCoordinate + 2 * i + 1];
                    fw += b[0] * w0 + b[3] * w1;
                    fx += b[1] * w0 + b[4] * w1;
                    fy += b[2] * w0 + b[5] * w1;
                }
                // The slides are the chain's last two entries (root first is the end of hostChain).
                const std::size_t chain = ring.hostChain.size();
                const auto& j = rw.hostJacobian[chain - 1 - r];
                px[r] += j.w * fw + j.x * fx + j.y * fy;
            }
        }
        return px;
    };
    const auto p0 = momentum(s);
    const double e0 = totalEnergy(*computeEnergy(model, s, ws));
    const auto h = sub0math::f64::Time::fromSi(2e-5);
    for (std::uint64_t t = 0; t < 2000; ++t)
    {
        REQUIRE(stepLinearlyImplicitMidpoint(model, {}, t, h, s, ws).has_value());
    }
    const auto p1 = momentum(s);
    CHECK(std::abs(p1[0] - p0[0]) < 1e-6 * (1.0 + std::abs(p0[0])));
    CHECK(std::abs(p1[1] - p0[1]) < 1e-6 * (1.0 + std::abs(p0[1])));
    const double e1 = totalEnergy(*computeEnergy(model, s, ws));
    CHECK(std::abs(e1 - e0) < 1e-4 * e0);
    MESSAGE("free flight, midpoint 20 us, 40 ms: momentum drift " << std::abs(p1[0] - p0[0]) << ", "
                                                                  << std::abs(p1[1] - p0[1]) << " N s; energy "
                                                                  << (e1 - e0) / e0 << " relative");
}

TEST_CASE("energy ledger closes with ring damping accounted")
{
    RingParameters p{.nodes = 16, .viscosity = 2.0};
    const Model model = ringmodels::ringRig(p, false);
    Workspace ws(model);
    // RK4 integrates the dissipated energy as a state: the ledger closes to the method's error.
    {
        State s = ringmodels::randomState(model, 5);
        const double e0 = totalEnergy(*computeEnergy(model, s, ws));
        Rk4Stepper rk(model, ws);
        const auto h = sub0math::f64::Time::fromSi(1e-5);
        for (std::uint64_t t = 0; t < 2000; ++t)
        {
            REQUIRE(rk.step(model, {}, t, h, s, ws).has_value());
        }
        const auto e = *computeEnergy(model, s, ws);
        const double residual = totalEnergy(e) - e0 + e.dissipated.si() - e.externalWork.si();
        CHECK(e.dissipated.si() > 1e-3);
        CHECK(std::abs(residual) < 1e-8 * e.dissipated.si());
        MESSAGE("RK4 10 us, 20 ms: dissipated " << e.dissipated.si() << " J, ledger residual " << residual << " J");
    }
    // Linearly implicit midpoint: dissipation integrated at the start of each step, so the residual is first order.
    std::array<double, 3> residuals{};
    const std::array<double, 3> steps{4e-5, 2e-5, 1e-5};
    for (std::size_t k = 0; k < steps.size(); ++k)
    {
        State s = ringmodels::randomState(model, 5);
        const double e0 = totalEnergy(*computeEnergy(model, s, ws));
        const auto h = sub0math::f64::Time::fromSi(steps[k]);
        const auto count = static_cast<std::uint64_t>(std::llround(0.02 / steps[k]));
        for (std::uint64_t t = 0; t < count; ++t)
        {
            REQUIRE(stepLinearlyImplicitMidpoint(model, {}, t, h, s, ws).has_value());
        }
        const auto e = *computeEnergy(model, s, ws);
        residuals[k] = totalEnergy(e) - e0 + e.dissipated.si() - e.externalWork.si();
        CHECK(std::abs(residuals[k]) < 0.05 * e.dissipated.si());
    }
    const double order = std::log2(std::abs(residuals[1] / residuals[2]));
    MESSAGE("midpoint ledger residuals at 40, 20, 10 us: " << residuals[0] << ", " << residuals[1] << ", "
                                                           << residuals[2] << " J; order " << order);
    CHECK(order > 0.8);
}

TEST_CASE("free-ring modal frequencies agree with the closed form")
{
    // A host of enormous inertia on its own hinge: the ring's j >= 2 modes do not move it. In each node's polar frame
    // the uniform ring is circulant; mode j has one foundation-only branch, w^2 = kf / m, and one with the edges,
    // w^2 = (kf + ke (s^2 |1 + z|^2 + c^2 |1 - z|^2)) / m, z = exp(i 2 pi j / N), s, c = sin, cos(pi / N).
    constexpr std::size_t N = 16;
    RingParameters p{.nodes = N, .viscosity = 0.0, .gravity = 0.0};
    ModelBuilder builder;
    builder.setGravity(testmodels::gravityAlong(0.0, 0.0, 0.0));
    const BodyId hub = builder.addBody({.name = "hub",
                                        .mass = sub0math::kilograms(1e9),
                                        .centreOfMass = {testmodels::metres(0.0), testmodels::metres(0.0),
                                                         testmodels::metres(0.0)},
                                        .inertia = testmodels::diagonalInertia(1e9, 1e9, 1e9),
                                        .joint = {.type = JointType::Revolute, .axis = {0.0, 1.0, 0.0}}});
    builder.addRing(*uniformRing(ringmodels::ringSpec("ring", hub, p)));
    const Model model = *builder.build();
    Workspace ws(model);
    const DeformableRing& ring = model.rings()[0];
    const double m = ring.mass[0];
    const double kf = ring.foundation[0].spec().stiffness.si();
    const double ke = ring.edgeLaws[0].spec().stiffness.si();
    const double s = std::sin(std::numbers::pi / N);
    const double c = std::cos(std::numbers::pi / N);
    const auto edgeBranch = [&](std::size_t j) {
        const double phi = 2.0 * std::numbers::pi * static_cast<double>(j) / N;
        return std::sqrt((kf + ke * (s * s * (2.0 + 2.0 * std::cos(phi)) + c * c * (2.0 - 2.0 * std::cos(phi)))) / m);
    };
    // Excite the radial-dominant mode j = 3 shape (a mix of both branches) and measure by projection over time: simulate
    // a pure eigenvector instead. Mode j = 3: in the polar frame (radial a, tangential b) the edge branch's vector is
    // u = (s (1 + z), c (z - 1)) (conjugated); a real standing mode uses cos/sin of j theta.
    const std::size_t j = 3;
    const double phi = 2.0 * std::numbers::pi * static_cast<double>(j) / N;
    // u in the polar frame for node offset 0: (s (1 + e^{-i phi}), c (e^{-i phi} - 1)) times e^{i j theta}.
    State st = State::atRest(model);
    const double amplitude = 1e-5;
    for (std::size_t i = 0; i < N; ++i)
    {
        const double theta = 2.0 * std::numbers::pi * static_cast<double>(i) / N;
        // Real part of (s (1 + e^{-i phi}) e^{i j theta}, c (e^{-i phi} - 1) e^{i j theta}).
        const double radial = s * (std::cos(j * theta) + std::cos(j * theta - phi));
        const double tangential = c * (std::cos(j * theta - phi) - std::cos(j * theta));
        // Polar frame at node i: radial (sin t, -cos t), tangential (cos t, sin t), in (a_0, a_1) = (x, z).
        st.q()[1 + 2 * i] = amplitude * (radial * std::sin(theta) + tangential * std::cos(theta));
        st.q()[2 + 2 * i] = amplitude * (-radial * std::cos(theta) + tangential * std::sin(theta));
    }
    const double omega = edgeBranch(j);
    const double h = 0.002 / omega; // w h = 0.002: the midpoint step's phase error is (w h)^2 / 12
    const auto step = sub0math::f64::Time::fromSi(h);
    // Track node 0's first coordinate through its zero crossings.
    double previous = st.q()[1];
    double firstCrossing = -1.0;
    double lastCrossing = -1.0;
    int crossings = 0;
    for (std::uint64_t t = 0; t < 40000; ++t)
    {
        REQUIRE(stepLinearlyImplicitMidpoint(model, {}, t, step, st, ws).has_value());
        const double now = st.q()[1];
        if ((previous < 0.0) != (now < 0.0))
        {
            const double at = (static_cast<double>(t) + previous / (previous - now)) * h;
            if (firstCrossing < 0.0)
            {
                firstCrossing = at;
            }
            lastCrossing = at;
            ++crossings;
        }
        previous = now;
    }
    REQUIRE(crossings > 10);
    const double measured = std::numbers::pi * (crossings - 1) / (lastCrossing - firstCrossing);
    MESSAGE("mode j = 3, edge branch: closed form " << omega << " rad/s, measured " << measured);
    CHECK(std::abs(measured - omega) < 1e-5 * omega);
    // The foundation-only branch.
    CHECK(std::sqrt(kf / m) < omega);
}

TEST_CASE("ring run digest (compare across compilers)")
{
    // Two 64-node rings on the tree, 2000 steps of each implicit scheme from a deformed, moving state. The digest of
    // every state bit and the energy flows is printed: the determinism gate compares it across GCC, GCC
    // -march=x86-64-v3 and Clang -march=x86-64-v3 builds (docs/PERFORMANCE.md).
    const Model model = ringmodels::ringRig({.nodes = 64}, false);
    Workspace ws(model);
    std::uint64_t digest = 1469598103934665603ULL;
    const auto mix = [&](double v) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &v, sizeof bits);
        for (int b = 0; b < 8; ++b)
        {
            digest = (digest ^ ((bits >> (8 * b)) & 0xFFU)) * 1099511628211ULL;
        }
    };
    const auto h = sub0math::f64::Time::fromSi(31.25e-6);
    for (int scheme = 0; scheme < 3; ++scheme)
    {
        State s = ringmodels::randomState(model, 21);
        for (std::uint64_t t = 0; t < 2000; ++t)
        {
            const auto r = scheme == 0   ? stepLinearlyImplicitEuler(model, {}, t, h, s, ws)
                           : scheme == 1 ? stepLinearlyImplicitMidpoint(model, {}, t, h, s, ws)
                                         : stepLinearlyImplicitEuler(model, {}, t, h, s, ws, ImplicitTerms::DampingOnly);
            REQUIRE(r.has_value());
        }
        for (const double v : s.values)
        {
            mix(v);
        }
        mix(s.dissipated.si());
    }
    char text[17];
    std::snprintf(text, sizeof text, "%016llx", static_cast<unsigned long long>(digest));
    MESSAGE("ring run digest " << text);
    // Pinned: identical on GCC, GCC -march=x86-64-v3 and Clang -march=x86-64-v3 (Windows), and expected on every
    // platform: the random state comes from the raw mt19937_64 output (ringmodels::uniformSigned, not
    // std::uniform_real_distribution), the revolute joints use sub0math::det::sincos (fdlibm, no C library), and
    // everything else is IEEE + - * / and sqrt under -ffp-contract=off. Before 2026-10-27 this was b96b9f5285379502
    // on Windows with libm sin/cos, abd8d2ab7e62fbc7 on Linux and f67fc6a4aa623f6d on macOS. A mismatch now is a bug.
    CHECK(std::strcmp(text, "f4b3d8f7baf6ed29") == 0);
}
