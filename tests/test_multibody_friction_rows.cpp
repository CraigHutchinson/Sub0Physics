// Friction rows on a deformable ring (docs/WHEEL.md section 3.3; step M3b): the box solve against enumeration, and the gate
// for exact Coulomb against the reference's elliptic cone on one fixture, a loaded wheel on an incline with its hub
// locked against rotation (BikeLab's casing values, 64 nodes, 31.25 us). Below the friction limit the exact rows hold with
// zero slip and the hub comes to rest; the reference's regularised cone creeps at a steady rate, measured here. Above the
// limit each contact slides at mu N and the hub accelerates at g (sin - mu cos). The friction work is never positive.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <span>
#include <vector>

#include <doctest.h>

#include <sub0math/units.hpp>

#include "sub0physics/geometry/convex_terrain.hpp"
#include "sub0physics/geometry/polyline_profile.hpp"
#include "sub0physics/multibody.hpp"
#include "test_ring_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using ringmodels::RingParameters;
namespace f64 = sub0math::f64;

namespace
{

// A locked hub on x and z slides (no spin coordinate) with a 64-node casing.
Model lockedWheel(double load, double foundationViscosity)
{
    using sub0math::kilograms;
    using testmodels::diagonalInertia;
    ModelBuilder builder;
    builder.setGravity(testmodels::gravityAlong(0.0, 0.0, -9.80665));
    BodyDescription travel{};
    travel.name = "travel";
    travel.joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}};
    const BodyId t = builder.addBody(travel);
    const BodyId hub = builder.addBody({.name = "hub",
                                        .mass = kilograms(load),
                                        .inertia = diagonalInertia(0.04, 0.04, 0.04),
                                        .joint = {.type = JointType::Prismatic, .parent = t, .axis = {0.0, 0.0, 1.0}}});
    RingParameters p;
    p.nodes = 64;
    UniformRingSpec spec = ringmodels::ringSpec("casing", hub, p);
    spec.foundationViscosity = foundationViscosity;
    const auto ring = uniformRing(spec);
    REQUIRE(ring.has_value());
    builder.addRing(*ring);
    auto model = builder.build();
    REQUIRE(model.has_value());
    return std::move(*model);
}

enum class Friction : std::uint8_t
{
    ReferenceCone,
    ExactCoulomb
};

// One reference-law step of the fixture, written as a consumer writes it (the bike's reference stepper).
struct Stepper
{
    Stepper(const Model& m, const ConvexTerrain& t, double mu_, Friction friction_)
        : model(m), terrain(t), mu(mu_), friction(friction_), n(m.dofCount()), ws(m), rest(m), solver(m, 2 * Capacity),
          points(64), found(Capacity + 1), inverseWeights(64), jn(Capacity * n), jt(Capacity * n), both(2 * Capacity * n),
          kinds(2 * Capacity), mus(2 * Capacity, mu_), reference(2 * Capacity), regularisation(2 * Capacity),
          force(2 * Capacity), supports(2 * Capacity), rates(Capacity), slips(Capacity), diag(Capacity), gaps(Capacity), free(n),
          generalized(n), tangential(n), tau(n)
    {
        for (std::size_t r = 0; r < kinds.size(); ++r)
        {
            kinds[r] = r % 2 == 0 ? SoftRowKind::ConeNormal : SoftRowKind::ConeTangent;
        }
        const State s = State::atRest(model);
        REQUIRE(prepareLinearlyImplicitStep(model, s, rest));
        REQUIRE(factorMassMatrix(model, rest));
        REQUIRE(ringNodeInverseWeights(model, rest, 0, inverseWeights));
    }

    static constexpr std::size_t Capacity = 64;

    [[nodiscard]] std::size_t contacts() const noexcept
    {
        return count;
    }
    [[nodiscard]] double normalForce(std::size_t c) const noexcept
    {
        return normal[c];
    }
    [[nodiscard]] double tangentForce(std::size_t c) const noexcept
    {
        return tangent[c];
    }

    // Advances one step; `worstStick`, `stuck`, `work` report the end of the step.
    void step(State& s, double h)
    {
        const f64::Time dt = f64::Time::fromSi(h);
        REQUIRE(prepareLinearlyImplicitStep(model, s, ws, ImplicitTerms::DampingOnlyMassPattern));
        const PlanarTransform& host = ws.planar.worldToBody[model.rings()[0].host];
        const auto xz = [&](std::size_t i) {
            std::array<double, 3> v{};
            const double lx = ws.rings[0].positionX[i];
            const double ly = ws.rings[0].positionY[i];
            v[model.planar().axisI] = host.x + host.c * lx - host.s * ly;
            v[model.planar().axisJ] = host.y + host.s * lx + host.c * ly;
            return std::array<double, 2>{v[0], v[2]};
        };
        for (std::size_t i = 0; i < 64; ++i)
        {
            const auto p = xz(i);
            points[i] = {f64::Length::fromSi(p[0]), f64::Length::fromSi(p[1])};
        }
        const double skin = model.rings()[0].skinRadius;
        const auto set = terrain.ringContacts(points, f64::Length::fromSi(skin), f64::Length{}, found);
        REQUIRE(!set.overflowed);
        const SoftConstraintParameters law{.timeConstant = 0.0005, .dampingRatio = 1.0, .impedance = {0.99, 0.99, 0.001, 0.5, 2.0}};
        count = 0;
        for (std::size_t m = 0; m < set.count; ++m)
        {
            const TerrainContact& t = found[m];
            if (!(t.gap.si() < 0.0))
            {
                continue;
            }
            const std::size_t e1 = (t.edge + 1) % 64;
            const std::array<double, 2> at{(t.witness.point.x.si() + t.capsulePoint.x.si()) * 0.5,
                                           (t.witness.point.y.si() + t.capsulePoint.y.si()) * 0.5};
            const auto a = xz(t.edge);
            const auto b = xz(e1);
            double w0 = 1.0 / std::hypot(at[0] - a[0], at[1] - a[1]);
            double w1 = 1.0 / std::hypot(at[0] - b[0], at[1] - b[1]);
            const double sum = 1.0 / (w0 + w1);
            w0 *= sum;
            w1 *= sum;
            const std::array<double, 2> nrm{t.witness.normal.x, t.witness.normal.y};
            RingContactPoint p{.ring = 0, .edge = t.edge, .weight0 = w0, .weight1 = w1,
                               .point = {at[0], 0.0, at[1]}, .normal = {nrm[0], 0.0, nrm[1]}};
            const std::span<double> rowN(jn.data() + count * n, n);
            REQUIRE(ringContactJacobian(model, ws, p, rowN));
            p.normal = {nrm[1], 0.0, -nrm[0]};
            const std::span<double> rowT(jt.data() + count * n, n);
            REQUIRE(ringContactJacobian(model, ws, p, rowT));
            double rate = 0.0;
            double slip = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                rate += rowN[i] * s.qd()[i];
                slip += rowT[i] * s.qd()[i];
            }
            rates[count] = rate;
            slips[count] = slip;
            diag[count] = inverseWeights[t.edge] * w0 + inverseWeights[e1] * w1;
            gaps[count] = t.gap.si();
            supports[count] = {.ring = 0, .node0 = t.edge, .node1 = e1};
            ++count;
        }
        const std::size_t k = count;
        for (std::size_t i = 0; i < n; ++i)
        {
            free[i] = ws.elementForce[i] - ws.bias[i];
        }
        normal.assign(k, 0.0);
        tangent.assign(k, 0.0);
        if (friction == Friction::ReferenceCone)
        {
            REQUIRE(factorMassMatrix(model, ws));
            REQUIRE(solveMassMatrix(model, ws, free));
            std::vector<SoftRowSupport> pairs(2 * k);
            for (std::size_t c = 0; c < k; ++c)
            {
                const SoftConstraintRow row = softConstraintRow(law, gaps[c], rates[c], diag[c], h);
                const SoftConstraintRow along = softConstraintTangentRow(row, slips[c]);
                reference[2 * c] = row.referenceAcceleration;
                regularisation[2 * c] = row.regularisation;
                reference[2 * c + 1] = along.referenceAcceleration;
                regularisation[2 * c + 1] = along.regularisation;
                std::copy_n(jn.data() + c * n, n, both.data() + 2 * c * n);
                std::copy_n(jt.data() + c * n, n, both.data() + (2 * c + 1) * n);
                pairs[2 * c] = supports[c];
                pairs[2 * c + 1] = supports[c];
            }
            REQUIRE(solver.solve(model, ws, std::span<const double>(both.data(), 2 * k * n),
                                 std::span<const double>(reference.data(), 2 * k),
                                 std::span<const double>(regularisation.data(), 2 * k), free,
                                 std::span<const SoftRowKind>(kinds.data(), 2 * k), std::span<double>(force.data(), 2 * k),
                                 generalized, pairs, DelassusForm::MassSolvePerRow,
                                 std::span<const double>(mus.data(), 2 * k)));
            for (std::size_t c = 0; c < k; ++c)
            {
                normal[c] = force[2 * c];
                tangent[c] = force[2 * c + 1];
            }
            tau = generalized;
        }
        else
        {
            // The normal rows with the step's matrix, then the Coulomb rows bounded by the normal forces.
            REQUIRE(factorPreparedStepMatrix(model, s, dt, ws));
            REQUIRE(solvePreparedStepMatrix(model, ws, free));
            for (std::size_t c = 0; c < k; ++c)
            {
                const SoftConstraintRow row = softConstraintRow(law, gaps[c], rates[c], diag[c], h);
                reference[c] = row.referenceAcceleration;
                regularisation[c] = row.regularisation;
            }
            const auto rowsN = std::span<const SoftRowSupport>(supports.data(), k);
            REQUIRE(solver.solve(model, ws, std::span<const double>(jn.data(), k * n), std::span<const double>(reference.data(), k),
                                 std::span<const double>(regularisation.data(), k), free, {},
                                 std::span<double>(force.data(), k), generalized, rowsN, DelassusForm::MassSolvePerRow, {},
                                 SoftRowSystem::PreparedStep));
            std::vector<double> bounds(k);
            for (std::size_t c = 0; c < k; ++c)
            {
                normal[c] = force[c];
                bounds[c] = mu * force[c];
                const SoftConstraintRow along = exactVelocityRow(slips[c], diag[c], h);
                reference[c] = along.referenceAcceleration;
                regularisation[c] = along.regularisation;
            }
            tau = generalized;
            REQUIRE(solvePreparedStepMatrix(model, ws, generalized));
            for (std::size_t i = 0; i < n; ++i)
            {
                free[i] += generalized[i];
            }
            REQUIRE(solver.solveBounded(model, ws, std::span<const double>(jt.data(), k * n),
                                        std::span<const double>(reference.data(), k),
                                        std::span<const double>(regularisation.data(), k), free, bounds,
                                        std::span<double>(force.data(), k), tangential, rowsN));
            for (std::size_t c = 0; c < k; ++c)
            {
                tangent[c] = force[c];
            }
            for (std::size_t i = 0; i < n; ++i)
            {
                tau[i] += tangential[i];
            }
        }
        REQUIRE(completeLinearlyImplicitStep(model, Forcing{.tau = tau}, s.tick, dt, s, ws));
        // The rows' own slip at the end of the step, and the friction work with it.
        stuck = 0;
        worstStick = 0.0;
        work = 0.0;
        workStart = 0.0;
        forceScale = 0.0;
        for (std::size_t c = 0; c < k; ++c)
        {
            double end = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                end += jt[c * n + i] * s.qd()[i];
            }
            work += h * tangent[c] * end;
            workStart += h * tangent[c] * slips[c];
            forceScale += std::abs(tangent[c]);
            if (normal[c] > 0.0 && std::abs(tangent[c]) < mu * normal[c] * (1.0 - 1e-9))
            {
                ++stuck;
                worstStick = std::max(worstStick, std::abs(end));
            }
        }
    }

    const Model& model;
    const ConvexTerrain& terrain;
    double mu;
    Friction friction;
    std::size_t n;
    Workspace ws;
    Workspace rest;
    SoftConstraintSolver solver;
    std::vector<sub0math::Point2<f64::Length>> points;
    std::vector<TerrainContact> found;
    std::vector<double> inverseWeights;
    std::vector<double> jn, jt, both;
    std::vector<SoftRowKind> kinds;
    std::vector<double> mus, reference, regularisation, force;
    std::vector<SoftRowSupport> supports;
    std::vector<double> rates, slips, diag, gaps, free, generalized, tangential, tau, normal, tangent;
    std::size_t count = 0;
    std::size_t stuck = 0;
    double worstStick = 0.0, work = 0.0, workStart = 0.0, forceScale = 0.0;
};

struct InclineResult
{
    double creep = 0.0;          // mean hub speed down the slope over the last window [m/s]
    double drift = 0.0;          // hub travel down the slope over the last window [m]
    double finalSpeed = 0.0;     // hub speed at the end [m/s]
    double acceleration = 0.0;   // mean hub acceleration down the slope over the last window [m/s^2]
    double worstStick = 0.0;     // largest end-of-step slip of a sticking row [m/s]
    double largestWork = -1e300; // largest step's friction work, end-of-step rates [J]
    double largestWorkStart = -1e300; // the same with start-of-step rates (the reference's power) [J]
    double dissipated = 0.0;     // total friction work [J]
    double coneRatio = 0.0;      // largest |f_t| / (mu f_n) [-]
    double slidingShare = 0.0;   // contacts on the cone over all pushing contacts, last window
    std::size_t contacts = 0;
    double clearance = 0.0;      // mean hub distance from the plane over the last window [m]
    std::uint64_t digest = 1469598103934665603ULL;
};

InclineResult runIncline(double slopeTangent, Friction friction, double duration, double windowLength)
{
    const double mu = 0.6;
    const double h = 3.125e-5;
    const Model model = lockedWheel(40.0, 2000.0);
    const double extent = 4.0;
    auto profile = PolylineProfile::create({{f64::Length::fromSi(-extent), f64::Length::fromSi(extent * slopeTangent)},
                                            {f64::Length::fromSi(extent), f64::Length::fromSi(-extent * slopeTangent)}});
    REQUIRE(profile.has_value());
    const ConvexTerrain terrain(*profile);
    Stepper stepper(model, terrain, mu, friction);
    State s = State::atRest(model);
    const double c = 1.0 / std::sqrt(1.0 + slopeTangent * slopeTangent);
    const double sn = slopeTangent * c;
    // Hub a little above touching, the casing at rest.
    s.q()[1] = (0.36 + 1e-4) / c;
    const std::size_t steps = static_cast<std::size_t>(duration / h);
    const std::size_t window = static_cast<std::size_t>(windowLength / h);
    double startAlong = 0.0;
    double startRate = 0.0;
    InclineResult r;
    std::size_t pushing = 0;
    std::size_t sliding = 0;
    const auto along = [&](const State& x) { return x.q()[0] * c - x.q()[1] * sn; };
    const auto alongRate = [&](const State& x) { return x.qd()[0] * c - x.qd()[1] * sn; };
    double clearanceSum = 0.0;
    for (std::size_t k = 0; k < steps; ++k)
    {
        if (k == steps - window)
        {
            startAlong = along(s);
            startRate = alongRate(s);
        }
        stepper.step(s, h);
        r.dissipated += stepper.work;
        r.largestWork = std::max(r.largestWork, stepper.work);
        r.largestWorkStart = std::max(r.largestWorkStart, stepper.workStart);
        r.worstStick = k >= steps / 2 ? std::max(r.worstStick, stepper.worstStick) : r.worstStick;
        for (std::size_t x = 0; x < stepper.contacts(); ++x)
        {
            if (stepper.normalForce(x) > 0.0)
            {
                r.coneRatio = std::max(r.coneRatio, std::abs(stepper.tangentForce(x)) / (mu * stepper.normalForce(x)));
                if (k >= steps - window)
                {
                    ++pushing;
                    // (Sampled per pushing contact, a contact-weighted mean.)
                    clearanceSum += s.q()[0] * sn + s.q()[1] * c;
                    if (std::abs(stepper.tangentForce(x)) >= mu * stepper.normalForce(x) * (1.0 - 1e-9))
                    {
                        ++sliding;
                    }
                }
            }
        }
    }
    r.contacts = stepper.contacts();
    r.clearance = pushing == 0 ? 0.0 : clearanceSum / static_cast<double>(pushing);
    r.drift = along(s) - startAlong;
    r.creep = r.drift / (static_cast<double>(window) * h);
    r.finalSpeed = alongRate(s);
    r.acceleration = (alongRate(s) - startRate) / (static_cast<double>(window) * h);
    r.slidingShare = pushing == 0 ? 0.0 : static_cast<double>(sliding) / static_cast<double>(pushing);
    for (const double v : s.values)
    {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &v, sizeof bits);
        for (int b = 0; b < 8; ++b)
        {
            r.digest = (r.digest ^ ((bits >> (8 * b)) & 0xFFU)) * 1099511628211ULL;
        }
    }
    return r;
}

} // namespace

TEST_CASE("bounded rows: the box solve against enumeration of every row's state")
{
    const Model model = ringmodels::ringRig({.nodes = 8}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    std::mt19937_64 rng(23);
    const auto u = [](std::mt19937_64& g) { return ringmodels::uniformSigned(g); };
    const State s = ringmodels::randomState(model, 3);
    REQUIRE(prepareLinearlyImplicitStep(model, s, ws, ImplicitTerms::DampingOnlyMassPattern));
    REQUIRE(factorMassMatrix(model, ws));
    REQUIRE(factorPreparedStepMatrix(model, s, f64::Time::fromSi(1e-4), ws));
    SoftConstraintSolver solver(model, 6);
    std::size_t atBound = 0;
    std::size_t inside = 0;
    for (int trial = 0; trial < 60; ++trial)
    {
        const std::size_t k = 1 + static_cast<std::size_t>(trial % 5);
        std::vector<double> jac(k * n, 0.0);
        std::vector<double> aref(k);
        std::vector<double> reg(k);
        std::vector<double> bounds(k);
        std::vector<SoftRowSupport> supports(k);
        for (std::size_t r = 0; r < k; ++r)
        {
            RingContactPoint p;
            p.edge = (static_cast<std::size_t>(trial) + 2 * r) % 8;
            p.weight0 = 0.5 + 0.4 * u(rng);
            p.weight1 = 1.0 - p.weight0;
            p.point = {0.01 * u(rng), 0.0, 0.01 * u(rng)};
            const double angle = 3.0 * u(rng);
            p.normal = {std::cos(angle), 0.0, std::sin(angle)};
            REQUIRE(ringContactJacobian(model, ws, p, std::span<double>(jac.data() + r * n, n)));
            supports[r] = {.ring = 0, .node0 = p.edge, .node1 = (p.edge + 1) % 8};
            aref[r] = 50.0 * u(rng);
            reg[r] = 1e-3 * (1.5 + u(rng));
            bounds[r] = trial % 7 == 3 && r == 0 ? 0.0 : 0.3 + 0.3 * u(rng);
        }
        std::vector<double> a0(n);
        for (double& v : a0)
        {
            v = 10.0 * u(rng);
        }
        for (const SoftRowSystem system : {SoftRowSystem::MassMatrix, SoftRowSystem::PreparedStep})
        {
            std::vector<double> f(k);
            std::vector<double> g(n);
            REQUIRE(solver.solveBounded(model, ws, jac, aref, reg, a0, bounds, f, g, {}, DelassusForm::MassSolvePerRow, system));
            // With supports: bit for bit. Condensed: to rounding.
            std::vector<double> fs(k);
            REQUIRE(solver.solveBounded(model, ws, jac, aref, reg, a0, bounds, fs, g, supports, DelassusForm::MassSolvePerRow, system));
            CHECK(fs == f);
            std::vector<double> fc(k);
            REQUIRE(solver.solveBounded(model, ws, jac, aref, reg, a0, bounds, fc, g, supports, DelassusForm::Condensed, system));
            for (std::size_t r = 0; r < k; ++r)
            {
                CHECK(std::abs(fc[r] - f[r]) <= 1e-9 * std::max(1.0, std::abs(f[r])));
            }
            // The matrix and vector of the QP: A = J S^-1 J^T + R, w = A f + q.
            const auto apply = [&](std::span<double> x) {
                REQUIRE((system == SoftRowSystem::MassMatrix ? solveMassMatrix(model, ws, x) : solvePreparedStepMatrix(model, ws, x)));
            };
            std::vector<double> y(n, 0.0);
            for (std::size_t r = 0; r < k; ++r)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    y[i] += jac[r * n + i] * f[r];
                }
            }
            apply(y);
            for (std::size_t r = 0; r < k; ++r)
            {
                double w = -aref[r] + reg[r] * f[r];
                for (std::size_t i = 0; i < n; ++i)
                {
                    w += jac[r * n + i] * (a0[i] + y[i]);
                }
                // Complementarity of the box: inside, the row's residual vanishes; at a bound it pushes outwards.
                CHECK(std::abs(f[r]) <= bounds[r]);
                if (bounds[r] == 0.0)
                {
                    CHECK(f[r] == 0.0);
                }
                else if (std::abs(f[r]) == bounds[r])
                {
                    ++atBound;
                    CHECK(f[r] * w <= 1e-9);
                }
                else
                {
                    ++inside;
                    CHECK(std::abs(w) <= 1e-8);
                }
            }
        }
    }
    MESSAGE("rows at a bound " << atBound << ", inside " << inside);
    CHECK(atBound > 0U);
    CHECK(inside > 0U);
    // Contract checks.
    std::vector<double> jac(n, 0.0);
    const std::vector<double> one{1.0};
    const std::vector<double> negative{-1.0};
    std::vector<double> a0(n, 0.0);
    std::vector<double> f(1);
    std::vector<double> g(n);
    CHECK_FALSE(solver.solveBounded(model, ws, jac, one, one, a0, negative, f, g).has_value());
    Workspace unfactored(model);
    CHECK_FALSE(solver.solveBounded(model, unfactored, jac, one, one, a0, one, f, g).has_value());
}

TEST_CASE("exact Coulomb on an incline below the friction limit: zero creep; the reference's cone creeps")
{
    // tan = 0.25 against mu = 0.6: well inside the cone. Foundation viscosity settles the hub on its casing (a description
    // value of this fixture, on both models), so the last window sees the steady state.
    const InclineResult exact = runIncline(0.25, Friction::ExactCoulomb, 4.0, 0.5);
    const InclineResult cone = runIncline(0.25, Friction::ReferenceCone, 4.0, 0.5);
    MESSAGE("exact rows: drift " << exact.drift << " m, mean speed " << exact.creep << " m/s, final speed " << exact.finalSpeed
                                 << " m/s, worst stick slip " << exact.worstStick << " m/s, largest step friction work "
                                 << exact.largestWork << " J, total " << exact.dissipated << " J, cone ratio " << exact.coneRatio
                                 << ", contacts " << exact.contacts);
    MESSAGE("reference cone: drift " << cone.drift << " m, creep " << cone.creep << " m/s, final speed " << cone.finalSpeed
                                     << " m/s, largest step friction work " << cone.largestWork << " J (start-of-step rates "
                                     << cone.largestWorkStart << " J), total " << cone.dissipated << " J, cone ratio " << cone.coneRatio
                                     << ", contacts " << cone.contacts);
    // Exact stick: every sticking row's slip at the end of the step is the row's compliance to rounding.
    CHECK(exact.worstStick <= 1e-11);
    CHECK(std::abs(exact.creep) <= 1e-9);
    CHECK(std::abs(exact.drift) <= 1e-10);
    CHECK(std::abs(exact.finalSpeed) <= 1e-8);
    CHECK(exact.coneRatio <= 1.0 + 1e-12);
    // The reference creeps at about R f / B: steady, down the slope, millimetres per second here.
    CHECK(cone.creep > 1e-4);
    CHECK(std::abs(cone.finalSpeed - cone.creep) <= 0.05 * cone.creep);
    // The work of the rows, end-of-step rates: never positive beyond rounding.
    CHECK(exact.largestWork <= 1e-14);
    CHECK(exact.dissipated < 0.0);
}

TEST_CASE("pivot rules: the same forces, bit for bit, on friction cones")
{
    const Model model = ringmodels::ringRig({.nodes = 8}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    std::mt19937_64 rng(31);
    const auto u = [](std::mt19937_64& g) { return ringmodels::uniformSigned(g); };
    const State s = ringmodels::randomState(model, 3);
    REQUIRE(prepareLinearlyImplicitStep(model, s, ws));
    REQUIRE(factorMassMatrix(model, ws));
    SoftConstraintSolver least(model, 16);
    SoftConstraintSolver block(model, 16);
    block.setPivotRule(SoftPivotRule::Block);
    std::size_t fewer = 0;
    for (int trial = 0; trial < 100; ++trial)
    {
        const std::size_t k = 2 * (1 + static_cast<std::size_t>(trial % 8));
        std::vector<double> jac(k * n, 0.0);
        std::vector<double> aref(k);
        std::vector<double> reg(k);
        std::vector<SoftRowKind> kinds(k);
        const std::vector<double> mu(k, 0.6);
        for (std::size_t r = 0; r < k; ++r)
        {
            RingContactPoint p;
            p.edge = (static_cast<std::size_t>(trial) + r / 2) % 8;
            p.point = {0.01 * u(rng), 0.0, 0.01 * u(rng)};
            const double angle = 3.0 * u(rng);
            p.normal = {std::cos(angle), 0.0, std::sin(angle)};
            REQUIRE(ringContactJacobian(model, ws, p, std::span<double>(jac.data() + r * n, n)));
            kinds[r] = r % 2 == 0 ? SoftRowKind::ConeNormal : SoftRowKind::ConeTangent;
            aref[r] = 50.0 * u(rng);
            reg[r] = 0.01 * (1.5 + u(rng));
        }
        std::vector<double> a0(n, 0.0);
        std::vector<double> f(k);
        std::vector<double> g(k);
        std::vector<double> gen(n);
        const auto a = least.solve(model, ws, jac, aref, reg, a0, kinds, f, gen, {}, DelassusForm::MassSolvePerRow, mu);
        const auto b = block.solve(model, ws, jac, aref, reg, a0, kinds, g, gen, {}, DelassusForm::MassSolvePerRow, mu);
        REQUIRE(a.has_value());
        REQUIRE(b.has_value());
        CHECK(f == g);
        fewer += b->pivots <= a->pivots ? 1U : 0U;
    }
    MESSAGE("block rule at or below least-index pivots on " << fewer << " of 100 sets");
    CHECK(fewer > 50U);
}

TEST_CASE("pivot rules on bounded rows: nearly dependent rows, few block pivots, the least-index forces")
{
    // A casing edge on two terrain pieces gives two tangent rows on one support with nearly one direction, and neighbouring
    // edges share nodes: the box's matrix is positive definite only by exactVelocityRow's rounding compliance. Least index
    // then needed 288 pivots for 12 rows in a private consumer's flow ride (its limit is 16 k + 64 = 256); the block rule, 12.
    const Model model = ringmodels::ringRig({.nodes = 8}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    std::mt19937_64 rng(41);
    const auto u = [](std::mt19937_64& g) { return ringmodels::uniformSigned(g); };
    const State s = ringmodels::randomState(model, 5);
    REQUIRE(prepareLinearlyImplicitStep(model, s, ws, ImplicitTerms::DampingOnlyMassPattern));
    REQUIRE(factorPreparedStepMatrix(model, s, f64::Time::fromSi(1.5625e-4), ws));
    SoftConstraintSolver least(model, 12);
    SoftConstraintSolver block(model, 12);
    block.setPivotRule(SoftPivotRule::Block);
    std::uint32_t mostLeast = 0;
    std::uint32_t mostBlock = 0;
    std::size_t leastFailed = 0;
    for (int trial = 0; trial < 200; ++trial)
    {
        const std::size_t k = 12;
        std::vector<double> jac(k * n, 0.0);
        std::vector<double> aref(k);
        std::vector<double> reg(k);
        std::vector<double> bounds(k);
        const double along = 0.2 * u(rng); // the patch's tangent, nearly the same for every row
        for (std::size_t r = 0; r < k; ++r)
        {
            RingContactPoint p;
            p.edge = (static_cast<std::size_t>(trial) + r / 2) % 8;
            p.weight0 = 0.5 + 0.45 * u(rng);
            p.weight1 = 1.0 - p.weight0;
            p.point = {0.01 * u(rng), 0.0, 0.01 * u(rng)};
            const double angle = along + 0.02 * u(rng);
            p.normal = {std::cos(angle), 0.0, std::sin(angle)};
            REQUIRE(ringContactJacobian(model, ws, p, std::span<double>(jac.data() + r * n, n)));
            const SoftConstraintRow row = exactVelocityRow(2.0 * u(rng), 1.0, 1.5625e-4);
            aref[r] = row.referenceAcceleration;
            reg[r] = row.regularisation;
            bounds[r] = 20.0 + 60.0 * (1.0 + u(rng));
        }
        std::vector<double> a0(n);
        for (double& v : a0)
        {
            v = 200.0 * u(rng);
        }
        std::vector<double> f(k);
        std::vector<double> g(k);
        std::vector<double> gen(n);
        const auto b = block.solveBounded(model, ws, jac, aref, reg, a0, bounds, g, gen);
        REQUIRE_MESSAGE(b.has_value(), "trial " << trial);
        mostBlock = std::max(mostBlock, b->pivots);
        for (std::size_t r = 0; r < k; ++r)
        {
            CHECK(std::abs(g[r]) <= bounds[r]);
        }
        const auto a = least.solveBounded(model, ws, jac, aref, reg, a0, bounds, f, gen);
        if (!a)
        {
            ++leastFailed;
            continue;
        }
        mostLeast = std::max(mostLeast, a->pivots);
        // One solution: the same forces to the matrix's conditioning (the rules' last factors match when their sets do).
        for (std::size_t r = 0; r < k; ++r)
        {
            CHECK(std::abs(g[r] - f[r]) <= 1e-6 * bounds[r]);
        }
    }
    MESSAGE("most pivots: least index " << mostLeast << " (" << leastFailed << " of 200 over its limit), block " << mostBlock);
    CHECK(mostBlock <= 4U * 12U);
    CHECK(mostBlock < mostLeast);
}

TEST_CASE("exact Coulomb above the friction limit: slides at mu N")
{
    // tan = 0.8 against mu = 0.6: hub acceleration g (sin - mu cos).
    const double tangent = 0.8;
    const double c = 1.0 / std::sqrt(1.0 + tangent * tangent);
    const double expected = 9.80665 * (tangent * c - 0.6 * c);
    const InclineResult exact = runIncline(tangent, Friction::ExactCoulomb, 1.5, 1.0);
    const InclineResult cone = runIncline(tangent, Friction::ReferenceCone, 1.5, 1.0);
    MESSAGE("expected " << expected << " m/s^2; exact rows " << exact.acceleration << " m/s^2, on the cone " << exact.slidingShare
                        << " of pushing contacts, largest step friction work " << exact.largestWork << " J, total "
                        << exact.dissipated << " J, cone ratio " << exact.coneRatio << ", hub clearance " << exact.clearance << " m");
    MESSAGE("reference cone " << cone.acceleration << " m/s^2, on the cone " << cone.slidingShare << " of pushing contacts, largest step friction work "
                              << cone.largestWork << " J, total " << cone.dissipated << " J, cone ratio " << cone.coneRatio
                              << ", hub clearance " << cone.clearance << " m");
    CHECK(std::abs(exact.acceleration - expected) <= 0.05 * expected);
    CHECK(exact.slidingShare >= 0.99);
    CHECK(exact.coneRatio <= 1.0 + 1e-12);
    CHECK(exact.largestWork <= 1e-14);
    CHECK(exact.dissipated < 0.0);
}
