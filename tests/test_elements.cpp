#include <algorithm>
#include <cmath>
#include <doctest.h>
#include <ostream> // doctest prints strings (MSVC needs it explicitly)
#include <sub0math/units.hpp>

#include "sub0physics/elements.hpp"

using namespace sub0physics;
using namespace sub0math;
using namespace sub0math::f64::literals;

namespace
{

AxialSpringDamper makeAxial()
{
    return *AxialSpringDamper::create({
        .preload = newtons(50.0),
        .stiffness = newtonsPerMetre(9000.0),
        .compressionDamping = newtonSecondsPerMetre(350.0),
        .reboundDamping = newtonSecondsPerMetre(650.0),
    });
}

UnilateralSpringDamper makeContact()
{
    return *UnilateralSpringDamper::create({.stiffness = newtonsPerMetre(90000.0), .damping = newtonSecondsPerMetre(400.0)});
}

} // namespace

TEST_CASE("Axial element uses direction-dependent damping")
{
    const AxialSpringDamper element = makeAxial();
    const AxialResponse compressing = element.respond(100_mm, metresPerSecond(0.2));
    CHECK(compressing.springForce == newtons(50.0) + newtons(900.0));
    CHECK(compressing.damperForce.in(newtons) == doctest::Approx(70.0));
    // The generalized force on the compression coordinate opposes compression.
    CHECK(compressing.force.in(newtons) == doctest::Approx(-1020.0));
    CHECK(compressing.force == -(compressing.springForce + compressing.damperForce));
    CHECK(compressing.dForceDPosition.si() == -9000.0);
    CHECK(compressing.dForceDVelocity.si() == -350.0);
    CHECK(compressing.potential.in(joules) == doctest::Approx(50.0 * 0.1 + 0.5 * 9000.0 * 0.01));
    CHECK(compressing.dissipationPower.in(watts) == doctest::Approx(350.0 * 0.04));

    const AxialResponse rebounding = element.respond(100_mm, metresPerSecond(-0.2));
    CHECK(rebounding.damperForce.in(newtons) == doctest::Approx(-130.0));
    CHECK(rebounding.dForceDVelocity.si() == -650.0);
    CHECK(rebounding.dissipationPower.in(watts) == doctest::Approx(650.0 * 0.04));
    // u = 0 takes the compression branch.
    CHECK(element.respond(100_mm, f64::Velocity{}).dForceDVelocity.si() == -350.0);
}

TEST_CASE("Unilateral contact pushes but never pulls")
{
    const UnilateralSpringDamper contact = makeContact();

    SUBCASE("separated")
    {
        const ContactResponse r = contact.respond(-10_mm, metresPerSecond(5.0));
        CHECK((r.force == f64::Force{} && r.compression == f64::Length{} && r.potential == f64::Energy{}));
        CHECK((r.dForceDPosition == f64::Stiffness{} && r.dForceDVelocity == f64::Damping{}));
    }
    SUBCASE("loading")
    {
        const ContactResponse r = contact.respond(10_mm, metresPerSecond(0.5));
        CHECK(r.normalForce.in(newtons) == doctest::Approx(900.0 + 200.0));
        CHECK(r.force == -r.normalForce);
        CHECK(r.dForceDPosition.si() == -90000.0);
        CHECK(r.dForceDVelocity.si() == -400.0);
        CHECK(r.dissipationPower.in(watts) == doctest::Approx(200.0 * 0.5));
    }
    SUBCASE("damper would pull: force clips and the spring's power is dissipated")
    {
        const ContactResponse r = contact.respond(1_mm, metresPerSecond(-5.0));
        CHECK(r.force == f64::Force{});
        CHECK(r.normalForce == f64::Force{});
        CHECK(r.dForceDPosition == f64::Stiffness{});
        CHECK(r.dForceDVelocity == f64::Damping{});
        CHECK(r.potential.in(joules) == doctest::Approx(0.5 * 90000.0 * 1e-6));
        CHECK(r.dissipationPower.in(watts) == doctest::Approx(-90.0 * -5.0));
        CHECK(r.dissipationPower >= f64::Power{});
    }
}

TEST_CASE("Law create refuses invalid parameters")
{
    const double nan = std::nan("");
    const auto axial = [](double p, double k, double cc, double cr) {
        return AxialSpringDamper::create({.preload = newtons(p),
                                          .stiffness = newtonsPerMetre(k),
                                          .compressionDamping = newtonSecondsPerMetre(cc),
                                          .reboundDamping = newtonSecondsPerMetre(cr)});
    };
    CHECK_FALSE(axial(nan, 1.0, 1.0, 1.0));
    CHECK_FALSE(axial(1.0, -1.0, 1.0, 1.0));
    CHECK_FALSE(axial(1.0, 1.0, -1.0, 1.0));
    CHECK_FALSE(axial(1.0, 1.0, 1.0, INFINITY));
    CHECK(axial(0.0, 0.0, 0.0, 0.0));
    const auto contact = [](double k, double c) {
        return UnilateralSpringDamper::create({.stiffness = newtonsPerMetre(k), .damping = newtonSecondsPerMetre(c)});
    };
    CHECK_FALSE(contact(-1.0, 1.0));
    CHECK_FALSE(contact(1.0, nan));
    CHECK(contact(0.0, 0.0));
    const auto refused = contact(INFINITY, 1.0);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::InvalidInput);
}

// ── Every law: the shared contract ────────────────────────────────────────────

namespace
{

// One fixture per law: the law, the coordinate range to sample (covering every branch) and the rates.
struct AxialFixture
{
    static AxialSpringDamper law()
    {
        return makeAxial();
    }
    static constexpr double qMin = -0.02, qMax = 0.12, qScale = 0.1, vScale = 1.5;
};

struct ContactFixture
{
    static UnilateralSpringDamper law()
    {
        return makeContact();
    }
    static constexpr double qMin = -0.01, qMax = 0.02, qScale = 0.01, vScale = 3.0;
};

struct TravelBufferFixture
{
    static TravelBuffer law()
    {
        return *TravelBuffer::create({.travel = metres(0.065),
                                      .zoneWidth = metres(0.008),
                                      .stiffness = newtonsPerMetre(1.0e6),
                                      .damping = newtonSecondsPerMetre(3000.0)});
    }
    static constexpr double qMin = -0.004, qMax = 0.069, qScale = 0.0365, vScale = 2.0;
};

template<typename L>
LawResponse<f64::Length> respondAt(const L& law, double q, double v)
{
    return law.respond(f64::Length::fromSi(q), f64::Velocity::fromSi(v));
}

// Checks a Jacobian entry against finite differences. Where the two one-sided differences agree, the
// entry must match the central one; at a kink it must match one of the two sides (the branch chosen).
void checkPartial(double exact, double backward, double forward, double scale)
{
    const double tolerance = 1e-5 * scale + 1e-6;
    if (std::abs(forward - backward) <= tolerance)
    {
        CHECK(exact == doctest::Approx(0.5 * (forward + backward)).epsilon(1e-5).scale(scale));
    }
    else
    {
        const bool matchesASide = std::abs(exact - backward) <= tolerance || std::abs(exact - forward) <= tolerance;
        CHECK(matchesASide);
    }
}

} // namespace

TEST_CASE_TEMPLATE("Law power balance holds along trajectories: force * qd = -dV/dt - dissipationPower", F, AxialFixture,
                   ContactFixture, TravelBufferFixture)
{
    const auto law = F::law();
    const double centre = 0.5 * (F::qMin + F::qMax);
    const double amplitude = 0.5 * (F::qMax - F::qMin);
    // Trajectories q(t) = centre + A sin(w t + phase) cross every branch in both directions.
    for (const double w : {3.0, 17.0, 41.0})
    {
        for (const double phase : {0.0, 0.7, 2.1})
        {
            const auto trajectory = [&](double t) { return centre + amplitude * std::sin(w * t + phase); };
            const double dt = 1e-7;
            for (double t = 0.0; t < 2.0 * 3.14159265358979 / w; t += 0.013 / w)
            {
                const double q = trajectory(t);
                const double qd = amplitude * w * std::cos(w * t + phase);
                const LawResponse<f64::Length> r = respondAt(law, q, qd);
                const double dVdt = (respondAt(law, trajectory(t + dt), qd).potential.si() -
                                     respondAt(law, trajectory(t - dt), qd).potential.si()) /
                                    (2.0 * dt);
                const double power = std::abs(r.force.si() * qd) + std::abs(dVdt) + 1.0;
                CHECK(r.force.si() * qd == doctest::Approx(-dVdt - r.dissipationPower.si()).epsilon(1e-6).scale(power));
                CHECK(r.dissipationPower.si() >= 0.0);
            }
        }
    }
}

TEST_CASE_TEMPLATE("Law Jacobians match finite differences", F, AxialFixture, ContactFixture, TravelBufferFixture)
{
    const auto law = F::law();
    const double hq = 1e-7 * F::qScale;
    const double hv = 1e-7 * F::vScale;
    for (double q = F::qMin; q <= F::qMax; q += (F::qMax - F::qMin) / 97.0)
    {
        for (double v = -F::vScale; v <= F::vScale; v += F::vScale / 7.3)
        {
            const LawResponse<f64::Length> r = respondAt(law, q, v);
            const double f = r.force.si();
            const double scaleQ = std::abs(r.dForceDPosition.si()) + 1.0;
            const double scaleV = std::abs(r.dForceDVelocity.si()) + 1.0;
            checkPartial(r.dForceDPosition.si(), (f - respondAt(law, q - hq, v).force.si()) / hq,
                         (respondAt(law, q + hq, v).force.si() - f) / hq, scaleQ);
            checkPartial(r.dForceDVelocity.si(), (f - respondAt(law, q, v - hv).force.si()) / hv,
                         (respondAt(law, q, v + hv).force.si() - f) / hv, scaleV);
        }
    }
}
