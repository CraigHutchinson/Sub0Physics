#include <cmath>
#include <doctest.h>
#include <ostream> // doctest prints strings (MSVC needs it explicitly)
#include <sub0math/units.hpp>

#include "sub0physics/elements.hpp"

using namespace sub0physics;
using namespace sub0math;

namespace
{

constexpr double Travel = 0.065;
constexpr double Width = 0.008;
constexpr double K = 1.0e6;
constexpr double C = 3000.0;

TravelBuffer makeBuffer()
{
    return *TravelBuffer::create({.travel = metres(Travel),
                                  .zoneWidth = metres(Width),
                                  .stiffness = newtonsPerMetre(K),
                                  .damping = newtonSecondsPerMetre(C)});
}

TravelBufferResponse at(const TravelBuffer& buffer, double q, double v)
{
    return buffer.respond(metres(q), metresPerSecond(v));
}

// BikeLab dist/travel-buffer.js (a9afcd5), evaluated in plain doubles.
struct Reference
{
    double force, potential, lossPower, topForce, bottomForce;
};

Reference bikeLab(double q, double velocity)
{
    const double top = std::max(0.0, Width - q), bottom = std::max(0.0, q - (Travel - Width));
    const double closingTop = std::max(0.0, -velocity), closingBottom = std::max(0.0, velocity);
    const double topForce = K * top * top / Width + C * closingTop * top / Width;
    const double bottomForce = K * bottom * bottom / Width + C * closingBottom * bottom / Width;
    return {topForce - bottomForce,
            K * (std::pow(top, 3) + std::pow(bottom, 3)) / (3 * Width),
            C * (std::pow(closingTop, 2) * top + std::pow(closingBottom, 2) * bottom) / Width,
            topForce,
            bottomForce};
}

double springPart(double q)
{
    const double top = std::max(0.0, Width - q), bottom = std::max(0.0, q - (Travel - Width));
    return K * (top * top - bottom * bottom) / Width;
}

} // namespace

TEST_CASE("Travel buffer is silent in the free zone")
{
    const TravelBuffer buffer = makeBuffer();
    for (const double q : {0.008, 0.02, 0.0325, 0.05, 0.057})
    {
        for (const double v : {-3.0, 0.0, 3.0})
        {
            const TravelBufferResponse r = at(buffer, q, v);
            CHECK(r.force.si() == 0.0);
            CHECK(r.potential.si() == 0.0);
            CHECK(r.lossPower.si() == 0.0);
            CHECK(r.topForce.si() == 0.0);
            CHECK(r.bottomForce.si() == 0.0);
            CHECK(r.dForceDPosition.si() == 0.0);
            CHECK(r.dForceDVelocity.si() == 0.0);
        }
    }
}

TEST_CASE("Travel buffer top zone pushes toward larger q")
{
    const TravelBuffer buffer = makeBuffer();
    // Half-way into the top buffer, entering at 0.5 m/s: depth 4 mm.
    const TravelBufferResponse r = at(buffer, 0.004, -0.5);
    CHECK(r.topForce.si() == doctest::Approx(K * 0.004 * 0.004 / Width + C * 0.5 * 0.004 / Width));
    CHECK(r.topForce.si() == doctest::Approx(2000.0 + 750.0));
    CHECK(r.bottomForce.si() == 0.0);
    CHECK(r.force.si() == r.topForce.si());
    CHECK(r.potential.si() == doctest::Approx(K * 0.004 * 0.004 * 0.004 / (3 * Width)));
    CHECK(r.lossPower.si() == doctest::Approx(C * 0.25 * 0.004 / Width));
    // Leaving the buffer is undamped.
    const TravelBufferResponse leaving = at(buffer, 0.004, 0.5);
    CHECK(leaving.topForce.si() == doctest::Approx(2000.0));
    CHECK(leaving.lossPower.si() == 0.0);
}

TEST_CASE("Travel buffer bottom zone pushes toward smaller q")
{
    const TravelBuffer buffer = makeBuffer();
    // 4 mm into the bottom buffer, entering at 0.5 m/s.
    const double q = Travel - Width + 0.004;
    const TravelBufferResponse r = at(buffer, q, 0.5);
    CHECK(r.bottomForce.si() == doctest::Approx(2000.0 + 750.0));
    CHECK(r.topForce.si() == 0.0);
    CHECK(r.force.si() == doctest::Approx(-2750.0));
    CHECK(r.potential.si() == doctest::Approx(K * 0.004 * 0.004 * 0.004 / (3 * Width)));
    CHECK(r.lossPower.si() == doctest::Approx(C * 0.25 * 0.004 / Width));
    const TravelBufferResponse leaving = at(buffer, q, -0.5);
    CHECK(leaving.force.si() == doctest::Approx(-2000.0));
    CHECK(leaving.lossPower.si() == 0.0);
}

TEST_CASE("Travel buffer at the kinks")
{
    const TravelBuffer buffer = makeBuffer();
    for (const double q : {Width, Travel - Width})
    {
        const TravelBufferResponse r = at(buffer, q, 0.7);
        CHECK(r.force.si() == 0.0);
        CHECK(r.potential.si() == 0.0);
        CHECK(r.dForceDPosition.si() == 0.0);
    }
    // Fully compressed: top zone fully engaged at q = 0 (depth = width), force k * width.
    const TravelBufferResponse full = at(buffer, 0.0, 0.0);
    CHECK(full.force.si() == doctest::Approx(K * Width));
    CHECK(full.potential.si() == doctest::Approx(K * Width * Width / 3.0));
    // Velocity kink: u = 0 selects the bottom branch for the velocity Jacobian.
    const double inBottom = Travel - Width + 0.002;
    CHECK(at(buffer, inBottom, 0.0).dForceDVelocity.si() == doctest::Approx(-C * 0.002 / Width));
    CHECK(at(buffer, 0.002, 0.0).dForceDVelocity.si() == 0.0);
    CHECK(at(buffer, 0.002, -1e-9).dForceDVelocity.si() == doctest::Approx(-C * (Width - 0.002) / Width));
}

TEST_CASE("Travel buffer potential is the integral of the spring force")
{
    const TravelBuffer buffer = makeBuffer();
    constexpr double h = 1e-7;
    for (double q = -0.002; q < Travel + 0.002; q += 0.0013)
    {
        const double dUdq = (at(buffer, q + h, 0.0).potential.si() - at(buffer, q - h, 0.0).potential.si()) / (2 * h);
        CHECK(-dUdq == doctest::Approx(springPart(q)).epsilon(1e-6).scale(1.0));
        // At zero velocity the response is the spring part alone.
        CHECK(at(buffer, q, 0.0).force.si() == doctest::Approx(springPart(q)).epsilon(1e-12));
    }
}

TEST_CASE("Travel buffer power balance: force * u = -dU/dt - loss")
{
    const TravelBuffer buffer = makeBuffer();
    constexpr double h = 1e-7;
    for (double q = -0.002; q < Travel + 0.002; q += 0.0017)
    {
        for (const double v : {-2.0, -0.3, 0.4, 1.5})
        {
            const TravelBufferResponse r = at(buffer, q, v);
            const double dUdq = (at(buffer, q + h, v).potential.si() - at(buffer, q - h, v).potential.si()) / (2 * h);
            CHECK(r.force.si() * v == doctest::Approx(-dUdq * v - r.lossPower.si()).epsilon(1e-6).scale(1.0));
        }
    }
}

TEST_CASE("Travel buffer loss power is never negative")
{
    const TravelBuffer buffer = makeBuffer();
    for (double q = -0.01; q < Travel + 0.01; q += 0.0007)
    {
        for (double v = -4.0; v <= 4.0; v += 0.5)
        {
            const TravelBufferResponse r = at(buffer, q, v);
            CHECK(r.lossPower.si() >= 0.0);
            CHECK(r.potential.si() >= 0.0);
            CHECK(r.topForce.si() >= 0.0);
            CHECK(r.bottomForce.si() >= 0.0);
        }
    }
}

TEST_CASE("Travel buffer Jacobians match central finite differences away from kinks")
{
    const TravelBuffer buffer = makeBuffer();
    constexpr double hq = 1e-7;
    constexpr double hv = 1e-5;
    int checked = 0;
    for (double q = -0.0031; q < Travel + 0.003; q += 0.00173)
    {
        // Skip samples whose stencil would straddle a q kink.
        const bool nearKink = std::abs(q - Width) < 2 * hq || std::abs(q - (Travel - Width)) < 2 * hq;
        if (nearKink)
        {
            continue;
        }
        for (const double v : {-2.0, -0.31, 0.27, 1.9})
        {
            const double dq = (at(buffer, q + hq, v).force.si() - at(buffer, q - hq, v).force.si()) / (2 * hq);
            const double dv = (at(buffer, q, v + hv).force.si() - at(buffer, q, v - hv).force.si()) / (2 * hv);
            const TravelBufferResponse r = at(buffer, q, v);
            CHECK(r.dForceDPosition.si() == doctest::Approx(dq).epsilon(1e-6).scale(1e4));
            CHECK(r.dForceDVelocity.si() == doctest::Approx(dv).epsilon(1e-6).scale(1.0));
            ++checked;
        }
    }
    CHECK(checked > 100);
}

TEST_CASE("Travel buffer matches BikeLab's formula on a grid")
{
    const TravelBuffer buffer = makeBuffer();
    const auto close = [](double actual, double expected) {
        return std::abs(actual - expected) <= 1e-12 * std::max(std::abs(expected), 1e-300);
    };
    int compared = 0;
    for (double q = -0.01; q < Travel + 0.01; q += 0.00093)
    {
        for (double v = -5.0; v <= 5.0; v += 0.37)
        {
            const Reference ref = bikeLab(q, v);
            const TravelBufferResponse r = at(buffer, q, v);
            CHECK(close(r.force.si(), ref.force));
            CHECK(close(r.potential.si(), ref.potential));
            CHECK(close(r.lossPower.si(), ref.lossPower));
            CHECK(close(r.topForce.si(), ref.topForce));
            CHECK(close(r.bottomForce.si(), ref.bottomForce));
            ++compared;
        }
    }
    CHECK(compared > 1000);
}

TEST_CASE("Travel buffer factory validates its parameters")
{
    const auto build = [](double travel, double width, double k, double c) {
        return TravelBuffer::create({.travel = metres(travel),
                                     .zoneWidth = metres(width),
                                     .stiffness = newtonsPerMetre(k),
                                     .damping = newtonSecondsPerMetre(c)});
    };
    CHECK(build(0.065, 0.008, 1e6, 3000.0).has_value());
    CHECK(build(0.0160001, 0.008, 0.0, 0.0).has_value()); // zones just apart
    CHECK_FALSE(build(0.016, 0.008, 0.0, 0.0).has_value()); // zones touch exactly: BikeLab rejects 2 * width >= travel
    CHECK_FALSE(build(0.0, 0.0, 1e6, 0.0).has_value());
    CHECK_FALSE(build(-0.065, 0.008, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.0, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.065, -0.008, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.0159, 0.008, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.008, -1.0, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.008, 1e6, -1.0).has_value());
    const double nan = std::nan("");
    const double inf = HUGE_VAL;
    CHECK_FALSE(build(nan, 0.008, 1e6, 0.0).has_value());
    CHECK_FALSE(build(inf, 0.008, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.065, nan, 1e6, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.008, nan, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.008, inf, 0.0).has_value());
    CHECK_FALSE(build(0.065, 0.008, 1e6, nan).has_value());
    CHECK_FALSE(build(0.065, 0.008, 1e6, inf).has_value());
    const auto bad = build(0.065, 0.0, 1e6, 0.0);
    CHECK(bad.error().kind == FaultKind::InvalidInput);
}
