#include <doctest.h>
#include <sub0math/units.hpp>

#include "sub0physics/elements.hpp"

using namespace sub0physics;
using namespace sub0math;
using namespace sub0math::f64::literals;

TEST_CASE("Axial element uses direction-dependent damping")
{
    constexpr AxialSpringDamper element{
        .preload = newtons(50.0),
        .stiffness = newtonsPerMetre(9000.0),
        .compressionDamping = newtonSecondsPerMetre(350.0),
        .reboundDamping = newtonSecondsPerMetre(650.0),
    };
    constexpr AxialResponse compressing = element.respond(100_mm, metresPerSecond(0.2));
    static_assert(compressing.springForce == newtons(50.0) + newtons(900.0));
    CHECK(compressing.damperForce.in(newtons) == doctest::Approx(70.0));
    CHECK(compressing.force.in(newtons) == doctest::Approx(1020.0));
    CHECK(compressing.storedEnergy.in(joules) == doctest::Approx(50.0 * 0.1 + 0.5 * 9000.0 * 0.01));
    CHECK(compressing.dissipationRate.in(watts) == doctest::Approx(350.0 * 0.04));

    const AxialResponse rebounding = element.respond(100_mm, metresPerSecond(-0.2));
    CHECK(rebounding.damperForce.in(newtons) == doctest::Approx(-130.0));
    CHECK(rebounding.dissipationRate.in(watts) == doctest::Approx(650.0 * 0.04));
}

TEST_CASE("Axial element stored energy is the integral of spring force")
{
    constexpr AxialSpringDamper element{newtons(40.0), newtonsPerMetre(22000.0), {}, {}};
    const f64::Length q0 = 10_mm;
    const f64::Length q1 = 20_mm;
    const f64::Length h = micrometres(1.0);
    f64::Energy work{};
    for (f64::Length q = q0; q < q1; q += h)
    {
        work += element.respond(q + 0.5 * h, f64::Velocity{}).springForce * h;
    }
    const f64::Energy stored = element.respond(q1, f64::Velocity{}).storedEnergy - element.respond(q0, f64::Velocity{}).storedEnergy;
    CHECK(work.si() == doctest::Approx(stored.si()).epsilon(1e-6));
}

TEST_CASE("Unilateral contact pushes but never pulls")
{
    constexpr UnilateralSpringDamper contact{newtonsPerMetre(90000.0), newtonSecondsPerMetre(400.0)};

    SUBCASE("separated")
    {
        constexpr ContactResponse r = contact.respond(-10_mm, metresPerSecond(5.0));
        static_assert(r.force == f64::Force{} && r.compression == f64::Length{} && r.storedEnergy == f64::Energy{});
    }
    SUBCASE("loading")
    {
        const ContactResponse r = contact.respond(10_mm, metresPerSecond(0.5));
        CHECK(r.force.in(newtons) == doctest::Approx(900.0 + 200.0));
        CHECK(r.dissipationRate.in(watts) == doctest::Approx(200.0 * 0.5));
    }
    SUBCASE("damper would pull: force clips and the spring's power is dissipated")
    {
        const ContactResponse r = contact.respond(1_mm, metresPerSecond(-5.0));
        CHECK(r.force == f64::Force{});
        CHECK(r.storedEnergy.in(joules) == doctest::Approx(0.5 * 90000.0 * 1e-6));
        CHECK(r.dissipationRate.in(watts) == doctest::Approx(-90.0 * -5.0));
        CHECK(r.dissipationRate >= f64::Power{});
    }
}
