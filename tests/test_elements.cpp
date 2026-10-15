#include <doctest.h>

#include "sub0physics/elements.hpp"

using namespace sub0physics;

TEST_CASE("Axial element uses direction-dependent damping")
{
    constexpr AxialSpringDamper element{50.0, 9000.0, 350.0, 650.0};
    constexpr AxialResponse compressing = element.respond(0.1, 0.2);
    static_assert(compressing.springForce == 50.0 + 900.0);
    CHECK(compressing.damperForce == doctest::Approx(70.0));
    CHECK(compressing.force == doctest::Approx(1020.0));
    CHECK(compressing.storedEnergy == doctest::Approx(50.0 * 0.1 + 0.5 * 9000.0 * 0.01));
    CHECK(compressing.dissipationRate == doctest::Approx(350.0 * 0.04));

    const AxialResponse rebounding = element.respond(0.1, -0.2);
    CHECK(rebounding.damperForce == doctest::Approx(-130.0));
    CHECK(rebounding.dissipationRate == doctest::Approx(650.0 * 0.04));
}

TEST_CASE("Axial element stored energy is the integral of spring force")
{
    constexpr AxialSpringDamper element{40.0, 22000.0, 0.0, 0.0};
    const double q0 = 0.01;
    const double q1 = 0.02;
    const double h = 1e-6;
    double work = 0.0;
    for (double q = q0; q < q1; q += h)
    {
        work += element.respond(q + 0.5 * h, 0.0).springForce * h;
    }
    CHECK(work == doctest::Approx(element.respond(q1, 0.0).storedEnergy - element.respond(q0, 0.0).storedEnergy).epsilon(1e-6));
}

TEST_CASE("Unilateral contact pushes but never pulls")
{
    constexpr UnilateralSpringDamper contact{90000.0, 400.0};

    SUBCASE("separated")
    {
        constexpr ContactResponse r = contact.respond(-0.01, 5.0);
        static_assert(r.force == 0.0 && r.compression == 0.0 && r.storedEnergy == 0.0);
    }
    SUBCASE("loading")
    {
        const ContactResponse r = contact.respond(0.01, 0.5);
        CHECK(r.force == doctest::Approx(900.0 + 200.0));
        CHECK(r.dissipationRate == doctest::Approx(200.0 * 0.5));
    }
    SUBCASE("damper would pull: force clips and the spring's power is dissipated")
    {
        const ContactResponse r = contact.respond(0.001, -5.0);
        CHECK(r.force == 0.0);
        CHECK(r.storedEnergy == doctest::Approx(0.5 * 90000.0 * 1e-6));
        CHECK(r.dissipationRate == doctest::Approx(-90.0 * -5.0));
        CHECK(r.dissipationRate >= 0.0);
    }
}
