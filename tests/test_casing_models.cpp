#include <doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <sub0math/number/det_math.hpp>

#include "sub0physics/multibody/casing_models.hpp"

// The cheaper casing models' pieces (docs/WHEEL.md, "Wheel-model study"): the tread's face laws and the nodal-to-rigid
// state transfer of a model switch.

using namespace sub0physics::multibody;

namespace
{

constexpr double K = 31566.0;
constexpr double C = 62.5;

// Central difference of the potential along (dp0, dp1).
double slope(auto law, double p0, double p1, double dp0, double dp1)
{
    constexpr double h = 1e-9;
    return (law(p0 + h * dp0, p1 + h * dp1).potential - law(p0 - h * dp0, p1 - h * dp1).potential) / (2.0 * h);
}

} // namespace

TEST_SUITE("casing models")
{
    TEST_CASE("smoothed tread: the witness law outside the band, continuous node forces across level")
    {
        const double band = 5e-4;
        const auto law = [&](double p0, double p1) { return smoothedTreadFace(K, 0.0, band, p0, p1, 0.0, 0.0); };
        // Outside the band: the whole force k max(p0, p1) at the deeper node.
        const auto far = law(0.003, 0.001);
        CHECK(far.force0 == doctest::Approx(K * 0.003).epsilon(1e-14));
        CHECK(far.force1 == 0.0);
        // Across level the shares move continuously; the total stays within the band's smoothing.
        double previous0 = law(0.002 - 120e-5, 0.002).force0;
        for (int i = -120; i <= 120; ++i)
        {
            const double d = 1e-5 * i;
            const auto r = law(0.002 + d, 0.002);
            CHECK(std::abs(r.force0 - previous0) < 0.03 * K * 0.002);
            CHECK(r.force0 + r.force1 <= K * (0.002 + std::max(d, 0.0)) * (1.0 + 1e-12));
            previous0 = r.force0;
        }
        // Forces are the potential's gradient (conservative), inside and outside the band.
        for (const auto& p : {std::array{0.002, 0.0021}, std::array{0.003, 0.001}, std::array{0.001, 0.0012}})
        {
            const auto r = law(p[0], p[1]);
            CHECK(slope(law, p[0], p[1], 1.0, 0.0) == doctest::Approx(r.force0).epsilon(1e-5));
            CHECK(slope(law, p[0], p[1], 0.0, 1.0) == doctest::Approx(r.force1).epsilon(1e-5));
        }
    }

    TEST_CASE("distributed tread: closed form against quadrature, conservative, energy identity with the damper")
    {
        const auto law = [](double p0, double p1) { return distributedTreadFace(K, 0.0, p0, p1, 0.0, 0.0); };
        for (const auto& p : {std::array{0.002, 0.003}, std::array{0.003, -0.002}, std::array{-0.001, 0.002}})
        {
            // Midpoint quadrature of k max(0, p(s)) (1 - s) and s.
            double f0 = 0.0, f1 = 0.0;
            constexpr int M = 200000;
            for (int i = 0; i < M; ++i)
            {
                const double s = (i + 0.5) / M;
                const double ps = std::max(0.0, (1.0 - s) * p[0] + s * p[1]);
                f0 += K * ps * (1.0 - s) / M;
                f1 += K * ps * s / M;
            }
            const auto r = law(p[0], p[1]);
            CHECK(r.force0 == doctest::Approx(f0).epsilon(1e-6));
            CHECK(r.force1 == doctest::Approx(f1).epsilon(1e-6));
            CHECK(slope(law, p[0], p[1], 1.0, 0.0) == doctest::Approx(r.force0).epsilon(1e-5));
            CHECK(slope(law, p[0], p[1], 0.0, 1.0) == doctest::Approx(r.force1).epsilon(1e-5));
        }
        // With the damper: the power taken equals dV/dt plus the dissipation, which is never negative.
        const auto r = distributedTreadFace(K, C, 0.002, 0.001, 0.3, -0.1);
        const auto e = distributedTreadFace(K, 0.0, 0.002, 0.001, 0.0, 0.0);
        const double dVdt = e.elastic0 * 0.3 + e.elastic1 * -0.1;
        CHECK(r.force0 * 0.3 + r.force1 * -0.1 == doctest::Approx(dVdt + r.dissipationPower).epsilon(1e-12));
        CHECK(r.dissipationPower >= 0.0);
        const auto s = smoothedTreadFace(K, C, 5e-4, 0.002, 0.0018, 0.3, -0.1);
        CHECK(s.dissipationPower >= 0.0);
    }

    TEST_CASE("nodal casing to rigid ring and back: momentum kept, energy never gained")
    {
        constexpr std::size_t N = 16;
        constexpr double R = 0.35;
        constexpr double m = 1.0 / N;
        constexpr double nodeSpin = N * m * 0.008 * 0.008 / 2.0;
        constexpr double hostInertia = 0.08;
        std::vector<double> rx(N), rz(N), x(N), z(N), vx(N), vz(N), mass(N, m);
        double ringInertia = nodeSpin;
        for (std::size_t i = 0; i < N; ++i)
        {
            const auto sc = sub0math::det::sincos(2.0 * 3.14159265358979323846 * static_cast<double>(i) / N);
            rx[i] = R * sc.sin;
            rz[i] = -R * sc.cos;
            ringInertia += m * R * R;
        }
        // A deformed, moving casing: the rest shape turned by 0.3 rad on a host at 0.2 rad, plus a k = 2 wobble.
        const double hostAngle = 0.2, hostRate = 3.0;
        const auto turn = sub0math::det::sincos(0.3);
        for (std::size_t i = 0; i < N; ++i)
        {
            const double wobble = 0.002 * (i % 2 == 0 ? 1.0 : -1.0);
            x[i] = 0.1 + turn.cos * rx[i] + turn.sin * rz[i] + wobble;
            z[i] = 0.4 + turn.cos * rz[i] - turn.sin * rx[i];
            vx[i] = 1.0 + 0.2 * rz[i] + 0.05 * wobble * 100.0;
            vz[i] = -0.5 - 0.2 * rx[i];
        }
        const NodalCasingView view{.x = x, .z = z, .vx = vx, .vz = vz, .mass = mass, .restX = rx, .restZ = rz,
                                   .hostAngle = hostAngle, .hostRate = hostRate, .nodeSpin = nodeSpin};
        const PlanarRigidState ring = fitRigidRing(view, ringInertia);
        // Linear momentum and angular momentum about the origin (host's own spin excluded: unchanged).
        double px = 0.0, pz = 0.0, l = nodeSpin * hostRate, kinetic = 0.5 * nodeSpin * hostRate * hostRate;
        for (std::size_t i = 0; i < N; ++i)
        {
            px += m * vx[i];
            pz += m * vz[i];
            l += m * (z[i] * vx[i] - x[i] * vz[i]);
            kinetic += 0.5 * m * (vx[i] * vx[i] + vz[i] * vz[i]);
        }
        const double ringMass = N * m;
        CHECK(ringMass * ring.vx == doctest::Approx(px).epsilon(1e-14));
        CHECK(ringMass * ring.vz == doctest::Approx(pz).epsilon(1e-14));
        const double lRing = ringInertia * ring.rate + ringMass * (ring.z * ring.vx - ring.x * ring.vz);
        CHECK(lRing == doctest::Approx(l).epsilon(1e-13));
        const double kineticRing = 0.5 * ringMass * (ring.vx * ring.vx + ring.vz * ring.vz) + 0.5 * ringInertia * ring.rate * ring.rate;
        CHECK(kineticRing <= kinetic);
        CHECK(ring.angle == doctest::Approx(0.3).epsilon(1e-12)); // the shape was turned by 0.3 rad in the world
        // Back to nodes: exact for the rigid ring (momentum kept, nodes on the ring).
        std::vector<double> x2(N), z2(N), vx2(N), vz2(N);
        const double rate = placeNodalCasing(ring, rx, rz, hostInertia, hostRate, nodeSpin, x2, z2, vx2, vz2);
        double px2 = 0.0, l2 = (hostInertia + nodeSpin) * rate;
        for (std::size_t i = 0; i < N; ++i)
        {
            px2 += m * vx2[i];
            l2 += m * (z2[i] * vx2[i] - x2[i] * vz2[i]);
        }
        CHECK(px2 == doctest::Approx(px).epsilon(1e-14));
        CHECK(l2 == doctest::Approx(hostInertia * hostRate + lRing).epsilon(1e-13));
        // A second fit of the placed casing is the same ring (the transfer back is exact).
        const NodalCasingView again{.x = x2, .z = z2, .vx = vx2, .vz = vz2, .mass = mass, .restX = rx, .restZ = rz,
                                    .hostAngle = hostAngle, .hostRate = ring.rate, .nodeSpin = 0.0};
        const PlanarRigidState refit = fitRigidRing(again, ringInertia - nodeSpin);
        CHECK(refit.x == doctest::Approx(ring.x).epsilon(1e-14));
        CHECK(refit.angle == doctest::Approx(ring.angle).epsilon(1e-12));
        CHECK(refit.rate == doctest::Approx(ring.rate).epsilon(1e-12));
    }
}
