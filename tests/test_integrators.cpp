#include <cmath>
#include <cstddef>
#include <cstdint>
#include <doctest.h>
#include <sub0math/units.hpp>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include "sub0physics/integrators.hpp"
#include "sub0physics/verification/convergence.hpp"

using namespace sub0physics;

namespace
{

// x'' = -w^2 x as [x, v]; exact solution x = cos(wt) from x(0) = 1, v(0) = 0.
struct HarmonicOscillator
{
    double omega = 2.0;
    double ceiling = std::numeric_limits<double>::infinity();
    std::size_t evaluations = 0U;
    double observedPosition = 0.0;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt)
    {
        ++evaluations;
        observedPosition = x[0];
        if (x[0] > ceiling)
        {
            return std::unexpected(Fault{FaultKind::OutsideDomain, "above ceiling"});
        }
        dxdt[0] = x[1];
        dxdt[1] = -omega * omega * x[0];
        return {};
    }
};

struct ConstantAcceleration
{
    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt)
    {
        dxdt[0] = x[1];
        dxdt[1] = -9.0;
        return {};
    }
};

struct Explosive
{
    std::expected<void, Fault> evaluate(double, const StateVector<1>&, StateVector<1>& dxdt)
    {
        dxdt[0] = std::numeric_limits<double>::infinity();
        return {};
    }
};

template<typename Integrator>
double finalError(double h)
{
    HarmonicOscillator system;
    Integrator integrator;
    StateVector<2> x{1.0, 0.0};
    const auto steps = static_cast<std::size_t>(std::llround(1.0 / h));
    for (std::size_t i = 0; i < steps; ++i)
    {
        REQUIRE(integrator.step(system, sub0math::seconds(static_cast<double>(i) * h), sub0math::seconds(h), x));
    }
    return std::abs(x[0] - std::cos(system.omega * 1.0));
}

} // namespace

TEST_CASE("Rk4 converges at fourth order")
{
    const double coarse = finalError<Rk4<2>>(0.02);
    const double fine = finalError<Rk4<2>>(0.01);
    CHECK(observedOrder(coarse, fine, 2.0) == doctest::Approx(4.0).epsilon(0.05));
    CHECK(fine < 1e-8);
}

TEST_CASE("Midpoint converges at second order")
{
    const double coarse = finalError<Midpoint<2>>(0.02);
    const double fine = finalError<Midpoint<2>>(0.01);
    CHECK(observedOrder(coarse, fine, 2.0) == doctest::Approx(2.0).epsilon(0.05));
}

TEST_CASE("Midpoint is exact under constant acceleration")
{
    ConstantAcceleration system;
    Midpoint<2> integrator;
    StateVector<2> x{10.0, 3.0};
    const double h = 0.1;
    for (int i = 0; i < 20; ++i)
    {
        REQUIRE(integrator.step(system, sub0math::seconds(i * h), sub0math::seconds(h), x));
    }
    CHECK(x[0] == doctest::Approx(10.0 + 3.0 * 2.0 - 0.5 * 9.0 * 4.0).epsilon(1e-12));
    CHECK(x[1] == doctest::Approx(3.0 - 9.0 * 2.0).epsilon(1e-12));
}

TEST_CASE("First stage is reused from the committed endpoint")
{
    HarmonicOscillator system;
    Rk4<2> rk4;
    StateVector<2> x{1.0, 0.0};
    REQUIRE(rk4.step(system, sub0math::seconds(0.0), sub0math::seconds(0.01), x));
    CHECK(system.evaluations == 5U);
    REQUIRE(rk4.step(system, sub0math::seconds(0.01), sub0math::seconds(0.01), x));
    CHECK(system.evaluations == 9U);

    SUBCASE("a changed state is re-evaluated")
    {
        x[1] += 1e-3;
        REQUIRE(rk4.step(system, sub0math::seconds(0.02), sub0math::seconds(0.01), x));
        CHECK(system.evaluations == 14U);
    }
    SUBCASE("invalidate forces a fresh first stage")
    {
        rk4.invalidate();
        REQUIRE(rk4.step(system, sub0math::seconds(0.02), sub0math::seconds(0.01), x));
        CHECK(system.evaluations == 14U);
    }

    HarmonicOscillator other;
    Midpoint<2> midpoint;
    StateVector<2> y{1.0, 0.0};
    REQUIRE(midpoint.step(other, sub0math::seconds(0.0), sub0math::seconds(0.01), y));
    REQUIRE(midpoint.step(other, sub0math::seconds(0.01), sub0math::seconds(0.01), y));
    CHECK(other.evaluations == 5U);
}

TEST_CASE("A refused step leaves state and observation unchanged")
{
    HarmonicOscillator system;
    system.ceiling = 1.0;
    Rk4<2> rk4;
    StateVector<2> x{0.999, 1.0};
    const StateVector<2> before = x;
    const auto result = rk4.step(system, sub0math::seconds(0.0), sub0math::seconds(0.01), x);
    REQUIRE_FALSE(result);
    CHECK(result.error().kind == FaultKind::OutsideDomain);
    CHECK(x == before);
    CHECK(system.observedPosition == before[0]);
}

TEST_CASE("Non-finite results are refused")
{
    Explosive system;
    Rk4<1> rk4;
    StateVector<1> x{0.0};
    const auto result = rk4.step(system, sub0math::seconds(0.0), sub0math::seconds(0.1), x);
    REQUIRE_FALSE(result);
    CHECK(result.error().kind == FaultKind::NonFiniteState);
    CHECK(x[0] == 0.0);
}

TEST_CASE("Invalid step sizes are refused")
{
    HarmonicOscillator system;
    Rk4<2> rk4;
    StateVector<2> x{1.0, 0.0};
    for (const double h : {0.0, -0.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
    {
        const auto result = rk4.step(system, sub0math::seconds(0.0), sub0math::seconds(h), x);
        REQUIRE_FALSE(result);
        CHECK(result.error().kind == FaultKind::InvalidInput);
    }
    CHECK(system.evaluations == 0U);
}

TEST_CASE("Observed order rejects meaningless input")
{
    CHECK(std::isnan(observedOrder(0.0, 1.0, 2.0)));
    CHECK(std::isnan(observedOrder(1.0, 1.0, 1.0)));
}

TEST_CASE("Tick stepping keeps exact times and always reuses the endpoint derivative")
{
    HarmonicOscillator timed;
    Rk4<2> byTick;
    StateVector<2> x{1.0, 0.0};
    const double h = 0.0005; // not representable: (n - 1) h + h often differs from n h
    const std::uint64_t ticks = 16000U;
    for (std::uint64_t tick = 0; tick < ticks; ++tick)
    {
        REQUIRE(byTick.stepTick(timed, tick, sub0math::seconds(h), x));
    }
    CHECK(timed.evaluations == 5U + 4U * (ticks - 1U));

    HarmonicOscillator accumulated;
    Rk4<2> byTime;
    StateVector<2> y{1.0, 0.0};
    for (std::uint64_t tick = 0; tick < ticks; ++tick)
    {
        REQUIRE(byTime.step(accumulated, sub0math::seconds(static_cast<double>(tick) * h), sub0math::seconds(h), y));
    }
    CHECK(accumulated.evaluations > timed.evaluations);
    CHECK(x[0] == doctest::Approx(y[0]).epsilon(1e-12));
}

TEST_CASE("Dynamic-size Runge-Kutta matches the fixed size exactly and checks the state size")
{
    struct SpanOscillator
    {
        HarmonicOscillator inner;
        std::expected<void, Fault> evaluate(double t, std::span<const double> x, std::span<double> dxdt)
        {
            StateVector<2> in{x[0], x[1]};
            StateVector<2> out{};
            auto r = inner.evaluate(t, in, out);
            dxdt[0] = out[0];
            dxdt[1] = out[1];
            return r;
        }
    };
    static_assert(OdeSystem<SpanOscillator, std::dynamic_extent>);
    HarmonicOscillator fixedSystem;
    SpanOscillator dynamicSystem;
    Rk4<2> fixed;
    Rk4<std::dynamic_extent> dynamic(2);
    CHECK(dynamic.size() == 2U);
    StateVector<2> x{1.0, 0.0};
    std::vector<double> y{1.0, 0.0};
    for (std::uint64_t tick = 0; tick < 1000U; ++tick)
    {
        REQUIRE(fixed.stepTick(fixedSystem, tick, sub0math::seconds(0.001), x));
        REQUIRE(dynamic.stepTick(dynamicSystem, tick, sub0math::seconds(0.001), y));
    }
    CHECK(x[0] == y[0]);
    CHECK(x[1] == y[1]);
    CHECK(fixedSystem.evaluations == dynamicSystem.inner.evaluations);

    std::vector<double> wrong{1.0, 0.0, 0.0};
    const auto refused = dynamic.stepTick(dynamicSystem, 1000U, sub0math::seconds(0.001), wrong);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::InvalidInput);
}
