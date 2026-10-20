#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <doctest.h>
#include <expected>
#include <limits>
#include <random>
#include <sub0math/units.hpp>
#include <vector>

#include "sub0physics/geometry/norm.hpp"
#include "sub0physics/integrators.hpp"
#include "sub0physics/solver/bounded_periodic_tridiagonal.hpp"
#include "sub0physics/verification/convergence.hpp"

// Fidelity levers (docs/LOD.md, "Levers"): norm arithmetic, warm-started and active-set chain solves, the
// embedded Runge-Kutta pairs, the linearly implicit Euler step and the adaptive advance.

using namespace sub0physics;

namespace
{

// x'' = -w^2 x - 2 zeta w x' as [x, v], with its Jacobian; from x(0) = 1, v(0) = 0.
struct Oscillator
{
    double omega = 2.0;
    double zeta = 0.0;
    std::uint64_t evaluations = 0U;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt)
    {
        ++evaluations;
        dxdt[0] = x[1];
        dxdt[1] = -omega * omega * x[0] - 2.0 * zeta * omega * x[1];
        return {};
    }

    std::expected<void, Fault> jacobian(double, const StateVector<2>&, JacobianMatrix<2>& j) const
    {
        j = {{{0.0, 1.0}, {-omega * omega, -2.0 * zeta * omega}}};
        return {};
    }

    [[nodiscard]] double exact(double t) const
    {
        return std::cos(omega * t);
    }
};

static_assert(OdeSystemWithJacobian<Oscillator, 2>);

template<typename Stepper>
double fixedStepError(int steps, double duration)
{
    Oscillator system;
    Stepper stepper;
    StateVector<2> x{1.0, 0.0};
    const double h = duration / steps;
    for (int i = 0; i < steps; ++i)
    {
        REQUIRE(stepper.attempt(system, i * h, h, x));
        stepper.accept(x);
    }
    return std::abs(x[0] - system.exact(duration));
}

template<typename Stepper>
double adaptiveError(double tolerance, std::uint64_t& accepted, StateVector<2>& end)
{
    Oscillator system;
    Stepper stepper;
    StateVector<2> x{1.0, 0.0};
    const ErrorTolerance<2> tol{{tolerance, tolerance * system.omega}, 0.0};
    const StepSizeController controller{.minStep = 1e-9, .maxStep = 0.5};
    AdaptiveCounters counters;
    double h = 0.01;
    double stopped = 0.0;
    // Twenty outer ticks of 0.25 s, as a caller's fixed grid.
    for (int tick = 0; tick < 20; ++tick)
    {
        REQUIRE(advanceAdaptively<2>(stepper, system, tick * 0.25, (tick + 1) * 0.25, x, h, tol, controller, counters,
                                     [](double, const StateVector<2>&) -> std::expected<void, Fault> { return {}; },
                                     stopped));
        CHECK(stopped == (tick + 1) * 0.25);
    }
    accepted = counters.accepted;
    end = x;
    return std::abs(x[0] - system.exact(5.0));
}

std::vector<double> randomBounds(std::mt19937& random, std::size_t n)
{
    std::uniform_real_distribution<double> depth(-0.01, 0.03);
    std::vector<double> lower(n);
    for (double& value : lower)
    {
        value = depth(random);
    }
    // A smooth contact patch over a quarter of the ring, as a footprint has.
    for (std::size_t i = 0; i < n / 4U; ++i)
    {
        lower[i] = 0.02 * std::sin(3.14159 * static_cast<double>(i) / static_cast<double>(n / 4U));
    }
    return lower;
}

} // namespace

TEST_CASE("Fast and compensated norms agree to rounding and survive extreme magnitudes")
{
    for (const auto& [a, b] : {std::pair{3.0, 4.0}, {0.1, 0.2}, {-1e-3, 7e-4}, {1e200, 1e200}, {1e-200, 3e-200}, {0.0, 0.0}})
    {
        const double reference = std::hypot(a, b);
        CHECK(compensatedNorm(a, b) == doctest::Approx(reference).epsilon(1e-15));
        CHECK(fastNorm(a, b) == doctest::Approx(reference).epsilon(1e-15));
        CHECK(norm(NormArithmetic::Compensated, a, b) == compensatedNorm(a, b));
        CHECK(norm(NormArithmetic::Fast, a, b) == fastNorm(a, b));
    }
    CHECK(fastNorm(3.0, 4.0) == 5.0);
    CHECK(std::isfinite(fastNorm(1e300, 1e300)));
    CHECK(fastNorm(1e-300, 0.0) == 1e-300);
}

TEST_CASE("A cold start is the default and a warm start from the answer needs one sweep")
{
    std::mt19937 random(7U);
    const std::vector<double> lower = randomBounds(random, 256);
    BoundedPeriodicSolution cold;
    BoundedPeriodicSolution again;
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 1300.0, 6800.0, cold));
    again.u.assign(256, 0.5); // ignored without warmStart
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 1300.0, 6800.0, again, {.warmStart = false}));
    CHECK(again.u == cold.u);
    CHECK(again.iterations == cold.iterations);

    BoundedPeriodicSolution warm = cold;
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 1300.0, 6800.0, warm, {.warmStart = true}));
    CHECK(warm.converged);
    CHECK(warm.iterations == 1);
    CHECK(cold.iterations > 10);
}

TEST_CASE("The active-set solve meets the KKT conditions to rounding and matches a tight Gauss-Seidel solve")
{
    std::mt19937 random(11U);
    for (const std::size_t n : {3U, 4U, 17U, 64U, 512U, 1024U})
    {
        for (const double c : {0.0, 30.0 / 0.0044, 200.0 / 0.0011})
        {
            CAPTURE(n);
            CAPTURE(c);
            const std::vector<double> lower = randomBounds(random, n);
            const double a = 300000.0 * 0.0044;
            BoundedPeriodicSolution direct;
            REQUIRE(solveBoundedPeriodicTridiagonalActiveSet(lower, a, c, direct));
            REQUIRE(direct.converged);
            CHECK(direct.iterations <= 40); // random rough bounds; smooth footprint patches take far fewer
            const double scale = a + 2.0 * c;
            CHECK(direct.residual <= 1e-12 * scale);
            for (std::size_t i = 0; i < n; ++i)
            {
                CHECK(direct.u[i] >= lower[i]);
                CHECK(direct.gradient[i] >= -1e-12 * scale);
                CHECK(std::abs(direct.gradient[i] * (direct.u[i] - lower[i])) <= 1e-12 * scale);
            }
            BoundedPeriodicSolution iterative;
            REQUIRE(solveBoundedPeriodicTridiagonal(lower, a, c, iterative, {.maxIterations = 2000000, .tolerance = 1e-9}));
            REQUIRE(iterative.converged);
            for (std::size_t i = 0; i < n; ++i)
            {
                CHECK(direct.u[i] == doctest::Approx(iterative.u[i]).scale(1e-3).epsilon(1e-9));
            }
            CHECK(direct.energy == doctest::Approx(iterative.energy).epsilon(1e-9));

            // From its own answer the active set is already right: one solve.
            BoundedPeriodicSolution warm = direct;
            REQUIRE(solveBoundedPeriodicTridiagonalActiveSet(lower, a, c, warm, {.warmStart = true}));
            CHECK(warm.iterations <= 2);
            CHECK(warm.u == direct.u);
        }
    }
}

TEST_CASE("With no bound active the active-set solve returns zero")
{
    const std::vector<double> lower(12, -0.5);
    BoundedPeriodicSolution solution;
    REQUIRE(solveBoundedPeriodicTridiagonalActiveSet(lower, 2.0, 1.0, solution));
    CHECK(solution.converged);
    for (const double u : solution.u)
    {
        CHECK(u == 0.0);
    }
    CHECK(!solveBoundedPeriodicTridiagonalActiveSet(std::vector<double>{1.0, 2.0}, 1.0, 1.0, solution));
}

TEST_CASE("Embedded Runge-Kutta pairs converge at their order")
{
    const double bs1 = fixedStepError<EmbeddedRungeKutta<2, BogackiShampine32>>(200, 5.0);
    const double bs2 = fixedStepError<EmbeddedRungeKutta<2, BogackiShampine32>>(400, 5.0);
    CHECK(observedOrder(bs1, bs2, 2.0) == doctest::Approx(3.0).epsilon(0.1));
    const double dp1 = fixedStepError<EmbeddedRungeKutta<2, DormandPrince54>>(50, 5.0);
    const double dp2 = fixedStepError<EmbeddedRungeKutta<2, DormandPrince54>>(100, 5.0);
    CHECK(observedOrder(dp1, dp2, 2.0) == doctest::Approx(5.0).epsilon(0.1));
    const double ie1 = fixedStepError<LinearlyImplicitEuler<2>>(2000, 5.0);
    const double ie2 = fixedStepError<LinearlyImplicitEuler<2>>(4000, 5.0);
    CHECK(observedOrder(ie1, ie2, 2.0) == doctest::Approx(1.0).epsilon(0.1));
}

TEST_CASE("Adaptive advance meets a global error that scales with the tolerance, deterministically")
{
    double previous = std::numeric_limits<double>::infinity();
    std::uint64_t previousSteps = 0U;
    for (const double tolerance : {1e-4, 1e-6, 1e-8})
    {
        CAPTURE(tolerance);
        std::uint64_t accepted = 0U;
        StateVector<2> end{};
        const double error = adaptiveError<EmbeddedRungeKutta<2, DormandPrince54>>(tolerance, accepted, end);
        // Local control over about 20 oscillation-scaled steps: the global error stays within a small multiple.
        CHECK(error < 50.0 * tolerance);
        CHECK(error < previous);
        CHECK(accepted > previousSteps);
        previous = error;
        previousSteps = accepted;
        std::uint64_t acceptedAgain = 0U;
        StateVector<2> endAgain{};
        CHECK(adaptiveError<EmbeddedRungeKutta<2, DormandPrince54>>(tolerance, acceptedAgain, endAgain) == error);
        CHECK(endAgain == end);
        CHECK(acceptedAgain == accepted);
    }
    std::uint64_t accepted = 0U;
    StateVector<2> end{};
    CHECK(adaptiveError<EmbeddedRungeKutta<2, BogackiShampine32>>(1e-6, accepted, end) < 1e-4);
    // First order: local control at 1e-4 over 5 s accumulates to a few percent.
    CHECK(adaptiveError<LinearlyImplicitEuler<2>>(1e-4, accepted, end) < 0.1);
}

TEST_CASE("Linearly implicit Euler is stable far beyond the explicit limit and never creates energy")
{
    Oscillator system{.omega = 1000.0};
    LinearlyImplicitEuler<2> stepper;
    StateVector<2> x{1.0, 0.0};
    const auto energy = [&](const StateVector<2>& s) { return 0.5 * s[1] * s[1] + 0.5 * system.omega * system.omega * s[0] * s[0]; };
    double last = energy(x);
    for (std::uint64_t tick = 0; tick < 200U; ++tick)
    {
        REQUIRE(stepper.stepTick(system, tick, sub0math::seconds(0.1), x)); // h w = 100
        const double now = energy(x);
        CHECK(now <= last);
        last = now;
    }
    CHECK(std::abs(x[0]) < 1e-6);
    CHECK(stepper.jacobians() == 200U);
}
