#include <algorithm>
#include <cmath>
#include <doctest.h>
#include <random>
#include <vector>

#include "sub0physics/solver/bounded_periodic_tridiagonal.hpp"

using namespace sub0physics;

namespace
{

// dE/du by the definition of the energy, not by the solver's rearrangement: g = (D u - c (previous + next)).
std::vector<double> gradientOf(const std::vector<double>& u, double a, double c)
{
    const std::size_t n = u.size();
    std::vector<double> g(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        g[i] = (a + 2.0 * c) * u[i] - c * (u[(i + n - 1U) % n] + u[(i + 1U) % n]);
    }
    return g;
}

double energyOf(const std::vector<double>& u, double a, double c)
{
    double e = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i)
    {
        const double step = u[(i + 1U) % u.size()] - u[i];
        e += 0.5 * a * u[i] * u[i] + 0.5 * c * step * step;
    }
    return e;
}

} // namespace

TEST_CASE("With every bound inactive the minimiser is zero")
{
    BoundedPeriodicSolution solution;
    const std::vector<double> lower(12, -0.5);
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 2.0, 1.0, solution, {.maxIterations = 30000, .tolerance = 1e-13}));
    CHECK(solution.converged);
    for (std::size_t i = 0; i < lower.size(); ++i)
    {
        CHECK(solution.u[i] == doctest::Approx(0.0).epsilon(1e-9).scale(1.0));
        CHECK(std::abs(solution.gradient[i]) < 1e-9);
    }
    CHECK(solution.energy == doctest::Approx(0.0).scale(1.0));
}

TEST_CASE("One active bound has the closed-form decaying solution")
{
    // Node 0 held at L: every other node satisfies c (u[j-1] + u[j+1]) = D u[j], so u[j] = L (rho^j + rho^(n-j)) / (1 + rho^n)
    // with rho the root of c rho^2 - D rho + c = 0 inside the unit interval.
    const double a = 1.0;
    const double c = 1.0;
    const double d = a + 2.0 * c;
    const double rho = (d - std::sqrt(d * d - 4.0 * c * c)) / (2.0 * c);
    for (const std::size_t n : {std::size_t{4}, std::size_t{9}, std::size_t{64}})
    {
        std::vector<double> lower(n, 0.0);
        const double bound = 0.37;
        lower[0] = bound;
        BoundedPeriodicSolution solution;
        REQUIRE(solveBoundedPeriodicTridiagonal(lower, a, c, solution, {.maxIterations = 30000, .tolerance = 1e-13}));
        REQUIRE(solution.converged);
        const double denominator = 1.0 + std::pow(rho, static_cast<double>(n));
        for (std::size_t j = 0; j < n; ++j)
        {
            const double expected = bound * (std::pow(rho, static_cast<double>(j)) + std::pow(rho, static_cast<double>(n - j))) / denominator;
            CHECK(solution.u[j] == doctest::Approx(expected).epsilon(1e-9).scale(1.0));
        }
        // Reaction at the held node, and the energy 1/2 L g0 (the bound is the only load).
        const double reaction = d * bound - c * 2.0 * bound * (rho + std::pow(rho, static_cast<double>(n - 1U))) / denominator;
        CHECK(solution.gradient[0] == doctest::Approx(reaction).epsilon(1e-9));
        CHECK(solution.energy == doctest::Approx(0.5 * bound * reaction).epsilon(1e-9));
        for (std::size_t j = 1; j < n; ++j)
        {
            CHECK(std::abs(solution.gradient[j]) < 1e-9);
        }
    }
}

TEST_CASE("Without coupling every node simply meets its bound")
{
    const std::vector<double> lower{0.0, 0.2, -1.0, 0.7, 0.0};
    BoundedPeriodicSolution solution;
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 3.0, 0.0, solution));
    CHECK(solution.converged);
    CHECK(solution.u[1] == 0.2);
    CHECK(solution.u[3] == 0.7);
    CHECK(solution.gradient[3] == doctest::Approx(3.0 * 0.7));
    CHECK(solution.u[2] == 0.0);
}

TEST_CASE("Random problems satisfy the KKT conditions")
{
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> bound(0.0, 1.0);
    std::bernoulli_distribution active(0.3);
    std::uniform_real_distribution<double> stiffness(0.5, 20.0);
    for (int round = 0; round < 60; ++round)
    {
        const std::size_t n = 4U + static_cast<std::size_t>(round) * 3U;
        std::vector<double> lower(n);
        for (double& v : lower)
        {
            v = active(rng) ? bound(rng) : -0.2;
        }
        const double a = stiffness(rng);
        const double c = stiffness(rng) * 0.2;
        BoundedPeriodicSolution solution;
        const BoundedPeriodicOptions options{.maxIterations = 30000, .tolerance = 1e-11};
        REQUIRE(solveBoundedPeriodicTridiagonal(lower, a, c, solution, options));
        REQUIRE(solution.converged);
        const std::vector<double> g = gradientOf(solution.u, a, c);
        for (std::size_t i = 0; i < n; ++i)
        {
            CHECK(solution.u[i] >= lower[i]);
            CHECK(g[i] >= -1e-9);
            CHECK(std::abs(g[i] * (solution.u[i] - lower[i])) < 1e-9);
            CHECK(solution.gradient[i] == doctest::Approx(g[i]).epsilon(1e-9).scale(1.0));
        }
        CHECK(solution.energy == doctest::Approx(energyOf(solution.u, a, c)).epsilon(1e-9));
        // Optimality against feasible perturbations: no random feasible point has lower energy.
        std::uniform_real_distribution<double> jitter(0.0, 0.05);
        for (int k = 0; k < 20; ++k)
        {
            std::vector<double> v = solution.u;
            for (std::size_t i = 0; i < n; ++i)
            {
                v[i] = std::max(lower[i], v[i] + (jitter(rng) - 0.025));
            }
            CHECK(energyOf(v, a, c) >= solution.energy - 1e-12);
        }
    }
}

TEST_CASE("The solver reports an unconverged problem and rejects invalid ones")
{
    BoundedPeriodicSolution solution;
    std::vector<double> lower(64, 0.0);
    lower[0] = 1.0;
    REQUIRE(solveBoundedPeriodicTridiagonal(lower, 1e-3, 1.0, solution, {.maxIterations = 3, .tolerance = 1e-12}));
    CHECK_FALSE(solution.converged);
    CHECK(solution.iterations == 3);
    CHECK(solution.residual > 1e-12);

    const std::vector<double> ok(8, 0.0);
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(std::vector<double>(2, 0.0), 1.0, 1.0, solution));
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(ok, 0.0, 1.0, solution));
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(ok, 1.0, -1.0, solution));
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(ok, std::nan(""), 1.0, solution));
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(ok, 1.0, 1.0, solution, {.maxIterations = 0, .tolerance = 1e-7}));
    std::vector<double> bad = ok;
    bad[3] = std::nan("");
    CHECK_FALSE(solveBoundedPeriodicTridiagonal(bad, 1.0, 1.0, solution));
}
