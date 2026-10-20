#include <algorithm>
#include <cmath>
#include <limits>

#include "sub0physics/solver/bounded_periodic_tridiagonal.hpp"

namespace sub0physics
{

namespace
{

bool validProblem(std::span<const double> lower, double a, double c) noexcept
{
    bool valid = lower.size() >= 3U && std::isfinite(a) && a > 0.0 && std::isfinite(c) && c >= 0.0;
    for (std::size_t i = 0; valid && i < lower.size(); ++i)
    {
        valid = std::isfinite(lower[i]);
    }
    return valid;
}

// E(u) = sum a/2 u^2 + c/2 (u[i+1] - u[i])^2, in the reference solve's operand order.
double chainEnergy(const double* u, std::size_t n, double a, double c) noexcept
{
    const std::size_t last = n - 1U;
    double energy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double step = u[i == last ? 0U : i + 1U] - u[i];
        energy += 0.5 * a * (u[i] * u[i]) + 0.5 * c * (step * step);
    }
    return energy;
}

// Starts u at the bounds, or at the caller's guess raised to them.
void startAt(std::span<const double> lower, bool warmStart, std::vector<double>& u)
{
    if (warmStart && u.size() == lower.size())
    {
        for (std::size_t i = 0; i < u.size(); ++i)
        {
            u[i] = u[i] > lower[i] ? u[i] : lower[i];
        }
        return;
    }
    u.assign(lower.begin(), lower.end());
}

} // namespace

std::expected<void, Fault> solveBoundedPeriodicTridiagonal(std::span<const double> lower, double nodeStiffness,
                                                           double couplingStiffness, BoundedPeriodicSolution& solution,
                                                           BoundedPeriodicOptions options)
{
    const std::size_t n = lower.size();
    const bool valid = validProblem(lower, nodeStiffness, couplingStiffness) && options.maxIterations >= 1
                       && std::isfinite(options.tolerance) && options.tolerance > 0.0;
    if (!valid)
    {
        return invalidInput("bounded periodic chain needs 3 or more finite bounds, a > 0, c >= 0 and a positive limit and tolerance");
    }
    const double a = nodeStiffness;
    const double c = couplingStiffness;
    const double diagonal = a + 2.0 * c;
    startAt(lower, options.warmStart, solution.u);
    solution.gradient.assign(n, 0.0);
    double* const u = solution.u.data();
    double* const gradient = solution.gradient.data();
    const double* const lo = lower.data();
    const std::size_t last = n - 1U;

    double residual = std::numeric_limits<double>::infinity();
    int iterations = 0;
    bool converged = false;
    for (; iterations < options.maxIterations; ++iterations)
    {
        // Cyclic sweep in node order: node i sees its already updated predecessor and the old successor.
        // PERF: peeled wrap nodes faster than modulo indices.
        const auto relax = [&](std::size_t i, double previous, double next) {
            const double relaxed = c * (previous + next) / diagonal;
            u[i] = relaxed > lo[i] ? relaxed : lo[i];
        };
        relax(0U, u[last], u[1]);
        for (std::size_t i = 1; i < last; ++i)
        {
            relax(i, u[i - 1U], u[i + 1U]);
        }
        relax(last, u[last - 1U], u[0]);

        residual = 0.0;
        const auto measure = [&](std::size_t i, double previous, double next) {
            gradient[i] = diagonal * u[i] - c * (previous + next);
            const double projected = u[i] - gradient[i] / diagonal;
            residual = std::max(residual, std::abs(u[i] - (projected > lo[i] ? projected : lo[i])) * diagonal);
        };
        measure(0U, u[last], u[1]);
        for (std::size_t i = 1; i < last; ++i)
        {
            measure(i, u[i - 1U], u[i + 1U]);
        }
        measure(last, u[last - 1U], u[0]);
        if (residual < options.tolerance)
        {
            converged = true;
            break;
        }
    }
    solution.energy = chainEnergy(u, n, a, c);
    solution.residual = residual;
    solution.iterations = converged ? iterations + 1 : options.maxIterations;
    solution.converged = converged;
    return {};
}

std::expected<void, Fault> solveBoundedPeriodicTridiagonalActiveSet(std::span<const double> lower, double nodeStiffness,
                                                                    double couplingStiffness,
                                                                    BoundedPeriodicSolution& solution,
                                                                    BoundedPeriodicActiveSetOptions options)
{
    const std::size_t n = lower.size();
    if (!validProblem(lower, nodeStiffness, couplingStiffness) || options.maxIterations < 1)
    {
        return invalidInput("bounded periodic chain needs 3 or more finite bounds, a > 0, c >= 0 and a positive limit");
    }
    const double a = nodeStiffness;
    const double c = couplingStiffness;
    const double diagonal = a + 2.0 * c;
    startAt(lower, options.warmStart, solution.u);
    solution.gradient.resize(n);
    solution.scratch.resize(n);
    solution.active.assign(n, 0U);
    double* const u = solution.u.data();
    double* const lambda = solution.gradient.data();
    double* const multiplier = solution.scratch.data();
    std::uint8_t* const active = solution.active.data();
    const double* const lo = lower.data();
    const std::size_t last = n - 1U;
    const auto previousOf = [last](std::size_t i) { return i == 0U ? last : i - 1U; };
    const auto nextOf = [last](std::size_t i) { return i == last ? 0U : i + 1U; };
    // Reaction of node i: dE/du = D u[i] - c (previous + next).
    const auto reaction = [&](std::size_t i) { return diagonal * u[i] - c * (u[previousOf(i)] + u[nextOf(i)]); };

    // The first prediction uses the reaction of the starting point everywhere.
    for (std::size_t i = 0; i < n; ++i)
    {
        lambda[i] = reaction(i);
    }
    int iterations = 0;
    bool converged = false;
    for (; iterations < options.maxIterations; ++iterations)
    {
        bool changed = iterations == 0;
        std::size_t anchor = n;
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::uint8_t predicted = lambda[i] - diagonal * (u[i] - lo[i]) > 0.0 ? 1U : 0U;
            changed = changed || predicted != active[i];
            active[i] = predicted;
            if (predicted != 0U && anchor == n)
            {
                anchor = i;
            }
        }
        if (!changed)
        {
            converged = true;
            break;
        }
        if (anchor == n)
        {
            // No active bound: H u = 0 with H positive definite, so u = 0.
            std::fill(u, u + n, 0.0);
            std::fill(lambda, lambda + n, 0.0);
            continue;
        }
        // Walk the ring once from an active node. Each run of free nodes is an open chain between active neighbours,
        // D u[k] - c (u[k-1] + u[k+1]) = 0 with the neighbours at their bounds: Thomas' algorithm, forward then back.
        std::size_t i = anchor;
        for (std::size_t visited = 0; visited < n; ++visited)
        {
            const std::size_t k = nextOf(i);
            if (active[i] != 0U)
            {
                u[i] = lo[i];
                if (active[k] == 0U)
                {
                    // Forward elimination of the free run starting at k; u holds the eliminated right-hand side.
                    double carriedMultiplier = 0.0;
                    double carriedValue = c * lo[i];
                    std::size_t j = k;
                    std::size_t length = 0;
                    while (active[j] == 0U)
                    {
                        const double pivot = diagonal - c * carriedMultiplier;
                        const double rhs = carriedValue + (active[nextOf(j)] != 0U ? c * lo[nextOf(j)] : 0.0);
                        multiplier[j] = c / pivot;
                        u[j] = rhs / pivot;
                        carriedMultiplier = multiplier[j];
                        carriedValue = c * u[j];
                        j = nextOf(j);
                        ++length;
                    }
                    // Back substitution from the run's last node.
                    std::size_t back = previousOf(j);
                    for (std::size_t m = 1; m < length; ++m)
                    {
                        const std::size_t before = previousOf(back);
                        u[before] += multiplier[before] * u[back];
                        back = before;
                    }
                }
            }
            i = k;
        }
        for (std::size_t m = 0; m < n; ++m)
        {
            lambda[m] = active[m] != 0U ? reaction(m) : 0.0;
        }
    }

    // Measured as the Gauss-Seidel solve measures: the gradient everywhere and the scaled projected residual.
    double residual = 0.0;
    for (std::size_t m = 0; m < n; ++m)
    {
        lambda[m] = reaction(m);
        const double projected = u[m] - lambda[m] / diagonal;
        residual = std::max(residual, std::abs(u[m] - (projected > lo[m] ? projected : lo[m])) * diagonal);
    }
    solution.energy = chainEnergy(u, n, a, c);
    solution.residual = residual;
    solution.iterations = converged ? iterations : options.maxIterations;
    solution.converged = converged;
    return {};
}

} // namespace sub0physics
