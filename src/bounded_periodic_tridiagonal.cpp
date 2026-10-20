#include <algorithm>
#include <cmath>
#include <limits>

#include "sub0physics/solver/bounded_periodic_tridiagonal.hpp"

namespace sub0physics
{

std::expected<void, Fault> solveBoundedPeriodicTridiagonal(std::span<const double> lower, double nodeStiffness,
                                                           double couplingStiffness, BoundedPeriodicSolution& solution,
                                                           BoundedPeriodicOptions options)
{
    const std::size_t n = lower.size();
    bool valid = n >= 3U && std::isfinite(nodeStiffness) && nodeStiffness > 0.0 && std::isfinite(couplingStiffness)
                 && couplingStiffness >= 0.0 && options.maxIterations >= 1 && std::isfinite(options.tolerance)
                 && options.tolerance > 0.0;
    for (std::size_t i = 0; valid && i < n; ++i)
    {
        valid = std::isfinite(lower[i]);
    }
    if (!valid)
    {
        return invalidInput("bounded periodic chain needs 3 or more finite bounds, a > 0, c >= 0 and a positive limit and tolerance");
    }
    const double a = nodeStiffness;
    const double c = couplingStiffness;
    const double diagonal = a + 2.0 * c;
    solution.u.assign(lower.begin(), lower.end());
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
    double energy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double step = u[i == last ? 0U : i + 1U] - u[i];
        energy += 0.5 * a * (u[i] * u[i]) + 0.5 * c * (step * step);
    }
    solution.energy = energy;
    solution.residual = residual;
    solution.iterations = converged ? iterations + 1 : options.maxIterations;
    solution.converged = converged;
    return {};
}

} // namespace sub0physics
