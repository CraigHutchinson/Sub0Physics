#pragma once

#include <cstddef>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Solution of a bound-constrained periodic chain problem. Reuse one object across solves: the solver
 * allocates only when the chain is longer than any it has held before.
 */
struct BoundedPeriodicSolution
{
    std::vector<double> u;        ///< The minimiser: at least its lower bound at every node.
    std::vector<double> gradient; ///< `dE/du` at the minimiser; positive where the bound is active (the reaction).
    double energy = 0.0;          ///< `E(u)`.
    double residual = 0.0;        ///< Largest diagonally scaled projected residual `|u - max(lo, u - g / D)| D`.
    int iterations = 0;           ///< Sweeps performed (the limit when not converged).
    bool converged = false;       ///< True when `residual` is below the tolerance.

    /** @param nodes Chain length to reserve for, so that the first solve does not allocate. */
    void reserve(std::size_t nodes)
    {
        u.reserve(nodes);
        gradient.reserve(nodes);
    }
};

/** Stopping rules of the chain solver. */
struct BoundedPeriodicOptions
{
    int maxIterations = 30000;   ///< Sweeps before a solve is declared unconverged.
    double tolerance = 1e-7;     ///< Converged when the scaled projected residual is below this (units of the gradient).
};

/**
 * Minimises a periodic chain energy under lower bounds.
 *
 * The energy is `E(u) = sum_i [ a/2 u[i]^2 + c/2 (u[i+1] - u[i])^2 ]` over a ring of nodes (the last neighbours
 * the first), minimised subject to `u >= lower`. With `a > 0` and `c >= 0` it is a strictly convex quadratic
 * whose Hessian is periodic tridiagonal: diagonal `D = a + 2c`, both neighbours `-c`. The solver is cyclic
 * exact coordinate minimisation (projected Gauss-Seidel) in node order, each node seeing its already updated
 * predecessor and the old successor, started at the lower bound; the wrap nodes are peeled out of the loop
 * rather than indexed modulo. The arithmetic is `relaxed = c (previous + next) / D`, `u = max(relaxed, lower)`,
 * so results are reproducible bit for bit wherever floating-point contraction is off.
 *
 * Convergence is the diagonally scaled projected residual: converged when it falls below
 * `options.tolerance`. At convergence the KKT conditions hold to that tolerance: `gradient >= 0`, `u >= lower`
 * and `gradient * (u - lower) = 0`.
 *
 * Once `solution` has capacity for the node count the solve does not allocate.
 *
 * @param lower Lower bounds, 3 or more finite values.
 * @param nodeStiffness `a`, the stiffness tying each node to zero; finite and positive.
 * @param couplingStiffness `c`, the stiffness between neighbours; finite and not negative.
 * @param solution Receives the result; its previous contents are replaced.
 * @param options Sweep limit and tolerance; the defaults are the ones the ported goldens use.
 * @return Nothing, or an InvalidInput fault for an invalid problem.
 */
[[nodiscard]] std::expected<void, Fault> solveBoundedPeriodicTridiagonal(std::span<const double> lower,
                                                                         double nodeStiffness,
                                                                         double couplingStiffness,
                                                                         BoundedPeriodicSolution& solution,
                                                                         BoundedPeriodicOptions options = {});

} // namespace sub0physics
