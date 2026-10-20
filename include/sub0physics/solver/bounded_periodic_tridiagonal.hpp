#pragma once

#include <cstddef>
#include <cstdint>
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
    int iterations = 0;           ///< Sweeps performed (the limit when not converged); active-set solves for the active-set method.
    bool converged = false;       ///< True when `residual` is below the tolerance (the active set is stable, for that method).
    std::vector<double> scratch;  ///< Active-set method only: the elimination's multipliers. Untouched by Gauss-Seidel.
    std::vector<std::uint8_t> active; ///< Active-set method only: one where the bound is active.

    /** @param nodes Chain length to reserve for, so that the first Gauss-Seidel solve does not allocate. */
    void reserve(std::size_t nodes)
    {
        u.reserve(nodes);
        gradient.reserve(nodes);
    }

    /** @param nodes Chain length to reserve for, so that the first active-set solve does not allocate either. */
    void reserveActiveSet(std::size_t nodes)
    {
        reserve(nodes);
        scratch.reserve(nodes);
        active.reserve(nodes);
    }
};

/** Stopping rules of the chain solver. */
struct BoundedPeriodicOptions
{
    int maxIterations = 30000;   ///< Sweeps before a solve is declared unconverged.
    double tolerance = 1e-7;     ///< Converged when the scaled projected residual is below this (units of the gradient).
    /**
     * Start from the caller's guess in `solution.u` (raised to the bounds) when it holds a vector of the chain's
     * length, instead of from the bounds. A lever (docs/LOD.md): a guess near the answer, such as the previous
     * solution of a slowly changing problem, needs far fewer sweeps. Off reproduces the cold start bit for bit.
     */
    bool warmStart = false;
};

/** Stopping rule of the active-set chain solver. */
struct BoundedPeriodicActiveSetOptions
{
    int maxIterations = 64; ///< Active-set changes before a solve is declared unconverged.
    /** Take the first active set from the caller's guess in `solution.u` (as BoundedPeriodicOptions::warmStart). */
    bool warmStart = false;
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

/**
 * Minimises the same periodic chain energy as solveBoundedPeriodicTridiagonal, directly, by a primal-dual
 * active-set iteration (Hintermueller, Ito and Kunisch 2002).
 *
 * Each iteration fixes the nodes whose bound is predicted active (`lambda - D (u - lower) > 0`, with `lambda`
 * the reaction) at their bound and solves the free nodes exactly. The free nodes form open chains between
 * active ones, each a symmetric tridiagonal system solved by Thomas' algorithm in O(n); when no node is
 * active the periodic system has no right-hand side, so its unique solution is zero and no periodic
 * (Sherman-Morrison) solve is needed. The Hessian is an M-matrix, so the iteration converges monotonically in
 * finitely many steps and stops when the active set repeats: the KKT conditions then hold to rounding, with no
 * iteration tolerance. Typical counts are 2 to 10 solves, against tens to thousands of Gauss-Seidel sweeps.
 *
 * `residual`, `gradient` and `energy` are measured exactly as the Gauss-Seidel solve measures them, so results
 * compare directly. Once `solution` has capacity for the node count (`reserveActiveSet`) the solve does not
 * allocate.
 *
 * @param lower Lower bounds, 3 or more finite values.
 * @param nodeStiffness `a`; finite and positive.
 * @param couplingStiffness `c`; finite and not negative.
 * @param solution Receives the result; with `options.warmStart`, `u` holds the guess.
 * @param options Iteration limit and warm start.
 * @return Nothing, or an InvalidInput fault for an invalid problem.
 */
[[nodiscard]] std::expected<void, Fault> solveBoundedPeriodicTridiagonalActiveSet(
    std::span<const double> lower, double nodeStiffness, double couplingStiffness, BoundedPeriodicSolution& solution,
    BoundedPeriodicActiveSetOptions options = {});

} // namespace sub0physics
