#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <span>

// Direct solve of a symmetric positive definite periodic block-tridiagonal system (docs/WHEEL.md section 2.3): the
// ring's block of a step matrix. Blocks are Dim x Dim (1 for a scalar chain, 2 for a planar ring, 3 for a spatial
// one), row-major. Block LDL^T in node order; the wrap coupling (node 0 to node N-1) is carried as one fill column
// into the last row, so the cost is linear in N and nothing is allocated: the caller owns every array.

namespace sub0physics
{

/** A Dim x Dim block, row-major. */
template<std::size_t Dim>
using Block = std::array<double, Dim * Dim>;

/** Dim rows of Columns right-hand sides (or solutions), row-major: one node's slice of a multi-column vector. */
template<std::size_t Dim, std::size_t Columns>
using BlockColumns = std::array<double, Dim * Columns>;

namespace detail
{

// c = a b, for Dim x Dim blocks.
template<std::size_t Dim>
[[nodiscard]] constexpr Block<Dim> multiply(const Block<Dim>& a, const Block<Dim>& b) noexcept
{
    Block<Dim> c{};
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t j = 0; j < Dim; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
            {
                sum += a[i * Dim + k] * b[k * Dim + j];
            }
            c[i * Dim + j] = sum;
        }
    }
    return c;
}

// c = a^T b.
template<std::size_t Dim>
[[nodiscard]] constexpr Block<Dim> multiplyTransposed(const Block<Dim>& a, const Block<Dim>& b) noexcept
{
    Block<Dim> c{};
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t j = 0; j < Dim; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
            {
                sum += a[k * Dim + i] * b[k * Dim + j];
            }
            c[i * Dim + j] = sum;
        }
    }
    return c;
}

// The inverse of a symmetric block, or false when it is not positive definite (leading minors) or not finite.
template<std::size_t Dim>
[[nodiscard]] bool invertSymmetricPositive(const Block<Dim>& a, Block<Dim>& inverse) noexcept
{
    static_assert(Dim >= 1 && Dim <= 3, "blocks of 1, 2 or 3 components");
    if constexpr (Dim == 1)
    {
        if (!(a[0] > 0.0) || !std::isfinite(a[0]))
        {
            return false;
        }
        inverse[0] = 1.0 / a[0];
        return true;
    }
    else if constexpr (Dim == 2)
    {
        // Symmetrised: the step matrix is symmetric, so a[1] and a[2] agree up to rounding of their sums.
        const double off = 0.5 * (a[1] + a[2]);
        const double det = a[0] * a[3] - off * off;
        if (!(a[0] > 0.0) || !(det > 0.0) || !std::isfinite(det))
        {
            return false;
        }
        const double r = 1.0 / det;
        inverse = {a[3] * r, -off * r, -off * r, a[0] * r};
        return true;
    }
    else
    {
        const double a01 = 0.5 * (a[1] + a[3]);
        const double a02 = 0.5 * (a[2] + a[6]);
        const double a12 = 0.5 * (a[5] + a[7]);
        const double c00 = a[4] * a[8] - a12 * a12;
        const double c01 = a02 * a12 - a01 * a[8];
        const double c02 = a01 * a12 - a02 * a[4];
        const double minor2 = a[0] * a[4] - a01 * a01;
        const double det = a[0] * c00 + a01 * c01 + a02 * c02;
        if (!(a[0] > 0.0) || !(minor2 > 0.0) || !(det > 0.0) || !std::isfinite(det))
        {
            return false;
        }
        const double r = 1.0 / det;
        const double c11 = a[0] * a[8] - a02 * a02;
        const double c12 = a01 * a02 - a[0] * a12;
        inverse = {c00 * r, c01 * r, c02 * r, c01 * r, c11 * r, c12 * r, c02 * r, c12 * r, minor2 * r};
        return true;
    }
}

// x -= l y for one node's columns.
template<std::size_t Dim, std::size_t Columns>
constexpr void subtractProduct(const Block<Dim>& l, const BlockColumns<Dim, Columns>& y,
                               BlockColumns<Dim, Columns>& x) noexcept
{
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t c = 0; c < Columns; ++c)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
            {
                sum += l[i * Dim + k] * y[k * Columns + c];
            }
            x[i * Columns + c] -= sum;
        }
    }
}

// x -= l^T y.
template<std::size_t Dim, std::size_t Columns>
constexpr void subtractTransposedProduct(const Block<Dim>& l, const BlockColumns<Dim, Columns>& y,
                                         BlockColumns<Dim, Columns>& x) noexcept
{
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t c = 0; c < Columns; ++c)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
            {
                sum += l[k * Dim + i] * y[k * Columns + c];
            }
            x[i * Columns + c] -= sum;
        }
    }
}

// a b for a block and one node's columns.
template<std::size_t Dim, std::size_t Columns>
[[nodiscard]] constexpr BlockColumns<Dim, Columns> applyBlock(const Block<Dim>& a,
                                                              const BlockColumns<Dim, Columns>& b) noexcept
{
    BlockColumns<Dim, Columns> c{};
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t j = 0; j < Columns; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
            {
                sum += a[i * Dim + k] * b[k * Columns + j];
            }
            c[i * Columns + j] = sum;
        }
    }
    return c;
}

} // namespace detail

/**
 * The factor of a periodic block-tridiagonal matrix, A = L D L^T, held in caller-owned arrays of N blocks each.
 *
 * `inversePivot[k]` is D_k^-1. `next[k]` is the multiplier L_{k+1,k} for k < N - 2, and `wrap[k]` the fill multiplier
 * L_{N-1,k} for k <= N - 2; the entries the layout does not use are left untouched.
 *
 * @tparam Dim Block size.
 */
template<std::size_t Dim>
struct PeriodicBlockFactor
{
    std::span<Block<Dim>> inversePivot; ///< N blocks: the inverse pivots.
    std::span<Block<Dim>> next;         ///< N blocks: multipliers onto the next node.
    std::span<Block<Dim>> wrap;         ///< N blocks: multipliers onto the last node (the wrap fill).
};

/**
 * Factorises a symmetric periodic block-tridiagonal matrix: A_kk = `diagonal[k]`, A_{k,k+1 mod N} = `upper[k]` and
 * A_{k+1 mod N,k} = `upper[k]`^T, every other block zero.
 *
 * Block LDL^T in node order with no pivoting, which is stable for a positive definite matrix. The coupling of node 0
 * to node N-1 is carried down the rows as one fill block, so the work is linear in N. Never allocates.
 *
 * @tparam Dim Block size (1, 2 or 3).
 * @param diagonal N symmetric diagonal blocks; N >= 3.
 * @param upper N blocks coupling node k to node k + 1 (the last one couples node N-1 to node 0).
 * @param factor Arrays of N blocks each, filled with the factor.
 * @return True on success; false when a pivot is not positive definite (the matrix is not) or the sizes disagree, in
 *         which case the factor's contents are unspecified.
 */
template<std::size_t Dim>
[[nodiscard]] bool factorPeriodicBlockTridiagonal(std::span<const Block<Dim>> diagonal,
                                                  std::span<const Block<Dim>> upper,
                                                  const PeriodicBlockFactor<Dim>& factor) noexcept
{
    const std::size_t n = diagonal.size();
    if (n < 3U || upper.size() != n || factor.inversePivot.size() != n || factor.next.size() != n ||
        factor.wrap.size() != n)
    {
        return false;
    }
    const std::size_t last = n - 1U;
    Block<Dim> pivot = diagonal[0];
    Block<Dim> lastPivot = diagonal[last];
    // The current row's block in the last column: A_{0,N-1} = upper[N-1]^T to start.
    Block<Dim> fill{};
    for (std::size_t i = 0; i < Dim; ++i)
    {
        for (std::size_t j = 0; j < Dim; ++j)
        {
            fill[i * Dim + j] = upper[last][j * Dim + i];
        }
    }
    for (std::size_t k = 0; k + 1U < last; ++k)
    {
        Block<Dim>& inverse = factor.inversePivot[k];
        if (!detail::invertSymmetricPositive<Dim>(pivot, inverse))
        {
            return false;
        }
        const Block<Dim>& e = upper[k];
        // L_{k+1,k} = E^T D^-1 and L_{N-1,k} = F^T D^-1 (D^-1 symmetric).
        const Block<Dim> nextMultiplier = detail::multiplyTransposed<Dim>(e, inverse);
        const Block<Dim> wrapMultiplier = detail::multiplyTransposed<Dim>(fill, inverse);
        factor.next[k] = nextMultiplier;
        factor.wrap[k] = wrapMultiplier;
        const Block<Dim> nextUpdate = detail::multiply<Dim>(nextMultiplier, e);
        const Block<Dim> lastUpdate = detail::multiply<Dim>(wrapMultiplier, fill);
        const Block<Dim> fillUpdate = detail::multiply<Dim>(nextMultiplier, fill);
        const bool beforeLast = k + 2U == last; // row k + 1 couples to N - 1 through upper[k + 1] as well
        for (std::size_t i = 0; i < Dim * Dim; ++i)
        {
            pivot[i] = diagonal[k + 1U][i] - nextUpdate[i];
            lastPivot[i] -= lastUpdate[i];
            fill[i] = (beforeLast ? upper[k + 1U][i] : 0.0) - fillUpdate[i];
        }
    }
    // Row N - 2: its only coupling ahead is to N - 1, through the fill.
    const std::size_t k = last - 1U;
    if (!detail::invertSymmetricPositive<Dim>(pivot, factor.inversePivot[k]))
    {
        return false;
    }
    factor.wrap[k] = detail::multiplyTransposed<Dim>(fill, factor.inversePivot[k]);
    const Block<Dim> lastUpdate = detail::multiply<Dim>(factor.wrap[k], fill);
    for (std::size_t i = 0; i < Dim * Dim; ++i)
    {
        lastPivot[i] -= lastUpdate[i];
    }
    return detail::invertSymmetricPositive<Dim>(lastPivot, factor.inversePivot[last]);
}

/**
 * Solves A X = B in place with the factor from factorPeriodicBlockTridiagonal, for `Columns` right-hand sides at once.
 * Never allocates.
 *
 * @tparam Dim Block size.
 * @tparam Columns Right-hand sides.
 * @param inversePivot The factor's inverse pivots (PeriodicBlockFactor::inversePivot), N blocks.
 * @param next The factor's next-node multipliers, N blocks.
 * @param wrap The factor's wrap multipliers, N blocks.
 * @param x N node slices: B on entry, X on return.
 */
template<std::size_t Dim, std::size_t Columns>
void solvePeriodicBlockTridiagonal(std::span<const Block<Dim>> inversePivot, std::span<const Block<Dim>> next,
                                   std::span<const Block<Dim>> wrap, std::span<BlockColumns<Dim, Columns>> x) noexcept
{
    struct
    {
        std::span<const Block<Dim>> inversePivot;
        std::span<const Block<Dim>> next;
        std::span<const Block<Dim>> wrap;
    } const factor{inversePivot, next, wrap};
    const std::size_t n = x.size();
    const std::size_t last = n - 1U;
    // L z = b.
    for (std::size_t k = 0; k < last; ++k)
    {
        if (k + 1U < last)
        {
            detail::subtractProduct<Dim, Columns>(factor.next[k], x[k], x[k + 1U]);
        }
        detail::subtractProduct<Dim, Columns>(factor.wrap[k], x[k], x[last]);
    }
    // D L^T x = z.
    x[last] = detail::applyBlock<Dim, Columns>(factor.inversePivot[last], x[last]);
    for (std::size_t k = last; k-- > 0;)
    {
        BlockColumns<Dim, Columns> v = detail::applyBlock<Dim, Columns>(factor.inversePivot[k], x[k]);
        if (k + 1U < last)
        {
            detail::subtractTransposedProduct<Dim, Columns>(factor.next[k], x[k + 1U], v);
        }
        detail::subtractTransposedProduct<Dim, Columns>(factor.wrap[k], x[last], v);
        x[k] = v;
    }
}

} // namespace sub0physics
