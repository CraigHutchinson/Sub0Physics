#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace sub0physics
{

/**
 * How a two-dimensional Euclidean norm is computed: a deterministic arithmetic lever (docs/LOD.md, "Levers").
 * Both are reproducible bit for bit wherever floating-point contraction is off; they differ from each other in
 * the last place.
 */
enum class NormArithmetic : std::uint8_t
{
    /**
     * Scaled by the larger magnitude, squares summed with Kahan compensation, then `sqrt(sum) * largest`: the
     * arithmetic of V8's `Math.hypot`. Independent of the C library, so ported upstream results reproduce bit for
     * bit. About 3x the cost of Fast (a division per argument).
     */
    Compensated,
    /**
     * `sqrt(a * a + b * b)`, falling back to Compensated only when the sum of squares is not a normal number
     * (overflow above about 1.3e154, or underflow below about 1.5e-154), so it is safe for every finite input.
     * Error under one unit in the last place plus the square root's half unit; no library call besides `sqrt`.
     */
    Fast
};

/**
 * Euclidean norm of (a, b) with NormArithmetic::Compensated.
 *
 * @param a First component.
 * @param b Second component.
 * @return The norm; 0 when both are 0. Not for NaN or infinite input.
 */
[[nodiscard]] inline double compensatedNorm(double a, double b) noexcept
{
    const double largest = std::max(std::abs(a), std::abs(b));
    if (largest == 0.0)
    {
        return 0.0;
    }
    double sum = 0.0;
    double compensation = 0.0;
    for (const double v : {a, b})
    {
        const double scaled = std::abs(v) / largest;
        const double summand = scaled * scaled - compensation;
        const double preliminary = sum + summand;
        compensation = (preliminary - sum) - summand;
        sum = preliminary;
    }
    return std::sqrt(sum) * largest;
}

/**
 * Euclidean norm of (a, b) with NormArithmetic::Fast.
 *
 * @param a First component.
 * @param b Second component.
 * @return The norm; 0 when both are 0. Not for NaN or infinite input.
 */
[[nodiscard]] inline double fastNorm(double a, double b) noexcept
{
    const double sum = a * a + b * b;
    // PERF: plain sqrt about 3x faster than the compensated norm.
    if (sum >= std::numeric_limits<double>::min() && sum <= std::numeric_limits<double>::max()) [[likely]]
    {
        return std::sqrt(sum);
    }
    return compensatedNorm(a, b);
}

/**
 * Euclidean norm of (a, b) in the given arithmetic.
 *
 * @param arithmetic Which arithmetic.
 * @param a First component.
 * @param b Second component.
 * @return The norm.
 */
[[nodiscard]] inline double norm(NormArithmetic arithmetic, double a, double b) noexcept
{
    return arithmetic == NormArithmetic::Fast ? fastNorm(a, b) : compensatedNorm(a, b);
}

} // namespace sub0physics
