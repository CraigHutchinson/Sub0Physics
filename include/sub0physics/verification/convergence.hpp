#pragma once

#include <cmath>

namespace sub0physics
{

/**
 * Estimates the observed order of accuracy from errors at two step sizes.
 *
 * For a method of order p, error ~ C h^p, so p = log(e_coarse / e_fine) / log(h_coarse / h_fine).
 * Compare the result with the method's nominal order to confirm a solver is in its asymptotic range.
 *
 * @param coarseError Error measured with the larger step; positive.
 * @param fineError Error measured with the smaller step; positive.
 * @param refinementRatio Larger step divided by smaller step; greater than one.
 * @return The observed order; NaN when an argument is outside its contract.
 */
[[nodiscard]] inline double observedOrder(double coarseError, double fineError, double refinementRatio) noexcept
{
    if (!(coarseError > 0.0) || !(fineError > 0.0) || !(refinementRatio > 1.0))
    {
        return std::nan("");
    }
    return std::log(coarseError / fineError) / std::log(refinementRatio);
}

} // namespace sub0physics
