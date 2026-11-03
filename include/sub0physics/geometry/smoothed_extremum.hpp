#pragma once

#include <cmath>

// Smoothed extrema of two values (docs/WHEEL.md, "M3d status"): the explicit ring's smoothed tread maximum
// (multibody/casing_models.hpp) and the continuous edge law's soft minimum (ConvexTerrain::continuousEdgeGap), each
// written once here. Plain arithmetic, no allocation, no transcendentals.

namespace sub0physics
{

/** A smoothed extremum of two values with its exact first and second derivatives. */
struct SmoothedExtremum
{
    double value = 0.0;     ///< The smoothed extremum.
    double weight0 = 0.5;   ///< d value / d first; in [0, 1].
    double weight1 = 0.5;   ///< d value / d second, `1 - weight0` (both 0.5 at a tie).
    double curvature = 0.0; ///< d2 value / d first^2 (= d2 / d second^2 = -d2 / d first d second): the weights' rate.
};

/**
 * The maximum of two values smoothed within a band: `m = (a + b) / 2 + h(|a - b| / 2)` with `h(x) = x` for
 * `x >= band` and `h(x) = x^2 (2 band - x) / band^2` within it. `h` and `h'` are continuous, so `m` is C1 (C2 at the
 * tie), equals `max(a, b)` exactly outside the band and at the tie, and is never above it. The weights are
 * `dm/da = (1 + h' s) / 2`, `dm/db = (1 - h' s) / 2` with `s` the sign of `a - b`: continuous, summing to one.
 *
 * @param a The first value.
 * @param b The second value.
 * @param band The band's half-width in `|a - b| / 2`; positive.
 * @return The smoothed maximum and its derivatives.
 */
[[nodiscard]] inline SmoothedExtremum smoothedMaximum(double a, double b, double band) noexcept
{
    const double half = 0.5 * (a - b);
    const double x = std::abs(half);
    const double sign = half >= 0.0 ? 1.0 : -1.0;
    double h = x;
    double slope = 1.0;
    double bend = 0.0;
    if (x < band)
    {
        h = x * x * (2.0 * band - x) / (band * band);
        slope = x * (4.0 * band - 3.0 * x) / (band * band);
        bend = (4.0 * band - 6.0 * x) / (band * band);
    }
    SmoothedExtremum out;
    out.value = 0.5 * (a + b) + h;
    out.weight0 = 0.5 + 0.5 * slope * sign;
    out.weight1 = 0.5 - 0.5 * slope * sign;
    out.curvature = 0.25 * bend;
    return out;
}

/**
 * The minimum of two values smoothed within a band, convex and monotone: `m = (a + b) / 2 - h(|a - b| / 2)` with
 * `h(x) = x` for `x >= band` and `h(x) = (x^2 + band^2) / (2 band)` within it. `m` is C1 (C-infinity inside the band),
 * equals `min(a, b)` exactly outside the band, and is never above it: `band / 2` below at the tie. The weights
 * `dm/da = (1 - h' s) / 2`, `dm/db = (1 + h' s) / 2` (`s` the sign of `a - b`, `h' = x / band`) stay in [0, 1] and
 * move linearly from one value to the other across the band.
 *
 * This is the continuous edge law's kernel (ConvexTerrain::continuousEdgeGap). `smoothedMaximum` keeps the extremum
 * exactly at the tie, which forces its weights outside [0, 1] within the band (a smoothing that is never beyond the
 * extremum, C1 and exact outside the band has `h(0) <= 0` and `h(band) = band`, so `h'` passes 1): fine for the tread's
 * force law, but a row built on it would pull a node towards the surface, so the rows use this one.
 *
 * @param a The first value.
 * @param b The second value.
 * @param band The band's half-width in `|a - b| / 2`; positive.
 * @return The smoothed minimum and its derivatives (weights `d/da`, `d/db`; curvature `d2/da2`, never positive).
 */
[[nodiscard]] inline SmoothedExtremum softMinimum(double a, double b, double band) noexcept
{
    const double half = 0.5 * (a - b);
    const double x = std::abs(half);
    const double sign = half >= 0.0 ? 1.0 : -1.0;
    double h = x;
    double slope = 1.0;
    double bend = 0.0;
    if (x < band)
    {
        h = (x * x + band * band) / (2.0 * band);
        slope = x / band;
        bend = 1.0 / band;
    }
    SmoothedExtremum out;
    out.value = 0.5 * (a + b) - h;
    out.weight0 = 0.5 - 0.5 * slope * sign;
    out.weight1 = 0.5 + 0.5 * slope * sign;
    out.curvature = -0.25 * bend;
    return out;
}

} // namespace sub0physics
