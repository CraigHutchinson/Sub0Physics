#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <type_traits>
#include <sub0math/linear/vec.hpp>
#include <sub0math/number/det_math.hpp>
#include <sub0math/units/quantities.hpp>
#include <variant>

#include "sub0physics/terrain/grid_spec.hpp"

/**
 * Tier C of the terrain test fixtures: one catalogue of analytic surfaces whose height and unit normal are
 * known exactly. Terrain components are checked against these where no frozen port is the oracle.
 *
 * Every surface is in the terrain frame (x forward, y left, z up, metres, plain `double`), is a value with
 * public parameters, and answers the same two queries. There is one name per surface.
 * - `height(x, y)`: the top surface.
 * - `normal(x, y)`: the unit normal of the top surface, pointing out of the solid. Where two pieces meet,
 *   pieces are half-open on their low side: the point belongs to the piece that starts there.
 *
 * Parameters are a contract, not validated: this is test support.
 */
namespace sub0physics::terrain::testing
{

/** An unlimited plane `z = z0 + gx x + gy y`. */
struct Plane
{
    static constexpr std::string_view Name = "Plane";

    double z0 = 0.0; ///< Height at the origin [m].
    double gx = 0.0; ///< Slope along x [m/m].
    double gy = 0.0; ///< Slope along y [m/m].

    /** @return Height at (x, y). */
    [[nodiscard]] double height(double x, double y) const noexcept
    {
        return z0 + gx * x + gy * y;
    }

    /** @return Unit normal at (x, y); the same everywhere. */
    [[nodiscard]] sub0math::Vec3d normal(double, double) const noexcept
    {
        const double n = std::sqrt(1.0 + gx * gx + gy * gy);
        return {-gx / n, -gy / n, 1.0 / n};
    }
};

/** Flat ground that rises at a constant slope over a length, then flat again at the top; independent of y. */
struct Ramp
{
    static constexpr std::string_view Name = "Ramp";

    double x0 = 1.0;     ///< Where the rise starts [m].
    double length = 2.0; ///< Horizontal length of the rise [m]; positive.
    double rise = 0.5;   ///< Height gained [m].

    /** @return Height at (x, y). */
    [[nodiscard]] double height(double x, double) const noexcept
    {
        return rise * std::clamp((x - x0) / length, 0.0, 1.0);
    }

    /** @return Unit normal at (x, y). */
    [[nodiscard]] sub0math::Vec3d normal(double x, double) const noexcept
    {
        if (!(x >= x0 && x < x0 + length))
        {
            return {0.0, 0.0, 1.0};
        }
        const double slope = rise / length;
        const double n = std::sqrt(1.0 + slope * slope);
        return {-slope / n, 0.0, 1.0 / n};
    }
};

/** A vertical riser at `x0`: ground at 0 before it, flat at `rise` from it on; independent of y. */
struct Step
{
    static constexpr std::string_view Name = "Step";

    double x0 = 1.0;   ///< Position of the riser [m].
    double rise = 0.1; ///< Height of the riser [m].

    /** @return Top height at (x, y): `rise` from the riser on, 0 before. */
    [[nodiscard]] double height(double x, double) const noexcept
    {
        return x >= x0 ? rise : 0.0;
    }

    /** @return Unit normal of the top surface; always up. */
    [[nodiscard]] sub0math::Vec3d normal(double, double) const noexcept
    {
        return {0.0, 0.0, 1.0};
    }

    /** @return Unit normal of the riser face (the set x = x0, 0 <= z <= rise): back towards -x. */
    [[nodiscard]] static constexpr sub0math::Vec3d riserNormal() noexcept
    {
        return {-1.0, 0.0, 0.0};
    }
};

/**
 * A block with a flat top, vertical sides and rounded top edges: the footprint is
 * [x0, x1] x [y0, y1], the top is at `top`, and each top edge is a quarter circle of radius `edgeRadius`
 * (a sphere at the corners). Below the rounding the sides are vertical, so the top surface drops by
 * `top - edgeRadius` at the footprint's edge.
 */
struct RoundedBlock
{
    static constexpr std::string_view Name = "RoundedBlock";

    double x0 = 1.0;         ///< Footprint minimum x [m].
    double x1 = 1.4;         ///< Footprint maximum x [m].
    double y0 = -0.2;        ///< Footprint minimum y [m].
    double y1 = 0.2;         ///< Footprint maximum y [m].
    double top = 0.15;       ///< Height of the flat top [m].
    double edgeRadius = 0.03; ///< Radius of the rounded edges [m]; positive, at most half the smaller footprint side and `top`.

    /** @return Height at (x, y). */
    [[nodiscard]] double height(double x, double y) const noexcept
    {
        const auto [qx, qy] = nearestInner(x, y);
        const double d2 = (x - qx) * (x - qx) + (y - qy) * (y - qy);
        if (d2 == 0.0)
        {
            return top;
        }
        return d2 <= edgeRadius * edgeRadius ? top - edgeRadius + std::sqrt(edgeRadius * edgeRadius - d2) : 0.0;
    }

    /** @return Unit normal at (x, y): up on the top and the ground, radial on the rounding. */
    [[nodiscard]] sub0math::Vec3d normal(double x, double y) const noexcept
    {
        const auto [qx, qy] = nearestInner(x, y);
        const double d2 = (x - qx) * (x - qx) + (y - qy) * (y - qy);
        if (d2 == 0.0 || d2 > edgeRadius * edgeRadius)
        {
            return {0.0, 0.0, 1.0};
        }
        return {(x - qx) / edgeRadius, (y - qy) / edgeRadius, std::sqrt(edgeRadius * edgeRadius - d2) / edgeRadius};
    }

private:
    // Nearest point of the flat top (the footprint shrunk by the radius).
    [[nodiscard]] std::array<double, 2> nearestInner(double x, double y) const noexcept
    {
        return {std::clamp(x, x0 + edgeRadius, x1 - edgeRadius), std::clamp(y, y0 + edgeRadius, y1 - edgeRadius)};
    }
};

/** A solid layer between two heights. */
struct Span
{
    double lower; ///< Underside [m].
    double upper; ///< Top [m].
};

/** The solid layers above one point, from the lowest. */
struct SpanList
{
    std::array<Span, 2> spans; ///< The layers; only the first `count` are meaningful.
    std::size_t count;         ///< How many layers (1 or 2).
};

/**
 * Ground at 0 with a deck over a gap: for x in [x0, x1) the top surface is the deck at `deckTop` and the
 * ground beneath is still there, with open space between. A height field sees only the deck; `spansAt`
 * gives both layers. Independent of y.
 */
struct Bridge
{
    static constexpr std::string_view Name = "Bridge";

    double x0 = 1.0;         ///< Start of the deck [m].
    double x1 = 2.0;         ///< End of the deck [m].
    double deckTop = 0.6;    ///< Top of the deck [m].
    double thickness = 0.1;  ///< Deck thickness [m]; the gap below is `deckTop - thickness`.

    /** @return Top height at (x, y). */
    [[nodiscard]] double height(double x, double) const noexcept
    {
        return onDeck(x) ? deckTop : 0.0;
    }

    /** @return Unit normal of the top surface; always up. */
    [[nodiscard]] sub0math::Vec3d normal(double, double) const noexcept
    {
        return {0.0, 0.0, 1.0};
    }

    /** @return Unit normal of the deck's underside: down. */
    [[nodiscard]] static constexpr sub0math::Vec3d undersideNormal() noexcept
    {
        return {0.0, 0.0, -1.0};
    }

    /** @return The solid layers above (x, y): the ground, and the deck where there is one. */
    [[nodiscard]] SpanList spansAt(double x, double) const noexcept
    {
        constexpr double Deep = -std::numeric_limits<double>::infinity();
        if (!onDeck(x))
        {
            return {{Span{Deep, 0.0}, Span{0.0, 0.0}}, 1U};
        }
        return {{Span{Deep, 0.0}, Span{deckTop - thickness, deckTop}}, 2U};
    }

private:
    [[nodiscard]] bool onDeck(double x) const noexcept
    {
        return x >= x0 && x < x1;
    }
};

/**
 * A banked turn: flat for y < y0, then a circular arc of radius `radius` that curves up to the angle
 * `maxAngle`, then a straight bank at that angle. The profile is smooth (C1); independent of x. The bank
 * rises towards +y.
 */
struct Berm
{
    static constexpr std::string_view Name = "Berm";

    double y0 = 0.0;         ///< Where the arc starts [m].
    double radius = 2.0;     ///< Radius of the arc [m]; positive.
    double maxAngle = 0.6;   ///< Bank angle at the end of the arc [rad]; between 0 and pi/2.

    /** @return Height at (x, y). */
    [[nodiscard]] double height(double, double y) const noexcept
    {
        const double u = y - y0;
        if (!(u > 0.0))
        {
            return 0.0;
        }
        const double arcEnd = radius * sub0math::det::sin(maxAngle);
        if (u <= arcEnd)
        {
            return radius - std::sqrt(radius * radius - u * u);
        }
        return radius * (1.0 - sub0math::det::cos(maxAngle)) + std::tan(maxAngle) * (u - arcEnd);
    }

    /** @return Unit normal at (x, y). */
    [[nodiscard]] sub0math::Vec3d normal(double, double y) const noexcept
    {
        const double u = y - y0;
        if (!(u > 0.0))
        {
            return {0.0, 0.0, 1.0};
        }
        const double arcEnd = radius * sub0math::det::sin(maxAngle);
        const double slope = u <= arcEnd ? u / std::sqrt(radius * radius - u * u) : std::tan(maxAngle);
        const double n = std::sqrt(1.0 + slope * slope);
        return {0.0, -slope / n, 1.0 / n};
    }
};

/** A parabolic trench along x, centred at `yc`: depth `depth` at the centre, flat ground beyond `halfWidth`. */
struct Rut
{
    static constexpr std::string_view Name = "Rut";

    double yc = 0.0;         ///< Centre line [m].
    double halfWidth = 0.1;  ///< Half width of the trench [m]; positive.
    double depth = 0.05;     ///< Depth at the centre [m].

    /** @return Height at (x, y): 0 outside the trench, down to `-depth` at its centre. */
    [[nodiscard]] double height(double, double y) const noexcept
    {
        const double u = (y - yc) / halfWidth;
        return std::abs(u) < 1.0 ? -depth * (1.0 - u * u) : 0.0;
    }

    /** @return Unit normal at (x, y). */
    [[nodiscard]] sub0math::Vec3d normal(double, double y) const noexcept
    {
        const double u = (y - yc) / halfWidth;
        if (!(std::abs(u) < 1.0))
        {
            return {0.0, 0.0, 1.0};
        }
        const double slope = 2.0 * depth * u / halfWidth;
        const double n = std::sqrt(1.0 + slope * slope);
        return {0.0, -slope / n, 1.0 / n};
    }
};

/** What every catalogue surface provides. */
template<typename S>
concept AnalyticSurface = requires(const S& s, double x, double y) {
    { S::Name } -> std::convertible_to<std::string_view>;
    { s.height(x, y) } -> std::same_as<double>;
    { s.normal(x, y) } -> std::same_as<sub0math::Vec3d>;
};

static_assert(AnalyticSurface<Plane> && AnalyticSurface<Ramp> && AnalyticSurface<Step> && AnalyticSurface<RoundedBlock> &&
              AnalyticSurface<Bridge> && AnalyticSurface<Berm> && AnalyticSurface<Rut>);

/** Any surface of the catalogue. */
using AnySurface = std::variant<Plane, Ramp, Step, RoundedBlock, Bridge, Berm, Rut>;

/** The names of the catalogue's surfaces, in `AnySurface` order. */
inline constexpr std::array<std::string_view, 7> SurfaceNames{Plane::Name,  Ramp::Name, Step::Name, RoundedBlock::Name,
                                                              Bridge::Name, Berm::Name, Rut::Name};

/** @return One instance of each surface with its default parameters, in `AnySurface` order. */
[[nodiscard]] inline std::array<AnySurface, 7> defaultSurfaces() noexcept
{
    return {Plane{0.3, 0.1, -0.05}, Ramp{}, Step{}, RoundedBlock{}, Bridge{}, Berm{}, Rut{}};
}

/**
 * Finds a surface's name.
 *
 * @param s Any surface.
 * @return Its catalogue name.
 */
[[nodiscard]] inline std::string_view nameOf(const AnySurface& s) noexcept
{
    return std::visit([](const auto& surface) { return std::remove_cvref_t<decltype(surface)>::Name; }, s);
}

/**
 * Finds a surface's height.
 *
 * @param s Any surface.
 * @param x Terrain x [m].
 * @param y Terrain y [m].
 * @return Its top height.
 */
[[nodiscard]] inline double heightOf(const AnySurface& s, double x, double y) noexcept
{
    return std::visit([=](const auto& surface) { return surface.height(x, y); }, s);
}

/**
 * Finds a surface's unit normal.
 *
 * @param s Any surface.
 * @param x Terrain x [m].
 * @param y Terrain y [m].
 * @return The unit normal of its top.
 */
[[nodiscard]] inline sub0math::Vec3d normalOf(const AnySurface& s, double x, double y) noexcept
{
    return std::visit([=](const auto& surface) { return surface.normal(x, y); }, s);
}

/**
 * Samples a surface at a page's sample positions (nodes or cell centres, as the grid is registered).
 *
 * @param grid The grid.
 * @param key The page.
 * @param surface The surface.
 * @param heights Output, `samplesPerSide()^2` heights in the terrain frame, row-major in j.
 * @param known Output, the same count: set to 1 everywhere (a catalogue surface is defined everywhere).
 */
template<AnalyticSurface S>
void sampleSurface(const GridSpec& grid, const TileKey& key, const S& surface, std::span<sub0math::f64::Length> heights,
                   std::span<std::uint8_t> known)
{
    const std::int32_t n = grid.samplesPerSide();
    for (std::int32_t j = 0; j < n; ++j)
    {
        for (std::int32_t i = 0; i < n; ++i)
        {
            const auto p = grid.samplePosition(key, i, j);
            const auto at = static_cast<std::size_t>(j) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
            heights[at] = sub0math::f64::Length::fromSi(surface.height(p.x.si(), p.y.si()));
            known[at] = 1U;
        }
    }
}

} // namespace sub0physics::terrain::testing
