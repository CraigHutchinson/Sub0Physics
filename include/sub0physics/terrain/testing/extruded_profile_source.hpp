#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "sub0physics/fault.hpp"
#include "sub0physics/geometry/circle_clearance.hpp"
#include "sub0physics/geometry/polyline_profile.hpp"
#include "sub0physics/terrain/feature.hpp"
#include "sub0physics/terrain/grid_spec.hpp"
#include "sub0physics/terrain/height_page.hpp"
#include "sub0physics/terrain/source_tile_cache.hpp"

/**
 * Tier B of the terrain test fixtures: adapters that feed an existing scenario into the generic components.
 * `ExtrudedProfileSource` turns a 2D ground profile (the form every wheel-terrain fixture takes) into a
 * terrain source: a surface that does not vary along y, served as height pages with its vertical faces as
 * exact features.
 */
namespace sub0physics::terrain::testing
{

/** A tile of `ExtrudedProfileSource`: a baked page, once loaded. */
struct BakedPage
{
    std::optional<HeightPage> page; ///< Empty until the source has built it.
};

/**
 * A `PolylineProfile` (x along the ground, height up, solid below) extruded along y over `[y0, y1]`.
 *
 * Three things hold, and the tests pin them:
 * - **The slice is the profile.** `clearanceSlice` at any y inside the extrusion is *exactly*
 *   `measureCircleClearance` on the profile: the same call, so the same bits.
 * - **The pages are samples of it.** `load` bakes a page at the grid's sample positions with
 *   `PolylineProfile::surfaceHeight`; a node outside the profile's x range or the extrusion is unknown, never
 *   extrapolated.
 * - **Vertical faces are features.** Each vertical profile segment is a `BoxFeature` of zero thickness (a
 *   face rectangle), because the pages' triangulation turns a riser into a ramp across one cell.
 *
 * Tests only: `load` allocates and may be called from `SourceTileCache::pump`.
 */
class ExtrudedProfileSource final : public TileSource<BakedPage>
{
public:
    /**
     * Builds a source.
     *
     * @param profile The ground profile.
     * @param y0 Lower y of the extrusion [m].
     * @param y1 Upper y of the extrusion [m]; above `y0`.
     * @param grid The grid the pages are baked on (256 cells, node-registered).
     * @return The source, or an InvalidInput fault for a bad extrusion or grid.
     */
    [[nodiscard]] static std::expected<ExtrudedProfileSource, Fault> create(PolylineProfile profile, sub0math::f64::Length y0,
                                                                            sub0math::f64::Length y1, const GridSpec& grid)
    {
        if (!(y0.si() < y1.si()) || !std::isfinite(y0.si()) || !std::isfinite(y1.si()))
        {
            return invalidInput("an extrusion needs a finite y0 below y1");
        }
        if (grid.cellsPerSide() != HeightPage::Cells || grid.registration() != Registration::NodeIsPoint)
        {
            return invalidInput("an extruded profile bakes onto 256-cell node-registered pages");
        }
        return ExtrudedProfileSource(std::move(profile), y0, y1, grid);
    }

    /** @return The profile. */
    [[nodiscard]] const PolylineProfile& profile() const noexcept
    {
        return profile_;
    }

    /** @return The grid the pages are baked on. */
    [[nodiscard]] const GridSpec& grid() const noexcept
    {
        return grid_;
    }

    /** @return The vertical faces as zero-thickness boxes, in profile order. */
    [[nodiscard]] std::span<const BoxFeature> features() const noexcept
    {
        return features_;
    }

    /**
     * Gives the exact surface height at a point.
     *
     * @param x Terrain x [m].
     * @param y Terrain y [m].
     * @return The height, or nothing outside the profile's x range or the extrusion.
     */
    [[nodiscard]] std::optional<double> height(double x, double y) const noexcept
    {
        if (!(x >= profile_.minX().si() && x <= profile_.maxX().si() && y >= y0_ && y <= y1_))
        {
            return std::nullopt;
        }
        return profile_.surfaceHeight(sub0math::f64::Length::fromSi(x)).si();
    }

    /**
     * Measures a circle's clearance in the plane y = const: the profile's own clearance.
     *
     * @param y The slice's y [m].
     * @param query The circle (centre x, height) and tie tolerance.
     * @param witnesses Caller's buffer, as `measureCircleClearance`.
     * @return The clearance; OutsideDomain when the slice is outside the extrusion.
     */
    [[nodiscard]] CircleClearance clearanceSlice(double y, const CircleQuery& query, std::span<Witness2> witnesses) const noexcept
    {
        if (!(y >= y0_ && y <= y1_))
        {
            return {ClearanceStatus::OutsideDomain, {}, {}, 0, false};
        }
        return measureCircleClearance(profile_, query, witnesses);
    }

    /**
     * Bakes a page: heights at the grid's sample positions, unknown outside the surface.
     *
     * @param key The page.
     * @param slot Receives the page.
     * @return Ready, or Failed for a key outside the grid.
     */
    [[nodiscard]] LoadResult load(const TileKey& key, BakedPage& slot) noexcept override
    {
        const auto n = static_cast<std::size_t>(HeightPage::Nodes);
        std::vector<sub0math::f64::Length> heights(n * n);
        std::vector<std::uint8_t> known(n * n, 0U);
        for (std::size_t j = 0; j < n; ++j)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                const auto p = grid_.samplePosition(key, static_cast<std::int32_t>(i), static_cast<std::int32_t>(j));
                if (const auto z = height(p.x.si(), p.y.si()))
                {
                    heights[j * n + i] = sub0math::f64::Length::fromSi(*z);
                    known[j * n + i] = 1U;
                }
            }
        }
        auto page = buildHeightPage(grid_, key, SourceRevision{1}, heights, known);
        if (!page)
        {
            return {LoadState::Failed, page.error().detail};
        }
        slot.page.emplace(std::move(*page));
        return {LoadState::Ready, {}};
    }

private:
    ExtrudedProfileSource(PolylineProfile profile, sub0math::f64::Length y0, sub0math::f64::Length y1, const GridSpec& grid)
        : profile_(std::move(profile)), y0_(y0.si()), y1_(y1.si()), grid_(grid)
    {
        const auto points = profile_.points();
        std::uint32_t id = 0;
        for (std::size_t i = 0; i + 1U < points.size(); ++i)
        {
            if (points[i].x.si() != points[i + 1U].x.si())
            {
                continue;
            }
            const double z0 = points[i].y.si(), z1 = points[i + 1U].y.si();
            const double x = points[i].x.si();
            // The solid is below the profile, so the air is on the side the profile comes from when it rises.
            const double nx = z1 > z0 ? -1.0 : 1.0;
            features_.push_back({FeatureId{id++},
                                 {{sub0math::f64::Length::fromSi(x), sub0math::f64::Length::fromSi(y0_), sub0math::f64::Length::fromSi(std::min(z0, z1))},
                                  {sub0math::f64::Length::fromSi(x), sub0math::f64::Length::fromSi(y1_), sub0math::f64::Length::fromSi(std::max(z0, z1))}},
                                 {nx, 0.0, 0.0}});
        }
    }

    PolylineProfile profile_;
    double y0_;
    double y1_;
    GridSpec grid_;
    std::vector<BoxFeature> features_;
};

/**
 * Rebuilds a profile from one row of a baked page plus the exact vertical faces, so the baked representation
 * can be asked the same clearance questions as the original.
 *
 * Authored nodes become vertices; a feature face inside a node interval inserts its two end points, so a riser
 * is vertical again. Nodes that are not authored are skipped, and so the rebuilt profile covers only the
 * authored x range of the row.
 *
 * @param page The baked page.
 * @param row Node row, 0 to Nodes - 1.
 * @param features The source's features.
 * @return The profile, or the profile rules' fault.
 */
[[nodiscard]] inline std::expected<PolylineProfile, Fault> profileFromPageRow(const GridSpec& grid, const HeightPage& page,
                                                                              std::int32_t row,
                                                                              std::span<const BoxFeature> features)
{
    std::vector<ProfilePoint> points;
    const auto add = [&](double x, double z) {
        if (!points.empty() && std::abs(points.back().x.si() - x) < 1e-9 && std::abs(points.back().y.si() - z) < 1e-9)
        {
            return;
        }
        points.push_back({sub0math::f64::Length::fromSi(x), sub0math::f64::Length::fromSi(z)});
    };
    const auto& frame = page.frame();
    const double spacing = grid.spacing(frame.key.level).si(); // the page holds it as a float
    const double x0 = frame.origin.x.si();
    const double z0 = frame.origin.z.si();
    std::optional<double> previousX;
    for (std::int32_t i = 0; i < HeightPage::Nodes; ++i)
    {
        if (!page.known(i, row))
        {
            previousX.reset();
            continue;
        }
        const double x = x0 + spacing * static_cast<double>(i);
        const double z = z0 + static_cast<double>(page.nodeHeight(i, row).si());
        bool onFace = false;
        if (previousX)
        {
            for (const BoxFeature& f : features)
            {
                const double fx = f.bounds.lo.x.si();
                if (fx > *previousX + 1e-9 && fx <= x + 1e-9)
                {
                    // Rising faces (normal -x) go from the lower end to the upper; falling ones the other way.
                    const bool rising = f.faceNormal.x < 0.0;
                    add(fx, rising ? f.bounds.lo.z.si() : f.bounds.hi.z.si());
                    add(fx, rising ? f.bounds.hi.z.si() : f.bounds.lo.z.si());
                    onFace = onFace || std::abs(fx - x) < 1e-9;
                }
            }
        }
        // A node exactly on a face holds the face's top; the face points above already say it.
        if (!onFace)
        {
            add(x, z);
        }
        previousX = x;
    }
    return PolylineProfile::create(std::move(points));
}

/**
 * Samples the Hausdorff distance between two profiles inside an x range: every segment of each is sampled at a
 * fixed step, and each sample's distance to the other profile (the nearest point of any of its segments) is
 * taken; the result is the largest. The distance from a point to a set moves by at most this much between the
 * two profiles, so it bounds how far any clearance query can differ.
 *
 * @param a One profile.
 * @param b The other.
 * @param xMin Samples left of this x are ignored [m].
 * @param xMax Samples right of this x are ignored [m].
 * @param step Sampling step along each segment [m].
 * @return The largest sampled distance [m].
 */
[[nodiscard]] inline double sampledHausdorff(const PolylineProfile& a, const PolylineProfile& b, double xMin, double xMax,
                                             double step = 0.00025)
{
    const auto distanceTo = [](const PolylineProfile& profile, double x, double z) {
        double best = std::numeric_limits<double>::infinity();
        const auto points = profile.points();
        for (std::size_t i = 0; i + 1U < points.size(); ++i)
        {
            const double ax = points[i].x.si(), az = points[i].y.si();
            const double dx = points[i + 1U].x.si() - ax, dz = points[i + 1U].y.si() - az;
            const double t = std::clamp(((x - ax) * dx + (z - az) * dz) / (dx * dx + dz * dz), 0.0, 1.0);
            best = std::min(best, std::hypot(x - (ax + t * dx), z - (az + t * dz)));
        }
        return best;
    };
    double worst = 0.0;
    for (const bool forward : {true, false})
    {
        const PolylineProfile& from = forward ? a : b;
        const PolylineProfile& to = forward ? b : a;
        const auto points = from.points();
        for (std::size_t i = 0; i + 1U < points.size(); ++i)
        {
            const double ax = points[i].x.si(), az = points[i].y.si();
            const double dx = points[i + 1U].x.si() - ax, dz = points[i + 1U].y.si() - az;
            const auto steps = static_cast<int>(std::max(1.0, std::hypot(dx, dz) / step));
            for (int k = 0; k <= steps; ++k)
            {
                const double x = ax + dx * k / steps;
                if (x >= xMin && x <= xMax)
                {
                    worst = std::max(worst, distanceTo(to, x, az + dz * k / steps));
                }
            }
        }
    }
    return worst;
}

} // namespace sub0physics::terrain::testing
