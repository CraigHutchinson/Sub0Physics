// The continuous, gradient-consistent edge law (docs/WHEEL.md, "M3d status"): the smoothed-extremum kernel against
// central differences, the edge gap against the exact capsule gap (never below it, equal outside the band) and its
// gradient against central differences on faces, vertices beyond an end and vertices facing the interior, its
// continuity as an edge rolls past level and over a corner, and the row (ringEdgeGradientRow) as the gap's rate along
// any motion of a ring on a moving host.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <span>
#include <vector>

#include <doctest.h>

#include <sub0math/geometry/capsule_contact.hpp>

#include "sub0physics/geometry/convex_terrain.hpp"
#include "sub0physics/geometry/smoothed_extremum.hpp"
#include "sub0physics/multibody.hpp"
#include "test_ring_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using sub0math::f64::Length;

namespace
{

PolylineProfile profileOf(const std::vector<std::array<double, 2>>& points)
{
    std::vector<ProfilePoint> v;
    for (const auto& p : points)
    {
        v.push_back({Length::fromSi(p[0]), Length::fromSi(p[1])});
    }
    auto p = PolylineProfile::create(std::move(v));
    REQUIRE(p);
    return *p;
}

CasingEdge edgeOf(double ax, double ay, double bx, double by, double skin)
{
    return CasingEdge{{{Length::fromSi(ax), Length::fromSi(ay)}, {Length::fromSi(bx), Length::fromSi(by)}},
                      Length::fromSi(skin)};
}

double exactGap(const ConvexTerrain& terrain, std::size_t piece, double ax, double ay, double bx, double by,
                double skin)
{
    const sub0math::Capsule<double, 2> e{{{ax, ay}, {bx, by}}, skin};
    return sub0math::contact(e, terrain.pieces()[piece]).gap;
}

double uniform(std::mt19937_64& rng)
{
    return ringmodels::uniformSigned(rng);
}

sub0math::Vec3d toWorld(const Model& model, double i, double j)
{
    sub0math::Vec3d v{};
    v[model.planar().axisI] = i;
    v[model.planar().axisJ] = j;
    return v;
}

std::array<double, 2> nodeWorld(const Model& model, const Workspace& ws, std::size_t r, std::size_t k)
{
    const PlanarTransform& t = ws.planar.worldToBody[model.rings()[r].host];
    const double lx = ws.rings[r].positionX[k];
    const double ly = ws.rings[r].positionY[k];
    return {t.x + t.c * lx - t.s * ly, t.y + t.s * lx + t.c * ly};
}

} // namespace

TEST_SUITE("continuous edge law")
{
    TEST_CASE("smoothed extremum: the exact extremum outside the band and at the tie, derivatives exact")
    {
        std::mt19937_64 rng(11);
        const double band = 2e-5;
        double worstWeight = 0.0;
        double worstCurvature = 0.0;
        for (int i = 0; i < 2000; ++i)
        {
            const double a = 1e-4 * uniform(rng);
            const double b = a + 3.0 * band * uniform(rng);
            const SmoothedExtremum m = softMinimum(a, b, band);
            CHECK(m.value <= std::min(a, b) + 1e-18); // rounding of (a + b) / 2 - |a - b| / 2 outside the band
            CHECK(m.value >= std::min(a, b) - 0.5 * band);
            CHECK(m.weight0 + m.weight1 == doctest::Approx(1.0).epsilon(1e-15));
            CHECK(m.weight0 >= 0.0);
            CHECK(m.weight1 >= 0.0);
            if (std::abs(a - b) >= 2.0 * band)
            {
                CHECK(m.value == doctest::Approx(std::min(a, b)).epsilon(1e-15).scale(1e-4));
            }
            const double h = 1e-10;
            const SmoothedExtremum plus = softMinimum(a + h, b, band);
            const SmoothedExtremum minus = softMinimum(a - h, b, band);
            worstWeight = std::max(worstWeight, std::abs((plus.value - minus.value) / (2 * h) - m.weight0));
            if (std::abs(std::abs(a - b) - 2.0 * band) > 4 * h)
            {
                worstCurvature = std::max(worstCurvature,
                                          std::abs((plus.weight0 - minus.weight0) / (2 * h) - m.curvature) * band);
            }
            // The tread's maximum: exact at the tie and outside the band, never beyond the maximum; its weights leave
            // [0, 1] within the band (by at most 1/6), why the rows use the soft minimum.
            const SmoothedExtremum mx = smoothedMaximum(a, b, band);
            CHECK(mx.value <= std::max(a, b) + 1e-18);
            CHECK(mx.weight0 >= -1.0 / 6.0 - 1e-12);
            const SmoothedExtremum mxPlus = smoothedMaximum(a + h, b, band);
            const SmoothedExtremum mxMinus = smoothedMaximum(a - h, b, band);
            worstWeight = std::max(worstWeight, std::abs((mxPlus.value - mxMinus.value) / (2 * h) - mx.weight0));
        }
        CHECK(worstWeight < 1e-5);
        CHECK(worstCurvature < 1e-5);
        const SmoothedExtremum tie = softMinimum(3e-4, 3e-4, band);
        CHECK(tie.value == doctest::Approx(3e-4 - 0.5 * band).epsilon(1e-14));
        CHECK(tie.weight0 == 0.5);
        CHECK(smoothedMaximum(3e-4, 3e-4, band).value == 3e-4);
        // The weights move from one end to the other linearly across the band, not in one jump.
        const SmoothedExtremum edge = softMinimum(0.0, 2.0 * band, band);
        CHECK(edge.weight0 == 1.0);
        CHECK(softMinimum(0.0, band, band).weight0 == doctest::Approx(0.75).epsilon(1e-14));
    }

    TEST_CASE("edge gap: never below the exact capsule gap, equal outside the band; gradient exact by central differences")
    {
        const ConvexTerrain terrain(profileOf({{-3.0, 0.0}, {0.6, 0.0}, {0.6, 0.04}, {0.9, 0.04}, {0.9, 0.0}, {1.2, 0.03},
                                               {1.5, 0.0}, {5.0, 0.0}}));
        const double skin = 0.008;
        const double band = 2e-5;
        std::mt19937_64 rng(7);
        int interior = 0;
        int smoothed = 0;
        int inBand = 0;
        double worstGradient = 0.0;
        double worstBelow = 0.0;
        double worstAbove = 0.0;
        double worstOutside = 0.0;
        const std::array<std::array<double, 2>, 5> anchors{{{0.3, 0.0}, {0.6, 0.04}, {0.9, 0.04}, {1.2, 0.03}, {0.75, 0.04}}};
        for (int i = 0; i < 4000; ++i)
        {
            const auto& c = anchors[static_cast<std::size_t>(i) % anchors.size()];
            const double angle = (i % 3 == 0 ? 0.002 : 0.8) * uniform(rng);
            const double length = 0.032;
            const double mx = c[0] + 0.02 * uniform(rng);
            const double my = c[1] + skin + 0.006 * uniform(rng);
            const double ax = mx - 0.5 * length * std::cos(angle);
            const double ay = my - 0.5 * length * std::sin(angle);
            const double bx = mx + 0.5 * length * std::cos(angle);
            const double by = my + 0.5 * length * std::sin(angle);
            for (std::size_t piece = 0; piece < terrain.pieces().size(); ++piece)
            {
                const double exact = exactGap(terrain, piece, ax, ay, bx, by, skin);
                if (exact > 0.01)
                {
                    continue;
                }
                const ContinuousEdgeGap g = terrain.continuousEdgeGap(piece, edgeOf(ax, ay, bx, by, skin), Length::fromSi(band));
                if (g.axisCrosses)
                {
                    continue;
                }
                worstAbove = std::max(worstAbove, g.gap - exact);
                worstBelow = std::max(worstBelow, exact - g.gap);
                CHECK(g.weight0 >= 0.0);
                CHECK(g.weight1 >= 0.0);
                interior += g.interior ? 1 : 0;
                smoothed += g.interior ? 0 : 1;
                const bool far = g.weight0 == 1.0 || g.weight1 == 1.0;
                inBand += !g.interior && !far ? 1 : 0;
                if (g.interior || far)
                {
                    worstOutside = std::max(worstOutside, std::abs(g.gap - exact));
                }
                // Central differences in each end's coordinates against weight * normal.
                const double h = 1e-9;
                const auto at = [&](double dax, double day, double dbx, double dby) {
                    return terrain
                        .continuousEdgeGap(piece, edgeOf(ax + dax, ay + day, bx + dbx, by + dby, skin), Length::fromSi(band))
                        .gap;
                };
                const std::array<double, 4> numeric{(at(h, 0, 0, 0) - at(-h, 0, 0, 0)) / (2 * h),
                                                    (at(0, h, 0, 0) - at(0, -h, 0, 0)) / (2 * h),
                                                    (at(0, 0, h, 0) - at(0, 0, -h, 0)) / (2 * h),
                                                    (at(0, 0, 0, h) - at(0, 0, 0, -h)) / (2 * h)};
                const std::array<double, 4> analytic{g.weight0 * g.normal0.x, g.weight0 * g.normal0.y,
                                                     g.weight1 * g.normal1.x, g.weight1 * g.normal1.y};
                for (std::size_t k = 0; k < 4; ++k)
                {
                    worstGradient = std::max(worstGradient, std::abs(numeric[k] - analytic[k]));
                }
            }
        }
        INFO("interior " << interior << ", smoothed " << smoothed << ", in band " << inBand);
        CHECK(interior > 100);
        CHECK(smoothed > 1000);
        CHECK(inBand > 100);
        CHECK(worstAbove <= 1e-15);
        CHECK(worstBelow <= 0.5 * band + 1e-15);
        CHECK(worstOutside <= 1e-15);
        CHECK(worstGradient < 1e-6);
    }

    TEST_CASE("edge gap: continuous with a continuous gradient as an edge rolls past level and over a corner")
    {
        const ConvexTerrain terrain(profileOf({{-3.0, 0.0}, {0.6, 0.0}, {0.6, 0.04}, {0.9, 0.04}, {0.9, 0.0}, {5.0, 0.0}}));
        const double skin = 0.008;
        const double band = 2e-5;
        const double length = 0.032;
        // Rolling past level on a face: the edge turns about its middle, 0.1 mm into the face.
        double previousGap = 0.0;
        double previousWeight = 0.0;
        double worstGapJump = 0.0;
        double worstWeightJump = 0.0;
        const int steps = 20001;
        const double turn = 4e-3;
        for (int k = 0; k < steps; ++k)
        {
            const double angle = turn * (2.0 * k / (steps - 1) - 1.0);
            const double mx = 0.3;
            const double my = skin - 1e-4;
            const auto g = terrain.continuousEdgeGap(0, edgeOf(mx - 0.5 * length * std::cos(angle), my - 0.5 * length * std::sin(angle),
                                                               mx + 0.5 * length * std::cos(angle), my + 0.5 * length * std::sin(angle), skin),
                                                     Length::fromSi(band));
            if (k > 0)
            {
                worstGapJump = std::max(worstGapJump, std::abs(g.gap - previousGap));
                worstWeightJump = std::max(worstWeightJump, std::abs(g.weight0 - previousWeight));
            }
            previousGap = g.gap;
            previousWeight = g.weight0;
        }
        const double dAngle = 2.0 * turn / (steps - 1);
        // Lipschitz bounds: the gap moves at most half an edge per radian, the weights at most |dg/d(g0 - g1)| / band.
        CHECK(worstGapJump <= 0.5 * length * dAngle * 1.000001);
        CHECK(worstWeightJump <= length * dAngle / band);
        CHECK(worstWeightJump > 0.0);

        // Translating over the rock's corner at (0.6, 0.04): a face end, then the vertex inside the edge, then the top.
        int regimes = 0;
        bool wasInterior = false;
        worstGapJump = 0.0;
        double worstNormalJump = 0.0;
        sub0math::Vec2d previousNormal{};
        const double y = 0.04 + skin - 2e-4;
        for (int k = 0; k < steps; ++k)
        {
            const double x = 0.6 - 0.06 + 0.12 * k / (steps - 1);
            // Tilted so the edge meets the corner with its interior.
            const double angle = 0.1;
            const double ax = x - 0.5 * length * std::cos(angle);
            const double ay = y - 0.5 * length * std::sin(angle) + 0.01;
            const double bx = x + 0.5 * length * std::cos(angle);
            const double by = y + 0.5 * length * std::sin(angle) + 0.01;
            // The rock's piece is piece 1 (the segment on its top).
            const auto g = terrain.continuousEdgeGap(1, edgeOf(ax, ay, bx, by, skin), Length::fromSi(band));
            const sub0math::Vec2d n{g.weight0 * g.normal0.x + g.weight1 * g.normal1.x,
                                    g.weight0 * g.normal0.y + g.weight1 * g.normal1.y};
            if (k > 0)
            {
                worstGapJump = std::max(worstGapJump, std::abs(g.gap - previousGap));
                worstNormalJump = std::max(worstNormalJump, std::hypot(n.x - previousNormal.x, n.y - previousNormal.y));
                regimes += g.interior != wasInterior ? 1 : 0;
            }
            wasInterior = g.interior;
            previousGap = g.gap;
            previousNormal = n;
        }
        const double dx = 0.12 / (steps - 1);
        CHECK(regimes >= 2);
        CHECK(worstGapJump <= dx * 1.000001);
        // The gradient turns with the vertex direction: at most the edge's travel over the skin, per step.
        CHECK(worstNormalJump <= 4.0 * dx / skin);
    }

    TEST_CASE("edge gradient row is the gap's rate along any motion of a ring on a moving host")
    {
        const Model model = ringmodels::ringRig({.nodes = 16}, false, false);
        Workspace ws(model);
        const std::size_t n = model.dofCount();
        const double skin = 0.008;
        const double band = 2e-5;
        double worst = 0.0;
        int checked = 0;
        for (unsigned seed = 1; seed <= 24; ++seed)
        {
            const State s = ringmodels::randomState(model, seed);
            REQUIRE(prepareLinearlyImplicitStep(model, s, ws).has_value());
            const std::size_t edge = seed % 16;
            const std::size_t next = (edge + 1) % 16;
            const auto a = nodeWorld(model, ws, 0, edge);
            const auto b = nodeWorld(model, ws, 0, next);
            // A face just under the lower end (both laws in play), or a rock corner under the edge's middle.
            const double lowest = std::min(a[1], b[1]) - skin + 3e-4;
            const double mid = 0.5 * (a[0] + b[0]);
            const ConvexTerrain terrain =
                seed % 2 == 0 ? ConvexTerrain(profileOf({{-10.0, lowest}, {10.0, lowest + 1e-5 * static_cast<double>(seed % 5)}}))
                              : ConvexTerrain(profileOf({{-10.0, lowest - 0.05}, {mid - 0.004, lowest - 0.05},
                                                         {mid - 0.004, 0.5 * (a[1] + b[1]) - skin + 1e-3},
                                                         {mid + 0.2, 0.5 * (a[1] + b[1]) - skin + 1e-3 - 0.1},
                                                         {10.0, lowest - 0.2}}));
            for (std::size_t piece = 0; piece < terrain.pieces().size(); ++piece)
            {
                const auto gapAt = [&](const State& t) {
                    REQUIRE(prepareLinearlyImplicitStep(model, t, ws).has_value());
                    const auto pa = nodeWorld(model, ws, 0, edge);
                    const auto pb = nodeWorld(model, ws, 0, next);
                    return terrain.continuousEdgeGap(piece, edgeOf(pa[0], pa[1], pb[0], pb[1], skin), Length::fromSi(band));
                };
                const ContinuousEdgeGap g = gapAt(s);
                if (g.axisCrosses || g.gap > 0.02)
                {
                    continue;
                }
                RingEdgeGradient gradient{.ring = 0,
                                          .edge = edge,
                                          .weight0 = g.weight0,
                                          .weight1 = g.weight1,
                                          .point0 = toWorld(model, g.point0.x, g.point0.y),
                                          .point1 = toWorld(model, g.point1.x, g.point1.y),
                                          .direction0 = toWorld(model, g.normal0.x, g.normal0.y),
                                          .direction1 = toWorld(model, g.normal1.x, g.normal1.y)};
                std::vector<double> row(n);
                REQUIRE(ringEdgeGradientRow(model, ws, gradient, row).has_value());
                double jqd = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    jqd += row[i] * s.qd()[i];
                }
                const double eps = 1e-7;
                State plus = s;
                State minus = s;
                for (std::size_t i = 0; i < n; ++i)
                {
                    plus.q()[i] += eps * s.qd()[i];
                    minus.q()[i] -= eps * s.qd()[i];
                }
                const double rate = (gapAt(plus).gap - gapAt(minus).gap) / (2 * eps);
                worst = std::max(worst, std::abs(jqd - rate) / std::max(1.0, std::abs(rate)));
                ++checked;
            }
        }
        CHECK(checked >= 12);
        CHECK(worst < 1e-6);
    }
}
