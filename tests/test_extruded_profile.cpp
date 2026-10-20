// Tier B: ExtrudedProfileSource. The slice is the profile (bit for bit), pages are samples of it, vertical faces
// are features, and the baked representation answers clearance questions within the bake's own error.
#include <algorithm>
#include <array>
#include <cmath>
#include <doctest.h>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "sub0physics/terrain.hpp"
#include "sub0physics/terrain/testing.hpp"

using namespace sub0physics;
using namespace sub0physics::terrain;
using namespace sub0physics::terrain::testing;
using sub0math::f64::Length;

namespace
{

ProfilePoint pt(double x, double z)
{
    return {Length::fromSi(x), Length::fromSi(z)};
}

PolylineProfile makeProfile(std::vector<ProfilePoint> points)
{
    auto p = PolylineProfile::create(std::move(points));
    REQUIRE(p);
    return std::move(*p);
}

PolylineProfile flat()
{
    return makeProfile({pt(-1.0, 0.0), pt(1.0, 0.0)});
}

PolylineProfile block()
{
    return makeProfile({pt(-1.0, 0.0), pt(-0.1, 0.0), pt(-0.1, 0.1), pt(0.1, 0.1), pt(0.1, 0.0), pt(1.0, 0.0)});
}

PolylineProfile rounded(int segments)
{
    std::vector<ProfilePoint> points{pt(-1.0, 0.0), pt(-0.1, 0.0)};
    for (int k = 1; k < segments; ++k)
    {
        const double a = std::numbers::pi * k / segments;
        points.push_back(pt(-0.1 * std::cos(a), 0.1 * std::sin(a)));
    }
    points.push_back(pt(0.1, 0.0));
    points.push_back(pt(1.0, 0.0));
    return makeProfile(std::move(points));
}

const GridSpec& bakeGrid()
{
    // 2.56 m pages at 1 cm, the profile (x in [-1, 1]) inside one page.
    static const GridSpec grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.01),
                                                    .anchor = {Length::fromSi(-1.28), Length::fromSi(-1.28)},
                                                    .pages = {1, 1}});
    return grid;
}

ExtrudedProfileSource makeSource(PolylineProfile profile)
{
    auto s = ExtrudedProfileSource::create(std::move(profile), Length::fromSi(-0.3), Length::fromSi(0.3), bakeGrid());
    REQUIRE(s);
    return std::move(*s);
}

BakedPage bake(ExtrudedProfileSource& source)
{
    BakedPage page;
    REQUIRE(source.load({0, 0, 0}, page).state == LoadState::Ready);
    return page;
}

struct Query
{
    double cx, cz, radius;
};

std::vector<Query> queries()
{
    std::vector<Query> out;
    for (const double radius : {0.05, 0.12, 0.2})
    {
        for (double cx = -0.7; cx <= 0.7; cx += 0.0371)
        {
            for (const double lift : {0.01, 0.06, 0.19})
            {
                out.push_back({cx, radius + lift, radius});
            }
        }
    }
    return out;
}

CircleQuery circle(const Query& q)
{
    return {{{Length::fromSi(q.cx), Length::fromSi(q.cz)}, Length::fromSi(q.radius)}, Length::fromSi(1e-10)};
}

} // namespace

TEST_CASE("Extruded profile: the slice is exactly measureCircleClearance, at every y inside the extrusion")
{
    for (const auto& profile : {flat(), block(), rounded(32)})
    {
        auto source = makeSource(profile);
        for (const Query& q : queries())
        {
            const CircleQuery query = circle(q);
            std::array<Witness2, DefaultWitnessCapacity> a{}, b{};
            const auto direct = measureCircleClearance(profile, query, a);
            for (const double y : {-0.3, -0.05, 0.0, 0.3})
            {
                const auto slice = source.clearanceSlice(y, query, b);
                REQUIRE(slice.status == direct.status);
                if (direct.status != ClearanceStatus::Available)
                {
                    continue;
                }
                CHECK(slice.distance.si() == direct.distance.si()); // bit for bit
                CHECK(slice.gap.si() == direct.gap.si());
                REQUIRE(slice.witnessCount == direct.witnessCount);
                for (std::size_t k = 0; k < slice.witnessCount; ++k)
                {
                    CHECK(b[k].index == a[k].index);
                    CHECK(b[k].kind == a[k].kind);
                    CHECK(b[k].normal.x == a[k].normal.x);
                }
            }
        }
        std::array<Witness2, 1> w{};
        CHECK(source.clearanceSlice(0.31, circle({0.0, 0.3, 0.1}), w).status == ClearanceStatus::OutsideDomain);
    }
}

TEST_CASE("Extruded profile: vertical faces are zero-thickness box features with outward normals")
{
    auto source = makeSource(block());
    REQUIRE(source.features().size() == 2U);
    const BoxFeature& rise = source.features()[0];
    const BoxFeature& fall = source.features()[1];
    CHECK(rise.bounds.lo.x.si() == -0.1);
    CHECK(rise.bounds.hi.x.si() == -0.1);
    CHECK(rise.bounds.lo.z.si() == 0.0);
    CHECK(rise.bounds.hi.z.si() == 0.1);
    CHECK(rise.bounds.lo.y.si() == -0.3);
    CHECK(rise.bounds.hi.y.si() == 0.3);
    CHECK(rise.faceNormal.x == -1.0); // air is towards -x of a rising face
    CHECK(fall.faceNormal.x == 1.0);
    CHECK(rise.id != fall.id);
    CHECK(makeSource(flat()).features().empty());
    CHECK(makeSource(rounded(32)).features().empty());
    std::array<std::size_t, 4> hits{};
    const auto found =
        featuresMeeting(source.features(), {{Length::fromSi(-0.2), Length::fromSi(0.0)}, {Length::fromSi(-0.05), Length::fromSi(0.1)}}, hits);
    CHECK(found.count == 1U);
    CHECK(hits[0] == 0U);
}

TEST_CASE("Extruded profile: pages are samples; nothing is extrapolated")
{
    auto source = makeSource(block());
    const BakedPage baked = bake(source);
    const HeightPage& page = *baked.page;
    for (const std::int32_t j : {0, 50, 128 + 30, 128 + 31, 256})
    {
        for (std::int32_t i = 0; i < HeightPage::Nodes; ++i)
        {
            const auto p = bakeGrid().samplePosition({0, 0, 0}, i, j);
            const auto exact = source.height(p.x.si(), p.y.si());
            REQUIRE(page.known(i, j) == exact.has_value());
            if (exact)
            {
                const double z = page.frame().origin.z.si() + static_cast<double>(page.nodeHeight(i, j).si());
                CHECK(std::abs(z - *exact) <= 6e-5); // 0.1 mm quantisation: half a step plus float
            }
        }
    }
    CHECK_FALSE(page.known(0, 128)); // x = -1.28, outside the profile
    CHECK_FALSE(page.known(128, 0)); // y = -1.28, outside the extrusion
    CHECK(page.known(128, 128));
}

TEST_CASE("Extruded profile: a riser is a ramp in the page and a face in the features")
{
    auto source = makeSource(block());
    const BakedPage baked = bake(source);
    const HeightPage& page = *baked.page;
    const double xf = -0.1;
    const auto i = static_cast<std::int32_t>(std::floor((xf - bakeGrid().anchor().x.si()) / 0.01 + 1e-9));
    const double z0 = page.frame().origin.z.si();
    const double left = z0 + static_cast<double>(page.nodeHeight(i - 1, 128).si());
    const double right = z0 + static_cast<double>(page.nodeHeight(i + 1, 128).si());
    CHECK(std::abs(left) < 1e-4);        // the node left of the face is on the ground
    CHECK(std::abs(right - 0.1) < 1e-4); // the node right of it is on the block
    // Just left of the face the page is already above the true ground (0): a ramp, which the feature replaces.
    const auto mid = page.heightAt(sub0math::Length::fromSi(static_cast<float>(xf - bakeGrid().anchor().x.si() + 0.005)),
                                   sub0math::Length::fromSi(1.28F));
    REQUIRE(mid.known);
    CHECK(static_cast<double>(mid.height.si()) + z0 > 0.0005);
}

TEST_CASE("Extruded profile: the baked page plus its faces answers clearance within the bake's own error")
{
    struct Case
    {
        const char* name;
        PolylineProfile profile;
        double maxDeviation; // loose upper bound the fixture must stay under
    };
    std::vector<Case> cases;
    cases.push_back({"flat", flat(), 1e-6});
    cases.push_back({"block", block(), 1e-6});
    cases.push_back({"rounded", rounded(32), 3e-3}); // measured 2.5 mm: 1 cm nodes against a 0.1 m radius
    for (auto& c : cases)
    {
        CAPTURE(c.name);
        auto source = makeSource(c.profile);
        const BakedPage baked = bake(source);
        auto rebuilt = profileFromPageRow(bakeGrid(), *baked.page, 128, source.features());
        if (!rebuilt) { MESSAGE(std::string(rebuilt.error().detail)); }
        REQUIRE(rebuilt);
        const double deviation = sampledHausdorff(c.profile, *rebuilt, -0.97, 0.97); // sampled Hausdorff distance of the bake
        CHECK(deviation <= c.maxDeviation);
        MESSAGE(std::string(c.name), ": Hausdorff distance ", deviation, " m");
        // The distance from a point to a set moves by at most the Hausdorff distance between the sets.
        double worst = 0.0;
        std::size_t compared = 0;
        for (const Query& q : queries())
        {
            std::array<Witness2, DefaultWitnessCapacity> wa{}, wb{};
            const auto exact = measureCircleClearance(c.profile, circle(q), wa);
            const auto bakedClearance = measureCircleClearance(*rebuilt, circle(q), wb);
            if (exact.status != ClearanceStatus::Available || bakedClearance.status != ClearanceStatus::Available)
            {
                continue;
            }
            ++compared;
            const double difference = std::abs(exact.distance.si() - bakedClearance.distance.si());
            worst = std::max(worst, difference);
            CHECK(difference <= deviation + 1e-9);
        }
        CHECK(compared > 100U);
        MESSAGE(std::string(c.name), ": worst clearance difference ", worst, " m over ", compared, " queries");
    }
}

TEST_CASE("Extruded profile: served through the cache, pages arrive on pump")
{
    auto source = makeSource(block());
    auto cache = SourceTileCache<BakedPage>::create({.residentTiles = 1, .pendingLoads = 2, .extraSlots = 1}, &source);
    REQUIRE(cache);
    CHECK(cache->access({0, 0, 0}).state == TileAccess::Started);
    std::int32_t slot = -1;
    CHECK(cache->pump([&](const LoadEvent& e) { slot = e.slot; }));
    REQUIRE(slot >= 0);
    CHECK(cache->tile(slot).page.has_value());
    CHECK_FALSE(ExtrudedProfileSource::create(flat(), Length::fromSi(1.0), Length::fromSi(0.0), bakeGrid()));
}
