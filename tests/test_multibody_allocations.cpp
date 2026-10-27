// Counts every heap allocation in this executable, to pin "a multibody step never allocates".
#include <array>
#include <atomic>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <doctest.h>
#include <ostream> // doctest prints strings (MSVC needs it explicitly)
#include <new>
#include <vector>

#include <sub0math/units.hpp>

#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"
#include "test_ring_models.hpp"

namespace
{

std::atomic<std::size_t> allocations{0};

} // namespace

void* operator new(std::size_t size)
{
    ++allocations;
    if (void* p = std::malloc(size == 0U ? 1U : size))
    {
        return p;
    }
    throw std::bad_alloc();
}

void operator delete(void* p) noexcept
{
    std::free(p);
}

void operator delete(void* p, std::size_t) noexcept
{
    std::free(p);
}

using namespace sub0physics::multibody;

TEST_CASE("Multibody dynamics and stepping (semi-implicit Euler, RK4) never allocate")
{
    const Model model = testmodels::planarTree();
    Workspace ws(model);
    State s = State::atRest(model);
    s.q()[2] = 0.4;
    std::vector<double> tau(model.dofCount(), 0.1);
    std::vector<double> out(model.dofCount());
    std::vector<double> jacobian(3 * model.dofCount());
    const SiteId tip = *model.findSite("tip");

    Rk4Stepper rk4(model, ws);
    std::vector<double> otherTau(model.dofCount(), -0.1);

    const std::size_t start = allocations.load();
    for (int k = 0; k < 2000; ++k)
    {
        REQUIRE(stepSemiImplicitEuler(model, {tau}, s.tick, sub0math::seconds(1e-3), s, ws));
    }
    for (int k = 0; k < 200; ++k)
    {
        REQUIRE(rk4.step(model, {k % 50 < 25 ? tau : otherTau}, s.tick, sub0math::seconds(1e-3), s, ws));
    }
    REQUIRE(inverseDynamics(model, s, ws.qdd, out, ws));
    REQUIRE(computeSiteJacobian(model, s, tip, ws, jacobian));
    const SiteKinematics site = computeSiteKinematics(model, s, tip, ws).value();
    const sub0physics::EnergyAccount energy = computeEnergy(model, s, ws).value();
    CHECK(allocations.load() - start == 0U);
    CHECK(std::isfinite(site.pose.translation.x.si() + energy.mechanical().si()));
}

TEST_CASE("Stepping with force elements (implicit Euler, semi-implicit Euler, RK4) never allocates")
{
    // Dense path: a tendon couples two branches.
    const Model model = testmodels::planarTreeWithElements();
    Workspace ws(model);
    State s = State::atRest(model);
    s.q()[2] = 0.4;
    s.q()[8] = -0.3;
    std::vector<double> tau(model.dofCount(), 0.1);
    Rk4Stepper rk4(model, ws);
    // Tree path: a tendon to a world anchor along one chain.
    const Model chain = [] {
        sub0physics::multibody::ModelBuilder b = testmodels::planarTreeBuilder();
        const SiteId anchor = b.addSite({.name = "anchor", .body = WorldBody, .frame = testmodels::offset(0.3, -1.5, 0.0)});
        b.addElement(TendonElementDescription<sub0physics::LinearSpringDamper<sub0math::f64::Length>>{
            .name = "cable",
            .first = SiteId{0},
            .second = anchor,
            .law = *sub0physics::LinearSpringDamper<sub0math::f64::Length>::create(
                {.stiffness = sub0math::f64::Stiffness::fromSi(5.0e4), .damping = sub0math::f64::Damping::fromSi(30.0)}),
            .reference = testmodels::metres(0.2)});
        return *b.build();
    }();
    REQUIRE(chain.elementsFollowTree());
    Workspace chainWs(chain);
    State c = State::atRest(chain);
    c.q()[8] = 0.2;

    const std::size_t start = allocations.load();
    for (int k = 0; k < 1000; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {tau}, s.tick, sub0math::seconds(1e-3), s, ws));
        REQUIRE(stepLinearlyImplicitEuler(chain, {}, c.tick, sub0math::seconds(1e-3), c, chainWs));
    }
    for (int k = 0; k < 200; ++k)
    {
        REQUIRE(stepSemiImplicitEuler(model, {tau}, s.tick, sub0math::seconds(1e-4), s, ws));
        REQUIRE(rk4.step(model, {tau}, s.tick, sub0math::seconds(1e-4), s, ws));
    }
    REQUIRE(evaluateElements(model, s, ws));
    const sub0physics::EnergyAccount energy = computeEnergy(model, s, ws).value();
    CHECK(allocations.load() - start == 0U);
    CHECK(std::isfinite(energy.balance().si()));
}

TEST_CASE("Constrained stepping (closures, limits, split steps), projection and reactions never allocate")
{
    // A slider-crank with a hard stop on the slider: closures, an LCP row and event location in one run.
    sub0physics::multibody::ModelBuilder builder = testmodels::sliderCrankBuilder();
    const State start = testmodels::sliderCrankState({}, 0.7, 6.0);
    builder.addConstraint(CoordinateLimitDescription{.name = "stop", .body = BodyId{2}, .upper = start.q()[2] + 0.002,
                                                     .restitution = 0.3});
    const Model model = *builder.build();
    const Model closed = *testmodels::sliderCrankBuilder().build();
    Workspace ws(model);
    Workspace closedWs(closed);
    Rk4Stepper rk4(closed, closedWs);
    State s = start;
    State c = start;
    std::vector<double> tau(model.dofCount(), 0.0);
    tau[0] = 0.3;
    std::size_t splits = 0;

    const std::size_t before = allocations.load();
    for (int k = 0; k < 500; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {tau}, s.tick, sub0math::seconds(1e-3), s, ws));
        splits += ws.constraintReport.substeps > 1U ? 1U : 0U;
        REQUIRE(stepSemiImplicitEuler(model, {tau}, s.tick, sub0math::seconds(1e-3), s, ws));
        splits += ws.constraintReport.substeps > 1U ? 1U : 0U;
        REQUIRE(rk4.step(closed, {tau}, c.tick, sub0math::seconds(1e-3), c, closedWs));
    }
    REQUIRE(projectConstraints(model, s, ws));
    REQUIRE(computeReactions(model, s, tau, ws));
    CHECK(allocations.load() - before == 0U);
    CHECK(splits > 0U);
}

TEST_CASE("M1d: site forces, forcing Jacobians, the midpoint step, the prepared step and reactions never allocate")
{
    for (const Formulation formulation : {Formulation::Automatic, Formulation::General})
    {
        ModelBuilder builder = testmodels::planarTreeWithElementsBuilder();
        builder.setFormulation(formulation);
        const Model model = *builder.build();
        const std::size_t n = model.dofCount();
        Workspace ws(model);
        Workspace aux(model);
        State s = State::atRest(model);
        s.q()[2] = 0.3;
        const SiteId tip = *model.findSite("tip");
        const std::array<SiteForce, 1> forces{SiteForce{tip, {0.0, 30.0, 0.0}, {}}};
        std::vector<double> stiffness(n * n, 0.0);
        std::vector<double> damping(n * n, 0.0);
        stiffness[0] = -100.0;
        damping[0] = -1.0;
        const Forcing forcing{{}, forces, stiffness, damping};
        const std::array<SiteId, 1> sites{tip};
        std::array<SiteKinematics, 1> kinematics{};
        std::array<SiteMotion, 1> motion{};
        sub0physics::EnergyAccount account{};
        const auto h = sub0math::f64::Time::fromSi(1e-3);

        const std::size_t start = allocations.load();
        for (int k = 0; k < 200; ++k)
        {
            REQUIRE(stepLinearlyImplicitEuler(model, forcing, s.tick, h, s, ws, ImplicitTerms::DampingOnly));
            REQUIRE(stepLinearlyImplicitMidpoint(model, forcing, s.tick, h, s, ws));
            REQUIRE(stepSemiImplicitEuler(model, forcing, s.tick, h, s, ws));
            REQUIRE(prepareLinearlyImplicitStep(model, s, ws));
            REQUIRE(computePreparedSiteKinematics(model, s, sites, kinematics, ws, &account));
            REQUIRE(computePreparedSiteMotion(model, s, sites, motion, ws, &account));
            REQUIRE(completeLinearlyImplicitStep(model, forcing, s.tick, h, s, ws, ImplicitScheme::Midpoint));
            REQUIRE(computeSiteKinematics(model, s, sites, kinematics, ws, &account));
            REQUIRE(computeReactionsWithForcing(model, s, forcing, aux));
            REQUIRE(computeAppliedForce(model, s, forcing, aux));
        }
        CHECK(allocations.load() - start == 0U);
        CHECK(std::isfinite(account.kinetic.si() + motion[0].position.x.si()));
    }
}

TEST_CASE("Ring substructures: steps, dynamics, elements, energy and reactions never allocate")
{
    const Model model = ringmodels::ringRig({.nodes = 64}, false);
    Workspace ws(model);
    Workspace aux(model);
    State s = ringmodels::randomState(model, 1);
    std::vector<double> tau(model.dofCount(), 0.1);
    std::vector<double> out(model.dofCount());
    Rk4Stepper rk4(model, ws);
    const auto h = sub0math::f64::Time::fromSi(1e-5);

    const std::size_t start = allocations.load();
    for (int k = 0; k < 50; ++k)
    {
        const Forcing forcing{tau};
        REQUIRE(stepLinearlyImplicitEuler(model, forcing, s.tick, h, s, ws));
        REQUIRE(stepLinearlyImplicitEuler(model, forcing, s.tick, h, s, ws, ImplicitTerms::DampingOnly));
        REQUIRE(stepLinearlyImplicitMidpoint(model, forcing, s.tick, h, s, ws));
        REQUIRE(stepSemiImplicitEuler(model, forcing, s.tick, h, s, ws));
        REQUIRE(rk4.step(model, forcing, s.tick, h, s, ws));
        REQUIRE(prepareLinearlyImplicitStep(model, s, ws));
        REQUIRE(completeLinearlyImplicitStep(model, forcing, s.tick, h, s, ws));
        REQUIRE(forwardDynamics(model, s, tau, aux));
        REQUIRE(solveMassMatrix(model, aux, out));
        REQUIRE(inverseDynamics(model, s, tau, out, aux));
        REQUIRE(evaluateElements(model, s, aux));
        REQUIRE(computeEnergy(model, s, aux));
        REQUIRE(computeReactions(model, s, tau, aux));
    }
    CHECK(allocations.load() - start == 0U);
    CHECK(std::isfinite(s.values[10]));
}
