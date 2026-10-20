// Counts every heap allocation in this executable, to pin "a multibody step never allocates".
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
