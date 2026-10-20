#include <doctest.h>
#include <sub0math/units.hpp>

#include "sub0physics/verification/energy_ledger.hpp"

using namespace sub0physics;
using sub0math::joules;

TEST_CASE("Energy account balance counts dissipation and external work")
{
    constexpr EnergyAccount account{.kinetic = joules(3.0), .potential = joules(4.0), .dissipated = joules(2.0), .externalWork = joules(5.0)};
    static_assert(account.mechanical() == joules(7.0));
    static_assert(account.balance() == joules(4.0));
}

TEST_CASE("Energy ledger tracks the peak balance residual")
{
    EnergyLedger ledger(EnergyAccount{.potential = joules(10.0)}, joules(10.0));
    CHECK(ledger.initialBalance() == joules(10.0));
    CHECK(ledger.residual(EnergyAccount{.potential = joules(10.0)}).si() == 0.0);
    ledger.record({.kinetic = joules(4.0), .potential = joules(5.0), .dissipated = joules(1.0)});
    CHECK(ledger.maxResidual().si() == 0.0);
    ledger.record({.kinetic = joules(12.0), .dissipated = joules(0.5), .externalWork = joules(3.0)});
    CHECK(ledger.maxResidual().si() == doctest::Approx(0.5));
    ledger.record({.potential = joules(10.0), .dissipated = joules(0.1)});
    CHECK(ledger.maxResidual().si() == doctest::Approx(0.5));
    ledger.includeScale(joules(-50.0));
    CHECK(ledger.scale() == joules(50.0));
    CHECK(ledger.relativeResidual() == doctest::Approx(0.01));
}

TEST_CASE("Relative residual never divides by a vanishing scale")
{
    EnergyLedger ledger(EnergyAccount{}, joules(0.0));
    ledger.record({.kinetic = joules(1e-9)});
    CHECK(ledger.relativeResidual() == doctest::Approx(1e-9 / EnergyLedger::MinimumScale.si()));
}
