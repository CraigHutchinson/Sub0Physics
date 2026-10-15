#include <doctest.h>

#include "sub0physics/verification/energy_ledger.hpp"

using namespace sub0physics;

TEST_CASE("Energy ledger tracks the peak balance residual")
{
    EnergyLedger ledger(10.0, 10.0);
    CHECK(ledger.residual(10.0, 0.0, 0.0) == 0.0);
    ledger.record(9.0, 0.0, 1.0);
    CHECK(ledger.maxResidual() == 0.0);
    ledger.record(12.0, 3.0, 0.5);
    CHECK(ledger.maxResidual() == doctest::Approx(0.5));
    ledger.record(10.0, 0.0, 0.1);
    CHECK(ledger.maxResidual() == doctest::Approx(0.5));
    ledger.includeScale(-50.0);
    CHECK(ledger.scale() == 50.0);
    CHECK(ledger.relativeResidual() == doctest::Approx(0.01));
}

TEST_CASE("Relative residual never divides by a vanishing scale")
{
    EnergyLedger ledger(0.0, 0.0);
    ledger.record(1e-9, 0.0, 0.0);
    CHECK(ledger.relativeResidual() == doctest::Approx(1e-9 / EnergyLedger::MinimumScale));
}
