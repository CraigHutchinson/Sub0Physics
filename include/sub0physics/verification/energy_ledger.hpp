#pragma once

#include <algorithm>
#include <cmath>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/verification/energy_account.hpp"

namespace sub0physics
{

/**
 * Tracks the energy balance of a run: `balance(t) = balance(0)`, where an account's balance is
 * kinetic + potential + dissipated - external work (EnergyAccount).
 *
 * The residual of that balance measures integration error. Dividing its peak by the largest energy
 * magnitude seen in the run gives a scale-free figure that can be compared with a tolerance.
 */
class EnergyLedger
{
public:
    /** Smallest scale used for the relative residual; keeps an idle run from dividing by zero. */
    static constexpr sub0math::f64::Energy MinimumScale = sub0math::f64::Energy::fromSi(1e-6);

    /**
     * Starts a ledger.
     *
     * @param initial The account at the start of the run.
     * @param initialScale Energy magnitude known at the start.
     */
    EnergyLedger(const EnergyAccount& initial, sub0math::f64::Energy initialScale) noexcept
        : initialBalance_(initial.balance().si())
        , scale_(std::abs(initialScale.si()))
    {
    }

    /**
     * Computes the balance residual of one account.
     *
     * @param account The account now.
     * @return Its balance minus the initial balance.
     */
    [[nodiscard]] sub0math::f64::Energy residual(const EnergyAccount& account) const noexcept
    {
        return sub0math::f64::Energy::fromSi(account.balance().si() - initialBalance_);
    }

    /**
     * Records one account's residual into the running peak.
     *
     * @param account The account now.
     */
    void record(const EnergyAccount& account) noexcept
    {
        maxResidual_ = std::max(maxResidual_, std::abs(residual(account).si()));
    }

    /**
     * Widens the run's energy scale.
     *
     * @param magnitude An energy magnitude seen in the run.
     */
    void includeScale(sub0math::f64::Energy magnitude) noexcept
    {
        scale_ = std::max(scale_, std::abs(magnitude.si()));
    }

    /** @return The balance of the initial account. */
    [[nodiscard]] sub0math::f64::Energy initialBalance() const noexcept
    {
        return sub0math::f64::Energy::fromSi(initialBalance_);
    }

    /** @return Largest absolute residual recorded. */
    [[nodiscard]] sub0math::f64::Energy maxResidual() const noexcept
    {
        return sub0math::f64::Energy::fromSi(maxResidual_);
    }

    /** @return Largest energy magnitude seen. */
    [[nodiscard]] sub0math::f64::Energy scale() const noexcept
    {
        return sub0math::f64::Energy::fromSi(scale_);
    }

    /** @return Peak residual divided by the energy scale (at least MinimumScale). */
    [[nodiscard]] double relativeResidual() const noexcept
    {
        return maxResidual_ / std::max(scale_, MinimumScale.si());
    }

private:
    double initialBalance_;
    double scale_;
    double maxResidual_ = 0.0;
};

} // namespace sub0physics
