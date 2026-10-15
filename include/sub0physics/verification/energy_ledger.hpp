#pragma once

#include <algorithm>
#include <cmath>

namespace sub0physics
{

/**
 * Tracks the energy balance of a run: E(t) + dissipated - external work = E(0).
 *
 * The residual of that balance measures integration error. Dividing its peak by the largest energy
 * magnitude seen in the run gives a scale-free figure that can be compared with a tolerance.
 */
class EnergyLedger
{
public:
    /** Smallest scale used for the relative residual, J; keeps an idle run from dividing by zero. */
    static constexpr double MinimumScale = 1e-6;

    /**
     * Starts a ledger.
     *
     * @param initialEnergy Mechanical energy at the start, J.
     * @param initialScale Energy magnitude known at the start, J.
     */
    EnergyLedger(double initialEnergy, double initialScale) noexcept
        : initialEnergy_(initialEnergy)
        , scale_(std::abs(initialScale))
    {
    }

    /**
     * Computes the balance residual for one state.
     *
     * @param energy Mechanical energy now, J.
     * @param work External work done on the system so far, J.
     * @param dissipated Energy dissipated so far, J.
     * @return energy + dissipated - work - initial energy, J.
     */
    [[nodiscard]] double residual(double energy, double work, double dissipated) const noexcept
    {
        return energy + dissipated - work - initialEnergy_;
    }

    /**
     * Records one state's residual into the running peak.
     *
     * @param energy Mechanical energy now, J.
     * @param work External work done on the system so far, J.
     * @param dissipated Energy dissipated so far, J.
     */
    void record(double energy, double work, double dissipated) noexcept
    {
        maxResidual_ = std::max(maxResidual_, std::abs(residual(energy, work, dissipated)));
    }

    /**
     * Widens the run's energy scale.
     *
     * @param magnitude An energy magnitude seen in the run, J.
     */
    void includeScale(double magnitude) noexcept
    {
        scale_ = std::max(scale_, std::abs(magnitude));
    }

    /** @return Mechanical energy at the start, J. */
    [[nodiscard]] double initialEnergy() const noexcept
    {
        return initialEnergy_;
    }

    /** @return Largest absolute residual recorded, J. */
    [[nodiscard]] double maxResidual() const noexcept
    {
        return maxResidual_;
    }

    /** @return Largest energy magnitude seen, J. */
    [[nodiscard]] double scale() const noexcept
    {
        return scale_;
    }

    /** @return Peak residual divided by the energy scale (at least MinimumScale). */
    [[nodiscard]] double relativeResidual() const noexcept
    {
        return maxResidual_ / std::max(scale_, MinimumScale);
    }

private:
    double initialEnergy_;
    double scale_;
    double maxResidual_ = 0.0;
};

} // namespace sub0physics
