#pragma once

#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/**
 * The energy of a run at one instant, by channel (docs/MULTIBODY.md, "M1b contracts").
 *
 * `kinetic` and `potential` describe the state: `potential` holds gravitational energy plus the
 * elements' stored energy. `dissipated` and `externalWork` are accumulated over the run: the energy
 * removed by dissipation and the work done on the system by external forcing. For an exact run the
 * `balance()` stays at its initial value.
 */
struct EnergyAccount
{
    sub0math::f64::Energy kinetic{};      ///< Kinetic energy now.
    sub0math::f64::Energy potential{};    ///< Gravitational plus elastic (element) potential energy now.
    sub0math::f64::Energy dissipated{};   ///< Energy dissipated so far; never decreases.
    sub0math::f64::Energy externalWork{}; ///< Work done on the system by external forcing so far.

    /** @return kinetic + potential. */
    [[nodiscard]] constexpr sub0math::f64::Energy mechanical() const noexcept
    {
        return kinetic + potential;
    }

    /** @return kinetic + potential + dissipated - externalWork: constant for an exact run. */
    [[nodiscard]] constexpr sub0math::f64::Energy balance() const noexcept
    {
        return kinetic + potential + dissipated - externalWork;
    }

    friend constexpr bool operator==(const EnergyAccount&, const EnergyAccount&) = default;
};

} // namespace sub0physics
