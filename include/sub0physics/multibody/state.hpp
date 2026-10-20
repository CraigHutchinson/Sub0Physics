#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <span>
#include <vector>

#include "sub0physics/multibody/model.hpp"

namespace sub0physics::multibody
{

/**
 * The evolving state of a multibody model, as a value: copy it to snapshot, assign it to restore.
 *
 * Coordinates and rates live in one contiguous buffer `values = [q; qd]`, so an integrator can treat
 * the whole state as one vector, and `q()` and `qd()` are its two halves. Entries are plain SI doubles
 * in the model's coordinate order: radians for revolute joints, metres for prismatic ones, and their
 * rates. Stepping updates the buffer in place and never resizes it.
 */
struct State
{
    std::vector<double> values; ///< `[q; qd]`: dofCount coordinates, then dofCount rates.
    std::uint64_t tick = 0;     ///< Index of the fixed step the state is at.

    /**
     * Makes the rest state of a model: every coordinate and rate zero, at tick zero.
     *
     * @param model The model.
     * @return A state sized for the model.
     */
    [[nodiscard]] static State atRest(const Model& model)
    {
        return State{std::vector<double>(2 * model.dofCount(), 0.0), 0};
    }

    /**
     * Makes a state from coordinates and rates of equal length (a convenience for setting up runs).
     *
     * @param q Generalized coordinates.
     * @param qd Generalized velocities; as many as `q`. A shorter list is padded with zeros, a longer one
     *        truncated, so the state is always `[q; qd]` of `2 q.size()` entries.
     * @param tick The tick it is at.
     * @return The state.
     */
    [[nodiscard]] static State fromCoordinates(std::span<const double> q, std::span<const double> qd,
                                               std::uint64_t tick = 0)
    {
        State s{std::vector<double>(2 * q.size(), 0.0), tick};
        std::copy(q.begin(), q.end(), s.values.begin());
        std::copy_n(qd.begin(), std::min(qd.size(), q.size()), s.values.begin() + static_cast<std::ptrdiff_t>(q.size()));
        return s;
    }

    /** @copydoc fromCoordinates */
    [[nodiscard]] static State fromCoordinates(std::initializer_list<double> q, std::initializer_list<double> qd,
                                               std::uint64_t tick = 0)
    {
        return fromCoordinates(std::span<const double>(q.begin(), q.size()),
                               std::span<const double>(qd.begin(), qd.size()), tick);
    }

    /** @return The number of coordinates: half the buffer. */
    [[nodiscard]] std::size_t dofCount() const noexcept
    {
        return values.size() / 2;
    }

    /** @return The generalized coordinates, the first half of `values`. */
    [[nodiscard]] std::span<double> q() noexcept
    {
        return {values.data(), dofCount()};
    }

    /** @return The generalized coordinates, the first half of `values`. */
    [[nodiscard]] std::span<const double> q() const noexcept
    {
        return {values.data(), dofCount()};
    }

    /** @return The generalized velocities, the second half of `values`. */
    [[nodiscard]] std::span<double> qd() noexcept
    {
        return {values.data() + dofCount(), dofCount()};
    }

    /** @return The generalized velocities, the second half of `values`. */
    [[nodiscard]] std::span<const double> qd() const noexcept
    {
        return {values.data() + dofCount(), dofCount()};
    }

    /**
     * Checks that the state has the size a model needs.
     *
     * @param model The model.
     * @return True when `values` holds one coordinate and one rate per degree of freedom.
     */
    [[nodiscard]] bool isSizedFor(const Model& model) const noexcept
    {
        return values.size() == 2 * model.dofCount();
    }

    friend bool operator==(const State&, const State&) = default;
};

} // namespace sub0physics::multibody
