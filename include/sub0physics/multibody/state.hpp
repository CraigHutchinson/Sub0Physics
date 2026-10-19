#pragma once

#include <cstdint>
#include <vector>

#include "sub0physics/multibody/model.hpp"

namespace sub0physics::multibody
{

/**
 * The evolving state of a multibody model, as a value: copy it to snapshot, assign it to restore.
 *
 * Coordinates are plain SI doubles in the model's coordinate order: radians for revolute joints, metres
 * for prismatic ones, and their rates. Stepping updates the vectors in place and never resizes them.
 */
struct State
{
    std::vector<double> q;  ///< Generalized coordinates.
    std::vector<double> qd; ///< Generalized velocities.
    std::uint64_t tick = 0; ///< Index of the fixed step the state is at.

    /**
     * Makes the rest state of a model: every coordinate and rate zero, at tick zero.
     *
     * @param model The model.
     * @return A state sized for the model.
     */
    [[nodiscard]] static State atRest(const Model& model)
    {
        return State{std::vector<double>(model.dofCount(), 0.0), std::vector<double>(model.dofCount(), 0.0), 0};
    }

    friend bool operator==(const State&, const State&) = default;
};

} // namespace sub0physics::multibody
