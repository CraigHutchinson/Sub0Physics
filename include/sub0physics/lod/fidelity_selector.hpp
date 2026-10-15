#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "sub0physics/fault.hpp"
#include "sub0physics/lod/model_descriptor.hpp"

namespace sub0physics
{

/**
 * The validity envelope of one model, as a band of a scalar demand.
 *
 * Demand is whatever the consumer chooses to measure the physical need for detail (contact activity,
 * strain rate, predicted error); it is not camera distance. A model stays valid while demand is at most
 * `leaveAbove`, and may be entered only when demand is at most `enterBelow`. The gap between the two
 * is hysteresis.
 */
struct FidelityBand
{
    ModelDescriptor model; ///< The model this band describes.
    double enterBelow;     ///< Highest demand at which the model may be selected.
    double leaveAbove;     ///< Highest demand at which the model stays valid; at least `enterBelow`.
};

/** Why a decision was made. */
enum class SwitchReason : std::uint8_t
{
    Stay,        ///< The active model remains valid and no cheaper one is due.
    Refine,      ///< The active model is no longer valid; a costlier valid model was chosen.
    Coarsen,     ///< A cheaper model is valid and the dwell time has elapsed.
    NoValidModel ///< Nothing covers the demand; the active model is kept and the caller must report it.
};

/** A proposed model change. Apply it with FidelitySelector::commit after transferring state. */
struct FidelityDecision
{
    std::size_t from;    ///< Band index active when decided.
    std::size_t to;      ///< Band index to run next.
    SwitchReason reason; ///< Why.

    /** @return True when the decision changes the active model. */
    [[nodiscard]] constexpr bool isSwitch() const noexcept
    {
        return from != to;
    }
};

/**
 * Chooses among alternative models of one mechanism from a physical demand, with hysteresis.
 *
 * Bands are ordered from cheapest to costliest. Refinement is immediate whenever the active model leaves
 * its envelope; coarsening waits for a minimum dwell, so the selection cannot chatter. Deciding and
 * committing are separate: the caller transfers state between models in between and commits only if
 * the transfer succeeded, so a failed transfer always leaves the previous model running.
 *
 * @note Not thread-safe. Decisions depend only on demand and tick, so replay reproduces them.
 */
class FidelitySelector
{
public:
    /**
     * Validates bands and creates a selector.
     *
     * @param bands Models ordered by strictly increasing relative cost; each with
     *        `enterBelow <= leaveAbove`.
     * @param minimumDwellTicks Ticks a model must stay active before a coarsening switch.
     * @param initial Index of the first active band.
     * @return The selector, or an InvalidInput fault.
     */
    [[nodiscard]] static std::expected<FidelitySelector, Fault> create(std::vector<FidelityBand> bands,
                                                                       std::uint64_t minimumDwellTicks,
                                                                       std::size_t initial);

    /**
     * Decides which model should run for the coming tick.
     *
     * @param demand Current demand; NaN is treated as exceeding every band.
     * @param tick Current simulation tick.
     * @return The decision; nothing changes until it is committed.
     */
    [[nodiscard]] FidelityDecision decide(double demand, std::uint64_t tick) const noexcept;

    /**
     * Applies a decision after its state transfer succeeded.
     *
     * @param decision A decision from `decide` made since the last commit.
     * @param tick Tick at which the new model takes over.
     * @return False, and nothing changed, when the decision is stale (its `from` is not active).
     */
    [[nodiscard]] bool commit(const FidelityDecision& decision, std::uint64_t tick) noexcept;

    /** @return Index of the active band. */
    [[nodiscard]] std::size_t active() const noexcept
    {
        return active_;
    }

    /** @return The bands, cheapest first. */
    [[nodiscard]] std::span<const FidelityBand> bands() const noexcept
    {
        return bands_;
    }

    /** @return Number of committed switches. */
    [[nodiscard]] std::uint64_t switchCount() const noexcept
    {
        return switchCount_;
    }

private:
    FidelitySelector(std::vector<FidelityBand> bands, std::uint64_t minimumDwellTicks, std::size_t initial) noexcept;

    std::vector<FidelityBand> bands_;
    std::uint64_t minimumDwellTicks_;
    std::size_t active_;
    std::uint64_t activeSince_ = 0U;
    std::uint64_t switchCount_ = 0U;
};

} // namespace sub0physics
