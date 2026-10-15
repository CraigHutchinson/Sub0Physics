#include <cmath>
#include <utility>

#include "sub0physics/lod/fidelity_selector.hpp"

namespace sub0physics
{

std::expected<FidelitySelector, Fault> FidelitySelector::create(std::vector<FidelityBand> bands,
                                                                std::uint64_t minimumDwellTicks, std::size_t initial)
{
    if (bands.empty() || initial >= bands.size())
    {
        return std::unexpected(Fault{FaultKind::InvalidInput, "selector needs bands and a valid initial band"});
    }
    for (std::size_t i = 0; i < bands.size(); ++i)
    {
        const FidelityBand& band = bands[i];
        if (std::isnan(band.enterBelow) || std::isnan(band.leaveAbove) || band.enterBelow > band.leaveAbove)
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "band needs enterBelow <= leaveAbove"});
        }
        if (!(band.model.relativeCost > 0.0) || (i > 0U && !(band.model.relativeCost > bands[i - 1U].model.relativeCost)))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "bands must be ordered by increasing positive cost"});
        }
    }
    return FidelitySelector(std::move(bands), minimumDwellTicks, initial);
}

FidelitySelector::FidelitySelector(std::vector<FidelityBand> bands, std::uint64_t minimumDwellTicks,
                                   std::size_t initial) noexcept
    : bands_(std::move(bands))
    , minimumDwellTicks_(minimumDwellTicks)
    , active_(initial)
{
}

FidelityDecision FidelitySelector::decide(double demand, std::uint64_t tick) const noexcept
{
    // Comparisons with NaN are false, so NaN demand is valid nowhere.
    if (!(demand <= bands_[active_].leaveAbove))
    {
        for (std::size_t i = active_ + 1U; i < bands_.size(); ++i)
        {
            if (demand <= bands_[i].leaveAbove)
            {
                return {active_, i, SwitchReason::Refine};
            }
        }
        return {active_, active_, SwitchReason::NoValidModel};
    }
    if (tick - activeSince_ >= minimumDwellTicks_)
    {
        for (std::size_t i = 0; i < active_; ++i)
        {
            if (demand <= bands_[i].enterBelow)
            {
                return {active_, i, SwitchReason::Coarsen};
            }
        }
    }
    return {active_, active_, SwitchReason::Stay};
}

bool FidelitySelector::commit(const FidelityDecision& decision, std::uint64_t tick) noexcept
{
    if (decision.from != active_ || decision.to >= bands_.size())
    {
        return false;
    }
    if (decision.isSwitch())
    {
        active_ = decision.to;
        activeSince_ = tick;
        ++switchCount_;
    }
    return true;
}

} // namespace sub0physics
