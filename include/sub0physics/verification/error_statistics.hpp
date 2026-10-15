#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sub0physics
{

/**
 * Peak and RMS of the differences between a candidate channel and its reference, against a limit.
 *
 * A non-finite difference poisons the statistics on purpose: `isFinite()` then reports false and
 * `passes()` can never succeed.
 */
class ErrorStatistics
{
public:
    /**
     * Creates empty statistics.
     *
     * @param limit Largest acceptable absolute difference, in the channel's unit.
     */
    explicit ErrorStatistics(double limit) noexcept
        : limit_(limit)
    {
    }

    /**
     * Adds one paired sample.
     *
     * @param difference Candidate minus reference.
     */
    void add(double difference) noexcept
    {
        const double magnitude = std::abs(difference);
        maxAbsolute_ = std::isnan(magnitude) ? magnitude : std::max(maxAbsolute_, magnitude);
        sumOfSquares_ += difference * difference;
        ++count_;
    }

    /** @return Largest absolute difference so far. */
    [[nodiscard]] double maxAbsolute() const noexcept
    {
        return maxAbsolute_;
    }

    /** @return Root mean square of the differences; zero when empty. */
    [[nodiscard]] double rms() const noexcept
    {
        return count_ == 0U ? 0.0 : std::sqrt(sumOfSquares_ / static_cast<double>(count_));
    }

    /** @return Number of samples added. */
    [[nodiscard]] std::size_t count() const noexcept
    {
        return count_;
    }

    /** @return The acceptance limit. */
    [[nodiscard]] double limit() const noexcept
    {
        return limit_;
    }

    /** @return True when every difference and their sum of squares were finite. */
    [[nodiscard]] bool isFinite() const noexcept
    {
        return std::isfinite(maxAbsolute_) && std::isfinite(sumOfSquares_);
    }

    /** @return True when at least one sample was added, all were finite and none exceeded the limit. */
    [[nodiscard]] bool passes() const noexcept
    {
        return count_ != 0U && isFinite() && maxAbsolute_ <= limit_;
    }

private:
    double limit_;
    double maxAbsolute_ = 0.0;
    double sumOfSquares_ = 0.0;
    std::size_t count_ = 0U;
};

} // namespace sub0physics
