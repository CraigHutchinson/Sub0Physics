#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace sub0physics
{

/**
 * Fixed-capacity recording of C channels, stored column by column.
 *
 * All storage is allocated by the constructor, so recording inside an integration loop never
 * allocates. Each column is contiguous, ready for plotting, comparison or export.
 *
 * @tparam C Number of channels.
 */
template<std::size_t C>
class ChannelTable
{
public:
    /** One recorded sample: a value per channel. */
    using Row = std::array<double, C>;

    /**
     * Allocates storage for `capacity` rows.
     *
     * @param capacity Maximum number of rows.
     */
    explicit ChannelTable(std::size_t capacity)
        : capacity_(capacity)
        , values_(C * capacity)
    {
    }

    /**
     * Appends a row.
     *
     * @param row Values in channel order.
     * @return False, and nothing recorded, when the table is full.
     */
    [[nodiscard]] bool append(const Row& row) noexcept
    {
        if (size_ == capacity_)
        {
            return false;
        }
        for (std::size_t c = 0; c < C; ++c)
        {
            values_[c * capacity_ + size_] = row[c];
        }
        ++size_;
        return true;
    }

    /**
     * Views one channel.
     *
     * @param channel Channel index, less than C.
     * @return The recorded values; valid until the table is destroyed or appended to past capacity.
     */
    [[nodiscard]] std::span<const double> column(std::size_t channel) const noexcept
    {
        return std::span<const double>(values_).subspan(channel * capacity_, size_);
    }

    /**
     * Reads one recorded row.
     *
     * @param index Row index, less than size().
     * @return The values of that row in channel order.
     */
    [[nodiscard]] Row row(std::size_t index) const noexcept
    {
        Row out{};
        for (std::size_t c = 0; c < C; ++c)
        {
            out[c] = values_[c * capacity_ + index];
        }
        return out;
    }

    /** @return Number of rows recorded. */
    [[nodiscard]] std::size_t size() const noexcept
    {
        return size_;
    }

    /** @return Maximum number of rows. */
    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return capacity_;
    }

    /** @return True when no rows are recorded. */
    [[nodiscard]] bool empty() const noexcept
    {
        return size_ == 0U;
    }

private:
    std::size_t capacity_;
    std::size_t size_ = 0U;
    std::vector<double> values_;
};

} // namespace sub0physics
