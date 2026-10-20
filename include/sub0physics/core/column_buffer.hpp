#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <type_traits>

namespace sub0physics
{

/**
 * Fixed-capacity array storage for one column of a structure-of-arrays result: allocated once by `reserve`,
 * then resized and written in place without allocating.
 *
 * Unlike `std::vector`, `resize` leaves new elements uninitialised (the caller is about to write all of
 * them), and nothing ever grows past the reserved capacity: a request that does not fit is refused with
 * `false`, never truncated, so a caller cannot deliver half a result without noticing. Indexing is unchecked.
 *
 * @tparam T Trivially copyable element type.
 */
template<typename T>
class ColumnBuffer
{
    static_assert(std::is_trivially_copyable_v<T>, "ColumnBuffer copies and leaves its elements uninitialised");

public:
    ColumnBuffer() = default;
    ColumnBuffer(ColumnBuffer&&) noexcept = default;
    ColumnBuffer& operator=(ColumnBuffer&&) noexcept = default;
    ~ColumnBuffer() = default;

    /**
     * Copies the elements in use into storage of the same capacity.
     *
     * @param other The buffer to copy.
     */
    ColumnBuffer(const ColumnBuffer& other)
        : data_(other.capacity_ != 0U ? std::make_unique_for_overwrite<T[]>(other.capacity_) : nullptr), size_(other.size_),
          capacity_(other.capacity_)
    {
        std::copy_n(other.data_.get(), size_, data_.get());
    }

    /**
     * Copies the elements in use.
     *
     * @param other The buffer to copy.
     * @return This buffer.
     */
    ColumnBuffer& operator=(const ColumnBuffer& other)
    {
        if (this != &other)
        {
            ColumnBuffer copy(other);
            *this = std::move(copy);
        }
        return *this;
    }

    /**
     * Allocates storage and empties the buffer. The only call that allocates.
     *
     * @param capacity The most elements ever held.
     */
    void reserve(std::size_t capacity)
    {
        data_ = capacity != 0U ? std::make_unique_for_overwrite<T[]>(capacity) : nullptr;
        capacity_ = capacity;
        size_ = 0;
    }

    /**
     * Sets the number of elements in use; new elements are uninitialised.
     *
     * @param size The count.
     * @return True, or false (and nothing changes) when `size` exceeds the capacity.
     */
    [[nodiscard]] bool resize(std::size_t size) noexcept
    {
        if (size > capacity_)
        {
            return false;
        }
        size_ = size;
        return true;
    }

    /** Empties the buffer; the capacity stays. */
    void clear() noexcept
    {
        size_ = 0;
    }

    /**
     * Appends an element.
     *
     * @param value The element.
     * @return True, or false (and nothing changes) when the buffer is full.
     */
    [[nodiscard]] bool pushBack(const T& value) noexcept
    {
        if (size_ == capacity_)
        {
            return false;
        }
        data_[size_++] = value;
        return true;
    }

    /** @return The elements in use. */
    [[nodiscard]] std::size_t size() const noexcept
    {
        return size_;
    }

    /** @return The most elements the buffer can hold. */
    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return capacity_;
    }

    /** @return True when no element is in use. */
    [[nodiscard]] bool empty() const noexcept
    {
        return size_ == 0U;
    }

    /**
     * Reads or writes an element.
     *
     * @param i An index below `size()`; not checked.
     * @return The element.
     */
    [[nodiscard]] T& operator[](std::size_t i) noexcept
    {
        return data_[i];
    }

    /**
     * Reads an element.
     *
     * @param i An index below `size()`; not checked.
     * @return The element.
     */
    [[nodiscard]] const T& operator[](std::size_t i) const noexcept
    {
        return data_[i];
    }

    /** @return The first element's address (null before `reserve`). */
    [[nodiscard]] T* data() noexcept
    {
        return data_.get();
    }

    /** @return The elements in use, read-only. */
    [[nodiscard]] std::span<const T> view() const noexcept
    {
        return {data_.get(), size_};
    }

    /** @return The elements in use, writable. */
    [[nodiscard]] std::span<T> mutableView() noexcept
    {
        return {data_.get(), size_};
    }

private:
    std::unique_ptr<T[]> data_;
    std::size_t size_ = 0;
    std::size_t capacity_ = 0;
};

} // namespace sub0physics
