#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "sub0physics/fault.hpp"
#include "sub0physics/terrain/grid_spec.hpp"

namespace sub0physics::terrain
{

/** Whether a tile load may finish now. */
enum class LoadState : std::uint8_t
{
    Ready,   ///< The tile was built into the slot it was given.
    Waiting, ///< Not yet available; the cache asks again at its next `pump`.
    Failed   ///< The load failed for good; whoever waits for the tile must be told.
};

/** A source's answer to one load attempt. */
struct LoadResult
{
    LoadState state;          ///< The decision.
    std::string_view failure; ///< Why, when Failed; static storage.
};

/**
 * Where tiles come from: a file, a pack, a procedural surface or a test fixture. The cache calls `load`
 * only inside `SourceTileCache::pump`, so *when* a tile becomes available is decided at a synchronisation
 * point, and a replay that reproduces the decisions reproduces the cache.
 *
 * @tparam Tile The immutable tile type; default constructible, because the cache owns a fixed pool of them.
 */
template<typename Tile>
class TileSource
{
public:
    virtual ~TileSource() = default;

    /**
     * Tries to build a tile.
     *
     * @param key The tile wanted.
     * @param slot Storage to build into, used only when the result is Ready. It holds whatever an earlier
     *             tile left there, so a source writes every field it later reads.
     * @return Ready (the slot holds the tile), Waiting (try again at the next pump) or Failed.
     */
    [[nodiscard]] virtual LoadResult load(const TileKey& key, Tile& slot) noexcept = 0;
};

/** Counters of a cache, as BikeLab's R0a `stats`. */
struct CacheStats
{
    std::uint64_t hits = 0;      ///< Accesses that found the tile resident or already loading.
    std::uint64_t misses = 0;    ///< Loads started.
    std::uint64_t loads = 0;     ///< Loads finished and made resident.
    std::uint64_t evictions = 0; ///< Tiles dropped from residency, least recently used first.

    friend constexpr bool operator==(const CacheStats&, const CacheStats&) = default;
};

/** Sizes of a cache; everything is allocated by `SourceTileCache::create`. */
struct CacheLimits
{
    std::uint32_t residentTiles = 32; ///< Tiles that may be resident, 1 to MaxTiles.
    std::uint32_t pendingLoads = 8;   ///< Loads that may be in flight at once, 1 to MaxTiles.
    /**
     * Storage beyond `residentTiles` for evicted tiles that readers still hold. A cache whose readers keep at
     * most `p` evicted tiles pinned never runs out with `extraSlots = p`. (One more slot, not counted here, is
     * always kept: a tile is built before the one it replaces is evicted.)
     */
    std::uint32_t extraSlots = 0;
    std::uint32_t tileBytes = 0; ///< Bytes one tile counts for in `bytes()`; the caller's own accounting.
};

/** What `SourceTileCache::access` found. */
enum class TileAccess : std::uint8_t
{
    Resident, ///< The tile is resident; `slot` is valid and the tile is now the most recently used.
    Loading,  ///< A load is in flight for it.
    Started,  ///< A new load was started.
    Refused   ///< The pending-load budget is full, so no load was started.
};

/** The answer of `SourceTileCache::access`. */
struct AccessResult
{
    TileAccess state;       ///< What happened.
    std::int32_t slot = -1; ///< The tile's storage, when Resident.
};

/** A load that finished inside `pump`, handed to the caller's callback at once. */
struct LoadEvent
{
    TileKey key;              ///< The tile.
    std::int32_t slot;        ///< Its storage, when it loaded; -1 when it failed.
    bool failed;              ///< The load failed.
    std::string_view failure; ///< Why, when failed; static storage.
};

/**
 * A bounded cache of immutable source tiles with deterministic least-recently-used eviction.
 *
 * - **De-duplicated loads.** Asking for a tile that is resident or already loading never starts a second
 *   load; both count as hits.
 * - **Bounded.** At most `residentTiles` tiles are resident and `pendingLoads` loads in flight. A request
 *   beyond the load budget is `Refused`: a visible outcome, never a silent truncation.
 * - **Pins.** A reader pins a slot while it uses the tile. An evicted tile's storage stays valid until its
 *   pins are released, and the slot is not reused before.
 * - **Deterministic.** Eviction takes the least recently *accessed* resident tile. Loads are decided oldest
 *   first, and only inside `pump`, so the sequence of cache states depends on the accesses and the source's
 *   answers, never on timing.
 * - **No allocation after `create`.** `access`, `pump`, `pin` and the rest allocate nothing. Whether a
 *   *load* allocates is up to the source: a tile that is a fixed-size array is built in place.
 *
 * Not thread-safe. The callback of `pump` must not call back into the cache.
 *
 * @tparam Tile The immutable tile type; default constructible.
 */
template<typename Tile>
class SourceTileCache
{
public:
    static constexpr std::uint32_t MaxTiles = 1U << 16; ///< Most resident tiles or pending loads.

    /**
     * A pin on one slot that releases itself: move-only. A lease must not outlive the cache or be held
     * across moving it.
     */
    class Lease
    {
    public:
        Lease() noexcept = default;

        /**
         * Pins a slot.
         *
         * @param cache The cache.
         * @param slot A slot returned by `access` or `LoadEvent`.
         */
        Lease(SourceTileCache& cache, std::int32_t slot) noexcept : cache_(&cache), slot_(slot)
        {
            cache.pin(slot);
        }

        Lease(Lease&& other) noexcept : cache_(std::exchange(other.cache_, nullptr)), slot_(std::exchange(other.slot_, -1))
        {
        }

        Lease& operator=(Lease&& other) noexcept
        {
            if (this != &other)
            {
                release();
                cache_ = std::exchange(other.cache_, nullptr);
                slot_ = std::exchange(other.slot_, -1);
            }
            return *this;
        }

        Lease(const Lease&) = delete;
        Lease& operator=(const Lease&) = delete;

        ~Lease()
        {
            release();
        }

        /** Unpins the slot, if any. */
        void release() noexcept
        {
            if (cache_ != nullptr)
            {
                cache_->unpin(slot_);
                cache_ = nullptr;
                slot_ = -1;
            }
        }

        /** @return True when a slot is held. */
        [[nodiscard]] explicit operator bool() const noexcept
        {
            return cache_ != nullptr;
        }

        /** @return The slot held, or -1. */
        [[nodiscard]] std::int32_t slot() const noexcept
        {
            return slot_;
        }

        /** @return The pinned tile; the lease must hold a slot. */
        [[nodiscard]] const Tile& tile() const noexcept
        {
            return cache_->tile(slot_);
        }

    private:
        SourceTileCache* cache_ = nullptr;
        std::int32_t slot_ = -1;
    };

    /**
     * Allocates a cache.
     *
     * @param limits Sizes; `residentTiles` and `pendingLoads` from 1 to MaxTiles, `extraSlots` up to MaxTiles.
     * @param source Where tiles come from; non-owning, must outlive the cache.
     * @return The cache, or an InvalidInput fault for a limit out of range or a null source.
     */
    [[nodiscard]] static std::expected<SourceTileCache, Fault> create(const CacheLimits& limits, TileSource<Tile>* source)
    {
        if (limits.residentTiles < 1U || limits.residentTiles > MaxTiles)
        {
            return invalidInput("a tile cache needs 1 to 65536 resident tiles");
        }
        if (limits.pendingLoads < 1U || limits.pendingLoads > MaxTiles)
        {
            return invalidInput("a tile cache needs 1 to 65536 pending loads");
        }
        if (limits.extraSlots > MaxTiles)
        {
            return invalidInput("a tile cache has at most 65536 extra slots");
        }
        if (source == nullptr)
        {
            return invalidInput("a tile cache needs a tile source");
        }
        return SourceTileCache(limits, source);
    }

    /** @return The counters. */
    [[nodiscard]] const CacheStats& stats() const noexcept
    {
        return stats_;
    }

    /** @return The limits. */
    [[nodiscard]] const CacheLimits& limits() const noexcept
    {
        return limits_;
    }

    /** @return Resident tiles. */
    [[nodiscard]] std::size_t residentCount() const noexcept
    {
        return resident_;
    }

    /** @return Loads in flight. */
    [[nodiscard]] std::size_t pendingCount() const noexcept
    {
        return loads_.size();
    }

    /** @return Resident tiles times `CacheLimits::tileBytes`. */
    [[nodiscard]] std::size_t bytes() const noexcept
    {
        return resident_ * limits_.tileBytes;
    }

    /**
     * Lists the resident tiles in key order (`TileKey` ordering), so the list does not depend on the order
     * of use.
     *
     * @param out Where to write; at most `out.size()` keys.
     * @return How many keys were written.
     */
    [[nodiscard]] std::size_t resident(std::span<TileKey> out) const noexcept
    {
        std::size_t n = 0;
        for (std::size_t slot = 0; slot < slotResident_.size() && n < out.size(); ++slot)
        {
            if (slotResident_[slot] != 0U)
            {
                out[n++] = keys_[slot];
            }
        }
        std::sort(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(n));
        return n;
    }

    /**
     * Reports whether a tile is resident, without touching it or the counters.
     *
     * @param key The tile.
     * @return The slot, or -1.
     */
    [[nodiscard]] std::int32_t residentSlot(const TileKey& key) const noexcept
    {
        for (std::size_t p = home(key);; p = (p + 1U) & mask_)
        {
            const std::int32_t slot = table_[p];
            if (slot < 0)
            {
                return -1;
            }
            if (keys_[static_cast<std::size_t>(slot)] == key)
            {
                return slot;
            }
        }
    }

    /**
     * Reports whether a load for a tile is in flight.
     *
     * @param key The tile.
     * @return True when one is.
     */
    [[nodiscard]] bool isLoading(const TileKey& key) const noexcept
    {
        return std::ranges::find(loads_, key) != loads_.end();
    }

    /**
     * Asks for a tile. A resident tile becomes the most recently used (a hit), a load in flight is joined
     * (a hit), otherwise a load starts if the budget allows (a miss) and is `Refused` if not.
     *
     * @param key The tile.
     * @return What happened and, when Resident, the slot. The slot is valid until the tile is evicted *and*
     *         unpinned; pin it to keep it.
     */
    [[nodiscard]] AccessResult access(const TileKey& key) noexcept
    {
        const std::int32_t slot = residentSlot(key);
        if (slot >= 0)
        {
            ++stats_.hits;
            touch(slot);
            return {TileAccess::Resident, slot};
        }
        if (isLoading(key))
        {
            ++stats_.hits;
            return {TileAccess::Loading, -1};
        }
        if (loads_.size() >= limits_.pendingLoads)
        {
            return {TileAccess::Refused, -1};
        }
        ++stats_.misses;
        loads_.push_back(key); // capacity reserved by the constructor
        return {TileAccess::Started, -1};
    }

    /**
     * Decides every load in flight once, oldest first. A load that finishes becomes resident (evicting the
     * least recently used tiles to make room) or is dropped, and `onEvent` is told at once, before the next
     * load is decided, so a reader can pin the tile before a later load evicts it.
     *
     * @param onEvent Called with a `LoadEvent` for each load that finished or failed; must not call back
     *                into the cache.
     * @return True when any load finished or failed.
     */
    template<typename F>
    bool pump(F&& onEvent) noexcept(noexcept(onEvent(std::declval<const LoadEvent&>())))
    {
        bool decided = false;
        for (std::size_t i = 0; i < loads_.size();)
        {
            const TileKey key = loads_[i];
            const std::int32_t slot = freeSlot();
            if (slot < 0)
            {
                // Unreachable when extraSlots covers every tile readers can pin; fails the load rather than hang.
                loads_.erase(loads_.begin() + static_cast<std::ptrdiff_t>(i));
                onEvent(LoadEvent{key, -1, true, "tile storage exhausted"});
                decided = true;
                continue;
            }
            const LoadResult result = source_->load(key, tiles_[static_cast<std::size_t>(slot)]);
            if (result.state == LoadState::Waiting)
            {
                ++i;
                continue;
            }
            loads_.erase(loads_.begin() + static_cast<std::ptrdiff_t>(i));
            decided = true;
            if (result.state == LoadState::Failed)
            {
                onEvent(LoadEvent{key, -1, true, result.failure});
                continue;
            }
            while (resident_ >= limits_.residentTiles)
            {
                evictOldest();
            }
            makeResident(key, slot);
            ++stats_.loads;
            onEvent(LoadEvent{key, slot, false, {}});
        }
        return decided;
    }

    /**
     * Pins a slot: its storage stays valid, even if the tile is evicted, until `unpin`.
     *
     * @param slot A slot from `access` or a `LoadEvent`.
     */
    void pin(std::int32_t slot) noexcept
    {
        ++pins_[static_cast<std::size_t>(slot)];
    }

    /**
     * Releases a pin.
     *
     * @param slot A pinned slot.
     */
    void unpin(std::int32_t slot) noexcept
    {
        --pins_[static_cast<std::size_t>(slot)];
    }

    /**
     * Reads a tile.
     *
     * @param slot A slot that is resident or pinned.
     * @return The tile.
     */
    [[nodiscard]] const Tile& tile(std::int32_t slot) const noexcept
    {
        return tiles_[static_cast<std::size_t>(slot)];
    }

private:
    SourceTileCache(const CacheLimits& limits, TileSource<Tile>* source)
        : limits_(limits), source_(source), tiles_(static_cast<std::size_t>(limits.residentTiles) + limits.extraSlots + 1U)
    {
        const std::size_t slots = tiles_.size();
        pins_.assign(slots, 0U);
        slotResident_.assign(slots, 0U);
        keys_.assign(slots, TileKey{0, 0, 0});
        prev_.assign(slots, -1);
        next_.assign(slots, -1);
        std::size_t capacity = 2;
        while (capacity < 2U * slots)
        {
            capacity *= 2U;
        }
        table_.assign(capacity, -1);
        mask_ = capacity - 1U;
        loads_.reserve(limits.pendingLoads);
    }

    // splitmix64 finaliser over the packed key.
    [[nodiscard]] std::size_t home(const TileKey& key) const noexcept
    {
        std::uint64_t x = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.i)) << 32U) |
                          static_cast<std::uint32_t>(key.j);
        x ^= static_cast<std::uint64_t>(key.level) * 0x9E3779B97F4A7C15ULL;
        x ^= x >> 30U;
        x *= 0xBF58476D1CE4E5B9ULL;
        x ^= x >> 27U;
        x *= 0x94D049BB133111EBULL;
        x ^= x >> 31U;
        return static_cast<std::size_t>(x) & mask_;
    }

    [[nodiscard]] std::int32_t freeSlot() const noexcept
    {
        for (std::size_t s = 0; s < tiles_.size(); ++s)
        {
            if (pins_[s] == 0U && slotResident_[s] == 0U)
            {
                return static_cast<std::int32_t>(s);
            }
        }
        return -1;
    }

    // Moves a resident slot to the newest end of the recency list.
    void touch(std::int32_t slot) noexcept
    {
        if (tail_ == slot)
        {
            return;
        }
        unlink(slot);
        link(slot);
    }

    void link(std::int32_t slot) noexcept
    {
        const auto s = static_cast<std::size_t>(slot);
        prev_[s] = tail_;
        next_[s] = -1;
        if (tail_ >= 0)
        {
            next_[static_cast<std::size_t>(tail_)] = slot;
        }
        else
        {
            head_ = slot;
        }
        tail_ = slot;
    }

    void unlink(std::int32_t slot) noexcept
    {
        const auto s = static_cast<std::size_t>(slot);
        if (prev_[s] >= 0)
        {
            next_[static_cast<std::size_t>(prev_[s])] = next_[s];
        }
        else
        {
            head_ = next_[s];
        }
        if (next_[s] >= 0)
        {
            prev_[static_cast<std::size_t>(next_[s])] = prev_[s];
        }
        else
        {
            tail_ = prev_[s];
        }
    }

    void makeResident(const TileKey& key, std::int32_t slot) noexcept
    {
        const auto s = static_cast<std::size_t>(slot);
        keys_[s] = key;
        slotResident_[s] = 1U;
        std::size_t p = home(key);
        while (table_[p] >= 0)
        {
            p = (p + 1U) & mask_;
        }
        table_[p] = slot;
        link(slot);
        ++resident_;
    }

    void evictOldest() noexcept
    {
        const std::int32_t slot = head_;
        unlink(slot);
        const auto s = static_cast<std::size_t>(slot);
        // Backward-shift deletion keeps probe sequences intact without tombstones.
        std::size_t i = home(keys_[s]);
        while (table_[i] != slot)
        {
            i = (i + 1U) & mask_;
        }
        table_[i] = -1;
        for (std::size_t j = (i + 1U) & mask_; table_[j] >= 0; j = (j + 1U) & mask_)
        {
            const std::size_t k = home(keys_[static_cast<std::size_t>(table_[j])]);
            const bool between = i <= j ? (i < k && k <= j) : (i < k || k <= j);
            if (!between)
            {
                table_[i] = table_[j];
                table_[j] = -1;
                i = j;
            }
        }
        slotResident_[s] = 0U; // the storage stays valid while pinned
        --resident_;
        ++stats_.evictions;
    }

    CacheLimits limits_;
    TileSource<Tile>* source_;
    CacheStats stats_;
    std::vector<Tile> tiles_;                // residentTiles + extraSlots
    std::vector<std::uint16_t> pins_;        // per slot
    std::vector<std::uint8_t> slotResident_; // per slot
    std::vector<TileKey> keys_;              // per slot; meaningful while resident
    std::vector<std::int32_t> prev_;         // recency list, towards the oldest
    std::vector<std::int32_t> next_;         // towards the newest
    std::vector<std::int32_t> table_;        // open-addressing map from key to slot; -1 empty
    std::size_t mask_ = 0;
    std::int32_t head_ = -1;                 // least recently used
    std::int32_t tail_ = -1;                 // most recently used
    std::size_t resident_ = 0;
    std::vector<TileKey> loads_;             // in flight, oldest first
};

} // namespace sub0physics::terrain
