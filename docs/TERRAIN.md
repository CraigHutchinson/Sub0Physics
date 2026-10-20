# Terrain foundations

What `sub0physics/terrain.hpp` and its neighbours provide today, and the rules they follow. The wider design
(representations, extraction, state, the request service) is the consumer's, and lives with the product that
drives it; this page describes only what is built here.

## Grids: `GridSpec`

A `GridSpec` says how pages tile the plane. It is built only by `GridSpec::create`, which validates a
`GridLayout`:

| Field | Meaning | Default |
|---|---|---|
| `finestSpacing` | Cell side at level 0 [m] | required |
| `levels` | Number of levels; each doubles the spacing, so a level-k page is `2^k` times the side | 1 |
| `cellsPerSide` | Cells along each side of every page (1 to 65536) | 256 |
| `registration` | `NodeIsPoint` or `CellCentre` (below) | `NodeIsPoint` |
| `anchor` | Terrain-frame position of page (0, 0)'s minimum corner [m] | (0, 0) |
| `pages` | Pages along x and y at level 0 | 65536 x 65536 |

```cpp
// An 8 m x 3.2 m source of 2 cm cells in 40-cell pages, sampled at cell centres: the shape of the R0a fixture.
const auto grid = GridSpec::create({.finestSpacing = metres(0.02),
                                    .cellsPerSide = 40,
                                    .registration = Registration::CellCentre,
                                    .anchor = {metres(0.0), metres(-1.6)},
                                    .pages = {10, 4}});
```

`GridSpec::create(spacing, levels)` is the short form: 256 cells, node registration and the default extent,
centred on the terrain origin.

**The grid is a finite rectangle.** Pages are indexed from the anchor towards +x and +y, so a `TileKey`
inside the grid has non-negative `i` and `j`. A point outside the rectangle has no page. (An authored
asset is finite, and the anchor is its placement; there is no unbounded world to index.) At level k the
rectangle has `ceil(pages / 2^k)` pages per axis, so a coarse page may extend past the authored edge, and a
level-k page is exactly the parent of the `2^k x 2^k` level-0 pages under it.

### Registration

The cell structure is the same for both; registration says where a source's samples sit.

- `NodeIsPoint` (GDAL "pixel is point"): a sample is the value at a grid node. A page has `cells + 1`
  samples per side, and the last row and column coincide with the neighbour's first, so pages share their
  edge nodes. `HeightPage` is this kind.
- `CellCentre` ("pixel is area"): a sample is the value over a cell, placed at its centre. A page has
  `cells` samples per side and shares nothing with its neighbours.

`samplesPerSide()` and `samplePosition(key, i, j)` give the count and the terrain-frame position, so a
source is sampled where its registration says.

### Ownership and coverage

Both operations are `sub0math::RegularGrid` on the cell grid of the level (`GridSpec::cellGrid(level)`), so
there is one set of boundary rules and no floor or clamp written here.

| Operation | Functions | Rule |
|---|---|---|
| Ownership | `pageOf`, `cellOf` | Each point belongs to exactly one cell, and so one page. Half-open: a point on a shared edge belongs to the cell on the positive side; the grid's far edge belongs to the last cell. A non-finite point, a level outside the grid or a point outside the rectangle is not known (a status aggregate, never a conversion to an integer). |
| Coverage | `cellsTouching`, `pagesTouching` | The cells (pages) a *closed* box touches, clipped to the grid. A box that reaches an edge exactly includes the cells on both sides. This is the rule of the R0a `cells(box)` function and of the separating-axis testers. |

Ownership answers "which cell holds this state"; coverage answers "which cells might this region need".
Using one for the other loses a cell at every boundary. A point exactly on a page boundary can land on
either side if the boundary's own computation rounded (the corner `anchor + side * k` is not always bit-equal
to `k * side` seen from the point); nothing relies on that, and a consumer that needs the page of a corner
should ask `pageCorner(key)` to be offset inward.

### Pages

`HeightPage` is a 256-cell, node-registered page with a fixed layout, so its `heightAt` has constant strides
and no page size to load. `buildHeightPage` refuses a grid with another cell count or cell-centre
registration, and a key outside the grid.

**Why not runtime-sized.** A runtime cell count costs the lookup: with the count read at run time instead of
as a constant, the hand-written kernel (`page size x1024` in `tests/bench_physics.cpp`) was 1.04x to 1.28x
slower over three runs (about 1.1x typical), on a path that is already measured against a 1.0x reference. A
source of a different shape (such as the 40-cell R0a fixture, or a 64-cell page) therefore gets its own page
type when it needs one, built from the same `GridSpec`; the grid describes any shape, the fast page serves
the default.

## Corridors: `CorridorBuilder`

Given a pose pair (position and heading, interpolated linearly) and a footprint (half length, half width,
and a lookbehind opposite the direction of travel), `CorridorBuilder::create` builds one oriented box per
sub-interval (1 to 64; 8 is typical), and `rasterise` marks the cells of any `RegularGrid` that the boxes
touch (separating-axis test of each candidate cell, by `AabbOverlapTester::overlapsCentred`). The output is a
bit mask in caller-owned words: no allocation, and a status for a corridor that misses the grid or a buffer
that is too small.

Each box is aligned with the heading at the middle of its sub-interval, and its half extents are the
footprint plus a bound for the rest of the motion: the heading's departure from the middle (using
`1 - cos e <= e^2/2` and `|sin e| <= |e|`, so it is rigorous for any rotation), half the sub-interval's
displacement, the lookbehind, and a nanometre of rounding slack. See the class comment for the bound.

**Containment is of the swept footprint, not of an axis-aligned box.** R0a's `sweptBox` is the axis-aligned
box of the corner hull, padded by `R |dyaw| / 2`; it contains the swept footprint, but it also contains the
empty corners of a diagonal or turning footprint, so no tight corridor can contain it. The property tests
therefore check, with oracles rebuilt from R0a's definitions: (1) every cell touched by the footprint at 1001
instants of random pose pairs (including the lookbehind copies, stationary pairs and any rotation) is in the
corridor's mask, and (2) R0a's swept box itself contains the footprint, while the corridor needs fewer cells
than the swept box. A mutation check (dropping the rotation term) fails (1).

## Identifiers: `Id<Tag, Rep>`

`sub0physics/core/id.hpp`. A class wrapping an unsigned integer, one type per tag. Construction from an
integer is explicit, and ids of different tags do not convert or compare. The maximum of `Rep` is the
invalid sentinel and a default-constructed id is invalid (so zero is a real id). Ordering is `<=>` on the
value, `std::hash` is provided, and `next()` returns the following id, saturating at the sentinel so a
counter never wraps onto a used value. `SourceRevision` is an `Id`.

## Witnesses: `Witness<N>`

`sub0physics/geometry/witness.hpp`. Every nearest-point or hit query returns a witness `{point, normal,
kind, index}`: `WitnessKind` is `Vertex`, `Edge` (the interior of a segment, `index` its first vertex) or
`Face` (3D). A query reports the feature itself and never averages the normals of several: at a square
corner the witness is the vertex with the true approach direction.

`measureCircleClearance(profile, query, std::span<Witness2> out)` writes the tied witnesses into the
caller's buffer, in profile order, each once. When more features tie than the buffer holds, the first ones
are kept, `overflowed` is set and `witnessCount` is the number written. A zero-length buffer asks for the
distance only. `DefaultWitnessCapacity` (8) is enough for corners and valleys. The query allocates nothing
(`Sub0Physics_AllocationTests`).

## Test surfaces: `terrain/testing/`

`sub0physics/terrain/testing.hpp` (not in the `sub0physics.hpp` umbrella) is the analytic-surface catalogue:
`Plane`, `Ramp`, `Step`, `RoundedBlock`, `Bridge`, `Berm`, `Rut`. Each has `height(x, y)` and an exact unit
`normal(x, y)`; pieces are half-open on their low side. `Step` and `Bridge` also give their vertical and
under faces (`Step::riserNormal`, `Bridge::spansAt`, `Bridge::undersideNormal`). `sampleSurface` fills a
page's samples at its registration's positions. These are the tier-C fixtures: the frozen ports stay the
oracles (tier A), and adapters feed the same scenarios into the generic components (tier B).
