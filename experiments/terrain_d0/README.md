# D0: terrain representation decision experiment

This program measures the candidate terrain representations against exact references before any
data format is fixed. It is measurement code, not library API; it has its own minimal vector type.
The plan it serves is BikeDynamics' terrain design (section 4.6). The raw numbers are in
[RESULTS.md](RESULTS.md).

```bash
cmake -S . -B build-d0 -G Ninja -DCMAKE_BUILD_TYPE=Release -DSUB0PHYSICS_BUILD_EXPERIMENTS=ON -DSUB0PHYSICS_BUILD_TESTING=OFF
cmake --build build-d0 --target terrain_d0
build-d0/experiments/terrain_d0/terrain_d0 RESULTS.md
```

## What is compared

| Label | Representation |
|---|---|
| B | Heightfield page: 16-bit heights (0.1 mm steps), fixed-diagonal triangulation, exact closest point on its triangles |
| B (all baked) | B with rocks and lips baked into the heights along vertical rays, as a heightmap-only engine would |
| B+E | B for the ground, plus rocks and lips as exact features (oriented boxes, convex meshes) |
| K | Sparse deformation delta on B: 64×64-node int16 chunks, allocated when first written |
| H | Sparse narrow-band SDF bricks: 8³ cells (9³ int16 samples), band max(4 voxels, 5 cm), trilinear value and gradient |

The fixtures, each with an exact analytic reference:

- a 70° berm, 0.6 m tall, with a 3 cm overhanging lip;
- two ruts, 8 cm wide and 6 cm deep;
- a half-buried rock garden: two oriented boxes and an 80-face mesh rock.

Accuracy is measured at 10,000–20,000 random points 2–30 mm outside the solid. Normal errors use
only points 2–5 mm from the surface. Farther out, the medial axis (points equidistant from two
surfaces) makes "the" nearest normal ambiguous.

The edits measured:

- applying a rut (a capsule sweep 3 m long, 3 cm radius);
- carving an 8 cm trench through the rock garden;
- carving a 10 m curved trail;
- painting a 2 m brush stroke.

## Findings

1. **Rocks, lips and faces need explicit features.**
   - Heights alone, with rocks baked in, err by 7–17 mm at p99 (up to 27 mm), and their normals by up
     to 100°.
   - After carving a trench, that becomes 39–58 mm: a heightfield cannot reveal a rock face it never
     stored.
   - B+E is exact to float precision.
2. **Heightfield pages are the economical base.**
   - Memory: 77 MiB per km of 4 m-wide trail at 1 cm, against 2.0 GiB for SDF bricks (26×).
   - Ruts at 1 cm: 1.6 mm p99 for B against 0.7 mm for H. Ruts at 2 cm: 3.5 mm against 2.0 mm.
   - Corridor extraction: 11–32 µs for B against 89–815 µs for H.
3. **Deformation deltas on heights are cheap and accurate.**
   - A 3 m rut takes 17 µs and adds 48 KiB at 1 cm.
   - The error after the rut is 2.6 mm p99. That is the chord error of a 1 cm triangulation on a
     3 cm-radius groove, so it shrinks with spacing.
   - H's CSG subtraction takes 1.2 ms and errs by 9.3 mm. Subtracting with `max()` gives a bound,
     not a distance, until the field is recomputed.
4. **One SDF for soil and rock cannot be carved.** Subtracting a trench from a single field removes
   the rock above the cut too: 47–77 mm error. Correct carving needs a separate field per material,
   or rebaking from separate sources. Rebaking took 1.5 s at 1 cm for the trench region, and still
   rounds rock edges by about the voxel size (2.6 mm p99).
5. **Editing on heights is interactive. On bricks it is not, yet.**

   | Edit | B + K | H |
   |---|---|---|
   | Carve 10 m of curved trail | 2.1 ms | 224 ms, allocating 5,436 bricks below the old band |
   | Paint a 2 m stroke | 50 µs | H has no surface grid to paint on |

6. **H answers point queries 10–40× faster.**
   - H: 19–29 ns per query.
   - B: 150–600 ns with an exact closest-point search over triangles.
   - B+E: up to 1.1 µs when brute-forcing an 80-face rock.

   B's query path needs acceleration (block min/max hierarchy, SIMD over cells, and feature
   bounding boxes with a BVH) before it carries contact.
7. **Seams between levels:**
   - An unstitched 1 cm/2 cm seam opens a 3.5 mm crack.
   - Stitching closes it by construction.
   - The stitch moves fine edge nodes by at most 3.5 mm, which is within the coarse page's own
     interpolation error of 6.9 mm.
8. **Normal accuracy near the surface is coarse for every representation at 1–2 cm.** On the
   8 cm ruts the error is 15–30° at p99. Contact on rut walls (R1) must test its sensitivity to
   spacing.

## Decision

The hybrid is confirmed:

- B pages for the ground;
- E for rocks, roots and lips;
- F spans where layers are authored;
- K deltas for ruts and compaction.

Two amendments follow from the results:

- **B's query path is a performance item for G3.** The budget is under 100 ns per exact query at
  1 cm. The target comes from H's 20–30 ns, allowing for B's 26× memory advantage.
- **H is not a storage format.** It stays a candidate in two roles, each to be measured when its
  trigger comes:
  1. **A transient contact cache:** bricks derived from B+E for the active wheel corridor only. This
     would give H's query speed at B's memory cost. The question is whether baking the corridor
     costs less than the queries it saves (R1).
  2. **Local heavy-deformation zones:** only with a separate field per material, and with
     re-distancing after edits (R2).
