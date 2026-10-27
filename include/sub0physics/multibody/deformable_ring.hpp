#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <string>
#include <string_view>
#include <vector>

#include <sub0math/geometry/point.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/id.hpp"
#include "sub0physics/element/linear_spring_damper.hpp"
#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/description.hpp"

// Substructures (docs/WHEEL.md sections 2.3 and 10, item 1): structured groups of generalized coordinates attached to
// one body of the tree, which the steppers eliminate before the tree factorisation. The first kind is the deformable
// ring: N material nodes on a host body, each with displacement coordinates along fixed axes of the host, held by a
// foundation spring on each coordinate and joined to its neighbours by edges (a law on the distance between two
// nodes). Planar first: two displacement coordinates per node, in the model's plane; the description and the block
// storage are written per component so a spatial ring adds a third axis rather than a rewrite.

namespace sub0physics::multibody
{

/** Identifies a ring within one ModelBuilder and the Model it builds; numbered as added. */
using RingId = Id<struct RingIdTag, std::uint32_t>;

/** Displacement coordinates per ring node in the planar formulation. */
inline constexpr std::size_t PlanarRingComponents = 2;

/**
 * A deformable ring attached to a host body (docs/WHEEL.md section 2.3).
 *
 * **Coordinates.** Node `i` has rest position `r_i` in the host's frame and displacement coordinates `w_i` along
 * `displacementAxes`: its position in the host frame is `r_i + w_i[0] a_0 + w_i[1] a_1`. The ring's coordinates follow
 * the tree's in `State`, in node order and per node in axis order, and never depend on contact. A node is a point
 * mass with its own rotational inertia about the plane's normal; it translates with its coordinates and rotates with
 * the host (the reference's node bodies on two slide joints).
 *
 * **Foundation.** `foundation[i]` acts on each displacement coordinate of node `i` independently, relaxed at zero:
 * energy `sum 1/2 k_i |w_i|^2` for a linear law. **Edges.** Edge `e` joins node `e` to node `e + 1 mod N`;
 * `edgeLaws[e]` acts on its extension `l_e - edgeRestLength[e]` (the law's coordinate grows with the length, like an
 * extension tendon), with the rate `dl_e/dt` as the law's rate. Coefficients are per node and per edge, so a
 * non-uniform ring needs nothing more.
 *
 * Built by ModelBuilder::addRing; validated at build (`std::expected`). The vectors are copied.
 */
struct DeformableRingDescription
{
    std::string_view name;                       ///< Unique among rings, non-empty; copied by the builder.
    BodyId host;                                 ///< The body the ring moves with (a hub).
    std::array<sub0math::Vec3d, PlanarRingComponents> displacementAxes{
        sub0math::Vec3d{1.0, 0.0, 0.0}, sub0math::Vec3d{0.0, 0.0, 1.0}}; ///< Orthonormal axes in the host frame, both in the model's plane.
    std::vector<sub0math::Point3<sub0math::f64::Length>> restPositions; ///< N >= 3 rest positions in the host frame, all at the same out-of-plane offset.
    std::vector<sub0math::f64::Mass> nodeMass;                    ///< N positive node masses.
    std::vector<sub0math::f64::MomentOfInertia> nodeInertia;      ///< N non-negative node inertias about the plane's normal, through the node.
    std::vector<LinearSpringDamper<sub0math::f64::Length>> foundation; ///< N foundation laws, one per node (on each of its coordinates).
    std::vector<LinearSpringDamper<sub0math::f64::Length>> edgeLaws;   ///< N edge laws; edge e joins node e to node e + 1 mod N.
    std::vector<sub0math::f64::Length> edgeRestLength;            ///< N positive rest lengths: where each edge law's coordinate is zero.
    sub0math::f64::Length skinRadius{};                           ///< The capsule radius of every edge (contact geometry, M2b); not negative.
    std::vector<std::uint32_t> nodeIds;                           ///< N unique stable material ids, or empty for 0..N-1.
    std::vector<std::uint32_t> edgeIds;                           ///< N unique stable material ids, or empty for 0..N-1.
    std::uint32_t topologyRevision = 0;                           ///< Revision of the ring's material topology; contact records refer to it.
};

/**
 * The sine and cosine of a canonical ring direction, `2 pi index / nodes`, from the binary64 table (no runtime
 * trigonometry: an ulp of `sin` differs between runtimes, and contact amplifies it).
 */
struct CanonicalRingDirection
{
    double sine;   ///< sin(2 pi index / nodes).
    double cosine; ///< cos(2 pi index / nodes).
};

/**
 * Gives a canonical ring direction from the quarter-circle table `sin(i pi / 64)`, i = 0..32, reflected for exact
 * symmetry. Exact table values, identical on every platform.
 *
 * @param nodes Node count: a divisor of 128 of at least 4 (4, 8, 16, 32, 64, 128).
 * @param index Node index, taken modulo `nodes`.
 * @return The direction, or InvalidInput for an unsupported node count.
 */
[[nodiscard]] std::expected<CanonicalRingDirection, Fault> canonicalRingDirection(std::size_t nodes,
                                                                                  std::size_t index) noexcept;

/**
 * Gives the canonical chord ratio `2 sin(pi / nodes)`: an edge's rest chord over the centreline radius.
 *
 * @param nodes Node count: a divisor of 128 of at least 4.
 * @return The ratio, or InvalidInput for an unsupported node count.
 */
[[nodiscard]] std::expected<double, Fault> canonicalChordRatio(std::size_t nodes) noexcept;

/**
 * A uniform ring in the reference's conventions (docs/WHEEL.md section 1): N nodes on the centreline radius R, a
 * foundation of `k ds` per coordinate (`ds = 2 pi R / N`), edges of stiffness `EA / L0` and viscosity
 * `c / L0` with `L0 = R chordRatio(N)`, node mass `m / N` and node inertia `(m / N) skin^2 / 2`. Node `i` rests at
 * `R (sin t_i a_0 - cos t_i a_1)`, `t_i = 2 pi i / N`: node 0 at `-R a_1`, advancing towards `+a_0`.
 */
struct UniformRingSpec
{
    std::string_view name;            ///< The ring's name.
    BodyId host;                      ///< The host body.
    std::size_t nodes = 64;           ///< Node count: a divisor of 128 of at least 4.
    sub0math::f64::Length centrelineRadius;        ///< R: the outer radius minus the skin; positive.
    sub0math::f64::Length skinRadius;              ///< The capsule radius; not negative.
    sub0math::f64::Mass ringMass;                  ///< Mass of all nodes together; positive.
    sub0math::f64::Pressure foundationStiffness;   ///< k, per unit centreline length [N/m^2]; not negative.
    sub0math::f64::Force extensionStiffness;       ///< EA [N]; not negative.
    sub0math::f64::Impulse edgeViscosity;          ///< c [N s]; not negative.
    std::array<sub0math::Vec3d, PlanarRingComponents> displacementAxes{
        sub0math::Vec3d{1.0, 0.0, 0.0}, sub0math::Vec3d{0.0, 0.0, 1.0}}; ///< As DeformableRingDescription.
};

/**
 * Builds the description of a uniform ring from canonical tables.
 *
 * Each edge's rest length is the distance between its rest nodes computed from the canonical positions; the
 * coefficients use the canonical chord `R chordRatio(N)` as BikeLab's reference does (the two agree to an ulp or so).
 *
 * @param spec The ring.
 * @return The description, or InvalidInput (node count, non-finite or negative values).
 */
[[nodiscard]] std::expected<DeformableRingDescription, Fault> uniformRing(const UniformRingSpec& spec);

/**
 * A compiled ring (Model::rings): flat per-node and per-edge arrays in the planar formulation's plane axes.
 */
struct DeformableRing
{
    RingId id;                     ///< The ring's id.
    std::string name;              ///< Its name.
    std::size_t host = 0;          ///< Host body index.
    std::size_t firstCoordinate = 0; ///< Index of node 0's first coordinate in q.
    std::size_t nodeCount = 0;     ///< N.
    std::vector<std::size_t> hostChain; ///< The host and its ancestors, host first: the bodies whose joints move it.
    std::array<double, 2> axis0{}; ///< a_0 in the host frame's plane axes (i, j).
    std::array<double, 2> axis1{}; ///< a_1 in the plane axes.
    double normalOffset = 0.0;     ///< The nodes' out-of-plane offset in the host frame [m].
    std::vector<double> restX;     ///< Rest position along the plane's i axis, host frame [m].
    std::vector<double> restY;     ///< Rest position along the plane's j axis [m].
    std::vector<double> rest0;     ///< Rest position along a_0 [m].
    std::vector<double> rest1;     ///< Rest position along a_1 [m].
    std::vector<double> mass;      ///< Node masses [kg].
    std::vector<double> inertia;   ///< Node inertias about the normal [kg m^2].
    std::vector<LinearSpringDamper<sub0math::f64::Length>> foundation; ///< Per node.
    std::vector<LinearSpringDamper<sub0math::f64::Length>> edgeLaws;   ///< Per edge.
    std::vector<double> edgeRestLength; ///< Per edge [m].
    double skinRadius = 0.0;       ///< [m].
    std::vector<std::uint32_t> nodeIds; ///< Material ids per node.
    std::vector<std::uint32_t> edgeIds; ///< Material ids per edge.
    std::uint32_t topologyRevision = 0; ///< As described.

    /** @return The ring's coordinate count: two per node. */
    [[nodiscard]] std::size_t coordinateCount() const noexcept
    {
        return PlanarRingComponents * nodeCount;
    }
};

} // namespace sub0physics::multibody
