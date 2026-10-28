#pragma once

#include <span>

// Pieces of the cheaper casing models (docs/WHEEL.md section 2.4 and "Wheel-model study"): the tread's contact laws on a
// terrain face between a casing edge's two nodes, and the state transfer between a nodal casing and a rigid ring that
// a model switch needs (docs/LOD.md, decide/transfer/commit). Planar, plain SI doubles, never allocating.

namespace sub0physics::multibody
{

/**
 * The response of a tread law on one casing edge against one terrain face: forces at the edge's two nodes along the
 * face normal (pushing the casing away from the face, never negative), the law's potential and dissipation.
 *
 * `force0 rate0 + force1 rate1` is the power the contact takes from the casing, and it equals
 * `dV/dt + dissipationPower` (rates are the end penetrations' rates).
 */
struct TreadFaceResponse
{
    double force0 = 0.0;           ///< Normal force at the edge's first node [N].
    double force1 = 0.0;           ///< At its second node [N].
    double elastic0 = 0.0;         ///< The spring's share of force0 (before the damper and the clipping) [N].
    double elastic1 = 0.0;         ///< Likewise for force1 [N].
    double potential = 0.0;        ///< The tread's elastic energy [J].
    double dissipationPower = 0.0; ///< Damper and clipping losses [W]; never negative.
    double penetration = 0.0;      ///< The law's penetration (its deflection measure, for the tread's range guard) [m].
};

/**
 * The witness law with its maximum smoothed (C1).
 *
 * A capsule edge against a face is penetrated by `max(p0, p1)`, the deeper of its two ends' penetrations below the face
 * line; the reference contact puts the whole force `k max(p0, p1)` at the deeper node, so as the edge rolls past level
 * the force jumps from one node to the other. Here the maximum is smoothed within a band:
 * `m = (p0 + p1) / 2 + h(|p0 - p1| / 2)`, `h(x) = x` for `x >= band` and `x^2 (2 band - x) / band^2` within it, so
 * `h <= x`, `h` and `h'` are continuous and the law is the witness law exactly outside the band. The potential is
 * `k m^2 / 2`, the force `max(0, k m + c m')` shared by `dm/dp0`, `dm/dp1` (continuous, summing to one).
 *
 * @param stiffness k per contact [N/m]; positive.
 * @param damping c per contact [N s/m]; not negative.
 * @param band The band's half-width in end-penetration difference [m]; positive.
 * @param p0 The first node's penetration below the face line (skin less its distance) [m].
 * @param p1 The second node's [m].
 * @param rate0 dp0/dt [m/s].
 * @param rate1 dp1/dt [m/s].
 * @return The response.
 */
[[nodiscard]] TreadFaceResponse smoothedTreadFace(double stiffness, double damping, double band, double p0, double p1,
                                                  double rate0, double rate1) noexcept;

/**
 * The tread spread along the edge: a spring k and a dashpot c per unit edge parameter over the penetrating part,
 * `p(s) = (1 - s) p0 + s p1`, integrated in closed form (a uniform overlap p gives k p, as the witness law on a level edge;
 * an edge with one end in is softer). Continuous and C1 in the end penetrations.
 *
 * @param stiffness k [N/m]; positive.
 * @param damping c [N s/m]; not negative.
 * @param p0 The first node's penetration below the face line [m].
 * @param p1 The second node's [m].
 * @param rate0 dp0/dt [m/s].
 * @param rate1 dp1/dt [m/s].
 * @return The response.
 */
[[nodiscard]] TreadFaceResponse distributedTreadFace(double stiffness, double damping, double p0, double p1,
                                                     double rate0, double rate1) noexcept;

/** A planar rigid body's pose and motion: position (x, z), angle about +y, and their rates. */
struct PlanarRigidState
{
    double x = 0.0;
    double z = 0.0;
    double angle = 0.0;
    double vx = 0.0;
    double vz = 0.0;
    double rate = 0.0;
};

/**
 * A nodal casing on a host: node world positions and velocities (x, z), node masses, rest positions in the host's frame,
 * and the host's angle and spin rate. The nodes' own rotational inertia (`nodeSpin`, their sum) turns with the host.
 */
struct NodalCasingView
{
    std::span<const double> x, z, vx, vz; ///< Node positions [m] and velocities [m/s], world.
    std::span<const double> mass;          ///< Node masses [kg].
    std::span<const double> restX, restZ;  ///< Rest positions in the host frame [m].
    double hostAngle = 0.0;                ///< [rad]; rotation about +y: x' = c x + s z, z' = -s x + c z.
    double hostRate = 0.0;                 ///< [rad/s].
    double nodeSpin = 0.0;                 ///< Sum of the nodes' own inertias [kg m^2].
};

/**
 * The rigid ring that carries a nodal casing's mean motion: its centre and velocity the nodes' (mass-weighted), its
 * angle the host's plus the best-fit rotation of the rest shape onto the nodes (Procrustes, `det::atan2`), and its spin
 * rate from the nodes' angular momentum about their centre plus their own spin (`nodeSpin` times the host's rate).
 *
 * Linear and angular momentum are kept exactly (to rounding) with a ring of mass `sum m` and spin inertia
 * `ringInertia` (`sum m |r|^2 + nodeSpin` for the rest shape); the kinetic energy of the nodes' deformation is the
 * transfer's loss, never a gain (the mass-weighted projection onto rigid motions). The host's own state is unchanged,
 * and its spin inertia loses `nodeSpin` to the ring.
 *
 * @param casing The nodal casing.
 * @param ringInertia The rigid ring's spin inertia about its centre [kg m^2]; positive.
 * @return The ring's state.
 */
[[nodiscard]] PlanarRigidState fitRigidRing(const NodalCasingView& casing, double ringInertia) noexcept;

/**
 * Places a nodal casing on a rigid ring exactly: node `i` at `centre + R(angle) rest_i`, moving with the ring. The
 * nodes' own spin, which turns with the host in the nodal casing, joins the host: the returned host rate keeps the
 * angular momentum, `(hostInertia hostRate + nodeSpin ringRate) / (hostInertia + nodeSpin)` (a merge, so the energy can
 * only fall, by the rate difference's share).
 *
 * @param ring The rigid ring.
 * @param restX Rest positions in the host frame [m].
 * @param restZ Likewise.
 * @param hostInertia The host's spin inertia without the nodes [kg m^2].
 * @param hostRate The host's spin rate [rad/s].
 * @param nodeSpin The nodes' own inertia, summed [kg m^2].
 * @param x Receive the node positions [m].
 * @param z Likewise.
 * @param vx Receive the node velocities [m/s].
 * @param vz Likewise.
 * @return The host's spin rate for the nodal casing [rad/s].
 */
[[nodiscard]] double placeNodalCasing(const PlanarRigidState& ring, std::span<const double> restX,
                                      std::span<const double> restZ, double hostInertia, double hostRate,
                                      double nodeSpin, std::span<double> x, std::span<double> z, std::span<double> vx,
                                      std::span<double> vz) noexcept;

} // namespace sub0physics::multibody
