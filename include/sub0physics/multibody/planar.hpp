#pragma once

#include <cstdint>

#include <sub0math/config.hpp>

// Planar spatial algebra for models whose motion stays in one coordinate plane (docs/MULTIBODY.md, M1d status:
// the planar formulation). A plane is named by its normal, one of the coordinate axes k; its in-plane axes are
// (i, j) = (k + 1, k + 2) mod 3, so (i, j, k) is right-handed: (x, y) for a normal along z, (z, x) along y and
// (y, z) along x. Vectors are three components instead of six (Featherstone, RBDA section 2.15): a motion is an
// angular rate about k and a linear velocity in the plane; a force is a moment about k and a force in the plane.
// Plain SI doubles, like the 6D vocabulary in Sub0Math.

namespace sub0physics::multibody
{

/** A planar motion vector: angular rate about the plane's normal and the velocity of the frame origin. */
struct PlanarMotion
{
    double w; ///< Angular rate about k [rad/s] (or its rate of change).
    double x; ///< Linear component along the plane's first axis i [m/s].
    double y; ///< Linear component along the plane's second axis j [m/s].
};

/** A planar force vector: moment about the plane's normal through the frame origin, and the force in the plane. */
struct PlanarForce
{
    double n; ///< Moment about k [N m].
    double x; ///< Force along i [N].
    double y; ///< Force along j [N].
};

/** A planar rigid-body inertia about a frame origin. */
struct PlanarInertia
{
    double m;  ///< Mass [kg].
    double hx; ///< First moment along i: mass times the centre of mass's i coordinate [kg m].
    double hy; ///< First moment along j [kg m].
    double i;  ///< Rotational inertia about k through the origin [kg m^2].
};

/**
 * A planar coordinate transform from frame A to frame B: B is A rotated by an angle about k (cosine `c`, sine `s`)
 * with its origin at (x, y) in A's axes. Applied to a motion it gives the same motion in B's frame.
 */
struct PlanarTransform
{
    double c; ///< Cosine of B's rotation relative to A.
    double s; ///< Sine of B's rotation relative to A.
    double x; ///< B's origin along A's i axis [m].
    double y; ///< B's origin along A's j axis [m].
};

/** @return The motion `m`, given in A's frame, in B's frame. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarMotion apply(const PlanarTransform& t, const PlanarMotion& m) noexcept
{
    // Velocity of B's origin (v + w k x r), then rotated into B's axes (by minus the angle).
    const double ux = m.x - t.y * m.w;
    const double uy = m.y + t.x * m.w;
    return {m.w, t.c * ux + t.s * uy, t.c * uy - t.s * ux};
}

/** @return The force `f`, given in B's frame, in A's frame (the transpose of `apply`, child to parent). */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarForce applyTranspose(const PlanarTransform& t, const PlanarForce& f) noexcept
{
    const double fx = t.c * f.x - t.s * f.y;
    const double fy = t.s * f.x + t.c * f.y;
    return {f.n + t.x * fy - t.y * fx, fx, fy};
}

/** @return The inertia `in`, given in B's frame, in A's frame (`X^T I X`, child to parent). */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarInertia applyTranspose(const PlanarTransform& t, const PlanarInertia& in) noexcept
{
    const double yx = t.c * in.hx - t.s * in.hy;
    const double yy = t.s * in.hx + t.c * in.hy;
    return {in.m, yx + in.m * t.x, yy + in.m * t.y,
            in.i + 2.0 * (t.x * yx + t.y * yy) + in.m * (t.x * t.x + t.y * t.y)};
}

/** @return `A_to_C = (B_to_C)(A_to_B)`: `outer` after `inner`. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarTransform compose(const PlanarTransform& outer, const PlanarTransform& inner) noexcept
{
    // C's origin in A: inner's origin plus outer's origin rotated into A's axes.
    return {inner.c * outer.c - inner.s * outer.s, inner.s * outer.c + inner.c * outer.s,
            inner.x + inner.c * outer.x - inner.s * outer.y, inner.y + inner.s * outer.x + inner.c * outer.y};
}

/** @return The sum of two motions. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarMotion operator+(const PlanarMotion& a, const PlanarMotion& b) noexcept
{
    return {a.w + b.w, a.x + b.x, a.y + b.y};
}

/** @return The sum of two forces. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarForce operator+(const PlanarForce& a, const PlanarForce& b) noexcept
{
    return {a.n + b.n, a.x + b.x, a.y + b.y};
}

/** @return The sum of two inertias about the same origin. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarInertia operator+(const PlanarInertia& a, const PlanarInertia& b) noexcept
{
    return {a.m + b.m, a.hx + b.hx, a.hy + b.hy, a.i + b.i};
}

/** @return The momentum `I m`. */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarForce operator*(const PlanarInertia& in, const PlanarMotion& m) noexcept
{
    return {in.i * m.w + in.hx * m.y - in.hy * m.x, in.m * m.x - in.hy * m.w, in.m * m.y + in.hx * m.w};
}

/** @return `v x m` for motions (crm). */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarMotion crossMotion(const PlanarMotion& v, const PlanarMotion& m) noexcept
{
    return {0.0, v.y * m.w - v.w * m.y, v.w * m.x - v.x * m.w};
}

/** @return `v x* f` for a motion and a force (crf). */
[[nodiscard]] SUB0MATH_FORCE_INLINE PlanarForce crossForce(const PlanarMotion& v, const PlanarForce& f) noexcept
{
    return {v.x * f.y - v.y * f.x, -v.w * f.y, v.w * f.x};
}

/** @return `m . f`: power of a force on a motion. */
[[nodiscard]] SUB0MATH_FORCE_INLINE double dot(const PlanarMotion& m, const PlanarForce& f) noexcept
{
    return m.w * f.n + m.x * f.x + m.y * f.y;
}

} // namespace sub0physics::multibody
