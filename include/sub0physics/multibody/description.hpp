#pragma once

#include <cstdint>
#include <limits>
#include <string_view>

#include <sub0math/geometry/point.hpp>
#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/id.hpp"

// The typed model description that ModelBuilder compiles into a Model (docs/MULTIBODY.md section 2).
// These aggregates are the asset that must not be locked in: they say what a mechanism is, in Sub0Math
// quantities, and nothing about how it is solved.

namespace sub0physics::multibody
{

/** Identifies a body within one ModelBuilder and the Model it builds; bodies are numbered as added. */
using BodyId = Id<struct BodyIdTag, std::uint32_t>;

/** Identifies a site within one ModelBuilder and the Model it builds. */
using SiteId = Id<struct SiteIdTag, std::uint32_t>;

/** The fixed world frame, as a joint's parent: the one sentinel, `BodyId::invalid()`. */
inline constexpr BodyId WorldBody = BodyId::invalid();

/** The kinds of joint the core supports; each has one generalized coordinate. */
enum class JointType : std::uint8_t
{
    Revolute, ///< Rotation about the axis; coordinate in radians, positive right-handed.
    Prismatic ///< Translation along the axis; coordinate in metres.
};

/**
 * The joint that connects a body to its parent.
 *
 * The body's frame is the joint's frame: at coordinate zero it sits at `placement` in the parent's
 * frame, and the coordinate rotates it about, or translates it along, `axis`.
 */
struct JointDescription
{
    JointType type = JointType::Revolute;                       ///< Kind of motion.
    BodyId parent = WorldBody;                                  ///< Parent body, or WorldBody.
    sub0math::RigidTransform<double, sub0math::f64::Length> placement{}; ///< Joint frame in the parent's frame at zero.
    sub0math::Vec3d axis{0.0, 0.0, 1.0};                        ///< Unit axis in the joint frame.
    double springReference = 0.0; ///< Coordinate at which a coordinate spring is relaxed (rad or m); used from M1b.
};

/** Rotational inertia about a body's centre of mass, in the body frame's axes. */
struct RotationalInertia
{
    sub0math::Vec3<sub0math::f64::MomentOfInertia> moments;    ///< Ixx, Iyy, Izz.
    sub0math::Vec3<sub0math::f64::MomentOfInertia> products{}; ///< Ixy, Ixz, Iyz (tensor entries, not negated).
};

/** A rigid body and the joint that attaches it to its parent. */
struct BodyDescription
{
    std::string_view name;                                 ///< Unique, non-empty; copied by the builder.
    sub0math::f64::Mass mass;                              ///< Positive.
    sub0math::Point3<sub0math::f64::Length> centreOfMass;  ///< In the body frame.
    RotationalInertia inertia;                             ///< About the centre of mass; positive definite, principal moments obeying the triangle inequality.
    JointDescription joint;                                ///< Connection to the parent.
};

/** A named frame fixed to a body: an attachment, an element endpoint or a load-reporting point. */
struct SiteDescription
{
    std::string_view name;                                           ///< Unique, non-empty; copied by the builder.
    BodyId body;                                                     ///< The body it moves with.
    sub0math::RigidTransform<double, sub0math::f64::Length> frame{}; ///< Site frame in the body frame.
};

} // namespace sub0physics::multibody
