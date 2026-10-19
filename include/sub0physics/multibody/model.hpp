#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/rigid_transform.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/multibody/description.hpp"

namespace sub0physics::multibody
{

/** Parent index of a body attached to the world. */
inline constexpr std::size_t NoParent = std::numeric_limits<std::size_t>::max();

/**
 * An immutable articulated-body model: flat arrays in topological order, built by ModelBuilder.
 *
 * Body `i` has parent `parent(i) < i` (or NoParent), one joint and therefore one generalized coordinate
 * at `coordinate(i)`. Dynamics quantities are plain SI doubles in each body's frame: the joint's motion
 * subspace, the spatial inertia, and the Plücker transform from the parent's frame to the joint frame at
 * coordinate zero. The model never changes after it is built, so one model can serve many states and
 * workspaces.
 */
class Model
{
public:
    /** @return The number of bodies. */
    [[nodiscard]] std::size_t bodyCount() const noexcept
    {
        return parent_.size();
    }

    /** @return The number of generalized coordinates (degrees of freedom). */
    [[nodiscard]] std::size_t dofCount() const noexcept
    {
        return parent_.size();
    }

    /** @return Parent index of each body, NoParent for the world; always less than the body's index. */
    [[nodiscard]] std::span<const std::size_t> parents() const noexcept
    {
        return parent_;
    }

    /** @return Index of each body's coordinate in q, qd and qdd. */
    [[nodiscard]] std::span<const std::size_t> coordinates() const noexcept
    {
        return coordinate_;
    }

    /** @return Joint type of each body. */
    [[nodiscard]] std::span<const JointType> jointTypes() const noexcept
    {
        return jointType_;
    }

    /** @return Unit joint axis of each body, in its joint frame. */
    [[nodiscard]] std::span<const sub0math::Vec3d> jointAxes() const noexcept
    {
        return axis_;
    }

    /** @return Each joint's frame in its parent's frame at coordinate zero, as a pose. */
    [[nodiscard]] std::span<const sub0math::RigidTransform<double, sub0math::f64::Length>> jointPlacements() const noexcept
    {
        return placement_;
    }

    /** @return The Plücker transform from each parent's frame to the joint frame at coordinate zero. */
    [[nodiscard]] std::span<const sub0math::PluckerTransformd> treeTransforms() const noexcept
    {
        return treeTransform_;
    }

    /** @return For each joint, whether its placement rotation is the identity (most planar links). */
    [[nodiscard]] std::span<const std::uint8_t> placementUnrotated() const noexcept
    {
        return placementUnrotated_;
    }

    /** @return Motion subspace S of each joint: the body's spatial velocity per unit coordinate rate. */
    [[nodiscard]] std::span<const sub0math::SpatialMotiond> motionSubspaces() const noexcept
    {
        return motionSubspace_;
    }

    /** @return Spatial inertia of each body in its own frame. */
    [[nodiscard]] std::span<const sub0math::SpatialInertiad> inertias() const noexcept
    {
        return inertia_;
    }

    /** @return Each joint's spring reference coordinate (rad or m). */
    [[nodiscard]] std::span<const double> springReferences() const noexcept
    {
        return springReference_;
    }

    /** @return Gravity, in the world frame. */
    [[nodiscard]] const sub0math::Vec3<sub0math::f64::Acceleration>& gravity() const noexcept
    {
        return gravity_;
    }

    /** @return The world frame's spatial acceleration that stands in for gravity, -g, in plain SI. */
    [[nodiscard]] const sub0math::SpatialMotiond& baseAcceleration() const noexcept
    {
        return baseAcceleration_;
    }

    /**
     * Gives a body's name.
     *
     * @param body A body of this model.
     * @return Its name.
     */
    [[nodiscard]] std::string_view bodyName(BodyId body) const noexcept
    {
        return bodyName_[body.value()];
    }

    /**
     * Finds a body by name.
     *
     * @param name The name given to the builder.
     * @return The body, or nothing if no body has that name.
     */
    [[nodiscard]] std::optional<BodyId> findBody(std::string_view name) const noexcept;

    /** @return The number of sites. */
    [[nodiscard]] std::size_t siteCount() const noexcept
    {
        return siteBody_.size();
    }

    /**
     * Finds a site by name.
     *
     * @param name The name given to the builder.
     * @return The site, or nothing if no site has that name.
     */
    [[nodiscard]] std::optional<SiteId> findSite(std::string_view name) const noexcept;

    /**
     * Gives the body a site is fixed to.
     *
     * @param site A site of this model.
     * @return Its body.
     */
    [[nodiscard]] BodyId siteBody(SiteId site) const noexcept
    {
        return BodyId{static_cast<std::uint32_t>(siteBody_[site.value()])};
    }

    /**
     * Gives a site's frame in its body's frame.
     *
     * @param site A site of this model.
     * @return The site's pose in the body frame.
     */
    [[nodiscard]] const sub0math::RigidTransform<double, sub0math::f64::Length>& siteFrame(SiteId site) const noexcept
    {
        return siteFrame_[site.value()];
    }

private:
    friend class ModelBuilder;

    Model() = default;

    std::vector<std::size_t> parent_;
    std::vector<std::size_t> coordinate_;
    std::vector<JointType> jointType_;
    std::vector<sub0math::Vec3d> axis_;
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> placement_;
    std::vector<sub0math::PluckerTransformd> treeTransform_;
    std::vector<std::uint8_t> placementUnrotated_;
    std::vector<sub0math::SpatialMotiond> motionSubspace_;
    std::vector<sub0math::SpatialInertiad> inertia_;
    std::vector<double> springReference_;
    std::vector<std::string> bodyName_;
    std::vector<std::size_t> siteBody_;
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> siteFrame_;
    std::vector<std::string> siteName_;
    sub0math::Vec3<sub0math::f64::Acceleration> gravity_{};
    sub0math::SpatialMotiond baseAcceleration_{};
};

} // namespace sub0physics::multibody
