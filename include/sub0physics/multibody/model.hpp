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

#include "sub0physics/multibody/constraint.hpp"
#include "sub0physics/multibody/deformable_ring.hpp"
#include "sub0physics/multibody/description.hpp"
#include "sub0physics/multibody/force_element.hpp"
#include "sub0physics/multibody/planar.hpp"

namespace sub0physics::multibody
{

/** Parent index of a body attached to the world. */
inline constexpr std::size_t NoParent = std::numeric_limits<std::size_t>::max();

/**
 * How a model's dynamics are evaluated (docs/MULTIBODY.md, M1d status). Both give the same results to rounding;
 * the planar one is cheaper.
 */
enum class Formulation : std::uint8_t
{
    Automatic, ///< A builder request only: Planar when the model qualifies, otherwise General (the default).
    General,   ///< The 3D spatial formulation, valid for every model.
    Planar     ///< Three-component planar algebra (planar.hpp), for a model whose motion stays in one coordinate plane.
};

/**
 * Which Jacobians a linearly implicit step treats implicitly: a discrete lever (docs/LOD.md) between the reference
 * stepper and MuJoCo's `implicitfast` structure. Used by stepLinearlyImplicitEuler and stepLinearlyImplicitMidpoint.
 */
enum class ImplicitTerms : std::uint8_t
{
    Full,              ///< D and K, including the tendons' geometric and rate-transport terms (default).
    MaterialStiffness, ///< D and the material part of K only (geometric stiffness off).
    DampingOnly        ///< D only: springs and stiffness explicit, as MuJoCo's `implicitfast`.
};

/**
 * The planar formulation's view of one body, in its plane's (i, j, k) axes: Model::planar.
 */
struct PlanarBody
{
    PlanarTransform placement; ///< The joint placement at coordinate zero: rotation about k and in-plane offset.
    bool unrotated;            ///< Whether the placement rotation is the identity.
    bool revolute;             ///< Joint type: revolute (about +-k) or prismatic (in the plane).
    double axisSign;           ///< Revolute: +1 for an axis along +k, -1 along -k.
    double axisX;              ///< Prismatic: the axis's i component, in the joint frame.
    double axisY;              ///< Prismatic: the axis's j component, in the joint frame.
    PlanarInertia inertia;     ///< In-plane mass properties about the body origin.
};

/**
 * What the planar formulation needs from a model (Formulation::Planar): the plane, and the in-plane view of every
 * body and site. A model qualifies when every revolute axis lies along one coordinate axis k (either sense), every
 * prismatic axis is perpendicular to it, every placement rotation is about k, and gravity has no k component; the
 * tests are exact, so a model built from in-plane numbers qualifies and one that is out of plane by rounding does
 * not. Out-of-plane offsets (placements, centres of mass, sites along k) are allowed: they never change the in-plane
 * motion, and tendon lengths include the constant k separation of their sites.
 */
struct PlanarModel
{
    std::size_t normal = 2;          ///< k: the axis normal to the plane (0 x, 1 y, 2 z).
    std::size_t axisI = 0;           ///< i: the plane's first axis.
    std::size_t axisJ = 1;           ///< j: the plane's second axis.
    double gravityX = 0.0;           ///< Gravity along i [m/s^2].
    double gravityY = 0.0;           ///< Gravity along j [m/s^2].
    std::vector<PlanarBody> bodies;  ///< One per body.
    std::vector<double> siteX;       ///< Each site's i offset in its body frame (world i for a world site) [m].
    std::vector<double> siteY;       ///< Each site's j offset [m].
    std::vector<double> siteNormal;  ///< Each site's world k coordinate, constant in planar motion [m].
};

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
    /**
     * Identifies the topology this model was built with. Every build gets a new id; copies share it.
     * A Workspace is bound to it so that cached structure (such as the zero entries of the mass matrix)
     * cannot be reused with a different tree.
     *
     * @return The id; never zero.
     */
    [[nodiscard]] std::uint64_t id() const noexcept
    {
        return id_;
    }

    /** @return The number of bodies. */
    [[nodiscard]] std::size_t bodyCount() const noexcept
    {
        return parent_.size();
    }

    /**
     * The number of generalized coordinates (degrees of freedom): the tree's, one per body, then every ring's
     * (docs/WHEEL.md section 2.3), in the order the rings were added. `State` holds `[q; qd]` of this size.
     *
     * @return treeDofCount() plus two per ring node.
     */
    [[nodiscard]] std::size_t dofCount() const noexcept
    {
        return parent_.size() + ringCoordinateCount_;
    }

    /**
     * The number of the tree's coordinates: one per body, the first ones in q. The dense matrices of the Workspace
     * (mass matrix, element Jacobians, step matrix) are this size; ring coordinates live in block storage.
     *
     * @return bodyCount().
     */
    [[nodiscard]] std::size_t treeDofCount() const noexcept
    {
        return parent_.size();
    }

    /** @return The compiled rings, in the order added (docs/WHEEL.md section 2.3). */
    [[nodiscard]] std::span<const DeformableRing> rings() const noexcept
    {
        return rings_;
    }

    /**
     * Finds a ring by name.
     *
     * @param name The name given to the builder.
     * @return The ring, or nothing if no ring has that name.
     */
    [[nodiscard]] std::optional<RingId> findRing(std::string_view name) const noexcept;

    /** @return Whether anything contributes element forces: a force element or a ring (its foundation and edges). */
    [[nodiscard]] bool hasForceElements() const noexcept
    {
        return !elementName_.empty() || !rings_.empty();
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

    /**
     * The joint placements: the Plücker transform from each parent's frame to the joint frame at
     * coordinate zero. It is the only stored form of a placement; `toRigidTransform` recovers the pose.
     *
     * @return One transform per body.
     */
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
     * @return Its body, or WorldBody for a site fixed in the world.
     */
    [[nodiscard]] BodyId siteBody(SiteId site) const noexcept
    {
        const std::size_t body = siteBody_[site.value()];
        return body == NoParent ? WorldBody : BodyId{static_cast<std::uint32_t>(body)};
    }

    /** @return Body index of each site, NoParent for a site fixed in the world. */
    [[nodiscard]] std::span<const std::size_t> siteBodies() const noexcept
    {
        return siteBody_;
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

    /** @return The number of force elements, of every kind. */
    [[nodiscard]] std::size_t elementCount() const noexcept
    {
        return elementName_.size();
    }

    /** @return The number of tendon elements; each reports a wrench at its two sites. */
    [[nodiscard]] std::size_t tendonCount() const noexcept
    {
        return tendonCount_;
    }

    /**
     * Whether the elements' Jacobians stay within the tree's ancestor pairs, so `M - h D - h^2 K` keeps the
     * mass matrix's sparsity and factorises without fill-in. True unless a tendon joins two branches.
     *
     * @return True when every tendon's ends lie on one chain from the root.
     */
    [[nodiscard]] bool elementsFollowTree() const noexcept
    {
        return elementsFollowTree_;
    }

    /**
     * Finds a force element by name.
     *
     * @param name The name given to the builder.
     * @return The element, or nothing if no element has that name.
     */
    [[nodiscard]] std::optional<ElementId> findElement(std::string_view name) const noexcept;

    /** @return The compiled coordinate elements, one list per law, in evaluation order. */
    [[nodiscard]] const CoordinateElementLists& coordinateElements() const noexcept
    {
        return coordinateElements_;
    }

    /** @return The compiled tendon elements, one list per law, in evaluation order. */
    [[nodiscard]] const TendonElementLists& tendonElements() const noexcept
    {
        return tendonElements_;
    }

    /** @return The number of constraints, of every kind, as described (a two-sided limit is one). */
    [[nodiscard]] std::size_t constraintCount() const noexcept
    {
        return constraintName_.size();
    }

    /**
     * Finds a constraint by name.
     *
     * @param name The name given to the builder.
     * @return The constraint, or nothing if no constraint has that name.
     */
    [[nodiscard]] std::optional<ConstraintId> findConstraint(std::string_view name) const noexcept;

    /** @return The compiled loop closures, in the order added. */
    [[nodiscard]] std::span<const LoopClosure> closures() const noexcept
    {
        return closures_;
    }

    /** @return The compiled one-sided limits, in the order added (a two-sided limit gives its lower one first). */
    [[nodiscard]] std::span<const UnilateralLimit> limits() const noexcept
    {
        return limits_;
    }

    /**
     * The constraint layout: every closure's rows (3 or 6 each, in order), then one row per one-sided limit. The
     * workspace's per-row results (constraintImpulse, constraintForce) follow it.
     *
     * @return The number of rows.
     */
    [[nodiscard]] std::size_t constraintRowCount() const noexcept
    {
        return closureRowCount_ + limits_.size();
    }

    /** @return The number of closure rows, the first rows of the layout. */
    [[nodiscard]] std::size_t closureRowCount() const noexcept
    {
        return closureRowCount_;
    }

    /** @return The number of constraint site wrenches: two per closure and per one-sided distance limit. */
    [[nodiscard]] std::size_t constraintWrenchCount() const noexcept
    {
        return constraintWrenchCount_;
    }

    /** @return How the dynamics are evaluated: General or Planar (never Automatic). */
    [[nodiscard]] Formulation formulation() const noexcept
    {
        return formulation_;
    }

    /** @return The planar formulation's view of the model; meaningful only when formulation() is Planar. */
    [[nodiscard]] const PlanarModel& planar() const noexcept
    {
        return planar_;
    }

    /** @return Whether some body is massless (it only carries a joint, like a body with several joints in MJCF). */
    [[nodiscard]] bool hasMasslessBodies() const noexcept
    {
        return hasMasslessBodies_;
    }

    /** @return The constraint solver's levers, set with ModelBuilder::setConstraintSettings. */
    [[nodiscard]] const ConstraintSettings& constraintSettings() const noexcept
    {
        return constraintSettings_;
    }

private:
    friend class ModelBuilder;

    Model() = default;

    std::uint64_t id_ = 0;
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> coordinate_;
    std::vector<JointType> jointType_;
    std::vector<sub0math::Vec3d> axis_;
    std::vector<sub0math::PluckerTransformd> treeTransform_;
    std::vector<std::uint8_t> placementUnrotated_;
    std::vector<sub0math::SpatialMotiond> motionSubspace_;
    std::vector<sub0math::SpatialInertiad> inertia_;
    std::vector<double> springReference_;
    std::vector<std::string> bodyName_;
    std::vector<std::size_t> siteBody_;
    std::vector<sub0math::RigidTransform<double, sub0math::f64::Length>> siteFrame_;
    std::vector<std::string> siteName_;
    CoordinateElementLists coordinateElements_;
    TendonElementLists tendonElements_;
    std::vector<std::string> elementName_;
    std::size_t tendonCount_ = 0;
    std::vector<LoopClosure> closures_;
    std::vector<UnilateralLimit> limits_;
    std::vector<std::string> constraintName_;
    std::size_t closureRowCount_ = 0;
    std::size_t constraintWrenchCount_ = 0;
    ConstraintSettings constraintSettings_{};
    std::vector<DeformableRing> rings_;
    std::size_t ringCoordinateCount_ = 0;
    bool elementsFollowTree_ = true;
    bool hasMasslessBodies_ = false;
    Formulation formulation_ = Formulation::General;
    PlanarModel planar_{};
    sub0math::Vec3<sub0math::f64::Acceleration> gravity_{};
    sub0math::SpatialMotiond baseAcceleration_{};
};

} // namespace sub0physics::multibody
