#pragma once

#include <cstdint>
#include <expected>
#include <string>
#include <string_view>
#include <tuple>
#include <variant>
#include <vector>

#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/constants.hpp"
#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/constraint.hpp"
#include "sub0physics/multibody/description.hpp"
#include "sub0physics/multibody/force_element.hpp"
#include "sub0physics/multibody/model.hpp"

namespace sub0physics::multibody
{

/**
 * Collects a typed model description and compiles it into an immutable Model.
 *
 * Bodies are numbered in the order they are added, and a joint's parent must already exist, so the
 * numbering is a topological order of the tree. Nothing is checked until `build`, which validates the
 * whole description at once: positive mass, positive definite inertia whose principal moments obey the
 * triangle inequality, unit axes and rotations, finite values, unique non-empty names, resolved
 * references, force elements whose law suits what it is bound to (a law on a length on a prismatic
 * joint or a tendon, a law on an angle on a revolute joint), and constraints with resolved sites on two
 * different bodies, ordered finite bounds and restitution in [0, 1].
 *
 * ```cpp
 * using namespace sub0math;
 * using namespace sub0math::f64::literals;
 * ModelBuilder builder;
 * const auto inertia = kilogramMetresSquared(0.01);
 * const BodyId arm = builder.addBody({.name = "arm", .mass = 1.0_kg, .centreOfMass = {0_m, -0.5_m, 0_m},
 *                                     .inertia = {{inertia, inertia, inertia}},
 *                                     .joint = {.type = JointType::Revolute, .axis = {0.0, 0.0, 1.0}}});
 * std::expected<Model, Fault> model = builder.build();
 * ```
 */
class ModelBuilder
{
public:
    /**
     * Sets gravity; the default is standard gravity along -z.
     *
     * @param gravity Gravitational acceleration in the world frame.
     */
    void setGravity(const sub0math::Vec3<sub0math::f64::Acceleration>& gravity) noexcept
    {
        gravity_ = gravity;
    }

    /**
     * Adds a body and the joint that attaches it to its parent.
     *
     * @param body The body; its joint's parent must be WorldBody or a body added earlier.
     * @return The new body's id, valid for this builder and the models it builds.
     */
    BodyId addBody(const BodyDescription& body);

    /**
     * Adds a site on a body.
     *
     * @param site The site.
     * @return The new site's id, valid for this builder and the models it builds.
     */
    SiteId addSite(const SiteDescription& site);

    /**
     * Adds a force element on a joint coordinate.
     *
     * @tparam L The law.
     * @param element The element.
     * @return The new element's id, valid for this builder and the models it builds.
     */
    template<ElementLaw L>
    ElementId addElement(const CoordinateElementDescription<L>& element)
    {
        const ElementId id = nextElement(element.name);
        auto& stored = std::get<std::vector<Stored<CoordinateElementDescription<L>>>>(coordinateElements_);
        stored.push_back({id, element});
        stored.back().description.name = {}; // the builder's name list owns the name
        return id;
    }

    /**
     * Adds a force element between two sites.
     *
     * @tparam L The law; on a length.
     * @param element The element.
     * @return The new element's id, valid for this builder and the models it builds.
     */
    template<TendonLaw L>
    ElementId addElement(const TendonElementDescription<L>& element)
    {
        const ElementId id = nextElement(element.name);
        auto& stored = std::get<std::vector<Stored<TendonElementDescription<L>>>>(tendonElements_);
        stored.push_back({id, element});
        stored.back().description.name = {};
        return id;
    }

    /**
     * Adds a loop closure between two sites (docs/MULTIBODY.md section 6).
     *
     * @param closure The closure.
     * @return The new constraint's id, valid for this builder and the models it builds.
     */
    ConstraintId addConstraint(const LoopClosureDescription& closure)
    {
        return addConstraintDescription(closure.name, closure);
    }

    /**
     * Adds a strict unilateral limit on a joint coordinate.
     *
     * @param limit The limit; at least one bound.
     * @return The new constraint's id, valid for this builder and the models it builds.
     */
    ConstraintId addConstraint(const CoordinateLimitDescription& limit)
    {
        return addConstraintDescription(limit.name, limit);
    }

    /**
     * Adds a strict unilateral limit on the distance between two sites.
     *
     * @param limit The limit; at least one bound.
     * @return The new constraint's id, valid for this builder and the models it builds.
     */
    ConstraintId addConstraint(const DistanceLimitDescription& limit)
    {
        return addConstraintDescription(limit.name, limit);
    }

    /**
     * Sets the constraint solver's levers (targets and guards); validated by `build`.
     *
     * @param settings The settings.
     */
    void setConstraintSettings(const ConstraintSettings& settings) noexcept
    {
        constraintSettings_ = settings;
    }

    /**
     * Validates the description and compiles it.
     *
     * @return The model, or the first fault found (FaultKind::InvalidInput) naming what is wrong.
     */
    [[nodiscard]] std::expected<Model, Fault> build() const;

private:
    struct StoredBody
    {
        std::string name;
        BodyDescription description;
    };

    struct StoredSite
    {
        std::string name;
        SiteDescription description;
    };

    template<typename Description>
    struct Stored
    {
        ElementId id;
        Description description;
    };

    template<typename L>
    using StoredCoordinateElement = Stored<CoordinateElementDescription<L>>;

    template<typename L>
    using StoredTendonElement = Stored<TendonElementDescription<L>>;

    using ConstraintDescription =
        std::variant<LoopClosureDescription, CoordinateLimitDescription, DistanceLimitDescription>;

    [[nodiscard]] std::expected<void, Fault> compileElements(Model& model) const;
    [[nodiscard]] std::expected<void, Fault> compileConstraints(Model& model) const;

    template<typename Description>
    ConstraintId addConstraintDescription(std::string_view name, Description description)
    {
        const ConstraintId id{static_cast<std::uint32_t>(constraints_.size())};
        constraintNames_.emplace_back(name);
        description.name = {}; // the builder's name list owns the name
        constraints_.emplace_back(description);
        return id;
    }

    ElementId nextElement(std::string_view name)
    {
        const ElementId id{static_cast<std::uint32_t>(elementNames_.size())};
        elementNames_.emplace_back(name);
        return id;
    }

    std::vector<StoredBody> bodies_;
    std::vector<StoredSite> sites_;
    ElementListsOf<StoredCoordinateElement, CoordinateLaws> coordinateElements_;
    ElementListsOf<StoredTendonElement, TendonLaws> tendonElements_;
    std::vector<std::string> elementNames_;
    std::vector<ConstraintDescription> constraints_;
    std::vector<std::string> constraintNames_;
    ConstraintSettings constraintSettings_{};
    sub0math::Vec3<sub0math::f64::Acceleration> gravity_{
        sub0math::f64::Acceleration{}, sub0math::f64::Acceleration{},
        sub0math::f64::Acceleration::fromSi(-StandardGravity)};
};

} // namespace sub0physics::multibody
