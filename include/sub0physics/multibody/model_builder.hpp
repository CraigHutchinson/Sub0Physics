#pragma once

#include <expected>
#include <string>
#include <vector>

#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/constants.hpp"
#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/description.hpp"
#include "sub0physics/multibody/model.hpp"

namespace sub0physics::multibody
{

/**
 * Collects a typed model description and compiles it into an immutable Model.
 *
 * Bodies are numbered in the order they are added, and a joint's parent must already exist, so the
 * numbering is a topological order of the tree. Nothing is checked until `build`, which validates the
 * whole description at once: positive mass, positive definite inertia, unit axes and rotations, finite
 * values, unique non-empty names and resolved references.
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

    std::vector<StoredBody> bodies_;
    std::vector<StoredSite> sites_;
    sub0math::Vec3<sub0math::f64::Acceleration> gravity_{
        sub0math::f64::Acceleration{}, sub0math::f64::Acceleration{},
        sub0math::f64::Acceleration::fromSi(-StandardGravity)};
};

} // namespace sub0physics::multibody
