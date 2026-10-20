#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <sub0math/linear/mat.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/quat.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/vec.hpp>

#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/model_builder.hpp"

namespace sub0physics::multibody
{

namespace
{

// Unit vectors and rotations must be unit to this tolerance; they are then renormalised exactly.
constexpr double UnitTolerance = 1e-9;

template<typename T>
bool isFinite(const sub0math::Vec3<T>& v)
{
    if constexpr (std::is_arithmetic_v<T>)
    {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    }
    else
    {
        return std::isfinite(v.x.si()) && std::isfinite(v.y.si()) && std::isfinite(v.z.si());
    }
}

bool isFinite(const sub0math::RigidTransform<double, sub0math::f64::Length>& t)
{
    const sub0math::Quatd& r = t.rotation;
    return std::isfinite(r.w) && std::isfinite(r.x) && std::isfinite(r.y) && std::isfinite(r.z) &&
           isFinite(t.translation);
}

bool isUnit(const sub0math::Quatd& r)
{
    return std::abs(r.w * r.w + r.x * r.x + r.y * r.y + r.z * r.z - 1.0) <= UnitTolerance;
}

sub0math::Mat3d inertiaMatrix(const RotationalInertia& inertia)
{
    sub0math::Mat3d m{};
    m(0, 0) = inertia.moments.x.si();
    m(1, 1) = inertia.moments.y.si();
    m(2, 2) = inertia.moments.z.si();
    m(0, 1) = m(1, 0) = inertia.products.x.si();
    m(0, 2) = m(2, 0) = inertia.products.y.si();
    m(1, 2) = m(2, 1) = inertia.products.z.si();
    return m;
}

// Sylvester's criterion on the symmetric matrix.
bool isPositiveDefinite(const sub0math::Mat3d& m)
{
    const double minor2 = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
    const double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
                       m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                       m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    return m(0, 0) > 0.0 && minor2 > 0.0 && det > 0.0;
}

std::optional<Fault> checkBody(const BodyDescription& body, std::uint32_t index, const std::vector<std::string>& names)
{
    if (body.name.empty())
    {
        return Fault{FaultKind::InvalidInput, "body name is empty"};
    }
    for (const std::string& other : names)
    {
        if (other == body.name)
        {
            return Fault{FaultKind::InvalidInput, "body name is not unique"};
        }
    }
    if (!std::isfinite(body.mass.si()) || !(body.mass.si() > 0.0))
    {
        return Fault{FaultKind::InvalidInput, "body mass must be finite and positive"};
    }
    const auto com = body.centreOfMass - sub0math::Point3<sub0math::f64::Length>::origin();
    if (!isFinite(com) || !isFinite(body.inertia.moments) || !isFinite(body.inertia.products))
    {
        return Fault{FaultKind::InvalidInput, "body centre of mass and inertia must be finite"};
    }
    if (!isPositiveDefinite(inertiaMatrix(body.inertia)))
    {
        return Fault{FaultKind::InvalidInput, "body inertia must be positive definite"};
    }
    const JointDescription& joint = body.joint;
    if (joint.parent != WorldBody && joint.parent.value() >= index)
    {
        return Fault{FaultKind::InvalidInput, "joint parent must be a body added earlier, or the world"};
    }
    if (joint.type != JointType::Revolute && joint.type != JointType::Prismatic)
    {
        return Fault{FaultKind::InvalidInput, "joint type is not supported"};
    }
    if (!isFinite(joint.placement) || !isUnit(joint.placement.rotation))
    {
        return Fault{FaultKind::InvalidInput, "joint placement must be finite with a unit rotation"};
    }
    if (!isFinite(joint.axis) || std::abs(sub0math::lengthSquared(joint.axis) - 1.0) > UnitTolerance)
    {
        return Fault{FaultKind::InvalidInput, "joint axis must be a finite unit vector"};
    }
    if (!std::isfinite(joint.springReference))
    {
        return Fault{FaultKind::InvalidInput, "joint spring reference must be finite"};
    }
    return std::nullopt;
}

} // namespace

BodyId ModelBuilder::addBody(const BodyDescription& body)
{
    const BodyId id{static_cast<std::uint32_t>(bodies_.size())};
    bodies_.push_back({std::string(body.name), body});
    bodies_.back().description.name = {}; // the stored copy owns the name
    return id;
}

SiteId ModelBuilder::addSite(const SiteDescription& site)
{
    const SiteId id{static_cast<std::uint32_t>(sites_.size())};
    sites_.push_back({std::string(site.name), site});
    sites_.back().description.name = {};
    return id;
}

std::expected<Model, Fault> ModelBuilder::build() const
{
    if (bodies_.empty())
    {
        return invalidInput("model has no bodies");
    }
    if (!isFinite(gravity_))
    {
        return invalidInput("gravity must be finite");
    }

    Model model;
    const std::size_t n = bodies_.size();
    model.parent_.reserve(n);
    model.coordinate_.reserve(n);
    model.jointType_.reserve(n);
    model.axis_.reserve(n);
    model.placement_.reserve(n);
    model.treeTransform_.reserve(n);
    model.placementUnrotated_.reserve(n);
    model.motionSubspace_.reserve(n);
    model.inertia_.reserve(n);
    model.springReference_.reserve(n);
    model.bodyName_.reserve(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        BodyDescription body = bodies_[i].description;
        body.name = bodies_[i].name;
        if (const auto fault = checkBody(body, static_cast<std::uint32_t>(i), model.bodyName_))
        {
            return std::unexpected(*fault);
        }
        const JointDescription& joint = body.joint;
        const sub0math::Vec3d axis = sub0math::normalized(joint.axis);
        auto placement = joint.placement;
        placement.rotation = sub0math::normalized(placement.rotation);

        model.parent_.push_back(joint.parent == WorldBody ? NoParent : joint.parent.value());
        model.coordinate_.push_back(i);
        model.jointType_.push_back(joint.type);
        model.axis_.push_back(axis);
        model.placement_.push_back(placement);
        model.treeTransform_.push_back(sub0math::PluckerTransformd::fromPose(placement));
        model.placementUnrotated_.push_back(model.treeTransform_.back().rotation == sub0math::Mat3d::identity() ? 1U : 0U);
        model.motionSubspace_.push_back(joint.type == JointType::Revolute
                                            ? sub0math::SpatialMotiond{axis, {}}
                                            : sub0math::SpatialMotiond{{}, axis});
        const auto com = body.centreOfMass - sub0math::Point3<sub0math::f64::Length>::origin();
        model.inertia_.push_back(sub0math::SpatialInertiad::fromCentreOfMass(
            body.mass.si(), {com.x.si(), com.y.si(), com.z.si()}, inertiaMatrix(body.inertia)));
        model.springReference_.push_back(joint.springReference);
        model.bodyName_.push_back(bodies_[i].name);
    }

    model.siteBody_.reserve(sites_.size());
    model.siteFrame_.reserve(sites_.size());
    model.siteName_.reserve(sites_.size());
    for (const StoredSite& stored : sites_)
    {
        const SiteDescription& site = stored.description;
        if (stored.name.empty())
        {
            return invalidInput("site name is empty");
        }
        for (const std::string& other : model.siteName_)
        {
            if (other == stored.name)
            {
                return invalidInput("site name is not unique");
            }
        }
        if (site.body.value() >= n)
        {
            return invalidInput("site body is not a body of this model");
        }
        if (!isFinite(site.frame) || !isUnit(site.frame.rotation))
        {
            return invalidInput("site frame must be finite with a unit rotation");
        }
        auto frame = site.frame;
        frame.rotation = sub0math::normalized(frame.rotation);
        model.siteBody_.push_back(site.body.value());
        model.siteFrame_.push_back(frame);
        model.siteName_.push_back(stored.name);
    }

    model.gravity_ = gravity_;
    model.baseAcceleration_ = {{}, {-gravity_.x.si(), -gravity_.y.si(), -gravity_.z.si()}};
    return model;
}

std::optional<BodyId> Model::findBody(std::string_view name) const noexcept
{
    for (std::size_t i = 0; i < bodyName_.size(); ++i)
    {
        if (bodyName_[i] == name)
        {
            return BodyId{static_cast<std::uint32_t>(i)};
        }
    }
    return std::nullopt;
}

std::optional<SiteId> Model::findSite(std::string_view name) const noexcept
{
    for (std::size_t i = 0; i < siteName_.size(); ++i)
    {
        if (siteName_[i] == name)
        {
            return SiteId{static_cast<std::uint32_t>(i)};
        }
    }
    return std::nullopt;
}

} // namespace sub0physics::multibody
