#include <atomic>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <sub0math/linear/mat.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/quat.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/vec.hpp>

#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/model_builder.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

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

// Eigenvalues of a symmetric 3x3 matrix in ascending order (the trigonometric solution of the
// characteristic cubic). For an inertia tensor these are the principal moments.
std::array<double, 3> principalMoments(const sub0math::Mat3d& m)
{
    const double p1 = m(0, 1) * m(0, 1) + m(0, 2) * m(0, 2) + m(1, 2) * m(1, 2);
    if (p1 == 0.0)
    {
        std::array<double, 3> diagonal{m(0, 0), m(1, 1), m(2, 2)};
        std::ranges::sort(diagonal);
        return diagonal;
    }
    const double q = (m(0, 0) + m(1, 1) + m(2, 2)) / 3.0;
    const double d0 = m(0, 0) - q;
    const double d1 = m(1, 1) - q;
    const double d2 = m(2, 2) - q;
    const double p = std::sqrt((d0 * d0 + d1 * d1 + d2 * d2 + 2.0 * p1) / 6.0);
    // det((m - q I) / p) / 2, expanded along the first row.
    const double det = d0 * (d1 * d2 - m(1, 2) * m(1, 2)) - m(0, 1) * (m(0, 1) * d2 - m(1, 2) * m(0, 2)) +
                       m(0, 2) * (m(0, 1) * m(1, 2) - d1 * m(0, 2));
    const double r = std::clamp(det / (2.0 * p * p * p), -1.0, 1.0);
    const double phi = std::acos(r) / 3.0;
    const double largest = q + 2.0 * p * std::cos(phi);
    const double smallest = q + 2.0 * p * std::cos(phi + 2.0 * std::numbers::pi / 3.0);
    return {smallest, 3.0 * q - largest - smallest, largest};
}

// A physical rigid body's principal moments obey I1 <= I2 + I3 for every permutation (the mass
// distribution is non-negative); equality is a planar lamina, so the test is relative.
bool satisfiesTriangleInequality(const sub0math::Mat3d& inertia)
{
    constexpr double RelativeTolerance = 1e-9;
    const std::array<double, 3> i = principalMoments(inertia);
    return i[2] <= (i[0] + i[1]) * (1.0 + RelativeTolerance);
}

// A massless body: mass and every inertia entry exactly zero.
bool isMassless(const BodyDescription& body)
{
    const auto& m = body.inertia.moments;
    const auto& p = body.inertia.products;
    return body.mass.si() == 0.0 && m.x.si() == 0.0 && m.y.si() == 0.0 && m.z.si() == 0.0 && p.x.si() == 0.0 &&
           p.y.si() == 0.0 && p.z.si() == 0.0;
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
    const auto com = body.centreOfMass - sub0math::Point3<sub0math::f64::Length>::origin();
    if (!isFinite(com) || !isFinite(body.inertia.moments) || !isFinite(body.inertia.products))
    {
        return Fault{FaultKind::InvalidInput, "body centre of mass and inertia must be finite"};
    }
    if (isMassless(body))
    {
        // A body that only carries a joint; build checks that the mass matrix stays positive definite.
    }
    else if (!std::isfinite(body.mass.si()) || !(body.mass.si() > 0.0))
    {
        return Fault{FaultKind::InvalidInput, "body mass must be finite and positive (or zero with zero inertia)"};
    }
    else if (!isPositiveDefinite(inertiaMatrix(body.inertia)))
    {
        return Fault{FaultKind::InvalidInput, "body inertia must be positive definite"};
    }
    else if (!satisfiesTriangleInequality(inertiaMatrix(body.inertia)))
    {
        return Fault{FaultKind::InvalidInput, "body inertia violates the triangle inequality (a principal moment exceeds the sum of the other two)"};
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

// Bodies whose joints move a point on `body`: itself and its ancestors.
std::vector<std::size_t> chainOf(const std::vector<std::size_t>& parents, std::size_t body)
{
    std::vector<std::size_t> chain;
    for (std::size_t j = body; j != NoParent; j = parents[j])
    {
        chain.push_back(j);
    }
    return chain;
}

template<typename Lists, typename F, std::size_t... I>
void forEachListImpl(Lists& lists, F& f, std::index_sequence<I...>)
{
    (f(std::get<I>(lists), std::integral_constant<std::size_t, I>{}), ...);
}

// Calls f(list, index) for each list of a tuple of lists, in order.
template<typename Lists, typename F>
void forEachList(Lists& lists, F&& f)
{
    forEachListImpl(lists, f, std::make_index_sequence<std::tuple_size_v<std::remove_cvref_t<Lists>>>{});
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

    // Ids only tell models apart (a Workspace checks it); they never influence a result.
    static std::atomic<std::uint64_t> nextId{1};
    Model model;
    model.id_ = nextId.fetch_add(1, std::memory_order_relaxed);
    const std::size_t n = bodies_.size();
    model.parent_.reserve(n);
    model.coordinate_.reserve(n);
    model.jointType_.reserve(n);
    model.axis_.reserve(n);
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
        model.hasMasslessBodies_ = model.hasMasslessBodies_ || isMassless(body);
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
        if (site.body != WorldBody && site.body.value() >= n)
        {
            return invalidInput("site body is not a body of this model, nor the world");
        }
        if (!isFinite(site.frame) || !isUnit(site.frame.rotation))
        {
            return invalidInput("site frame must be finite with a unit rotation");
        }
        auto frame = site.frame;
        frame.rotation = sub0math::normalized(frame.rotation);
        model.siteBody_.push_back(site.body == WorldBody ? NoParent : site.body.value());
        model.siteFrame_.push_back(frame);
        model.siteName_.push_back(stored.name);
    }

    if (auto compiled = compileElements(model); !compiled)
    {
        return std::unexpected(compiled.error());
    }
    if (auto compiled = compileConstraints(model); !compiled)
    {
        return std::unexpected(compiled.error());
    }

    model.gravity_ = gravity_;
    model.baseAcceleration_ = {{}, {-gravity_.x.si(), -gravity_.y.si(), -gravity_.z.si()}};
    if (!rings_.empty())
    {
        if (formulation_ == Formulation::General)
        {
            return invalidInput("a model with rings needs the planar formulation (planar first, docs/WHEEL.md)");
        }
        if (!constraints_.empty())
        {
            return invalidInput("constraints are not supported with ring substructures yet");
        }
    }
    if (auto compiled = compileFormulation(model); !compiled)
    {
        return std::unexpected(compiled.error());
    }
    if (auto compiled = compileRings(model); !compiled)
    {
        return std::unexpected(compiled.error());
    }
    if (model.hasMasslessBodies_)
    {
        // A massless body is valid only when every joint still moves some mass: M positive definite. Checked at the
        // zero configuration (a configuration-dependent singularity, such as a massless revolute body whose
        // children's mass lies on its axis there, faults the first step that meets it instead).
        Workspace workspace(model);
        const State zero = State::atRest(model);
        if (auto mass = computeMassMatrix(model, zero, workspace); !mass)
        {
            return std::unexpected(mass.error());
        }
        if (!factorMassMatrix(model, workspace))
        {
            return invalidInput("a massless body leaves a joint that moves no mass: the mass matrix is singular");
        }
    }
    return model;
}

std::expected<void, Fault> ModelBuilder::compileFormulation(Model& model) const
{
    model.formulation_ = Formulation::General;
    if (formulation_ == Formulation::General)
    {
        return {};
    }
    // The plane's normal k: the first coordinate axis (z, then y, then x) for which every test holds exactly.
    const std::size_t n = model.parent_.size();
    const std::array<double, 3> gravity{model.gravity_.x.si(), model.gravity_.y.si(), model.gravity_.z.si()};
    for (const std::size_t k : {std::size_t{2}, std::size_t{1}, std::size_t{0}})
    {
        const std::size_t i = (k + 1) % 3;
        const std::size_t j = (k + 2) % 3;
        if (gravity[k] != 0.0)
        {
            continue;
        }
        PlanarModel planar{.normal = k, .axisI = i, .axisJ = j, .gravityX = gravity[i], .gravityY = gravity[j],
                           .bodies = {}, .siteX = {}, .siteY = {}, .siteNormal = {}};
        planar.bodies.reserve(n);
        std::vector<double> normalOffset(n, 0.0); // each body origin's world k coordinate
        bool qualifies = true;
        for (std::size_t b = 0; b < n; ++b)
        {
            const sub0math::PluckerTransformd& tree = model.treeTransform_[b];
            const sub0math::Mat3d& e = tree.rotation;
            const sub0math::Vec3d& axis = model.axis_[b];
            const bool revolute = model.jointType_[b] == JointType::Revolute;
            if (e(k, i) != 0.0 || e(k, j) != 0.0 || e(i, k) != 0.0 || e(j, k) != 0.0 || e(k, k) != 1.0 ||
                (revolute ? (axis[i] != 0.0 || axis[j] != 0.0) : axis[k] != 0.0))
            {
                qualifies = false;
                break;
            }
            const sub0math::SpatialInertiad& inertia = model.inertia_[b];
            // E = R^T for a frame rotated by theta about k: E(i, i) = cos, E(i, j) = sin.
            planar.bodies.push_back({.placement = {e(i, i), e(i, j), tree.translation[i], tree.translation[j]},
                                     .unrotated = model.placementUnrotated_[b] != 0U,
                                     .revolute = revolute,
                                     .axisSign = revolute ? axis[k] : 0.0,
                                     .axisX = revolute ? 0.0 : axis[i],
                                     .axisY = revolute ? 0.0 : axis[j],
                                     .inertia = {inertia.mass, inertia.firstMoment[i], inertia.firstMoment[j],
                                                 inertia.rotational(k, k)}});
            const std::size_t p = model.parent_[b];
            normalOffset[b] = (p == NoParent ? 0.0 : normalOffset[p]) + tree.translation[k];
        }
        if (!qualifies)
        {
            continue;
        }
        const std::size_t sites = model.siteBody_.size();
        planar.siteX.reserve(sites);
        planar.siteY.reserve(sites);
        planar.siteNormal.reserve(sites);
        for (std::size_t s = 0; s < sites; ++s)
        {
            const auto& t = model.siteFrame_[s].translation;
            const std::array<double, 3> offset{t.x.si(), t.y.si(), t.z.si()};
            const std::size_t body = model.siteBody_[s];
            planar.siteX.push_back(offset[i]);
            planar.siteY.push_back(offset[j]);
            planar.siteNormal.push_back((body == NoParent ? 0.0 : normalOffset[body]) + offset[k]);
        }
        model.formulation_ = Formulation::Planar;
        model.planar_ = std::move(planar);
        return {};
    }
    if (formulation_ == Formulation::Planar)
    {
        return invalidInput("the planar formulation needs every revolute axis along one coordinate axis, prismatic "
                            "axes and placement rotations in its plane, and no gravity along it");
    }
    return {};
}

std::expected<void, Fault> ModelBuilder::compileRings(Model& model) const
{
    if (rings_.empty())
    {
        return {};
    }
    if (model.formulation_ != Formulation::Planar)
    {
        return invalidInput("a model with rings needs the planar formulation, and this one does not qualify");
    }
    const PlanarModel& pm = model.planar_;
    const std::size_t i = pm.axisI;
    const std::size_t j = pm.axisJ;
    const std::size_t k = pm.normal;
    std::size_t coordinate = model.parent_.size();
    model.rings_.reserve(rings_.size());
    for (std::size_t r = 0; r < rings_.size(); ++r)
    {
        const DeformableRingDescription& d = rings_[r].description;
        const std::string& name = rings_[r].name;
        if (name.empty())
        {
            return invalidInput("ring name is empty");
        }
        for (std::size_t s = 0; s < r; ++s)
        {
            if (rings_[s].name == name)
            {
                return invalidInput("ring name is not unique");
            }
        }
        if (d.host == WorldBody || d.host.value() >= model.parent_.size())
        {
            return invalidInput("ring host is not a body of this model");
        }
        const std::size_t n = d.restPositions.size();
        if (n < 3U || d.nodeMass.size() != n || d.nodeInertia.size() != n || d.foundation.size() != n ||
            d.edgeLaws.size() != n || d.edgeRestLength.size() != n || (!d.nodeIds.empty() && d.nodeIds.size() != n) ||
            (!d.edgeIds.empty() && d.edgeIds.size() != n))
        {
            return invalidInput("ring needs 3 or more nodes and one mass, inertia, foundation, edge law, rest length "
                                "(and id, when given) per node or edge");
        }
        const sub0math::Vec3d& a0 = d.displacementAxes[0];
        const sub0math::Vec3d& a1 = d.displacementAxes[1];
        if (!isFinite(a0) || !isFinite(a1) || std::abs(dot(a0, a0) - 1.0) > UnitTolerance ||
            std::abs(dot(a1, a1) - 1.0) > UnitTolerance || std::abs(dot(a0, a1)) > UnitTolerance || a0[k] != 0.0 ||
            a1[k] != 0.0)
        {
            return invalidInput("ring displacement axes must be orthonormal and in the model's plane");
        }
        const double skin = d.skinRadius.si();
        if (!std::isfinite(skin) || !(skin >= 0.0))
        {
            return invalidInput("ring skin radius must be finite and not negative");
        }
        DeformableRing ring;
        ring.id = RingId{static_cast<std::uint32_t>(r)};
        ring.name = name;
        ring.host = d.host.value();
        ring.firstCoordinate = coordinate;
        ring.nodeCount = n;
        for (std::size_t b = ring.host; b != NoParent; b = model.parent_[b])
        {
            ring.hostChain.push_back(b);
        }
        // Exact unit axes, as the builder renormalises joint axes.
        const sub0math::Vec3d u0 = sub0math::normalized(a0);
        const sub0math::Vec3d u1 = sub0math::normalized(a1);
        ring.axis0 = {u0[i], u0[j]};
        ring.axis1 = {u1[i], u1[j]};
        ring.skinRadius = skin;
        ring.topologyRevision = d.topologyRevision;
        ring.foundation = d.foundation;
        ring.edgeLaws = d.edgeLaws;
        for (std::size_t e = 0; e < n; ++e)
        {
            const auto& p = d.restPositions[e];
            const std::array<double, 3> rest{p.x.si(), p.y.si(), p.z.si()};
            const double mass = d.nodeMass[e].si();
            const double inertia = d.nodeInertia[e].si();
            const double length = d.edgeRestLength[e].si();
            if (!std::isfinite(rest[0]) || !std::isfinite(rest[1]) || !std::isfinite(rest[2]) ||
                !std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(inertia) || !(inertia >= 0.0) ||
                !std::isfinite(length) || !(length > 0.0))
            {
                return invalidInput("ring nodes need finite rest positions, positive masses, non-negative inertias "
                                    "and positive edge rest lengths");
            }
            if (e == 0U)
            {
                ring.normalOffset = rest[k];
            }
            else if (rest[k] != ring.normalOffset)
            {
                return invalidInput("ring nodes must share one out-of-plane offset");
            }
            ring.restX.push_back(rest[i]);
            ring.restY.push_back(rest[j]);
            ring.rest0.push_back(ring.axis0[0] * rest[i] + ring.axis0[1] * rest[j]);
            ring.rest1.push_back(ring.axis1[0] * rest[i] + ring.axis1[1] * rest[j]);
            ring.mass.push_back(mass);
            ring.inertia.push_back(inertia);
            ring.edgeRestLength.push_back(length);
            ring.nodeIds.push_back(d.nodeIds.empty() ? static_cast<std::uint32_t>(e) : d.nodeIds[e]);
            ring.edgeIds.push_back(d.edgeIds.empty() ? static_cast<std::uint32_t>(e) : d.edgeIds[e]);
        }
        for (const std::vector<std::uint32_t>* ids : {&ring.nodeIds, &ring.edgeIds})
        {
            std::vector<std::uint32_t> sorted = *ids;
            std::ranges::sort(sorted);
            if (std::ranges::adjacent_find(sorted) != sorted.end())
            {
                return invalidInput("ring material ids must be unique");
            }
        }
        coordinate += ring.coordinateCount();
        model.rings_.push_back(std::move(ring));
    }
    model.ringCoordinateCount_ = coordinate - model.parent_.size();
    return {};
}

std::expected<void, Fault> ModelBuilder::compileElements(Model& model) const
{
    for (std::size_t i = 0; i < elementNames_.size(); ++i)
    {
        if (elementNames_[i].empty())
        {
            return invalidInput("element name is empty");
        }
        for (std::size_t j = 0; j < i; ++j)
        {
            if (elementNames_[j] == elementNames_[i])
            {
                return invalidInput("element name is not unique");
            }
        }
    }
    model.elementName_ = elementNames_;

    const std::size_t n = model.parent_.size();
    std::optional<Fault> fault;
    forEachList(coordinateElements_, [&](const auto& stored, auto index) {
        auto& compiled = std::get<decltype(index)::value>(model.coordinateElements_);
        using Element = typename std::remove_cvref_t<decltype(compiled)>::value_type;
        using L = decltype(Element::law);
        compiled.reserve(stored.size());
        for (const auto& [id, element] : stored)
        {
            if (fault)
            {
                return;
            }
            if (element.body == WorldBody || element.body.value() >= n)
            {
                fault = Fault{FaultKind::InvalidInput, "coordinate element body is not a body of this model"};
                return;
            }
            const std::size_t body = element.body.value();
            constexpr bool linear = std::same_as<LawCoordinate<L>, sub0math::f64::Length>;
            if ((model.jointType_[body] == JointType::Prismatic) != linear)
            {
                fault = Fault{FaultKind::InvalidInput, "coordinate element law does not suit the joint: a law on a "
                                                       "length needs a prismatic joint, a law on an angle a revolute one"};
                return;
            }
            const double reference = element.reference.value_or(model.springReference_[body]);
            if (!std::isfinite(reference))
            {
                fault = Fault{FaultKind::InvalidInput, "coordinate element reference must be finite"};
                return;
            }
            compiled.push_back(Element{.id = id,
                                       .law = element.law,
                                       .coordinate = model.coordinate_[body],
                                       .sign = element.sense == LawSense::Increasing ? 1.0 : -1.0,
                                       .reference = reference});
        }
    });
    if (fault)
    {
        return std::unexpected(*fault);
    }

    std::size_t tendons = 0;
    forEachList(tendonElements_, [&](const auto& stored, auto index) {
        auto& compiled = std::get<decltype(index)::value>(model.tendonElements_);
        using Element = typename std::remove_cvref_t<decltype(compiled)>::value_type;
        compiled.reserve(stored.size());
        for (const auto& [id, element] : stored)
        {
            if (fault)
            {
                return;
            }
            const std::size_t sites = model.siteBody_.size();
            if (element.first.value() >= sites || element.second.value() >= sites)
            {
                fault = Fault{FaultKind::InvalidInput, "tendon site is not a site of this model"};
                return;
            }
            const std::size_t firstBody = model.siteBody_[element.first.value()];
            const std::size_t secondBody = model.siteBody_[element.second.value()];
            if (firstBody == secondBody)
            {
                fault = Fault{FaultKind::InvalidInput, "tendon sites must be on different bodies (or one on the world)"};
                return;
            }
            const double reference = element.reference.si();
            if (!std::isfinite(reference) || reference < 0.0)
            {
                fault = Fault{FaultKind::InvalidInput, "tendon reference length must be finite and not negative"};
                return;
            }
            const std::vector<std::size_t> first =
                firstBody == NoParent ? std::vector<std::size_t>{} : chainOf(model.parent_, firstBody);
            const std::vector<std::size_t> second =
                secondBody == NoParent ? std::vector<std::size_t>{} : chainOf(model.parent_, secondBody);
            Element compiledElement{.id = id,
                                    .law = element.law,
                                    .first = element.first,
                                    .second = element.second,
                                    .sign = element.sense == LawSense::Increasing ? 1.0 : -1.0,
                                    .reference = reference,
                                    .wrenchSlot = 2 * tendons,
                                    .support = {},
                                    .supportCoordinates = {},
                                    .onFirst = {},
                                    .onSecond = {}};
            for (std::size_t body = 0; body < n; ++body)
            {
                const bool onFirst = std::ranges::find(first, body) != first.end();
                const bool onSecond = std::ranges::find(second, body) != second.end();
                // A prismatic joint that carries both ends translates them together: it changes neither the length
                // nor any other entry's Jacobian, so it is left out.
                const bool inert = onFirst && onSecond && model.jointType_[body] == JointType::Prismatic;
                if ((onFirst || onSecond) && !inert)
                {
                    compiledElement.support.push_back(body);
                    compiledElement.supportCoordinates.push_back(model.coordinate_[body]);
                    compiledElement.onFirst.push_back(onFirst ? 1U : 0U);
                    compiledElement.onSecond.push_back(onSecond ? 1U : 0U);
                }
            }
            const bool oneChain = std::ranges::all_of(compiledElement.onFirst, [](std::uint8_t on) { return on != 0U; }) ||
                                  std::ranges::all_of(compiledElement.onSecond, [](std::uint8_t on) { return on != 0U; });
            model.elementsFollowTree_ = model.elementsFollowTree_ && oneChain;
            compiled.push_back(std::move(compiledElement));
            ++tendons;
        }
    });
    if (fault)
    {
        return std::unexpected(*fault);
    }
    model.tendonCount_ = tendons;
    return {};
}

std::expected<void, Fault> ModelBuilder::compileConstraints(Model& model) const
{
    const ConstraintSettings& settings = constraintSettings_;
    const auto positiveFinite = [](double v) { return std::isfinite(v) && v > 0.0; };
    if (!positiveFinite(settings.closureTolerance) || !positiveFinite(settings.eventTolerance) ||
        !positiveFinite(settings.solverTolerance) || !(settings.solverTolerance < 1e-3) ||
        settings.maxEventIterations == 0U)
    {
        return invalidInput("constraint settings need positive finite tolerances (solver tolerance below 1e-3) and "
                            "at least one event iteration");
    }
    model.constraintSettings_ = settings;
    for (std::size_t i = 0; i < constraintNames_.size(); ++i)
    {
        if (constraintNames_[i].empty())
        {
            return invalidInput("constraint name is empty");
        }
        for (std::size_t j = 0; j < i; ++j)
        {
            if (constraintNames_[j] == constraintNames_[i])
            {
                return invalidInput("constraint name is not unique");
            }
        }
    }
    model.constraintName_ = constraintNames_;

    const std::size_t sites = model.siteBody_.size();
    // Two sites that the constraint can act between: resolved, and on different bodies (one may be the world).
    const auto checkSites = [&](SiteId first, SiteId second) -> std::optional<Fault> {
        if (first.value() >= sites || second.value() >= sites)
        {
            return Fault{FaultKind::InvalidInput, "constraint site is not a site of this model"};
        }
        if (model.siteBody_[first.value()] == model.siteBody_[second.value()])
        {
            return Fault{FaultKind::InvalidInput, "constraint sites must be on different bodies (or one on the world)"};
        }
        return std::nullopt;
    };
    const auto validRestitution = [](double e) { return std::isfinite(e) && e >= 0.0 && e <= 1.0; };

    std::size_t rows = 0;
    std::size_t wrenches = 0;
    for (std::size_t i = 0; i < constraints_.size(); ++i)
    {
        const auto* closure = std::get_if<LoopClosureDescription>(&constraints_[i]);
        if (closure == nullptr)
        {
            continue;
        }
        if (const auto fault = checkSites(closure->first, closure->second))
        {
            return std::unexpected(*fault);
        }
        if (closure->kind != ClosureKind::Point && closure->kind != ClosureKind::Frame)
        {
            return invalidInput("closure kind is not supported");
        }
        model.closures_.push_back({.id = ConstraintId{static_cast<std::uint32_t>(i)},
                                   .first = closure->first,
                                   .second = closure->second,
                                   .kind = closure->kind,
                                   .firstRow = rows,
                                   .wrenchSlot = wrenches});
        rows += model.closures_.back().rows();
        wrenches += 2;
    }
    model.closureRowCount_ = rows;

    for (std::size_t i = 0; i < constraints_.size(); ++i)
    {
        const ConstraintId id{static_cast<std::uint32_t>(i)};
        if (const auto* limit = std::get_if<CoordinateLimitDescription>(&constraints_[i]))
        {
            const std::size_t n = model.parent_.size();
            if (limit->body == WorldBody || limit->body.value() >= n)
            {
                return invalidInput("coordinate limit body is not a body of this model");
            }
            if (!limit->lower && !limit->upper)
            {
                return invalidInput("limit needs at least one bound");
            }
            if ((limit->lower && !std::isfinite(*limit->lower)) || (limit->upper && !std::isfinite(*limit->upper)) ||
                (limit->lower && limit->upper && !(*limit->lower < *limit->upper)))
            {
                return invalidInput("limit bounds must be finite, with lower below upper");
            }
            if (!validRestitution(limit->restitution))
            {
                return invalidInput("limit restitution must be in [0, 1]");
            }
            const std::size_t coordinate = model.coordinate_[limit->body.value()];
            const auto add = [&](double bound, double sign) {
                model.limits_.push_back({.id = id,
                                         .kind = LimitKind::Coordinate,
                                         .coordinate = coordinate,
                                         .first = {},
                                         .second = {},
                                         .bound = bound,
                                         .sign = sign,
                                         .restitution = limit->restitution,
                                         .row = rows + model.limits_.size(),
                                         .wrenchSlot = 0});
            };
            if (limit->lower)
            {
                add(*limit->lower, 1.0);
            }
            if (limit->upper)
            {
                add(*limit->upper, -1.0);
            }
        }
        else if (const auto* distance = std::get_if<DistanceLimitDescription>(&constraints_[i]))
        {
            if (const auto fault = checkSites(distance->first, distance->second))
            {
                return std::unexpected(*fault);
            }
            if (!distance->minimum && !distance->maximum)
            {
                return invalidInput("limit needs at least one bound");
            }
            const double minimum = distance->minimum ? distance->minimum->si() : 0.0;
            const double maximum = distance->maximum ? distance->maximum->si() : 0.0;
            if ((distance->minimum && (!std::isfinite(minimum) || !(minimum > 0.0))) ||
                (distance->maximum && (!std::isfinite(maximum) || !(maximum > 0.0))) ||
                (distance->minimum && distance->maximum && !(minimum < maximum)))
            {
                return invalidInput("distance limit bounds must be finite and positive, with minimum below maximum");
            }
            if (!validRestitution(distance->restitution))
            {
                return invalidInput("limit restitution must be in [0, 1]");
            }
            const auto add = [&](double bound, double sign) {
                model.limits_.push_back({.id = id,
                                         .kind = LimitKind::Distance,
                                         .coordinate = 0,
                                         .first = distance->first,
                                         .second = distance->second,
                                         .bound = bound,
                                         .sign = sign,
                                         .restitution = distance->restitution,
                                         .row = rows + model.limits_.size(),
                                         .wrenchSlot = wrenches});
                wrenches += 2;
            };
            if (distance->minimum)
            {
                add(minimum, 1.0);
            }
            if (distance->maximum)
            {
                add(maximum, -1.0);
            }
        }
    }
    model.constraintWrenchCount_ = wrenches;
    return {};
}

std::optional<ConstraintId> Model::findConstraint(std::string_view name) const noexcept
{
    for (std::size_t i = 0; i < constraintName_.size(); ++i)
    {
        if (constraintName_[i] == name)
        {
            return ConstraintId{static_cast<std::uint32_t>(i)};
        }
    }
    return std::nullopt;
}

std::optional<ElementId> Model::findElement(std::string_view name) const noexcept
{
    for (std::size_t i = 0; i < elementName_.size(); ++i)
    {
        if (elementName_[i] == name)
        {
            return ElementId{static_cast<std::uint32_t>(i)};
        }
    }
    return std::nullopt;
}

std::optional<RingId> Model::findRing(std::string_view name) const noexcept
{
    for (const DeformableRing& ring : rings_)
    {
        if (ring.name == name)
        {
            return ring.id;
        }
    }
    return std::nullopt;
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
