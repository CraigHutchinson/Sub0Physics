#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/id.hpp"
#include "sub0physics/element/axial_spring_damper.hpp"
#include "sub0physics/element/law_response.hpp"
#include "sub0physics/element/linear_spring_damper.hpp"
#include "sub0physics/element/travel_buffer.hpp"
#include "sub0physics/element/unilateral_spring_damper.hpp"
#include "sub0physics/multibody/description.hpp"

// Force elements (docs/MULTIBODY.md section 4 and "M1b contracts"): an element binds a law (element/) to the
// model, either to one joint coordinate or to the distance between two sites. The law's formula exists once,
// in the law; an element only maps the model's coordinates to the law's coordinate and the law's force and
// Jacobians back. Elements are described here, compiled by ModelBuilder and evaluated by evaluateElements
// (multibody/dynamics.hpp) into the Workspace, without allocating.

namespace sub0physics::multibody
{

class Model;
struct State;
struct Workspace;

/** Identifies a force element within one ModelBuilder and the Model it builds; numbered as added. */
using ElementId = Id<struct ElementIdTag, std::uint32_t>;

/** How a law's coordinate follows the model quantity it is bound to (a joint coordinate or a length). */
enum class LawSense : std::uint8_t
{
    Increasing, ///< The law's coordinate is `x - reference`: for a tendon, its extension.
    Decreasing  ///< The law's coordinate is `reference - x`: for a tendon, its compression.
};

/**
 * A compile-time list of law types.
 *
 * @tparam Laws The laws.
 */
template<typename... Laws>
struct LawList
{
};

/**
 * The laws an element can bind on a coordinate, a closed list so evaluation dispatches at compile time (no
 * virtual calls). A new law joins by being added here (and to TendonLaws if it is on a length); the explicit
 * instantiations in `multibody_elements.cpp` follow these lists.
 */
using CoordinateLaws = LawList<LinearSpringDamper<sub0math::f64::Length>, LinearSpringDamper<sub0math::f64::Angle>,
                               AxialSpringDamper, UnilateralSpringDamper, TravelBuffer>;

/** The laws a tendon can bind: those on a length. */
using TendonLaws = LawList<LinearSpringDamper<sub0math::f64::Length>, AxialSpringDamper, UnilateralSpringDamper,
                           TravelBuffer>;

namespace detail
{

template<typename T, typename List>
inline constexpr bool InLawList = false;

template<typename T, typename... Laws>
inline constexpr bool InLawList<T, LawList<Laws...>> = (std::same_as<T, Laws> || ...);

template<template<typename> class Element, typename List>
struct ListsOf;

template<template<typename> class Element, typename... Laws>
struct ListsOf<Element, LawList<Laws...>>
{
    using type = std::tuple<std::vector<Element<Laws>>...>;
};

template<typename L>
consteval auto lawCoordinateOf()
{
    if constexpr (Law<L, sub0math::f64::Length>)
    {
        return sub0math::f64::Length{};
    }
    else
    {
        return sub0math::f64::Angle{};
    }
}

} // namespace detail

/**
 * A law that an element can bind: one of CoordinateLaws.
 *
 * @tparam L A law type.
 */
template<typename L>
concept ElementLaw = detail::InLawList<L, CoordinateLaws>;

/**
 * The laws a tendon can bind: those on a linear coordinate.
 *
 * @tparam L A law type.
 */
template<typename L>
concept TendonLaw = detail::InLawList<L, TendonLaws>;

/**
 * The coordinate quantity a law works on: `f64::Length` or `f64::Angle`.
 *
 * @tparam L An element law.
 */
template<ElementLaw L>
using LawCoordinate = decltype(detail::lawCoordinateOf<L>());

/**
 * A law on one joint coordinate: a coordinate spring-damper, a travel buffer on a slider, a torsion spring.
 *
 * The law's coordinate is `q - reference` (or `reference - q` for LawSense::Decreasing). A law on a length
 * binds a prismatic joint and a law on an angle a revolute one; ModelBuilder refuses any other pairing.
 *
 * @tparam L The law.
 */
template<ElementLaw L>
struct CoordinateElementDescription
{
    std::string_view name;                ///< Unique among elements, non-empty; copied by the builder.
    BodyId body;                          ///< The body whose joint coordinate the law acts on.
    L law;                                ///< The law, built by its validating `create`.
    LawSense sense = LawSense::Increasing; ///< Direction of the law's coordinate.
    std::optional<double> reference{};    ///< Joint coordinate (rad or m) where the law's coordinate is zero; the joint's `springReference` when empty.
};

/**
 * A law on the straight-line distance between two sites: a strut, a cable or a shock between two bodies.
 *
 * With `l` the distance, the law's coordinate is `l - reference` (extension) or `reference - l`
 * (compression, LawSense::Decreasing). The force acts along the line between the sites, and a site may be
 * fixed to the world (a SiteDescription on WorldBody).
 *
 * @tparam L The law; on a length.
 */
template<TendonLaw L>
struct TendonElementDescription
{
    std::string_view name;                ///< Unique among elements, non-empty; copied by the builder.
    SiteId first;                         ///< One end.
    SiteId second;                        ///< The other end; a different site.
    L law;                                ///< The law, built by its validating `create`.
    LawSense sense = LawSense::Increasing; ///< Extension (Increasing) or compression (Decreasing).
    sub0math::f64::Length reference{};    ///< Length at which the law's coordinate is zero; finite, not negative.
};

/**
 * The wrench an element applies to the body carrying one of its sites, for reactions and load reporting.
 * Plain SI in world axes, like the workspace it is written to.
 */
struct ElementSiteWrench
{
    ElementId element;      ///< The element.
    SiteId site;            ///< The site; the wrench acts on the site's body.
    sub0math::Vec3d force;  ///< Force on the body [N], world axes, through the site's origin.
    sub0math::Vec3d moment; ///< Moment about the site's origin [N m], world axes; zero for a tendon.
};

/** What an element reads: the model, the state and the kinematics already in the workspace. */
struct ElementInputs
{
    const Model& model;         ///< The model the element belongs to.
    const State& state;         ///< Coordinates and rates.
    const Workspace& workspace; ///< `worldToBody` holds the world transforms at `state`.
};

/**
 * Where an element writes, without allocating: sums over every element, and per-site wrenches.
 *
 * Matrices are row-major `dofCount x dofCount`. An empty `stiffness` asks for the force, energy and wrenches
 * only. Elements add to the sums; they never clear them.
 */
struct ElementOutputs
{
    std::span<double> force;     ///< Generalized force tau [N or N m], dofCount entries.
    std::span<double> stiffness; ///< dtau/dq, or empty when Jacobians are not wanted.
    std::span<double> damping;   ///< dtau/dqd, or empty when Jacobians are not wanted.
    bool geometricStiffness = true; ///< Whether `stiffness` includes the tendon terms from d2l/dq2 (geometric and rate transport).
    double potential = 0.0;        ///< Stored energy [J].
    double dissipationPower = 0.0; ///< Power removed [W]; never negative.
    std::span<ElementSiteWrench> siteWrenches; ///< One slot per tendon end, in the model's order.
    std::span<sub0math::Vec3d> vectors;        ///< Scratch: 3 dofCount vectors.
    std::span<double> scalars;                 ///< Scratch: dofCount (dofCount + 2) values.
};

/**
 * A force element (docs/MULTIBODY.md section 4): evaluated in place into ElementOutputs, never allocating,
 * and inlined in the evaluation loop because the element set is a compile-time list (no virtual calls).
 * `evaluate` returns false when the state is outside the element's domain (a tendon of zero length).
 *
 * @tparam E The element type.
 */
template<typename E>
concept ForceElement = requires(const E e, const ElementInputs& in, ElementOutputs& out) {
    { e.evaluate(in, out) } noexcept -> std::same_as<bool>;
    { e.hasJacobians() } -> std::same_as<bool>;
};

/**
 * A compiled coordinate element: the law on joint coordinate `coordinate`.
 *
 * tau on the coordinate is `sign * F(s, sd)` with `s = sign (q - reference)`, `sd = sign qd`, so the
 * Jacobian entries are the law's own: `dtau/dq = dF/ds`, `dtau/dqd = dF/dsd` (sign squared is one).
 *
 * @tparam L The law.
 */
template<ElementLaw L>
struct CoordinateElement
{
    ElementId id;           ///< The element's id.
    L law;                  ///< The law.
    std::size_t coordinate; ///< Index of the coordinate in q.
    double sign;            ///< +1 (Increasing) or -1 (Decreasing).
    double reference;       ///< Coordinate where the law's coordinate is zero.

    /**
     * Adds the element's force, Jacobians, energy and power.
     *
     * @param in Model, state and workspace.
     * @param out Sums to add to.
     * @return Always true: every coordinate is in a coordinate element's domain.
     */
    [[nodiscard]] bool evaluate(const ElementInputs& in, ElementOutputs& out) const noexcept;

    /** @return True: the laws supply exact Jacobians. */
    [[nodiscard]] static constexpr bool hasJacobians() noexcept
    {
        return true;
    }
};

/**
 * A compiled tendon element: the law on the distance between two sites.
 *
 * With `d = p1 - p2`, `l = |d|`, `u = d / l` and `J_l = u^T (J_1 - J_2)`, tau is `J_l^T F` where
 * `F = sign * law.force` is the force conjugate to `l`. The stiffness is exact:
 * `dtau/dq = J_l^T (dF/dl) J_l + F d2l/dq2 + J_l^T (dF/dl') (d2l/dq2 qd)^T`: the material term, the
 * geometric term (which keeps a preloaded tendon's implicit step stable) and the transport of the rate
 * `l' = J_l qd` with the configuration (non-symmetric; zero without damping). `dtau/dqd = J_l^T (dF/dl') J_l`.
 * `d2l/dq2` follows from the site Hessians: for joints j, k on a site's chain, with j the ancestor,
 * `d J_k / d q_j = w_j x J_k` for a revolute j (zero for a prismatic one), where `w_j` is its world axis.
 *
 * @tparam L The law; on a length.
 */
template<TendonLaw L>
struct TendonElement
{
    ElementId id;           ///< The element's id.
    L law;                  ///< The law.
    SiteId first;           ///< First end.
    SiteId second;          ///< Second end.
    double sign;            ///< +1 (extension) or -1 (compression).
    double reference;       ///< Length at which the law's coordinate is zero [m].
    std::size_t wrenchSlot; ///< Index of the first end's wrench in ElementOutputs::siteWrenches; the second follows.
    std::vector<std::size_t> support;   ///< Bodies whose joints move either end, ascending (topological order).
    std::vector<std::uint8_t> onFirst;  ///< Per support entry: whether it moves the first end.
    std::vector<std::uint8_t> onSecond; ///< Per support entry: whether it moves the second end.

    /**
     * Adds the element's force, Jacobians, energy, power and site wrenches.
     *
     * @param in Model, state and workspace with world transforms.
     * @param out Sums to add to.
     * @return False when the two sites coincide (the direction is undefined); nothing is written then.
     */
    [[nodiscard]] bool evaluate(const ElementInputs& in, ElementOutputs& out) const noexcept;

    /** @return True: the laws supply exact Jacobians and the geometry is differentiated exactly. */
    [[nodiscard]] static constexpr bool hasJacobians() noexcept
    {
        return true;
    }
};

/**
 * One list per law for an element template, in the law list's order (also the evaluation order).
 *
 * @tparam Element An element template, such as CoordinateElement.
 * @tparam List A LawList.
 */
template<template<typename> class Element, typename List>
using ElementListsOf = typename detail::ListsOf<Element, List>::type;

/** The compiled coordinate elements of a model, by law. */
using CoordinateElementLists = ElementListsOf<CoordinateElement, CoordinateLaws>;

/** The compiled tendon elements of a model, by law. */
using TendonElementLists = ElementListsOf<TendonElement, TendonLaws>;

static_assert(ForceElement<CoordinateElement<TravelBuffer>>);
static_assert(ForceElement<CoordinateElement<LinearSpringDamper<sub0math::f64::Angle>>>);
static_assert(ForceElement<TendonElement<AxialSpringDamper>>);

} // namespace sub0physics::multibody
