#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include <sub0math/linear/vec.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/soft_constraint.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// Contact on a deformable ring's casing (docs/WHEEL.md sections 3.2 and 10, item 3; step M2c).
//
// A contact acts at a point of one casing edge, along a unit normal: its signed gap g (negative in overlap) moves with
// the generalized coordinates through the host's motion and the two nodes of the edge. Two laws use the same gap
// kinematics:
//
// - **The reference soft-constraint law** (the parity model, soft_constraint.hpp): the contact's row
//   (ringContactJacobian) with the reference engine's regularised law, solved at the start of the step with the mass
//   matrix and held over the step.
// - **The contact element** (the improved model): a scalar law on the gap whose force, stiffness and damping are
//   exact, including the rotation of the normal at a vertex witness and the sliding of the casing's witness along the
//   edge (evaluateRingContactElement), applied through Forcing::ringLoads so the linearly implicit steps treat it
//   implicitly through the ring's block storage.
//
// Planar formulation only (a ring needs it). World vectors are 3D at the interface and must lie in the model's plane.

namespace sub0physics::multibody
{

/**
 * A point of a ring's casing edge where a contact force acts, and its direction.
 *
 * The point is attached to the edge's two nodes with weights `weight0` (node `edge`) and `weight1` (node `edge + 1
 * mod N`): its velocity is the host's motion at the point plus the weighted node rates along the ring's axes. That
 * is the reference engine's attachment of a flex contact (its point moves with each node body, weighted), and with
 * barycentric weights and the point on the edge's axis it is the material point of the edge.
 */
struct RingContactPoint
{
    std::size_t ring = 0;      ///< Index in Model::rings().
    std::size_t edge = 0;      ///< Edge e: nodes e and e + 1 mod N.
    double weight0 = 0.5;      ///< Share of node `edge`; with weight1, non-negative and summing to one.
    double weight1 = 0.5;      ///< Share of node `edge + 1`.
    sub0math::Vec3d point{};   ///< World position of the point [m], in the model's plane.
    sub0math::Vec3d normal{};  ///< Unit world direction of the force on the ring (from the obstacle towards the ring).
};

/**
 * The contact's row: the rate of the gap along the normal per unit generalized velocity, `J = n^T dp/dq` (a dense row
 * of Model::dofCount entries; zero off the host's chain and the edge's four node coordinates). With the normal held
 * fixed this is the gradient of the gap (the envelope theorem: the witnesses' own motion does not enter it).
 *
 * Needs the planar kinematics of a prepared step (prepareLinearlyImplicitStep) or any query that wrote
 * `planar.worldToBody` and the ring's node positions at the state.
 *
 * @param model A planar model with rings.
 * @param workspace The workspace with the kinematics at the state.
 * @param contact The contact; its ring and edge must exist.
 * @param row Receives dofCount entries.
 * @return Nothing, or InvalidInput (not planar, an index, the row's size).
 */
[[nodiscard]] std::expected<void, Fault> ringContactJacobian(const Model& model, const Workspace& workspace,
                                                             const RingContactPoint& contact,
                                                             std::span<double> row) noexcept;

/**
 * The reference engine's average translational inverse inertia of a ring node: the mean over the three world axes of
 * the diagonal of `J M^-1 J^T` for the node's translation, where the axis normal to the plane contributes zero. The
 * reference computes it once at the model's reference configuration (MuJoCo `body_invweight0`,
 * engine_setconst.c:1250-1257) and uses it to scale its contact regularisation.
 *
 * Setup work: it allocates a scratch vector.
 *
 * @param model A planar model with rings.
 * @param workspace A workspace whose mass matrix is factorised (factorMassMatrix) at the configuration wanted.
 * @param ring Ring index.
 * @param out Receives one weight per node of the ring [1/kg].
 * @return Nothing, or InvalidInput.
 */
[[nodiscard]] std::expected<void, Fault> ringNodeInverseWeights(const Model& model, const Workspace& workspace,
                                                                std::size_t ring, std::span<double> out);

/** What a contact element is attached to: the kind of terrain feature that carries its witness. */
enum class ContactWitness : std::uint8_t
{
    Face,  ///< A face of fixed normal: the normal does not turn with the casing.
    Vertex ///< A vertex: the normal is along the casing's nearest axis point minus the vertex, and turns with it.
};

/**
 * The geometry of a contact element at a state: the casing edge's axis point nearest the obstacle (at barycentric
 * parameter `axisT`), and the obstacle's witness.
 */
struct RingContactGeometry
{
    std::size_t ring = 0;      ///< Index in Model::rings().
    std::size_t edge = 0;      ///< Edge e: nodes e and e + 1 mod N.
    double axisT = 0.0;        ///< The casing witness on the edge's axis: 0 at node e, 1 at node e + 1.
    ContactWitness witness = ContactWitness::Face; ///< The obstacle's feature.
    sub0math::Vec3d witnessPoint{}; ///< World position of the obstacle's witness [m]: on the face, or the vertex.
    sub0math::Vec3d normal{};  ///< Unit world normal from the obstacle towards the casing (fixed for a face).
    double skin = 0.0;         ///< Radius of the casing edge [m]: the gap is measured from the axis less this.
};

/**
 * A contact element's contribution, in the element's local coordinates `(host motion w, x, y in the host frame; node e
 * along a_0, a_1; node e + 1 along a_0, a_1)`: the generalized force and its exact Jacobians. The steps map it to the
 * tree's coordinates through the host's chain and into the ring's blocks (Forcing::ringLoads).
 */
struct RingContactLoad
{
    /** Local coordinates: the host's motion (w, x, y) and the edge's two nodes (two each). */
    static constexpr std::size_t Coordinates = 7;

    std::size_t ring = 0;               ///< Index in Model::rings().
    std::size_t edge = 0;               ///< Edge e.
    std::array<double, 7> force{};      ///< Generalized force on the local coordinates [N m, N]: the host-frame wrench (moment, force along i, j), then the nodes'.
    std::array<double, 49> stiffness{}; ///< d force / d local position, row-major 7 x 7, in the chart `p = o + R (v + Rot(theta) l)` (exact; not symmetric).
    std::array<double, 49> damping{};   ///< d force / d local velocity, row-major 7 x 7.
    std::array<double, 7> gradient{};   ///< The gap's gradient in the local coordinates (for the host chain's curvature terms).
    double normalForce = 0.0;           ///< F along the normal on the ring [N].
    double dForceDGapRate = 0.0;        ///< dF/dgdot [N s/m].
};

/** A contact element's state at evaluation. */
struct RingContactResponse
{
    double gap = 0.0;          ///< Signed gap [m].
    double gapRate = 0.0;      ///< Its rate [m/s].
    double normalForce = 0.0;  ///< Force along the normal on the ring [N]; never negative for a unilateral law.
    double potential = 0.0;    ///< Stored energy [J].
    double dissipationPower = 0.0; ///< Power removed [W].
};

/** A contact law's value at the measured gap and rate (from UnilateralSpringDamper::respond or any scalar law). */
struct ContactLawValue
{
    double force = 0.0;          ///< F along the normal on the ring [N].
    double dForceDGap = 0.0;     ///< dF/dg [N/m].
    double dForceDGapRate = 0.0; ///< dF/dgdot [N s/m].
};

/**
 * Measures a contact element's gap and gap rate at a state, the inputs of its law.
 *
 * @param model A planar model with rings.
 * @param state The state.
 * @param workspace Workspace with the state's kinematics (a prepared step).
 * @param geometry The contact's geometry.
 * @param gap Receives g [m].
 * @param gapRate Receives dg/dt [m/s].
 * @return Nothing, or InvalidInput.
 */
[[nodiscard]] std::expected<void, Fault> measureRingContact(const Model& model, const State& state,
                                                            const Workspace& workspace,
                                                            const RingContactGeometry& geometry, double& gap,
                                                            double& gapRate) noexcept;

/**
 * Forms a contact element's load from its law's value: a unilateral law on the signed gap along the witness normal,
 * applied at the casing's axis witness, with exact Jacobians in the local coordinates:
 * `K = J^T (dF/dg) J + F H + (dF/dgdot) J^T (H v)^T` and `D = J^T (dF/dgdot) J`, where H is the gap's Hessian. H
 * contains the host's rotation (the point's second derivative), and for a vertex witness the rotation of the normal
 * and the witness's sliding along the edge (`H = P^T (I - n n^T) P / rho + n^T P'' - phi_qt phi_tq / phi_tt`, the
 * envelope of the distance over the axis parameter). Past the axis (deep overlap at a vertex) the normal's terms are
 * omitted. See docs/WHEEL.md, "M2c status".
 *
 * @param model A planar model with rings.
 * @param state The state.
 * @param workspace Workspace with the state's kinematics (a prepared step).
 * @param geometry The contact's geometry.
 * @param law The law's force and derivatives at the measured gap and rate.
 * @param load Receives the load.
 * @return Nothing, or InvalidInput.
 */
[[nodiscard]] std::expected<void, Fault> ringContactLoad(const Model& model, const State& state,
                                                         const Workspace& workspace,
                                                         const RingContactGeometry& geometry,
                                                         const ContactLawValue& law, RingContactLoad& load) noexcept;

/**
 * Maps ring contact loads to generalized coordinates: the force `tau` and, when given, the exact dense Jacobians
 * (row-major dofCount x dofCount) including the host chain's own curvature. A test and oracle path: dense, so for
 * small models; the steps use the structured path.
 *
 * @param model A planar model with rings.
 * @param workspace Workspace with the state's kinematics.
 * @param loads The loads.
 * @param state The state (its rates enter the chain's rate transport).
 * @param tau Receives dofCount entries (added to).
 * @param stiffness Empty, or dofCount^2 entries (added to).
 * @param damping Empty, or dofCount^2 entries (added to).
 * @return Nothing, or InvalidInput.
 */
[[nodiscard]] std::expected<void, Fault> ringLoadsToGeneralized(const Model& model, const State& state,
                                                                const Workspace& workspace,
                                                                std::span<const RingContactLoad> loads,
                                                                std::span<double> tau, std::span<double> stiffness,
                                                                std::span<double> damping) noexcept;

} // namespace sub0physics::multibody
