#pragma once

#include <cstddef>

#include <array>
#include <span>

#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/ring_contact.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// The deformable ring substructure in the planar formulation (docs/WHEEL.md section 2.3): the hooks the planar
// algorithms call for a model with rings, the ring's elements and the condensation of its block of a step matrix
// onto its host's motion. Every function loops over the model's rings and never allocates.

namespace sub0physics::multibody::detail
{

// Which condensation a solve uses: the mass matrix's (forward dynamics) or the step matrix's.
enum class RingSystem : unsigned char
{
    Mass,
    Step
};

// Node positions (host frame) and, with velocities, node velocities; needs planar.velocity of the hosts.
void ringTransforms(const Model& model, const State& state, Workspace& ws, bool withVelocity) noexcept;

// Newton-Euler over the nodes, between the tree's forward and backward passes: node accelerations from their host's
// (planar.acceleration) and `qdd` (null for zero), tau for the ring coordinates, and the nodes' force added to their
// host's planar.force (also kept in RingWorkspace::hostForce). Needs ringTransforms with velocities.
void ringNewtonEuler(const Model& model, const State& state, const double* qdd, double* tau, Workspace& ws) noexcept;

// The ring's rigid inertia at its current shape added to its host's planar.composite, before the composite pass, and
// the border (mass coupling to the host's motion). Needs ringTransforms.
void ringCompositeInertia(const Model& model, Workspace& ws) noexcept;

// The host's motion per unit rate of each joint on its chain (RingWorkspace::hostJacobian); needs
// planar.parentToBody.
void ringHostJacobians(const Model& model, Workspace& ws) noexcept;

// The ring elements (foundation and edges): forces into elementForce's ring entries, potential and dissipation added,
// Jacobian blocks into the RingWorkspace when asked. False when an edge's nodes coincide or are not finite.
[[nodiscard]] bool ringElements(const Model& model, const State& state, Workspace& ws, bool jacobians,
                                bool geometricStiffness, double& potential, double& dissipationPower) noexcept;

// The nodes' kinetic and gravitational energy; needs ringTransforms with velocities and planar.worldToBody.
void ringMotionEnergy(const Model& model, const Workspace& ws, double& kinetic, double& potential) noexcept;

// Forms each ring's block of `M - hd D - h2 K~` (without the neighbour blocks of D when `neighbourDamping` is false: the
// reference engine's mass-matrix pattern, ImplicitTerms::DampingOnlyMassPattern) (the element blocks are read only when hd or h2 is not zero),
// factorises it, and removes its condensation from the tree's matrix `tree` (lower triangle, treeDofCount stride).
// False when a ring block is not positive definite. Needs ringCompositeInertia and ringHostJacobians.
[[nodiscard]] bool condenseRings(const Model& model, Workspace& ws, RingSystem system, double hd,
                                 bool neighbourDamping, double h2,
                                 double* tree) noexcept;

// The ring rows of a linearly implicit step's right-hand side: b_w = dt (tau + tau_e - c + thetaStep K~ wd), with K~
// from the element blocks (skipped when thetaStep is zero). `tau` may be empty; `elementForce` may be null.
void ringStepRhs(const Model& model, const State& state, const Workspace& ws, const double* tau,
                 const double* elementForce, double dt, double thetaStep, double* b) noexcept;

// A ring load in generalized coordinates over the coordinates it touches (the host chain, then the edge's four node
// coordinates): force, stiffness and damping (row-major, Capacity stride). `structured`: the rate transport through the
// chain's curvature only on the chain's rows (the rest does not factor through the host's motion; the steps drop it).
struct RingLoadBlocks
{
    static constexpr std::size_t MaxChain = 8;
    static constexpr std::size_t Capacity = MaxChain + 4;
    static constexpr std::size_t Size = Capacity * Capacity;
    std::size_t chain = 0;
    std::size_t count = 0;
    std::array<std::size_t, Capacity> index{};
    std::array<double, Capacity> force{};
    std::array<double, Size> stiffness{};
    std::array<double, Size> damping{};
};

// Needs the host Jacobians (ringHostJacobians) and the transforms at the state; the chain at most MaxChain long.
void ringLoadBlocks(const Model& model, const State& state, const Workspace& ws, const RingContactLoad& load,
                    bool structured, RingLoadBlocks& out) noexcept;

// The ring loads' Jacobians in a linearly implicit step `-hd D - h2 K~` (Forcing::ringLoads; their force is in the
// applied force): the host chain's block into the tree's matrix `tree` (lower triangle), the nodes' blocks and the
// border into each touched ring's contact blocks, which condenseRings adds for the step system, and
// `dt thetaStep K~ qd` into the right-hand side `b` (when h2 is not zero); without `neighbourDamping` (the mass pattern)
// the damping between the edge's two nodes stays explicit. Call after ringStepRhs.
void ringLoadStepTerms(const Model& model, const State& state, Workspace& ws, std::span<const RingContactLoad> loads,
                       double hd, bool neighbourDamping, double h2, double dt, double thetaStep, double* tree,
                       double* b) noexcept;

// Before the tree solve: eliminates the ring rows of x (dofCount entries) onto the tree rows.
void ringForward(const Model& model, const Workspace& ws, RingSystem system, double* x) noexcept;

// After the tree solve: the ring rows of x from the tree rows.
void ringBackward(const Model& model, const Workspace& ws, RingSystem system, double* x) noexcept;

} // namespace sub0physics::multibody::detail
