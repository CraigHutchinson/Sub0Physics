#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string_view>

#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>

#include "multibody_joint.hpp"
#include "multibody_planar.hpp"
#include "multibody_ring.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/implicit_euler.hpp"
#include "sub0physics/multibody/stepper.hpp"

namespace sub0physics::multibody
{

namespace
{

// Recursive Newton-Euler over transforms and velocities already in the workspace. With accelerations,
// tau = M qdd + c; without, tau = c.
// With external wrenches (body frame, one per body), they are subtracted from each body's force, so bodyForce ends
// as the wrench each joint transmits to its body.
template<bool WithAcceleration, bool WithExternal = false>
void recursiveNewtonEuler(const Model& model, const State& state, const double* qdd, double* tau, Workspace& ws,
                          const sub0math::SpatialForced* external = nullptr) noexcept
{
    const auto parents = model.parents();
    const auto subspace = model.motionSubspaces();
    const auto inertias = model.inertias();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const std::size_t n = parents.size();
    const double* qd = state.values.data() + model.dofCount();
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::size_t p = parents[i];
        const sub0math::SpatialMotiond& v = ws.bodyVelocity[i];
        const sub0math::SpatialMotiond& parentAcceleration =
            p == NoParent ? model.baseAcceleration() : ws.bodyAcceleration[p];
        sub0math::SpatialMotiond a = ws.parentToBody[i].apply(parentAcceleration) + cross(v, subspace[i] * qd[i]);
        if constexpr (WithAcceleration)
        {
            a += subspace[i] * qdd[i];
        }
        ws.bodyAcceleration[i] = a;
        ws.bodyForce[i] = inertias[i] * a + cross(v, inertias[i] * v);
        if constexpr (WithExternal)
        {
            ws.bodyForce[i] -= external[i];
        }
    }
    for (std::size_t i = n; i-- > 0;)
    {
        tau[i] = detail::subspaceDot(types[i], axes[i], ws.bodyForce[i]);
        if (const std::size_t p = parents[i]; p != NoParent)
        {
            ws.bodyForce[p] += ws.parentToBody[i].applyInverse(ws.bodyForce[i]);
        }
    }
}

// Composite-rigid-body algorithm over transforms already in the workspace.
void compositeRigidBody(const Model& model, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const auto inertias = model.inertias();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const std::size_t n = parents.size();
    double* h = ws.massMatrix.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        ws.compositeInertia[i] = inertias[i];
    }
    for (std::size_t i = n; i-- > 0;)
    {
        if (const std::size_t p = parents[i]; p != NoParent)
        {
            ws.compositeInertia[p] += ws.parentToBody[i].applyInverse(ws.compositeInertia[i]);
        }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        // PERF: sparse S products faster than general spatial ones.
        sub0math::SpatialForced f = detail::inertiaTimesSubspace(types[i], axes[i], ws.compositeInertia[i]);
        h[i * n + i] = detail::subspaceDot(types[i], axes[i], f);
        // Only ancestor pairs are written; the rest stays zero from the workspace's construction.
        for (std::size_t j = i; parents[j] != NoParent;)
        {
            f = ws.parentToBody[j].applyInverse(f);
            j = parents[j];
            const double value = detail::subspaceDot(types[j], axes[j], f);
            h[i * n + j] = value;
            h[j * n + i] = value;
        }
    }
}

// Featherstone, RBDA section 6.5: eliminate from the leaves; only ancestor entries are touched. Factorises
// the lower triangle of `h` in place: L below the diagonal, D on it. False if a pivot is not positive.
bool factorTree(const Model& model, double* h) noexcept
{
    const auto parents = model.parents();
    const std::size_t n = parents.size();
    for (std::size_t k = n; k-- > 0;)
    {
        const double pivot = h[k * n + k];
        if (!(pivot > 0.0) || !std::isfinite(pivot))
        {
            return false;
        }
        for (std::size_t i = parents[k]; i != NoParent; i = parents[i])
        {
            const double a = h[k * n + i] / pivot;
            for (std::size_t j = i; j != NoParent; j = parents[j])
            {
                h[i * n + j] -= a * h[k * n + j];
            }
            h[k * n + i] = a;
        }
    }
    return true;
}

std::expected<void, Fault> factorUnchecked(const Model& model, Workspace& workspace) noexcept
{
    const std::size_t n = model.treeDofCount();
    std::copy_n(workspace.massMatrix.data(), n * n, workspace.massFactor.data());
    // A ring's block is eliminated onto its host first (its border and host Jacobians come with the mass matrix).
    if (!model.rings().empty() &&
        !detail::condenseRings(model, workspace, detail::RingSystem::Mass, 0.0, true, 0.0, workspace.massFactor.data()))
    {
        return invalidInput("a ring's mass block is not positive definite");
    }
    if (!factorTree(model, workspace.massFactor.data()))
    {
        return invalidInput("mass matrix is not positive definite");
    }
    return {};
}

void solveTree(const Model& model, const double* h, double* x) noexcept
{
    const auto parents = model.parents();
    const std::size_t n = parents.size();
    // M = L^T D L: solve L^T y = b, then D z = y, then L x = z.
    for (std::size_t i = n; i-- > 0;)
    {
        for (std::size_t j = parents[i]; j != NoParent; j = parents[j])
        {
            x[j] -= h[i * n + j] * x[i];
        }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] /= h[i * n + i];
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = parents[i]; j != NoParent; j = parents[j])
        {
            x[i] -= h[i * n + j] * x[j];
        }
    }
}

void solveUnchecked(const Model& model, const Workspace& workspace, std::span<double> x) noexcept
{
    const bool rings = !model.rings().empty();
    if (rings)
    {
        detail::ringForward(model, workspace, detail::RingSystem::Mass, x.data());
    }
    solveTree(model, workspace.massFactor.data(), x.data());
    if (rings)
    {
        detail::ringBackward(model, workspace, detail::RingSystem::Mass, x.data());
    }
}

// Dense A = L D L^T of the lower triangle of a row-major n x n matrix, in place (L below the diagonal, D on
// it), for a matrix that the elements have coupled across branches. `scratch` holds n values. False if a
// pivot is not positive: the matrix is not positive definite.
bool factorDense(std::size_t n, double* a, double* scratch) noexcept
{
    for (std::size_t j = 0; j < n; ++j)
    {
        double* rowJ = a + j * n;
        double pivot = rowJ[j];
        for (std::size_t k = 0; k < j; ++k)
        {
            scratch[k] = rowJ[k] * a[k * n + k]; // L_jk d_k
            pivot -= rowJ[k] * scratch[k];
        }
        if (!(pivot > 0.0) || !std::isfinite(pivot))
        {
            return false;
        }
        rowJ[j] = pivot;
        const double inverse = 1.0 / pivot;
        for (std::size_t i = j + 1; i < n; ++i)
        {
            double* rowI = a + i * n;
            double sum = rowI[j];
            for (std::size_t k = 0; k < j; ++k)
            {
                sum -= rowI[k] * scratch[k];
            }
            rowI[j] = sum * inverse;
        }
    }
    return true;
}

// Solves L D L^T x = b in place with the factor from factorDense.
void solveDense(std::size_t n, const double* a, double* x) noexcept
{
    for (std::size_t i = 0; i < n; ++i)
    {
        double sum = x[i];
        for (std::size_t k = 0; k < i; ++k)
        {
            sum -= a[i * n + k] * x[k];
        }
        x[i] = sum;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] /= a[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;)
    {
        double sum = x[i];
        for (std::size_t k = i + 1; k < n; ++k)
        {
            sum -= a[k * n + i] * x[k];
        }
        x[i] = sum;
    }
}

// Elements at a state whose parent-to-body transforms are in the workspace (planar ones for a planar model).
std::expected<void, Fault> addElements(const Model& model, const State& state, Workspace& ws, bool jacobians,
                                       bool geometricStiffness = true) noexcept
{
    const bool planar = detail::isPlanar(model);
    if (model.tendonCount() > 0)
    {
        if (planar)
        {
            detail::planarWorldTransforms(model, ws);
        }
        else
        {
            detail::composeWorldTransforms(model, ws);
        }
    }
    return detail::evaluateElementsUnchecked(model, state, ws, jacobians, geometricStiffness, planar);
}

[[nodiscard]] std::expected<void, Fault> checkStep(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                   double dt, const State& state, Workspace& ws) noexcept
{
    if (!std::isfinite(dt) || !(dt > 0.0))
    {
        return invalidInput("step needs a finite positive step size");
    }
    if (tick != state.tick)
    {
        return invalidInput("tick differs from the state's tick");
    }
    if (auto checked = detail::checkInputs(model, state, ws); !checked)
    {
        return checked;
    }
    return detail::checkForcing(model, forcing);
}

// The start of an implicit (or constrained) step: transforms, velocities, bias forces, mass matrix and the elements
// (with their Jacobians when asked).
std::expected<void, Fault> prepareStep(const Model& model, const State& state, Workspace& ws, bool jacobians,
                                       bool geometricStiffness = true) noexcept
{
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, ws);
        detail::planarNewtonEuler<false>(model, state, nullptr, ws.bias.data(), ws);
        detail::planarCompositeRigidBody(model, ws);
    }
    else
    {
        detail::updateTransforms<true>(model, state, ws);
        recursiveNewtonEuler<false>(model, state, nullptr, ws.bias.data(), ws);
        compositeRigidBody(model, ws);
    }
    if (model.hasForceElements())
    {
        return addElements(model, state, ws, jacobians, geometricStiffness);
    }
    return {};
}

// Forms and factorises the step matrix in `implicitMatrix` and its right-hand side in `implicitRhs`:
// implicit: (M - h D - h^2 K~) dqd = h (tau + tau_e - c + h K~ qd), with K~ the symmetric part of K (a valid
// linearization for a W-method: the step stays first order with any Jacobian; only stability depends on it);
// explicit: M dqd = h (tau + tau_e - c). D and K are the elements' plus the forcing's; DampingOnly drops K. Only the
// lower triangle is formed; both factorisations read only it. Needs prepareStep at `state`; `elementForce` is tau_e
// (null without elements). Returns whether the factor is dense (else along the tree). The geometric fallback
// re-evaluates the elements, so it needs the start's transforms (`allowFallback`).
std::expected<bool, Fault> formStepSystem(const Model& model, const State& state, const Forcing& applied,
                                          const double* elementForce, Workspace& ws, double dt, bool implicit,
                                          ImplicitTerms terms, bool allowFallback, double theta = 1.0) noexcept
{
    const std::size_t n = model.treeDofCount();
    const double* qd = state.values.data() + model.dofCount();
    const double* m = ws.massMatrix.data();
    const double* k = ws.elementStiffness.data();
    const double* d = ws.elementDamping.data();
    const double* fk = applied.stiffness.empty() ? nullptr : applied.stiffness.data();
    const double* fd = applied.damping.empty() ? nullptr : applied.damping.data();
    const std::span<const double> tau = applied.tau;
    double* a = ws.implicitMatrix.data();
    double* b = ws.implicitRhs.data();
    const bool hasElements = model.hasForceElements();
    const bool rings = !model.rings().empty();
    // theta = 1: backward Euler linearised; theta = 1/2: the linearised trapezoid (implicit midpoint).
    const double thetaStep = theta * dt;
    const double hd = implicit ? thetaStep : 0.0;
    // DampingOnlyMassPattern: keep D only between coordinates whose bodies are on one chain (the mass matrix's pattern).
    const bool massPattern = terms == ImplicitTerms::DampingOnlyMassPattern && !model.elementsFollowTree();
    const auto related = [&](std::size_t i, std::size_t j) {
        const auto parents = model.parents();
        for (std::size_t body = i; body != NoParent; body = parents[body])
        {
            if (body == j)
            {
                return true;
            }
        }
        return false;
    };
    const auto assemble = [&](bool withStiffness) {
        const double h2 = withStiffness ? thetaStep * thetaStep : 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            double stiffnessRate = 0.0;
            for (std::size_t j = 0; j <= i; ++j)
            {
                double damping = d[i * n + j];
                if (massPattern && !related(i, j))
                {
                    damping = 0.0; // outside the mass matrix's sparsity: explicit, as the reference engine
                }
                if (fd != nullptr)
                {
                    damping += 0.5 * (fd[i * n + j] + fd[j * n + i]);
                }
                double entry = m[i * n + j] - hd * damping;
                if (withStiffness)
                {
                    // PERF: K~ only when stiffness is implicit (not for DampingOnly or explicit steps).
                    double symmetric = 0.5 * (k[i * n + j] + k[j * n + i]);
                    if (fk != nullptr)
                    {
                        symmetric += 0.5 * (fk[i * n + j] + fk[j * n + i]);
                    }
                    entry -= h2 * symmetric;
                }
                a[i * n + j] = entry;
            }
            if (withStiffness)
            {
                for (std::size_t j = 0; j < n; ++j)
                {
                    double symmetric = 0.5 * (k[i * n + j] + k[j * n + i]);
                    if (fk != nullptr)
                    {
                        symmetric += 0.5 * (fk[i * n + j] + fk[j * n + i]);
                    }
                    stiffnessRate += symmetric * qd[j];
                }
            }
            const double force = tau.empty() ? 0.0 : tau[i];
            const double element = elementForce == nullptr ? 0.0 : elementForce[i];
            b[i] = dt * (force + element - ws.bias[i] + thetaStep * stiffnessRate);
        }
        if (rings)
        {
            // The ring rows, and the rings condensed onto their hosts (block storage, docs/WHEEL.md section 2.3).
            detail::ringStepRhs(model, state, ws, tau.empty() ? nullptr : tau.data(), elementForce, dt,
                                withStiffness ? thetaStep : 0.0, b);
            // Contact loads' Jacobians (their force is already in tau), or none: clears last step's terms.
            detail::ringLoadStepTerms(model, state, ws, implicit ? applied.ringLoads : std::span<const RingContactLoad>{},
                                      hd, terms != ImplicitTerms::DampingOnlyMassPattern, h2, dt, thetaStep, a, b);
            return detail::condenseRings(model, ws, detail::RingSystem::Step, hd,
                                         terms != ImplicitTerms::DampingOnlyMassPattern,
                                         withStiffness ? thetaStep * thetaStep : 0.0, a);
        }
        return true;
    };
    // Forcing Jacobians may couple any coordinates: they take the dense factorisation.
    const bool dense = implicit && (!model.elementsFollowTree() || fk != nullptr || fd != nullptr);
    const auto factor = [&]() { return dense ? factorDense(n, a, ws.elementScalars.data()) : factorTree(model, a); };
    // Without elements, K and D are zero from the workspace's construction and never written.
    const bool stiffness = implicit && (hasElements || fk != nullptr) && terms != ImplicitTerms::DampingOnly &&
                           terms != ImplicitTerms::DampingOnlyMassPattern;
    const bool condensed = assemble(stiffness);
    ws.geometricStiffnessDropped = false;
    if (!condensed || !factor())
    {
        if (!stiffness || !hasElements || !allowFallback || terms != ImplicitTerms::Full)
        {
            return invalidInput("step matrix is not positive definite");
        }
        // K made the matrix indefinite: a tendon's geometric term (a compressed strut, a softening that grows
        // faster than h can follow). Fall back to the material stiffness, which is negative semidefinite for
        // passive laws, so M - h D - h^2 K_material is positive definite and stiff elements stay implicit.
        ws.geometricStiffnessDropped = true;
        if (auto material = detail::evaluateElementsUnchecked(model, state, ws, true, false, detail::isPlanar(model));
            !material)
        {
            return std::unexpected(material.error());
        }
        if (!assemble(true) || !factor())
        {
            return invalidInput("implicit step matrix is not positive definite, even without geometric stiffness");
        }
    }
    return dense;
}

// Solves the factorised step matrix in place.
void solveStepSystem(const Model& model, const Workspace& ws, bool dense, double* x) noexcept
{
    const bool rings = !model.rings().empty();
    if (rings)
    {
        detail::ringForward(model, ws, detail::RingSystem::Step, x);
    }
    if (dense)
    {
        solveDense(model.treeDofCount(), ws.implicitMatrix.data(), x);
    }
    else
    {
        solveTree(model, ws.implicitMatrix.data(), x);
    }
    if (rings)
    {
        detail::ringBackward(model, ws, detail::RingSystem::Step, x);
    }
}

// The free (unconstrained) accelerations M^-1 (tau + tau_e - c) at `state`, the dynamics written on the way.
std::expected<void, Fault> freeAccelerations(const Model& model, const State& state, std::span<const double> tau,
                                             Workspace& ws) noexcept
{
    const std::size_t n = model.dofCount();
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, ws);
        detail::planarNewtonEuler<false>(model, state, nullptr, ws.bias.data(), ws);
        detail::planarCompositeRigidBody(model, ws);
    }
    else
    {
        detail::updateTransforms<true>(model, state, ws);
        recursiveNewtonEuler<false>(model, state, nullptr, ws.bias.data(), ws);
        compositeRigidBody(model, ws);
    }
    if (auto factored = factorUnchecked(model, ws); !factored)
    {
        return factored;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        ws.qdd[i] = (tau.empty() ? 0.0 : tau[i]) - ws.bias[i];
    }
    if (model.hasForceElements())
    {
        if (auto added = addElements(model, state, ws, false); !added)
        {
            return added;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            ws.qdd[i] += ws.elementForce[i];
        }
    }
    solveUnchecked(model, ws, ws.qdd);
    return {};
}

// Forward dynamics after the input checks; `tau` is the applied generalized force.
std::expected<void, Fault> forwardUnchecked(const Model& model, const State& state, std::span<const double> tau,
                                            Workspace& ws) noexcept
{
    if (auto free = freeAccelerations(model, state, tau, ws); !free)
    {
        return free;
    }
    if (model.constraintCount() > 0)
    {
        if (detail::isPlanar(model))
        {
            detail::expandPlanarKinematics(model, ws, true);
        }
        return detail::constrainAccelerations(model, state, ws);
    }
    return {};
}

} // namespace

std::expected<void, Fault> computeMassMatrix(const Model& model, const State& state, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<false>(model, state, workspace);
        detail::planarCompositeRigidBody(model, workspace);
        return {};
    }
    detail::updateTransforms<false>(model, state, workspace);
    compositeRigidBody(model, workspace);
    return {};
}

std::expected<void, Fault> computeBiasForces(const Model& model, const State& state, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, workspace);
        detail::planarNewtonEuler<false>(model, state, nullptr, workspace.bias.data(), workspace);
        return {};
    }
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<false>(model, state, nullptr, workspace.bias.data(), workspace);
    return {};
}

std::expected<void, Fault> inverseDynamics(const Model& model, const State& state, std::span<const double> qdd,
                                           std::span<double> tau, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (qdd.size() != model.dofCount() || tau.size() != model.dofCount())
    {
        return invalidInput("inverse dynamics needs qdd and tau of dofCount entries");
    }
    if (detail::isPlanar(model))
    {
        detail::planarTransforms<true>(model, state, workspace);
        detail::planarNewtonEuler<true>(model, state, qdd.data(), tau.data(), workspace);
    }
    else
    {
        detail::updateTransforms<true>(model, state, workspace);
        recursiveNewtonEuler<true>(model, state, qdd.data(), tau.data(), workspace);
    }
    if (model.hasForceElements())
    {
        if (auto added = addElements(model, state, workspace, false); !added)
        {
            return added;
        }
        for (std::size_t i = 0; i < tau.size(); ++i)
        {
            tau[i] -= workspace.elementForce[i];
        }
    }
    return {};
}

std::expected<void, Fault> factorMassMatrix(const Model& model, Workspace& workspace) noexcept
{
    if (!workspace.isBuiltFor(model))
    {
        return invalidInput("workspace was not built for this model");
    }
    return factorUnchecked(model, workspace);
}

std::expected<void, Fault> solveMassMatrix(const Model& model, const Workspace& workspace, std::span<double> x) noexcept
{
    if (!workspace.isBuiltFor(model))
    {
        return invalidInput("workspace was not built for this model");
    }
    if (x.size() != model.dofCount())
    {
        return invalidInput("mass matrix solve needs x of dofCount entries");
    }
    solveUnchecked(model, workspace, x);
    return {};
}

std::expected<void, Fault> forwardDynamics(const Model& model, const State& state, std::span<const double> tau,
                                           Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (!tau.empty() && tau.size() != model.dofCount())
    {
        return invalidInput("forward dynamics needs tau of dofCount entries, or none");
    }
    return forwardUnchecked(model, state, tau, workspace);
}

std::expected<void, Fault> forwardDynamicsWithForcing(const Model& model, const State& state, const Forcing& forcing,
                                           Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (auto checked = detail::checkForcing(model, forcing); !checked)
    {
        return checked;
    }
    const Forcing applied = detail::applyForcing(model, state, forcing, workspace);
    return forwardUnchecked(model, state, applied.tau, workspace);
}

std::expected<void, Fault> computeAppliedForce(const Model& model, const State& state, const Forcing& forcing,
                                               Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    if (auto checked = detail::checkForcing(model, forcing); !checked)
    {
        return checked;
    }
    const Forcing applied = detail::applyForcing(model, state, forcing, workspace);
    if (applied.tau.data() != workspace.appliedForce.data())
    {
        for (std::size_t i = 0; i < model.dofCount(); ++i)
        {
            workspace.appliedForce[i] = applied.tau.empty() ? 0.0 : applied.tau[i];
        }
    }
    return {};
}

std::expected<void, Fault> stepSemiImplicitEuler(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                 sub0math::f64::Time h, State& state, Workspace& workspace) noexcept
{
    const double dt = h.si();
    if (auto checked = checkStep(model, forcing, tick, dt, state, workspace); !checked)
    {
        return checked;
    }
    const Forcing applied = detail::applyForcing(model, state, forcing, workspace);
    if (model.constraintCount() > 0)
    {
        return detail::stepConstrained(model, applied, dt, state, workspace, false);
    }
    if (auto accelerated = forwardUnchecked(model, state, applied.tau, workspace); !accelerated)
    {
        return accelerated;
    }
    const std::size_t n = model.dofCount();
    double* q = state.values.data();
    double* qd = q + n;
    // Validate the whole endpoint before committing any of it.
    for (std::size_t i = 0; i < n; ++i)
    {
        const double rate = qd[i] + dt * workspace.qdd[i];
        if (!std::isfinite(rate) || !std::isfinite(q[i] + dt * rate))
        {
            return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
        }
    }
    double work = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        qd[i] += dt * workspace.qdd[i];
        const double step = dt * qd[i];
        q[i] += step;
        work += (applied.tau.empty() ? 0.0 : applied.tau[i]) * step;
    }
    if (model.hasForceElements())
    {
        state.dissipated += sub0math::f64::Energy::fromSi(dt * workspace.elementDissipationPower);
    }
    state.externalWork += sub0math::f64::Energy::fromSi(work);
    state.tick = tick + 1U;
    return {};
}

// The rest of a linearly implicit theta step after prepareStep at `state`:
// (M - theta h D - theta^2 h^2 K~) dqd = h (tau + tau_e - c + theta h K~ qd), then q += h (qd + theta dqd),
// qd += dqd. theta = 1 is linearly implicit Euler (bit for bit its earlier form).
std::expected<void, Fault> finishLinearlyImplicit(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                  double dt, State& state, Workspace& ws, ImplicitTerms terms,
                                                  double theta) noexcept
{
    const std::size_t n = model.dofCount();
    // PERF: site forces mapped on the step's own transforms, not a second pass.
    const Forcing applied = detail::applyForcing(model, state, forcing, ws, true);
    const bool hasElements = model.hasForceElements();
    const auto dense = formStepSystem(model, state, applied, hasElements ? ws.elementForce.data() : nullptr, ws, dt,
                                      true, terms, true, theta);
    if (!dense)
    {
        return std::unexpected(dense.error());
    }
    double* b = ws.implicitRhs.data();
    solveStepSystem(model, ws, *dense, b);

    // Validate the whole endpoint before committing any of it.
    const double* q = state.values.data();
    const double* qd = q + n;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double rate = qd[i] + b[i];
        if (!std::isfinite(rate) || !std::isfinite(q[i] + dt * (qd[i] + theta * b[i])))
        {
            return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
        }
    }
    double* qOut = state.values.data();
    double* qdOut = qOut + n;
    double work = 0.0;
    const double inverseStep = 1.0 / dt;
    for (std::size_t i = 0; i < n; ++i)
    {
        ws.qdd[i] = b[i] * inverseStep;
        const double step = dt * (qdOut[i] + theta * b[i]);
        qdOut[i] += b[i];
        qOut[i] += step;
        work += (applied.tau.empty() ? 0.0 : applied.tau[i]) * step;
    }
    if (hasElements)
    {
        state.dissipated += sub0math::f64::Energy::fromSi(dt * ws.elementDissipationPower);
    }
    state.externalWork += sub0math::f64::Energy::fromSi(work);
    state.tick = tick + 1U;
    return {};
}

std::expected<void, Fault> stepLinearlyImplicit(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                double dt, State& state, Workspace& ws, ImplicitTerms terms,
                                                double theta) noexcept
{
    if (auto checked = checkStep(model, forcing, tick, dt, state, ws); !checked)
    {
        return checked;
    }
    if (model.constraintCount() > 0)
    {
        if (theta != 1.0)
        {
            return invalidInput("the linearly implicit midpoint step does not support constraints yet");
        }
        return detail::stepConstrained(model, detail::applyForcing(model, state, forcing, ws), dt, state, ws, true,
                                       terms);
    }
    if (auto prepared = prepareStep(model, state, ws, true, terms == ImplicitTerms::Full); !prepared)
    {
        return prepared;
    }
    return finishLinearlyImplicit(model, forcing, tick, dt, state, ws, terms, theta);
}

std::expected<void, Fault> stepLinearlyImplicitEuler(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                     sub0math::f64::Time h, State& state, Workspace& ws,
                                                     ImplicitTerms terms) noexcept
{
    return stepLinearlyImplicit(model, forcing, tick, h.si(), state, ws, terms, 1.0);
}

std::expected<void, Fault> stepLinearlyImplicitMidpoint(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                        sub0math::f64::Time h, State& state, Workspace& ws,
                                                        ImplicitTerms terms) noexcept
{
    return stepLinearlyImplicit(model, forcing, tick, h.si(), state, ws, terms, 0.5);
}

std::expected<void, Fault> prepareLinearlyImplicitStep(const Model& model, const State& state, Workspace& ws,
                                                       ImplicitTerms terms) noexcept
{
    if (auto checked = detail::checkInputs(model, state, ws); !checked)
    {
        return checked;
    }
    if (model.constraintCount() > 0)
    {
        return invalidInput("a prepared step does not support constraints; use stepLinearlyImplicitEuler");
    }
    if (auto prepared = prepareStep(model, state, ws, true, terms == ImplicitTerms::Full); !prepared)
    {
        return prepared;
    }
    // The world transforms the site queries read (the elements composed them only when the model has tendons).
    if (detail::isPlanar(model))
    {
        detail::planarWorldTransforms(model, ws);
    }
    else
    {
        detail::composeWorldTransforms(model, ws);
        for (std::size_t i = 0; i < model.bodyCount(); ++i)
        {
            ws.bodyPose[i] = sub0math::toRigidTransform<sub0math::f64::Length>(ws.worldToBody[i]);
        }
    }
    std::ranges::copy(state.values, ws.preparedValues.begin());
    ws.preparedTerms = terms;
    ws.stepPrepared = true;
    return {};
}

std::expected<void, Fault> completeLinearlyImplicitStep(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                        sub0math::f64::Time h, State& state, Workspace& ws,
                                                        ImplicitScheme scheme) noexcept
{
    if (!ws.isBuiltFor(model) || !state.isSizedFor(model) || !ws.stepPrepared ||
        !std::ranges::equal(state.values, ws.preparedValues))
    {
        return invalidInput("the workspace holds no step prepared at this state");
    }
    const double dt = h.si();
    if (!std::isfinite(dt) || !(dt > 0.0))
    {
        return invalidInput("step needs a finite positive step size");
    }
    if (tick != state.tick)
    {
        return invalidInput("tick differs from the state's tick");
    }
    if (auto checked = detail::checkForcing(model, forcing); !checked)
    {
        return checked;
    }
    ws.stepPrepared = false;
    return finishLinearlyImplicit(model, forcing, tick, dt, state, ws, ws.preparedTerms,
                                  scheme == ImplicitScheme::Midpoint ? 0.5 : 1.0);
}

namespace detail
{

std::expected<void, Fault> checkForcing(const Model& model, const Forcing& forcing) noexcept
{
    const std::size_t n = model.dofCount();
    if (!forcing.tau.empty() && forcing.tau.size() != n)
    {
        return invalidInput("forcing needs tau of dofCount entries, or none");
    }
    if ((!forcing.stiffness.empty() && forcing.stiffness.size() != n * n) ||
        (!forcing.damping.empty() && forcing.damping.size() != n * n))
    {
        return invalidInput("forcing Jacobians need dofCount x dofCount entries, or none");
    }
    if (!model.rings().empty() && (!forcing.stiffness.empty() || !forcing.damping.empty()))
    {
        return invalidInput("dense forcing Jacobians are not supported with ring substructures: use Forcing::ringLoads");
    }
    for (const RingContactLoad& load : forcing.ringLoads)
    {
        if (!isPlanar(model) || load.ring >= model.rings().size() || load.edge >= model.rings()[load.ring].nodeCount ||
            model.rings()[load.ring].hostChain.size() > RingLoadBlocks::MaxChain)
        {
            return invalidInput("a ring load names a ring or edge that is not in the model, or its host chain is too long");
        }
    }
    for (const SiteForce& f : forcing.siteForces)
    {
        if (f.site.value() >= model.siteCount())
        {
            return invalidInput("site force names a site that is not in the model");
        }
    }
    return {};
}

Forcing applyForcing(const Model& model, const State& state, const Forcing& forcing, Workspace& ws,
                     bool transformsCurrent) noexcept
{
    if (forcing.siteForces.empty() && forcing.ringLoads.empty())
    {
        return forcing; // nothing to map: the caller's tau runs as given
    }
    const std::size_t n = model.dofCount();
    double* out = ws.appliedForce.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        out[i] = forcing.tau.empty() ? 0.0 : forcing.tau[i];
    }
    const auto parents = model.parents();
    if (forcing.siteForces.empty())
    {
        // Ring loads only (a planar model, checked): through the host's motion.
        if (!transformsCurrent)
        {
            planarTransforms<false>(model, state, ws);
            ringHostJacobians(model, ws);
        }
    }
    else if (isPlanar(model))
    {
        if (!transformsCurrent)
        {
            planarTransforms<false>(model, state, ws);
        }
        planarWorldTransforms(model, ws);
        const PlanarModel& pm = model.planar();
        for (const SiteForce& f : forcing.siteForces)
        {
            const std::size_t body = model.siteBodies()[f.site.value()];
            if (body == NoParent)
            {
                continue;
            }
            const PlanarPoint p = planarSitePoint(model, ws, f.site.value());
            const double fx = f.force[pm.axisI];
            const double fy = f.force[pm.axisJ];
            const double mk = f.moment[pm.normal];
            for (std::size_t j = body; j != NoParent; j = parents[j])
            {
                const PlanarBody& b = pm.bodies[j];
                const PlanarTransform& x = ws.planar.worldToBody[j];
                if (b.revolute)
                {
                    out[j] += b.axisSign * (mk + (p.x - x.x) * fy - (p.y - x.y) * fx);
                }
                else
                {
                    out[j] += (x.c * b.axisX - x.s * b.axisY) * fx + (x.s * b.axisX + x.c * b.axisY) * fy;
                }
            }
        }
    }
    else
    {
        if (!transformsCurrent)
        {
            updateTransforms<false>(model, state, ws);
        }
        composeWorldTransforms(model, ws);
        const auto subspace = model.motionSubspaces();
        for (const SiteForce& f : forcing.siteForces)
        {
            const std::size_t body = model.siteBodies()[f.site.value()];
            if (body == NoParent)
            {
                continue;
            }
            const sub0math::PluckerTransformd& xb = ws.worldToBody[body];
            const auto& offset = model.siteFrame(f.site).translation;
            const sub0math::Vec3d p =
                xb.translation + bodyToWorld(xb.rotation, sub0math::Vec3d{offset.x.si(), offset.y.si(), offset.z.si()});
            for (std::size_t j = body; j != NoParent; j = parents[j])
            {
                const sub0math::PluckerTransformd& x = ws.worldToBody[j];
                const sub0math::Vec3d omega = bodyToWorld(x.rotation, subspace[j].angular);
                const sub0math::Vec3d linear = bodyToWorld(x.rotation, subspace[j].linear);
                out[j] += dot(linear + cross(omega, p - x.translation), f.force) + dot(omega, f.moment);
            }
        }
    }
    if (!forcing.ringLoads.empty())
    {
        if (!transformsCurrent && !forcing.siteForces.empty())
        {
            ringHostJacobians(model, ws);
        }
        (void)multibody::ringLoadsToGeneralized(model, state, ws, forcing.ringLoads, ws.appliedForce, {}, {});
    }
    return Forcing{ws.appliedForce, {}, forcing.stiffness, forcing.damping, forcing.ringLoads};
}

std::expected<void, Fault> prepareStepDynamics(const Model& model, const State& state, Workspace& ws,
                                               bool elementJacobians, ImplicitTerms terms) noexcept
{
    if (auto prepared = prepareStep(model, state, ws, elementJacobians, terms == ImplicitTerms::Full); !prepared)
    {
        return prepared;
    }
    if (isPlanar(model))
    {
        expandPlanarKinematics(model, ws, true);
        return {};
    }
    composeWorldTransforms(model, ws);
    return {};
}

std::expected<bool, Fault> formStepMatrix(const Model& model, const State& state, const Forcing& applied,
                                          const double* elementForce, Workspace& ws, double dt, bool implicit,
                                          ImplicitTerms terms, bool allowFallback) noexcept
{
    return formStepSystem(model, state, applied, elementForce, ws, dt, implicit, terms, allowFallback);
}

void solveStepMatrix(const Model& model, const Workspace& ws, bool dense, double* x) noexcept
{
    solveStepSystem(model, ws, dense, x);
}

void solveMassFactor(const Model& model, const Workspace& ws, double* x) noexcept
{
    solveTree(model, ws.massFactor.data(), x);
}

void transmittedForces(const Model& model, const State& state, const double* qdd, double* tau,
                       const sub0math::SpatialForced* external, Workspace& ws) noexcept
{
    updateTransforms<true>(model, state, ws);
    composeWorldTransforms(model, ws);
    recursiveNewtonEuler<true, true>(model, state, qdd, tau, ws, external);
}

} // namespace detail

} // namespace sub0physics::multibody
