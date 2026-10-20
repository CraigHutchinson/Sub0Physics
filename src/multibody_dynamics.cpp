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
    const double* qd = state.values.data() + n;
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
    const std::size_t n = model.dofCount();
    std::copy_n(workspace.massMatrix.data(), n * n, workspace.massFactor.data());
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
    solveTree(model, workspace.massFactor.data(), x.data());
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

// Elements at a state whose parent-to-body transforms are in the workspace.
std::expected<void, Fault> addElements(const Model& model, const State& state, Workspace& ws, bool jacobians) noexcept
{
    if (model.tendonCount() > 0)
    {
        detail::composeWorldTransforms(model, ws);
    }
    return detail::evaluateElementsUnchecked(model, state, ws, jacobians);
}

[[nodiscard]] std::expected<void, Fault> checkStep(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                   double dt, const State& state, const Workspace& ws) noexcept
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
    if (!forcing.tau.empty() && forcing.tau.size() != model.dofCount())
    {
        return invalidInput("forcing needs tau of dofCount entries, or none");
    }
    return {};
}

// The start of an implicit (or constrained) step: transforms, velocities, bias forces, mass matrix and the elements
// (with their Jacobians when asked).
std::expected<void, Fault> prepareStep(const Model& model, const State& state, Workspace& ws, bool jacobians) noexcept
{
    detail::updateTransforms<true>(model, state, ws);
    recursiveNewtonEuler<false>(model, state, nullptr, ws.bias.data(), ws);
    compositeRigidBody(model, ws);
    if (model.elementCount() > 0)
    {
        return addElements(model, state, ws, jacobians);
    }
    return {};
}

// Forms and factorises the step matrix in `implicitMatrix` and its right-hand side in `implicitRhs`:
// implicit: (M - h D - h^2 K~) dqd = h (tau + tau_e - c + h K~ qd), with K~ the symmetric part of K (a valid
// linearization for a W-method: the step stays first order with any Jacobian; only stability depends on it);
// explicit: M dqd = h (tau + tau_e - c). Only the lower triangle is formed; both factorisations read only it.
// Needs prepareStep at `state`; `elementForce` is tau_e (null without elements). Returns whether the factor is
// dense (else along the tree). The geometric fallback re-evaluates the elements, so it needs the start's transforms
// (`allowFallback`).
std::expected<bool, Fault> formStepSystem(const Model& model, const State& state, std::span<const double> tau,
                                          const double* elementForce, Workspace& ws, double dt, bool implicit,
                                          bool allowFallback) noexcept
{
    const std::size_t n = model.dofCount();
    const double* qd = state.values.data() + n;
    const double* m = ws.massMatrix.data();
    const double* k = ws.elementStiffness.data();
    const double* d = ws.elementDamping.data();
    double* a = ws.implicitMatrix.data();
    double* b = ws.implicitRhs.data();
    const bool hasElements = model.elementCount() > 0;
    const double hd = implicit ? dt : 0.0;
    const auto assemble = [&](bool withStiffness) {
        const double h2 = withStiffness ? dt * dt : 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            double stiffnessRate = 0.0;
            for (std::size_t j = 0; j <= i; ++j)
            {
                const double symmetric = 0.5 * (k[i * n + j] + k[j * n + i]);
                a[i * n + j] = m[i * n + j] - hd * d[i * n + j] - h2 * symmetric;
            }
            if (withStiffness)
            {
                for (std::size_t j = 0; j < n; ++j)
                {
                    stiffnessRate += 0.5 * (k[i * n + j] + k[j * n + i]) * qd[j];
                }
            }
            const double applied = tau.empty() ? 0.0 : tau[i];
            const double element = elementForce == nullptr ? 0.0 : elementForce[i];
            b[i] = dt * (applied + element - ws.bias[i] + dt * stiffnessRate);
        }
    };
    const bool dense = implicit && !model.elementsFollowTree();
    const auto factor = [&]() { return dense ? factorDense(n, a, ws.elementScalars.data()) : factorTree(model, a); };
    // Without elements, K and D are zero from the workspace's construction and never written.
    assemble(implicit && hasElements);
    ws.geometricStiffnessDropped = false;
    if (!factor())
    {
        if (!implicit || !hasElements || !allowFallback)
        {
            return invalidInput("step matrix is not positive definite");
        }
        // K made the matrix indefinite: a tendon's geometric term (a compressed strut, a softening that grows
        // faster than h can follow). Fall back to the material stiffness, which is negative semidefinite for
        // passive laws, so M - h D - h^2 K_material is positive definite and stiff elements stay implicit.
        ws.geometricStiffnessDropped = true;
        if (auto material = detail::evaluateElementsUnchecked(model, state, ws, true, false); !material)
        {
            return std::unexpected(material.error());
        }
        assemble(true);
        if (!factor())
        {
            return invalidInput("implicit step matrix is not positive definite, even without geometric stiffness");
        }
    }
    return dense;
}

// Solves the factorised step matrix in place.
void solveStepSystem(const Model& model, const Workspace& ws, bool dense, double* x) noexcept
{
    if (dense)
    {
        solveDense(model.dofCount(), ws.implicitMatrix.data(), x);
    }
    else
    {
        solveTree(model, ws.implicitMatrix.data(), x);
    }
}

} // namespace

std::expected<void, Fault> computeMassMatrix(const Model& model, const State& state, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
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
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<true>(model, state, qdd.data(), tau.data(), workspace);
    if (model.elementCount() > 0)
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
    const std::size_t n = model.dofCount();
    if (!tau.empty() && tau.size() != n)
    {
        return invalidInput("forward dynamics needs tau of dofCount entries, or none");
    }
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<false>(model, state, nullptr, workspace.bias.data(), workspace);
    compositeRigidBody(model, workspace);
    if (auto factored = factorUnchecked(model, workspace); !factored)
    {
        return factored;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        workspace.qdd[i] = (tau.empty() ? 0.0 : tau[i]) - workspace.bias[i];
    }
    if (model.elementCount() > 0)
    {
        if (auto added = addElements(model, state, workspace, false); !added)
        {
            return added;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            workspace.qdd[i] += workspace.elementForce[i];
        }
    }
    solveUnchecked(model, workspace, workspace.qdd);
    if (model.constraintCount() > 0)
    {
        return detail::constrainAccelerations(model, state, workspace);
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
    if (model.constraintCount() > 0)
    {
        return detail::stepConstrained(model, forcing, dt, state, workspace, false);
    }
    if (auto accelerated = forwardDynamics(model, state, forcing.tau, workspace); !accelerated)
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
        work += (forcing.tau.empty() ? 0.0 : forcing.tau[i]) * step;
    }
    if (model.elementCount() > 0)
    {
        state.dissipated += sub0math::f64::Energy::fromSi(dt * workspace.elementDissipationPower);
    }
    state.externalWork += sub0math::f64::Energy::fromSi(work);
    state.tick = tick + 1U;
    return {};
}

std::expected<void, Fault> stepLinearlyImplicitEuler(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                     sub0math::f64::Time h, State& state, Workspace& ws) noexcept
{
    const double dt = h.si();
    if (auto checked = checkStep(model, forcing, tick, dt, state, ws); !checked)
    {
        return checked;
    }
    if (model.constraintCount() > 0)
    {
        return detail::stepConstrained(model, forcing, dt, state, ws, true);
    }
    const std::size_t n = model.dofCount();
    if (auto prepared = prepareStep(model, state, ws, true); !prepared)
    {
        return prepared;
    }
    const bool hasElements = model.elementCount() > 0;
    const auto dense = formStepSystem(model, state, forcing.tau, hasElements ? ws.elementForce.data() : nullptr, ws,
                                      dt, true, true);
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
        if (!std::isfinite(rate) || !std::isfinite(q[i] + dt * rate))
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
        qdOut[i] += b[i];
        const double step = dt * qdOut[i];
        qOut[i] += step;
        work += (forcing.tau.empty() ? 0.0 : forcing.tau[i]) * step;
    }
    if (hasElements)
    {
        state.dissipated += sub0math::f64::Energy::fromSi(dt * ws.elementDissipationPower);
    }
    state.externalWork += sub0math::f64::Energy::fromSi(work);
    state.tick = tick + 1U;
    return {};
}

namespace detail
{

std::expected<void, Fault> prepareStepDynamics(const Model& model, const State& state, Workspace& ws,
                                               bool elementJacobians) noexcept
{
    if (auto prepared = prepareStep(model, state, ws, elementJacobians); !prepared)
    {
        return prepared;
    }
    composeWorldTransforms(model, ws);
    return {};
}

std::expected<bool, Fault> formStepMatrix(const Model& model, const State& state, std::span<const double> tau,
                                          const double* elementForce, Workspace& ws, double dt, bool implicit,
                                          bool allowFallback) noexcept
{
    return formStepSystem(model, state, tau, elementForce, ws, dt, implicit, allowFallback);
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
