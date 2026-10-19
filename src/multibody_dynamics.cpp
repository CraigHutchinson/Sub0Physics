#include <cmath>
#include <cstddef>
#include <expected>
#include <span>
#include <string_view>

#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>

#include "multibody_joint.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/symplectic_euler.hpp"

namespace sub0physics::multibody
{

namespace
{

std::unexpected<Fault> invalid(std::string_view detail) noexcept
{
    return std::unexpected(Fault{FaultKind::InvalidInput, detail});
}

// Recursive Newton-Euler over transforms and velocities already in the workspace. With accelerations,
// tau = M qdd + c; without, tau = c.
template<bool WithAcceleration>
void recursiveNewtonEuler(const Model& model, const State& state, const double* qdd, double* tau,
                          Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const auto subspace = model.motionSubspaces();
    const auto inertias = model.inertias();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const std::size_t n = parents.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::size_t p = parents[i];
        const sub0math::SpatialMotiond& v = ws.bodyVelocity[i];
        const sub0math::SpatialMotiond& parentAcceleration =
            p == NoParent ? model.baseAcceleration() : ws.bodyAcceleration[p];
        sub0math::SpatialMotiond a = ws.parentToBody[i].apply(parentAcceleration) + cross(v, subspace[i] * state.qd[i]);
        if constexpr (WithAcceleration)
        {
            a += subspace[i] * qdd[i];
        }
        ws.bodyAcceleration[i] = a;
        ws.bodyForce[i] = inertias[i] * a + cross(v, inertias[i] * v);
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

} // namespace

void computeMassMatrix(const Model& model, const State& state, Workspace& workspace) noexcept
{
    detail::updateTransforms<false>(model, state, workspace);
    compositeRigidBody(model, workspace);
}

void computeBiasForces(const Model& model, const State& state, Workspace& workspace) noexcept
{
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<false>(model, state, nullptr, workspace.bias.data(), workspace);
}

std::expected<void, Fault> inverseDynamics(const Model& model, const State& state, std::span<const double> qdd,
                                           std::span<double> tau, Workspace& workspace) noexcept
{
    if (qdd.size() != model.dofCount() || tau.size() != model.dofCount())
    {
        return invalid("inverse dynamics needs qdd and tau of dofCount entries");
    }
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<true>(model, state, qdd.data(), tau.data(), workspace);
    return {};
}

std::expected<void, Fault> factorMassMatrix(const Model& model, Workspace& workspace) noexcept
{
    // Featherstone, RBDA section 6.5: eliminate from the leaves; only ancestor entries are touched.
    const auto parents = model.parents();
    const std::size_t n = parents.size();
    double* h = workspace.massFactor.data();
    const double* m = workspace.massMatrix.data();
    for (std::size_t k = 0; k < n * n; ++k)
    {
        h[k] = m[k];
    }
    for (std::size_t k = n; k-- > 0;)
    {
        const double pivot = h[k * n + k];
        if (!(pivot > 0.0) || !std::isfinite(pivot))
        {
            return invalid("mass matrix is not positive definite");
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
    return {};
}

void solveMassMatrix(const Model& model, const Workspace& workspace, std::span<double> x) noexcept
{
    const auto parents = model.parents();
    const std::size_t n = parents.size();
    const double* h = workspace.massFactor.data();
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

std::expected<void, Fault> forwardDynamics(const Model& model, const State& state, std::span<const double> tau,
                                           Workspace& workspace) noexcept
{
    const std::size_t n = model.dofCount();
    if (!tau.empty() && tau.size() != n)
    {
        return invalid("forward dynamics needs tau of dofCount entries, or none");
    }
    detail::updateTransforms<true>(model, state, workspace);
    recursiveNewtonEuler<false>(model, state, nullptr, workspace.bias.data(), workspace);
    compositeRigidBody(model, workspace);
    if (auto factored = factorMassMatrix(model, workspace); !factored)
    {
        return factored;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        workspace.qdd[i] = (tau.empty() ? 0.0 : tau[i]) - workspace.bias[i];
    }
    solveMassMatrix(model, workspace, workspace.qdd);
    return {};
}

std::expected<void, Fault> stepSymplecticEuler(const Model& model, sub0math::f64::Time h, std::span<const double> tau,
                                               State& state, Workspace& workspace) noexcept
{
    const double dt = h.si();
    if (!std::isfinite(dt) || !(dt > 0.0))
    {
        return invalid("step needs a finite positive step size");
    }
    if (state.q.size() != model.dofCount() || state.qd.size() != model.dofCount())
    {
        return invalid("state is not sized for the model");
    }
    if (auto accelerated = forwardDynamics(model, state, tau, workspace); !accelerated)
    {
        return accelerated;
    }
    const std::size_t n = model.dofCount();
    // Validate the whole endpoint before committing any of it.
    for (std::size_t i = 0; i < n; ++i)
    {
        const double qd = state.qd[i] + dt * workspace.qdd[i];
        if (!std::isfinite(qd) || !std::isfinite(state.q[i] + dt * qd))
        {
            return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
        }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        state.qd[i] += dt * workspace.qdd[i];
        state.q[i] += dt * state.qd[i];
    }
    ++state.tick;
    return {};
}

} // namespace sub0physics::multibody
