#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/explicit_runge_kutta.hpp"
#include "sub0physics/integrator/ode_system.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/stepper.hpp"
#include "sub0physics/multibody/workspace.hpp"

namespace sub0physics::multibody
{

/**
 * A model as a first-order system for the dynamic-size integrators: `x = [q; qd]` (State::values) and
 * `dx/dt = [qd; qdd]`, with `qdd` from forwardDynamics under constant applied forces.
 *
 * It refers to a model and a workspace it does not own; both must outlive it, or be rebound with
 * `bind`. Its scratch is allocated once at construction, so evaluation never allocates.
 */
class ModelOde
{
public:
    /**
     * Builds the adapter for a model, with no applied forces.
     *
     * @param model The model; its dofCount fixes every size.
     * @param workspace Workspace built for the model; receives the dynamics of each evaluation.
     */
    ModelOde(const Model& model, Workspace& workspace)
        : model_(&model), workspace_(&workspace), scratch_(State::atRest(model)), tau_(model.dofCount(), 0.0)
    {
    }

    /**
     * Rebinds the model and workspace (same dofCount as the one it was built for).
     *
     * @param model The model.
     * @param workspace Workspace built for the model.
     */
    void bind(const Model& model, Workspace& workspace) noexcept
    {
        model_ = &model;
        workspace_ = &workspace;
    }

    /**
     * Sets the applied generalized forces, held constant until changed.
     *
     * @param tau dofCount entries, or empty for none.
     * @return Whether the forces changed, or InvalidInput for a wrong size.
     */
    [[nodiscard]] std::expected<bool, Fault> setForcing(std::span<const double> tau) noexcept
    {
        if (!tau.empty() && tau.size() != tau_.size())
        {
            return invalidInput("forcing needs tau of dofCount entries, or none");
        }
        bool changed = false;
        for (std::size_t i = 0; i < tau_.size(); ++i)
        {
            const double value = tau.empty() ? 0.0 : tau[i];
            changed = changed || value != tau_[i];
            tau_[i] = value;
        }
        return changed;
    }

    /**
     * Evaluates `[qd; qdd]` at `[q; qd]`.
     *
     * @param t Time [s]; the forcing is constant, so unused.
     * @param x `[q; qd]`, 2 dofCount entries.
     * @param dxdt Receives `[qd; qdd]`.
     * @return Nothing, or the forward-dynamics fault; InvalidInput for wrong sizes.
     */
    [[nodiscard]] std::expected<void, Fault> evaluate(double t, std::span<const double> x,
                                                      std::span<double> dxdt) noexcept
    {
        (void)t;
        const std::size_t n = tau_.size();
        if (x.size() != 2 * n || dxdt.size() != 2 * n)
        {
            return invalidInput("model ODE needs [q; qd] of 2 dofCount entries");
        }
        std::copy(x.begin(), x.end(), scratch_.values.begin());
        if (auto accelerated = forwardDynamics(*model_, scratch_, tau_, *workspace_); !accelerated)
        {
            return accelerated;
        }
        std::copy(x.begin() + static_cast<std::ptrdiff_t>(n), x.end(), dxdt.begin());
        std::copy(workspace_->qdd.begin(), workspace_->qdd.end(), dxdt.begin() + static_cast<std::ptrdiff_t>(n));
        return {};
    }

private:
    const Model* model_;
    Workspace* workspace_;
    State scratch_;
    std::vector<double> tau_;
};

static_assert(OdeSystem<ModelOde, std::dynamic_extent>);

/**
 * An explicit Runge-Kutta Stepper for a model: the non-stiff reference run (RK4: energy conserved to
 * the method's order, unlike SemiImplicitEuler).
 *
 * It integrates `[q; qd]` with the dynamic-size ExplicitRungeKutta through ModelOde, reusing the
 * endpoint derivative of one tick as the first stage of the next while the forcing is unchanged
 * (4 forward-dynamics evaluations per RK4 tick). It is bound to the model it was built for, like a
 * Workspace, and allocates only at construction.
 *
 * @tparam Method The tableau.
 */
template<RungeKuttaMethod Method>
class RungeKuttaStepper
{
public:
    /** Integrator order. */
    static constexpr int Order = ExplicitRungeKutta<std::dynamic_extent, Method>::Order;

    /**
     * Builds a stepper for a model.
     *
     * @param model The model.
     * @param workspace A workspace built for the model (the one passed to `step` is used when stepping).
     */
    RungeKuttaStepper(const Model& model, Workspace& workspace)
        : modelId_(model.id()), integrator_(2 * model.dofCount()), ode_(model, workspace)
    {
    }

    /**
     * Advances one tick; see Stepper.
     *
     * @param model The model it was built for.
     * @param forcing Applied generalized forces, constant over the step.
     * @param tick The tick the state is at; must equal `state.tick`.
     * @param h Step size; finite and positive.
     * @param state State at tick; replaced by the state at tick + 1 only on success.
     * @param workspace Workspace built for the model; holds the dynamics at the new state on success.
     * @return Nothing, or the fault: InvalidInput (another model, sizes, tick, step size) or a dynamics
     *         fault, with the state unchanged.
     */
    [[nodiscard]] std::expected<void, Fault> step(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                  sub0math::f64::Time h, State& state, Workspace& workspace) noexcept
    {
        if (model.id() != modelId_ || !workspace.isBuiltFor(model) || !state.isSizedFor(model))
        {
            return invalidInput("stepper, workspace and state must be built for the model");
        }
        if (tick != state.tick)
        {
            return invalidInput("tick differs from the state's tick");
        }
        const auto changed = ode_.setForcing(forcing.tau);
        if (!changed)
        {
            return std::unexpected(changed.error());
        }
        if (*changed)
        {
            integrator_.invalidate();
        }
        ode_.bind(model, workspace);
        if (auto stepped = integrator_.stepTick(ode_, tick, h, state.values); !stepped)
        {
            return stepped;
        }
        state.tick = tick + 1U;
        return {};
    }

private:
    std::uint64_t modelId_;
    ExplicitRungeKutta<std::dynamic_extent, Method> integrator_;
    ModelOde ode_;
};

/** Classic RK4 on a model: the reference stepper. */
using Rk4Stepper = RungeKuttaStepper<RungeKuttaMethod::ClassicRk4>;

/** Explicit midpoint on a model. */
using MidpointStepper = RungeKuttaStepper<RungeKuttaMethod::Midpoint>;

static_assert(Stepper<Rk4Stepper>);
static_assert(Stepper<MidpointStepper>);

} // namespace sub0physics::multibody
