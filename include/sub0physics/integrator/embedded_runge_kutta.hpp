#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/explicit_runge_kutta.hpp"
#include "sub0physics/integrator/ode_system.hpp"

namespace sub0physics
{

/**
 * Bogacki-Shampine 3(2): third-order step, second-order error estimate, four stages of which the last is the
 * next step's first (FSAL), so three evaluations per accepted step. Cheap at loose tolerances.
 */
struct BogackiShampine32
{
    static constexpr std::size_t Stages = 4; ///< Stages, the last at the proposal.
    static constexpr int Order = 3;          ///< Order of the propagated solution.
    static constexpr int ErrorOrder = 2;     ///< Order of the embedded estimate; the controller's exponent is 1/(ErrorOrder + 1).
    static constexpr std::array<double, Stages> C{0.0, 0.5, 0.75, 1.0}; ///< Stage times as fractions of the step.
    /** Stage coefficients `a[i][j]`, j < i. */
    static constexpr std::array<std::array<double, Stages>, Stages> A{{
        {0.0, 0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0, 0.0},
        {0.0, 0.75, 0.0, 0.0},
        {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0},
    }};
    /** Weights of the propagated solution (the last row of A: first same as last). */
    static constexpr std::array<double, Stages> B{2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0};
    /** Weights of the error estimate: B minus the embedded solution's weights. */
    static constexpr std::array<double, Stages> E{2.0 / 9.0 - 7.0 / 24.0, 1.0 / 3.0 - 0.25, 4.0 / 9.0 - 1.0 / 3.0,
                                                  -0.125};
};

/**
 * Dormand-Prince 5(4): fifth-order step with a fourth-order error estimate, seven stages with FSAL, so six
 * evaluations per accepted step. The usual choice at tight tolerances.
 */
struct DormandPrince54
{
    static constexpr std::size_t Stages = 7; ///< Stages, the last at the proposal.
    static constexpr int Order = 5;          ///< Order of the propagated solution.
    static constexpr int ErrorOrder = 4;     ///< Order of the embedded estimate.
    static constexpr std::array<double, Stages> C{0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0}; ///< Stage times.
    /** Stage coefficients `a[i][j]`, j < i. */
    static constexpr std::array<std::array<double, Stages>, Stages> A{{
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
        {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
        {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
    }};
    /** Weights of the propagated solution. */
    static constexpr std::array<double, Stages> B{35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
                                                  -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
    /** Weights of the error estimate: B minus the embedded fourth-order weights. */
    static constexpr std::array<double, Stages> E{
        35.0 / 384.0 - 5179.0 / 57600.0,     0.0,
        500.0 / 1113.0 - 7571.0 / 16695.0,   125.0 / 192.0 - 393.0 / 640.0,
        -2187.0 / 6784.0 + 92097.0 / 339200.0, 11.0 / 84.0 - 187.0 / 2100.0,
        -1.0 / 40.0};
};

/**
 * Explicit Runge-Kutta pair with an embedded error estimate, fixed size, never allocating: a stepper for
 * `advanceAdaptively` (step_size_controller.hpp).
 *
 * `attempt` computes a proposal and its local error estimate without touching the caller's state; `accept`
 * commits the proposal. The last stage is evaluated at the proposal, so after an accepted attempt the system
 * observes the committed state, and that derivative is the next step's first stage (first same as last).
 *
 * @tparam N Number of state variables.
 * @tparam Pair The tableau, BogackiShampine32 or DormandPrince54.
 * @note Not thread-safe; one instance per integrated system.
 */
template<std::size_t N, typename Pair>
class EmbeddedRungeKutta
{
public:
    /** Order of the propagated solution. */
    static constexpr int Order = Pair::Order;
    /** Order of the error estimate. */
    static constexpr int ErrorOrder = Pair::ErrorOrder;

    /**
     * Computes a step from `(t, x)` of size `h`, and its error estimate.
     *
     * @tparam S The system type.
     * @param system The system; it may refuse a stage state.
     * @param t Time at the start, s.
     * @param h Step size, s; finite and positive.
     * @param x State at `t`; not modified.
     * @return Nothing, or the system's fault (the proposal is then invalid), or NonFiniteState.
     */
    template<OdeSystem<N> S>
    [[nodiscard]] SUB0PHYSICS_FLATTEN std::expected<void, Fault> attempt(S& system, double t, double h,
                                                                         const StateVector<N>& x)
    {
        if (!cacheValid_ || cachedTime_ != t || cachedState_ != x)
        {
            ++evaluations_;
            if (auto r = system.evaluate(t, x, k_[0]); !r)
            {
                return r;
            }
            cachedTime_ = t;
            cachedState_ = x;
            cacheValid_ = true;
        }
        for (std::size_t s = 1; s < Pair::Stages; ++s)
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                double sum = 0.0;
                for (std::size_t j = 0; j < s; ++j)
                {
                    sum += Pair::A[s][j] * k_[j][i];
                }
                trial_[i] = x[i] + h * sum;
            }
            if (s + 1U == Pair::Stages)
            {
                for (const double value : trial_)
                {
                    if (!std::isfinite(value))
                    {
                        return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
                    }
                }
            }
            ++evaluations_;
            if (auto r = system.evaluate(t + Pair::C[s] * h, trial_, k_[s]); !r)
            {
                return r;
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            double sum = 0.0;
            for (std::size_t j = 0; j < Pair::Stages; ++j)
            {
                sum += Pair::E[j] * k_[j][i];
            }
            error_[i] = h * sum;
        }
        proposalTime_ = t + h;
        return {};
    }

    /** @return The proposed state of the last successful attempt (the last stage's state). */
    [[nodiscard]] const StateVector<N>& proposal() const noexcept
    {
        return trial_;
    }

    /** @return The local error estimate of the last successful attempt, per component. */
    [[nodiscard]] const StateVector<N>& error() const noexcept
    {
        return error_;
    }

    /**
     * Commits the last successful attempt.
     *
     * @param x Receives the proposal.
     */
    void accept(StateVector<N>& x) noexcept
    {
        x = trial_;
        k_[0] = k_[Pair::Stages - 1U];
        cachedTime_ = proposalTime_;
        cachedState_ = x;
        cacheValid_ = true;
    }

    /** @return Derivative evaluations so far (the cost counter). */
    [[nodiscard]] std::uint64_t evaluations() const noexcept
    {
        return evaluations_;
    }

    /** @return Jacobian evaluations so far; always zero for an explicit pair. */
    [[nodiscard]] static constexpr std::uint64_t jacobians() noexcept
    {
        return 0U;
    }

private:
    std::array<StateVector<N>, Pair::Stages> k_{};
    StateVector<N> trial_{};
    StateVector<N> error_{};
    StateVector<N> cachedState_{};
    double cachedTime_ = 0.0;
    double proposalTime_ = 0.0;
    std::uint64_t evaluations_ = 0U;
    bool cacheValid_ = false;
};

} // namespace sub0physics
