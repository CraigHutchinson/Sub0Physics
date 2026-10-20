#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <utility>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/explicit_runge_kutta.hpp"
#include "sub0physics/integrator/ode_system.hpp"

namespace sub0physics
{

/**
 * Linearly implicit Euler for a small fixed-size system with a Jacobian: backward Euler linearised once at the
 * start of the step, never allocating. The fixed-size counterpart of `multibody::LinearlyImplicitEuler`, for
 * systems written as `dx/dt = f(t, x)` (station models).
 *
 * One step solves `(I - h J) dx = h f(t + h, x)` with `J = df/dx (t + h, x)` by Gaussian elimination with partial
 * pivoting, and sets `x += dx`: exact backward Euler when `f` is linear in `x`, first order otherwise, and stable
 * for stiff elements because they enter through their Jacobians (`dForceDPosition`, `dForceDVelocity`). It
 * damps numerically: an energy ledger shows that as a negative residual of order h.
 *
 * The error estimate is `h/2 (f(t + h, x1) - f(t, x0))`, Euler's local truncation error `h^2/2 x''` by a
 * difference, so `advanceAdaptively` can control the step. Both derivatives are evaluated anyway (the first is
 * the previous step's endpoint), so a step costs two evaluations and one Jacobian.
 *
 * **Quadratures.** Slots marked with `setQuadrature` (integrals nothing depends on, such as accumulated work and
 * dissipation) are integrated by the trapezoid rule from the step's two endpoint derivatives, which the step has
 * anyway, instead of the linearised update: second order and consistent with the committed trajectory, so an
 * energy ledger built on them measures the method's own energy error rather than the bookkeeping's lag.
 *
 * `attempt` leaves the caller's state alone; `accept` commits. `stepTick` does both on a fixed tick grid.
 *
 * @tparam N Number of state variables.
 * @note Not thread-safe; one instance per integrated system.
 */
template<std::size_t N>
class LinearlyImplicitEuler
{
public:
    /** Order of the method. */
    static constexpr int Order = 1;
    /** Order of the error estimate's method (the controller's exponent is 1/(ErrorOrder + 1)). */
    static constexpr int ErrorOrder = 1;

    /**
     * Marks a slot as a quadrature: an integral that no derivative depends on, integrated by the trapezoid rule.
     *
     * @param index The slot; less than N.
     * @param quadrature True to integrate it by the trapezoid rule, false for the linearised update (the default).
     */
    void setQuadrature(std::size_t index, bool quadrature = true) noexcept
    {
        quadrature_[index] = quadrature;
    }

    /**
     * Computes a step from `(t, x)` of size `h`, and its error estimate.
     *
     * @tparam S The system type.
     * @param system The system; it may refuse a state.
     * @param t Time at the start, s.
     * @param h Step size, s; finite and positive.
     * @param x State at `t`; not modified.
     * @return Nothing; or the system's fault; NonFiniteState; or InvalidInput when `I - h J` is singular.
     */
    template<OdeSystemWithJacobian<N> S>
    [[nodiscard]] std::expected<void, Fault> attempt(S& system, double t, double h, const StateVector<N>& x)
    {
        if (!cacheValid_ || cachedTime_ != t || cachedState_ != x)
        {
            ++evaluations_;
            if (auto r = system.evaluate(t, x, start_); !r)
            {
                return r;
            }
            cachedTime_ = t;
            cachedState_ = x;
            cacheValid_ = true;
        }
        const double end = t + h;
        ++evaluations_;
        if (auto r = system.evaluate(end, x, rhs_); !r)
        {
            return r;
        }
        ++jacobians_;
        if (auto r = system.jacobian(end, x, matrix_); !r)
        {
            return r;
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t j = 0; j < N; ++j)
            {
                matrix_[i][j] = (i == j ? 1.0 : 0.0) - h * matrix_[i][j];
            }
            rhs_[i] *= h;
        }
        if (!solveInPlace())
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "linearly implicit step matrix is singular"});
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            proposal_[i] = x[i] + rhs_[i];
            if (!std::isfinite(proposal_[i]))
            {
                return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
            }
        }
        ++evaluations_;
        if (auto r = system.evaluate(end, proposal_, endDerivative_); !r)
        {
            return r;
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            error_[i] = 0.5 * h * (endDerivative_[i] - start_[i]);
            if (quadrature_[i])
            {
                proposal_[i] = x[i] + 0.5 * h * (start_[i] + endDerivative_[i]);
            }
        }
        proposalTime_ = end;
        return {};
    }

    /** @return The proposed state of the last successful attempt. */
    [[nodiscard]] const StateVector<N>& proposal() const noexcept
    {
        return proposal_;
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
        x = proposal_;
        start_ = endDerivative_;
        cachedTime_ = proposalTime_;
        cachedState_ = x;
        cacheValid_ = true;
    }

    /**
     * Advances over fixed-step tick `tick`, from `tick * h` to `(tick + 1) * h`; times come from the tick.
     *
     * @tparam S The system type.
     * @param system The system.
     * @param tick Index of the step.
     * @param h Step size.
     * @param state State at `tick * h`; replaced only on success.
     * @return Nothing, or the first fault with `state` unchanged and the system observing it.
     */
    template<OdeSystemWithJacobian<N> S>
    [[nodiscard]] std::expected<void, Fault> stepTick(S& system, std::uint64_t tick, sub0math::f64::Time h,
                                                      StateVector<N>& state)
    {
        const double dt = h.si();
        const double t = static_cast<double>(tick) * dt;
        if (!std::isfinite(dt) || !(dt > 0.0))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "step needs a finite time and positive step"});
        }
        // Step from the tick's own time; the end comes from the next tick, so it is exact and the cache hits.
        auto result = attempt(system, t, static_cast<double>(tick + 1U) * dt - t, state);
        if (!result)
        {
            restore(system, t, state);
            return result;
        }
        accept(state);
        return {};
    }

    /**
     * Re-evaluates `(t, x)` so that the system observes it again, after a failed attempt.
     *
     * @tparam S The system type.
     * @param system The system.
     * @param t Time, s.
     * @param x The state to observe.
     */
    template<OdeSystem<N> S>
    void restore(S& system, double t, const StateVector<N>& x)
    {
        ++evaluations_;
        cacheValid_ = system.evaluate(t, x, start_).has_value();
        cachedTime_ = t;
        cachedState_ = x;
    }

    /** @return Derivative evaluations so far. */
    [[nodiscard]] std::uint64_t evaluations() const noexcept
    {
        return evaluations_;
    }

    /** @return Jacobian evaluations (and factorisations) so far. */
    [[nodiscard]] std::uint64_t jacobians() const noexcept
    {
        return jacobians_;
    }

private:
    // Solves matrix_ * y = rhs_ into rhs_ (Gaussian elimination, partial pivoting); false when a pivot is zero.
    bool solveInPlace() noexcept
    {
        for (std::size_t col = 0; col < N; ++col)
        {
            std::size_t pivot = col;
            for (std::size_t row = col + 1U; row < N; ++row)
            {
                if (std::abs(matrix_[row][col]) > std::abs(matrix_[pivot][col]))
                {
                    pivot = row;
                }
            }
            if (!(matrix_[pivot][col] != 0.0) || !std::isfinite(matrix_[pivot][col]))
            {
                return false;
            }
            if (pivot != col)
            {
                std::swap(matrix_[pivot], matrix_[col]);
                std::swap(rhs_[pivot], rhs_[col]);
            }
            for (std::size_t row = col + 1U; row < N; ++row)
            {
                const double factor = matrix_[row][col] / matrix_[col][col];
                if (factor == 0.0)
                {
                    continue;
                }
                for (std::size_t k = col; k < N; ++k)
                {
                    matrix_[row][k] -= factor * matrix_[col][k];
                }
                rhs_[row] -= factor * rhs_[col];
            }
        }
        for (std::size_t i = N; i-- > 0U;)
        {
            double sum = rhs_[i];
            for (std::size_t k = i + 1U; k < N; ++k)
            {
                sum -= matrix_[i][k] * rhs_[k];
            }
            rhs_[i] = sum / matrix_[i][i];
        }
        return true;
    }

    JacobianMatrix<N> matrix_{};
    std::array<bool, N> quadrature_{};
    StateVector<N> start_{};
    StateVector<N> rhs_{};
    StateVector<N> proposal_{};
    StateVector<N> endDerivative_{};
    StateVector<N> error_{};
    StateVector<N> cachedState_{};
    double cachedTime_ = 0.0;
    double proposalTime_ = 0.0;
    std::uint64_t evaluations_ = 0U;
    std::uint64_t jacobians_ = 0U;
    bool cacheValid_ = false;
};

} // namespace sub0physics
