#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>

#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/soft_constraint.hpp"

namespace sub0physics::multibody
{

SoftConstraintRow softConstraintRow(const SoftConstraintParameters& parameters, double gap, double rate,
                                    double diagApprox, double step) noexcept
{
    // MuJoCo 3.16.0 engine_core_constraint.c, at commit a892750 (docs/WHEEL.md, "M2c status").
    constexpr double MinVal = 1e-15; // mjMINVAL
    constexpr double MinImp = 0.0001; // mjMINIMP
    constexpr double MaxImp = 0.9999; // mjMAXIMP
    // getsolparam (:1777-1780): integrator safety, then the solimp bounds (:1794-1798).
    double solref[2] = {parameters.timeConstant, parameters.dampingRatio};
    solref[0] = std::max(solref[0], 2 * step);
    double solimp[5] = {parameters.impedance[0], parameters.impedance[1], parameters.impedance[2],
                        parameters.impedance[3], parameters.impedance[4]};
    solimp[0] = std::min(MaxImp, std::max(MinImp, solimp[0]));
    solimp[1] = std::min(MaxImp, std::max(MinImp, solimp[1]));
    solimp[2] = std::max(0.0, solimp[2]);
    solimp[3] = std::min(MaxImp, std::max(MinImp, solimp[3]));
    solimp[4] = std::max(1.0, solimp[4]);
    // getimpedance (:1851-1900), margin zero.
    double imp = 0.0;
    if (solimp[0] == solimp[1] || solimp[2] <= MinVal)
    {
        imp = 0.5 * (solimp[0] + solimp[1]);
    }
    else
    {
        double x = gap / solimp[2];
        if (x < 0)
        {
            x = -x;
        }
        if (x >= 1 || x <= 0)
        {
            imp = x >= 1 ? solimp[1] : solimp[0];
        }
        else
        {
            double y = 0.0;
            if (solimp[4] == 1)
            {
                y = x;
            }
            else if (x <= solimp[3])
            {
                const double a = 1 / std::pow(solimp[3], solimp[4] - 1);
                y = a * std::pow(x, solimp[4]);
            }
            else
            {
                const double b = 1 / std::pow(1 - solimp[3], solimp[4] - 1);
                y = 1 - b * std::pow(1 - x, solimp[4]);
            }
            imp = solimp[0] + y * (solimp[1] - solimp[0]);
        }
    }
    // mj_makeImpedance (:1930, :1947, :1956): R, then K and B in the standard format.
    SoftConstraintRow row;
    row.impedance = imp;
    row.regularisation = std::max(MinVal, (1 - imp) * diagApprox / imp);
    row.stiffness = 1 / std::max(MinVal, solimp[1] * solimp[1] * solref[0] * solref[0] * solref[1] * solref[1]);
    row.damping = 2 / std::max(MinVal, solimp[1] * solref[0]);
    // mj_referenceConstraint (:3194-3196): aref = -B v - K I (pos - margin).
    row.referenceAcceleration = -row.damping * rate - row.stiffness * row.impedance * gap;
    return row;
}

SoftConstraintSolver::SoftConstraintSolver(const Model& model, std::size_t capacity)
    : capacity_(capacity), dof_(model.dofCount()), solved_(capacity * model.dofCount(), 0.0),
      delassus_(capacity * capacity, 0.0), factor_(capacity * capacity, 0.0), rhs_(capacity, 0.0),
      slack_(capacity, 0.0), work_(capacity, 0.0), subset_(capacity, 0), active_(capacity, 0)
{
}

std::expected<SoftConstraintSolution, Fault>
SoftConstraintSolver::solve(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                         std::span<const double> reference, std::span<const double> regularisation,
                         std::span<const double> freeAcceleration, std::span<const SoftRowKind> kinds,
                         std::span<double> force, std::span<double> generalized) noexcept
{
    const std::size_t n = dof_;
    const std::size_t k = reference.size();
    if (model.dofCount() != n || k > capacity_ || jacobian.size() != k * n || regularisation.size() != k ||
        freeAcceleration.size() != n || force.size() != k || generalized.size() != n ||
        (!kinds.empty() && kinds.size() != k))
    {
        return invalidInput("soft contact solve needs consistent sizes within the solver's capacity");
    }
    std::ranges::fill(generalized, 0.0);
    SoftConstraintSolution solution;
    if (k == 0)
    {
        return solution;
    }
    // Y = M^-1 J^T, A = J Y + R, q = J a0 - aref.
    double scale = 0.0;
    for (std::size_t r = 0; r < k; ++r)
    {
        if (!(regularisation[r] > 0.0))
        {
            return invalidInput("soft contact rows need a positive regularisation");
        }
        const std::span<double> y(solved_.data() + r * n, n);
        std::copy_n(jacobian.data() + r * n, n, y.begin());
        if (auto solved = solveMassMatrix(model, workspace, y); !solved)
        {
            return std::unexpected(solved.error());
        }
        double q = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            q += jacobian[r * n + i] * freeAcceleration[i];
        }
        rhs_[r] = q - reference[r];
        scale = std::max(scale, std::abs(rhs_[r]));
    }
    for (std::size_t r = 0; r < k; ++r)
    {
        for (std::size_t c = 0; c <= r; ++c)
        {
            double sum = 0.0;
            const double* jr = jacobian.data() + r * n;
            const double* yc = solved_.data() + c * n;
            for (std::size_t i = 0; i < n; ++i)
            {
                sum += jr[i] * yc[i];
            }
            delassus_[r * k + c] = sum;
            delassus_[c * k + r] = sum;
        }
        delassus_[r * k + r] += regularisation[r];
    }
    // Murty's least-index principal pivoting, from the rows the free motion violates.
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    const auto equality = [&](std::size_t r) { return !kinds.empty() && kinds[r] == SoftRowKind::Equality; };
    for (std::size_t r = 0; r < k; ++r)
    {
        active_[r] = equality(r) || rhs_[r] < 0.0 ? 1U : 0U;
    }
    const std::uint32_t maxPivots = static_cast<std::uint32_t>(16U * k + 64U);
    for (;;)
    {
        // Solve A_SS f_S = -q_S by Cholesky.
        std::size_t m = 0;
        for (std::size_t r = 0; r < k; ++r)
        {
            if (active_[r] != 0U)
            {
                subset_[m++] = r;
            }
        }
        for (std::size_t a = 0; a < m; ++a)
        {
            for (std::size_t c = 0; c <= a; ++c)
            {
                double sum = delassus_[subset_[a] * k + subset_[c]];
                for (std::size_t l = 0; l < c; ++l)
                {
                    sum -= factor_[a * k + l] * factor_[c * k + l];
                }
                if (a == c)
                {
                    if (!(sum > 0.0))
                    {
                        return invalidInput("soft contact matrix is not positive definite");
                    }
                    factor_[a * k + a] = std::sqrt(sum);
                }
                else
                {
                    factor_[a * k + c] = sum / factor_[c * k + c];
                }
            }
        }
        for (std::size_t a = 0; a < m; ++a)
        {
            double sum = -rhs_[subset_[a]];
            for (std::size_t l = 0; l < a; ++l)
            {
                sum -= factor_[a * k + l] * work_[l];
            }
            work_[a] = sum / factor_[a * k + a];
        }
        for (std::size_t a = m; a-- > 0;)
        {
            double sum = work_[a];
            for (std::size_t l = a + 1; l < m; ++l)
            {
                sum -= factor_[l * k + a] * work_[l];
            }
            work_[a] = sum / factor_[a * k + a];
        }
        std::ranges::fill(force, 0.0);
        for (std::size_t a = 0; a < m; ++a)
        {
            force[subset_[a]] = work_[a];
        }
        // w = A f + q, and the first row that breaks complementarity.
        std::size_t flip = k;
        for (std::size_t r = 0; r < k && flip == k; ++r)
        {
            if (active_[r] != 0U)
            {
                if (force[r] < 0.0 && !equality(r))
                {
                    flip = r;
                }
                continue;
            }
            double w = rhs_[r];
            for (std::size_t c = 0; c < k; ++c)
            {
                w += delassus_[r * k + c] * force[c];
            }
            slack_[r] = w;
            if (w < -tolerance)
            {
                flip = r;
            }
        }
        if (flip == k)
        {
            solution.activeRows = m;
            break;
        }
        if (++solution.pivots > maxPivots)
        {
            return invalidInput("soft contact pivoting did not terminate");
        }
        active_[flip] = active_[flip] != 0U ? 0U : 1U;
    }
    for (std::size_t r = 0; r < k; ++r)
    {
        const double f = force[r];
        if (!std::isfinite(f))
        {
            return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite soft contact force"});
        }
        if (f == 0.0)
        {
            continue;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            generalized[i] += jacobian[r * n + i] * f;
        }
    }
    return solution;
}

} // namespace sub0physics::multibody
