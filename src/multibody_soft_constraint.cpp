#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include <sub0math/number/det_math.hpp>

#include "multibody_mass.hpp"
#include "sub0physics/multibody/implicit_euler.hpp"
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
                const double a = 1 / sub0math::det::pow(solimp[3], solimp[4] - 1);
                y = a * sub0math::det::pow(x, solimp[4]);
            }
            else
            {
                const double b = 1 / sub0math::det::pow(1 - solimp[3], solimp[4] - 1);
                y = 1 - b * sub0math::det::pow(1 - x, solimp[4]);
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

SoftConstraintRow softConstraintTangentRow(const SoftConstraintRow& normal, double rate, double impratio) noexcept
{
    constexpr double MinVal = 1e-15; // mjMINVAL
    // mj_makeImpedance (:1940-1943, :1954-1957 with solreffriction zero, :2012), mj_referenceConstraint (:3195).
    SoftConstraintRow row;
    row.impedance = normal.impedance;
    row.regularisation = normal.regularisation / std::max(MinVal, impratio);
    row.damping = normal.damping;
    row.referenceAcceleration = -row.damping * rate;
    return row;
}

SoftConstraintRow softConstraintFrictionLossRow(const SoftConstraintParameters& parameters, double rate,
                                                double diagApprox, double step) noexcept
{
    // The law at zero position (getimpedance at x = 0 gives solimp[0]); friction rows have K = 0 (mj_makeImpedance
    // :1941-1943), so aref is the damping term alone (:3194-3196).
    SoftConstraintRow row = softConstraintRow(parameters, 0.0, rate, diagApprox, step);
    row.stiffness = 0.0;
    row.referenceAcceleration = -row.damping * rate;
    return row;
}

std::expected<double, Fault> dofInverseWeight(const Model& model, const Workspace& workspace, std::size_t dof)
{
    if (!workspace.isBuiltFor(model) || dof >= model.dofCount())
    {
        return invalidInput("a dof's inverse weight needs a workspace for the model and one of its coordinates");
    }
    // engine_setconst.c :1288-1300: A = J M^-1 J^T with J the dof's unit row; one component, so A[0].
    std::vector<double> x(model.dofCount(), 0.0);
    x[dof] = 1.0;
    if (auto solved = solveMassMatrix(model, workspace, x); !solved)
    {
        return std::unexpected(solved.error());
    }
    return x[dof];
}

SoftConstraintRow exactVelocityRow(double rate, double diagApprox, double step) noexcept
{
    constexpr double Compliance = 1e-14; // relative to the row's scale: rounding, only so dependent rows factor
    SoftConstraintRow row;
    row.regularisation = Compliance * diagApprox;
    row.referenceAcceleration = -rate / step;
    return row;
}

namespace
{

// A row's support: the tree's coordinates and the edge's four node coordinates.
constexpr std::size_t RingSupportSize = 4;

// A/B uniform pass 3: the batched solve's ring elimination skips nodes no row touches (bit for bit).
constexpr bool SkipZeroRingRows = true;

} // namespace

SoftConstraintSolver::SoftConstraintSolver(const Model& model, std::size_t capacity)
    : capacity_(capacity), dof_(model.dofCount()), solved_(capacity * model.dofCount(), 0.0),
      index_(capacity * (model.treeDofCount() + RingSupportSize), 0),
      structured_(capacity * (2 * model.treeDofCount() + 2 * RingSupportSize), 0.0), sweep_(6 * capacity, 0.0),
      key_(2 * capacity, 0), touched_(model.dofCount(), 0),
      delassus_(capacity * capacity, 0.0),
      factor_(capacity * capacity, 0.0), rhs_(capacity, 0.0), slack_(capacity, 0.0), work_(capacity, 0.0),
      subset_(capacity, 0), prior_(capacity, 0), flagged_(capacity, 0), active_(capacity, 0), bound_(capacity, 0)
{
}

void SoftConstraintSolver::formUniform(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                                       std::span<const double> reference, std::span<const double> freeAcceleration,
                                       std::span<const SoftRowSupport> supports, SoftRowSystem system,
                                       std::size_t k) noexcept
{
    const std::size_t n = dof_;
    // Y = M^-1 J^T, A = J Y + R, q = J a0 - aref.
    if (supports.empty())
    {
        for (std::size_t r = 0; r < k; ++r)
        {
            const std::span<double> y(solved_.data() + r * n, n);
            std::copy_n(jacobian.data() + r * n, n, y.begin());
            // Sizes checked by the caller.
            if (system == SoftRowSystem::PreparedStep)
            {
                (void)solvePreparedStepMatrix(model, workspace, y);
            }
            else
            {
                (void)solveMassMatrix(model, workspace, y);
            }
        }
        for (std::size_t r = 0; r < k; ++r)
        {
            double q = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                q += jacobian[r * n + i] * freeAcceleration[i];
            }
            rhs_[r] = q - reference[r];
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
        }
        return;
    }
    // PERF (A/B uniform pass 2): all rows solved in one sweep, Y laid out coordinate-major (Y[i k + r]) so the rows
    // share SIMD lanes; each column bit for bit its own solve.
    // PERF (A/B uniform pass 4): only the tree and the touched nodes are laid out; no zero scan.
    const std::size_t nt = model.treeDofCount();
    std::ranges::fill(touched_, 0U);
    for (std::size_t r = 0; r < k; ++r)
    {
        const double* jr = jacobian.data() + r * n;
        for (std::size_t i = 0; i < nt; ++i)
        {
            solved_[i * k + r] = jr[i];
        }
        const std::size_t* index = index_.data() + r * (nt + RingSupportSize);
        for (const std::size_t at : {index[nt], index[nt + 2]})
        {
            if (touched_[at] == 0U)
            {
                touched_[at] = 1U;
                std::fill_n(solved_.data() + at * k, 2 * k, 0.0);
            }
        }
    }
    for (std::size_t r = 0; r < k; ++r)
    {
        const double* jr = jacobian.data() + r * n;
        const std::size_t* index = index_.data() + r * (nt + RingSupportSize);
        for (std::size_t s = nt; s < nt + RingSupportSize; ++s)
        {
            solved_[index[s] * k + r] = jr[index[s]];
        }
    }
    detail::solveMassMatrixColumns(model, workspace, solved_.data(), k, sweep_.data(), SkipZeroRingRows,
                                   touched_.data(),
                                   system == SoftRowSystem::PreparedStep ? detail::RingSystem::Step : detail::RingSystem::Mass);
    // PERF (A/B uniform pass 1): with supports, the products read only each row's support, in ascending coordinate
    // order: the terms skipped are exact zeros, so the sums are bit for bit those over every coordinate.
    const std::size_t width = model.treeDofCount() + RingSupportSize;
    for (std::size_t r = 0; r < k; ++r)
    {
        const double* jr = jacobian.data() + r * n;
        const std::size_t* index = index_.data() + r * width;
        double q = 0.0;
        for (std::size_t s = 0; s < width; ++s)
        {
            q += jr[index[s]] * freeAcceleration[index[s]];
        }
        rhs_[r] = q - reference[r];
        for (std::size_t c = 0; c <= r; ++c)
        {
            double sum = 0.0;
            for (std::size_t s = 0; s < width; ++s)
            {
                sum += jr[index[s]] * solved_[index[s] * k + c];
            }
            delassus_[r * k + c] = sum;
            delassus_[c * k + r] = sum;
        }
    }
}

void SoftConstraintSolver::formCondensed(const Model& model, const Workspace& workspace,
                                         std::span<const double> jacobian, std::span<const double> reference,
                                         std::span<const double> freeAcceleration,
                                         std::span<const SoftRowSupport> supports, SoftRowSystem system,
                                         std::size_t k) noexcept
{
    const bool step = system == SoftRowSystem::PreparedStep;
    const double* factor = step ? workspace.implicitMatrix.data() : workspace.massFactor.data();
    const std::size_t nt = model.treeDofCount();
    const std::size_t stride = 2 * nt + 2 * RingSupportSize;
    const auto coordinates = model.coordinates();
    // Per row: w = L^-T t (tree), D^-1 w, the node values j_w and y = A^-1 j_w at node0 then node1.
    for (std::size_t r = 0; r < k; ++r)
    {
        const SoftRowSupport& sup = supports[r];
        const DeformableRing& ring = model.rings()[sup.ring];
        const RingWorkspace& rw = workspace.rings[sup.ring];
        const RingCondensation& cond = step ? rw.step : rw.mass; // block-diagonal: always for the mass matrix, checked for the step
        const double* jr = jacobian.data() + r * dof_;
        double* row = structured_.data() + r * stride;
        double* t = row;
        double* v = row + nt;
        double* jw = row + 2 * nt;
        double* y = jw + RingSupportSize;
        double g0 = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        const std::array<std::size_t, 2> nodes{sup.node0, sup.node1};
        for (std::size_t a = 0; a < 2; ++a)
        {
            const std::size_t node = nodes[a];
            const std::size_t at = ring.firstCoordinate + PlanarRingComponents * node;
            const BlockColumns<PlanarRingComponents, 1> jn{jr[at], jr[at + 1]};
            const auto yn = sub0physics::detail::applyBlock<PlanarRingComponents, 1>(cond.inversePivot[node], jn);
            jw[2 * a] = jn[0];
            jw[2 * a + 1] = jn[1];
            y[2 * a] = yn[0];
            y[2 * a + 1] = yn[1];
            const BlockColumns<PlanarRingComponents, 3>& b = cond.border[node];
            g0 += b[0] * yn[0] + b[3] * yn[1];
            g1 += b[1] * yn[0] + b[4] * yn[1];
            g2 += b[2] * yn[0] + b[5] * yn[1];
        }
        std::copy_n(jr, nt, t);
        for (std::size_t c = 0; c < ring.hostChain.size(); ++c)
        {
            const PlanarMotion& jh = rw.hostJacobian[c];
            t[coordinates[ring.hostChain[c]]] -= jh.w * g0 + jh.x * g1 + jh.y * g2;
        }
        // PERF (A/B structured pass 1): half the tree solve. With S = L^T D L, t_r^T S^-1 t_c = w_r^T D^-1 w_c for
        // w = L^-T t (one sweep), so each row keeps w and D^-1 w and a pair costs a dot of tree length.
        detail::sweepMassTree(model, factor, t);
        for (std::size_t i = 0; i < nt; ++i)
        {
            v[i] = t[i] / factor[i * nt + i];
        }
        // PERF (A/B structured pass 2): J a0 straight from the tree's coordinates and the two nodes, no index lists.
        double q = 0.0;
        for (std::size_t i = 0; i < nt; ++i)
        {
            q += jr[i] * freeAcceleration[i];
        }
        for (std::size_t a = 0; a < 2; ++a)
        {
            const std::size_t at = ring.firstCoordinate + PlanarRingComponents * nodes[a];
            q += jr[at] * freeAcceleration[at] + jr[at + 1] * freeAcceleration[at + 1];
        }
        rhs_[r] = q - reference[r];
        key_[2 * r] = (sup.ring << 32U) | sup.node0;
        key_[2 * r + 1] = (sup.ring << 32U) | sup.node1;
    }
    for (std::size_t r = 0; r < k; ++r)
    {
        const double* rr = structured_.data() + r * stride;
        for (std::size_t c = 0; c <= r; ++c)
        {
            const double* rc = structured_.data() + c * stride;
            double sum = 0.0;
            for (std::size_t i = 0; i < nt; ++i)
            {
                sum += rr[i] * rc[nt + i];
            }
            // Shared nodes: j_r^T A^-1 j_c = j_r . y_c. PERF (A/B structured pass 3): one integer key per node
            // (ring, node), compared directly.
            for (std::size_t a = 0; a < 2; ++a)
            {
                for (std::size_t b = 0; b < 2; ++b)
                {
                    if (key_[2 * r + a] == key_[2 * c + b])
                    {
                        sum += rr[2 * nt + 2 * a] * rc[2 * nt + RingSupportSize + 2 * b] +
                               rr[2 * nt + 2 * a + 1] * rc[2 * nt + RingSupportSize + 2 * b + 1];
                    }
                }
            }
            delassus_[r * k + c] = sum;
            delassus_[c * k + r] = sum;
        }
    }
}

std::expected<void, Fault> SoftConstraintSolver::checkRows(const Model& model, const Workspace& workspace,
                                                           std::span<const double> jacobian, std::size_t k,
                                                           std::span<const double> regularisation,
                                                           std::span<const double> freeAcceleration,
                                                           std::span<const double> force,
                                                           std::span<const double> generalized,
                                                           std::span<const SoftRowSupport> supports, DelassusForm form,
                                                           SoftRowSystem system) const noexcept
{
    const std::size_t n = dof_;
    if (model.dofCount() != n || k > capacity_ || jacobian.size() != k * n || regularisation.size() != k ||
        freeAcceleration.size() != n || force.size() != k || generalized.size() != n ||
        (!supports.empty() && supports.size() != k) || !workspace.isBuiltFor(model))
    {
        return invalidInput("soft contact solve needs consistent sizes within the solver's capacity");
    }
    for (const SoftRowSupport& sup : supports)
    {
        if (sup.ring >= model.rings().size() || sup.node0 >= model.rings()[sup.ring].nodeCount ||
            sup.node1 >= model.rings()[sup.ring].nodeCount || sup.node0 == sup.node1)
        {
            return invalidInput("soft contact row support names two distinct nodes of a ring of the model");
        }
    }
    if (form == DelassusForm::Condensed && supports.empty() && k > 0)
    {
        return invalidInput("the condensed Delassus form needs each row's ring support");
    }
    if (system == SoftRowSystem::PreparedStep && !workspace.stepFactored)
    {
        return invalidInput("a solve with the step's matrix needs factorPreparedStepMatrix");
    }
    return {};
}

std::expected<double, Fault> SoftConstraintSolver::formRows(const Model& model, const Workspace& workspace,
                                                            std::span<const double> jacobian,
                                                            std::span<const double> reference,
                                                            std::span<const double> regularisation,
                                                            std::span<const double> freeAcceleration,
                                                            std::span<const SoftRowSupport> supports, DelassusForm form,
                                                            SoftRowSystem system, std::size_t k) noexcept
{
    for (std::size_t r = 0; r < k; ++r)
    {
        if (!(regularisation[r] > 0.0))
        {
            return invalidInput("soft contact rows need a positive regularisation");
        }
    }
    // The structure the supports' short forms read: always for the mass matrix, for a step when nothing coupled the nodes.
    const auto structured = [&] {
        if (system != SoftRowSystem::PreparedStep)
        {
            return true;
        }
        return !workspace.stepFactorDense &&
               std::ranges::all_of(workspace.rings, [](const RingWorkspace& rw) { return rw.step.blockDiagonal; });
    };
    if (form == DelassusForm::Condensed && !structured())
    {
        return invalidInput("the condensed Delassus form on a step needs a block-diagonal ring block and a tree factor");
    }
    if (!supports.empty() && !structured())
    {
        supports = {}; // the same numbers through the full solves
    }
    if (!supports.empty() && form == DelassusForm::MassSolvePerRow)
    {
        // Each row's support in ascending coordinate order: the tree's coordinates, then the edge's two nodes.
        const std::size_t nt = model.treeDofCount();
        const std::size_t width = nt + RingSupportSize;
        for (std::size_t r = 0; r < k; ++r)
        {
            std::size_t* index = index_.data() + r * width;
            for (std::size_t i = 0; i < nt; ++i)
            {
                index[i] = i;
            }
            const DeformableRing& ring = model.rings()[supports[r].ring];
            const std::size_t lo = std::min(supports[r].node0, supports[r].node1);
            const std::size_t hi = std::max(supports[r].node0, supports[r].node1);
            index[nt] = ring.firstCoordinate + PlanarRingComponents * lo;
            index[nt + 1] = index[nt] + 1;
            index[nt + 2] = ring.firstCoordinate + PlanarRingComponents * hi;
            index[nt + 3] = index[nt + 2] + 1;
        }
    }
    if (form == DelassusForm::Condensed)
    {
        formCondensed(model, workspace, jacobian, reference, freeAcceleration, supports, system, k);
    }
    else
    {
        formUniform(model, workspace, jacobian, reference, freeAcceleration, supports, system, k);
    }
    double scale = 0.0;
    for (std::size_t r = 0; r < k; ++r)
    {
        scale = std::max(scale, std::abs(rhs_[r]));
        delassus_[r * k + r] += regularisation[r];
    }
    return scale;
}

bool SoftConstraintSolver::factorSubset(std::size_t m, std::size_t k, std::size_t& priorCount) noexcept
{
    // PERF: factor rows before the first changed row kept across pivots.
    std::size_t kept = 0;
    while (kept < m && kept < priorCount && prior_[kept] == subset_[kept])
    {
        ++kept;
    }
    std::copy_n(subset_.begin(), m, prior_.begin());
    priorCount = m;
    for (std::size_t a = kept; a < m; ++a)
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
                    return false;
                }
                factor_[a * k + a] = std::sqrt(sum);
            }
            else
            {
                factor_[a * k + c] = sum / factor_[c * k + c];
            }
        }
    }
    return true;
}

void SoftConstraintSolver::solveSubset(std::size_t m, std::size_t k) noexcept
{
    for (std::size_t a = 0; a < m; ++a)
    {
        double sum = work_[a];
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
}

std::expected<void, Fault> SoftConstraintSolver::mapForces(const Model& model, std::span<const double> jacobian,
                                                           std::span<const double> force, std::span<double> generalized,
                                                           std::span<const SoftRowSupport> supports,
                                                           std::size_t k) noexcept
{
    const std::size_t n = dof_;
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
        if (supports.empty())
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                generalized[i] += jacobian[r * n + i] * f;
            }
            continue;
        }
        // Over the support only: the other entries would add exact zeros.
        const std::size_t nt = model.treeDofCount();
        for (std::size_t i = 0; i < nt; ++i)
        {
            generalized[i] += jacobian[r * n + i] * f;
        }
        const std::size_t first = model.rings()[supports[r].ring].firstCoordinate;
        const std::size_t lo = first + PlanarRingComponents * std::min(supports[r].node0, supports[r].node1);
        const std::size_t hi = first + PlanarRingComponents * std::max(supports[r].node0, supports[r].node1);
        for (const std::size_t at : {lo, lo + 1, hi, hi + 1})
        {
            generalized[at] += jacobian[r * n + at] * f;
        }
    }
    return {};
}

std::expected<SoftConstraintSolution, Fault>
SoftConstraintSolver::solve(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                            std::span<const double> reference, std::span<const double> regularisation,
                            std::span<const double> freeAcceleration, std::span<const SoftRowKind> kinds,
                            std::span<double> force, std::span<double> generalized,
                            std::span<const SoftRowSupport> supports, DelassusForm form,
                            std::span<const double> friction, SoftRowSystem system) noexcept
{
    const std::size_t k = reference.size();
    if (auto checked = checkRows(model, workspace, jacobian, k, regularisation, freeAcceleration, force, generalized,
                                 supports, form, system);
        !checked)
    {
        return std::unexpected(checked.error());
    }
    if (!kinds.empty() && kinds.size() != k)
    {
        return invalidInput("soft contact solve needs consistent sizes within the solver's capacity");
    }
    const auto isKind = [&](std::size_t r, SoftRowKind kind) { return !kinds.empty() && kinds[r] == kind; };
    bool cones = false;
    bool boxes = false;
    for (std::size_t r = 0; r < k; ++r)
    {
        if (isKind(r, SoftRowKind::FrictionLoss))
        {
            if (friction.size() != k || !std::isfinite(friction[r]) || friction[r] < 0.0)
            {
                return invalidInput("a friction-loss row needs a finite bound that is not negative");
            }
            boxes = true;
            continue;
        }
        if (isKind(r, SoftRowKind::ConeTangent))
        {
            return invalidInput("a friction cone is a normal row, then its tangent row");
        }
        if (isKind(r, SoftRowKind::ConeNormal))
        {
            if (r + 1 == k || !isKind(r + 1, SoftRowKind::ConeTangent) || friction.size() != k ||
                !std::isfinite(friction[r + 1]) || !(friction[r + 1] > 0.0))
            {
                return invalidInput("a friction cone is a normal row, then its tangent row, with a positive mu");
            }
            cones = true;
            ++r;
        }
    }
    std::ranges::fill(generalized, 0.0);
    SoftConstraintSolution solution;
    if (k == 0)
    {
        return solution;
    }
    const auto formed = formRows(model, workspace, jacobian, reference, regularisation, freeAcceleration, supports, form,
                                 system, k);
    if (!formed)
    {
        return std::unexpected(formed.error());
    }
    double scale = *formed;
    if (cones)
    {
        // Each cone's rows to its edges' coordinates, in place: A <- E^T A E, q <- E^T q (the class comment).
        scale = 0.0;
        for (std::size_t r = 0; r < k; ++r)
        {
            if (isKind(r, SoftRowKind::ConeNormal))
            {
                const double mu = friction[r + 1];
                const double qn = rhs_[r];
                const double qt = rhs_[r + 1];
                rhs_[r] = qn + mu * qt;
                rhs_[r + 1] = qn - mu * qt;
                for (std::size_t c = 0; c < k; ++c)
                {
                    const double an = delassus_[r * k + c];
                    const double at = delassus_[(r + 1) * k + c];
                    delassus_[r * k + c] = an + mu * at;
                    delassus_[(r + 1) * k + c] = an - mu * at;
                }
                for (std::size_t c = 0; c < k; ++c)
                {
                    const double an = delassus_[c * k + r];
                    const double at = delassus_[c * k + r + 1];
                    delassus_[c * k + r] = an + mu * at;
                    delassus_[c * k + r + 1] = an - mu * at;
                }
            }
            scale = std::max(scale, std::abs(rhs_[r]));
        }
    }
    // Murty's least-index principal pivoting, from the rows the free motion violates. A row is inactive (zero force),
    // active (free), or, for friction loss, at its upper or its lower bound.
    constexpr std::uint8_t Inactive = 0U;
    constexpr std::uint8_t Active = 1U;
    constexpr std::uint8_t AtUpper = 2U;
    constexpr std::uint8_t AtLower = 3U;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    const auto equality = [&](std::size_t r) { return isKind(r, SoftRowKind::Equality); };
    const auto box = [&](std::size_t r) { return boxes && isKind(r, SoftRowKind::FrictionLoss); };
    for (std::size_t r = 0; r < k; ++r)
    {
        // A friction-loss row starts free; a zero bound keeps it inactive, its force pinned at zero.
        active_[r] = box(r) ? (friction[r] > 0.0 ? Active : Inactive) : equality(r) || rhs_[r] < 0.0 ? Active : Inactive;
    }
    const std::uint32_t maxPivots = static_cast<std::uint32_t>(16U * k + 64U);
    std::size_t priorCount = 0;
    // The block rule: flip every violating row while their count falls or within a few stalls, then one at a time.
    bool block = rule_ == SoftPivotRule::Block;
    std::size_t fewest = k + 1;
    int stalls = 3;
    for (;;)
    {
        // Solve A_SS f_S = -(q_S + A_SB f_B) by Cholesky, the rows at a bound (B) held there.
        std::size_t m = 0;
        std::size_t held = 0;
        std::ranges::fill(force, 0.0);
        for (std::size_t r = 0; r < k; ++r)
        {
            if (active_[r] == Active)
            {
                subset_[m++] = r;
            }
            else if (active_[r] != Inactive)
            {
                bound_[held++] = r;
                force[r] = active_[r] == AtUpper ? friction[r] : -friction[r];
            }
        }
        if (!factorSubset(m, k, priorCount))
        {
            return invalidInput("soft contact matrix is not positive definite");
        }
        for (std::size_t a = 0; a < m; ++a)
        {
            double sum = -rhs_[subset_[a]];
            for (std::size_t b = 0; b < held; ++b)
            {
                sum -= delassus_[subset_[a] * k + bound_[b]] * force[bound_[b]];
            }
            work_[a] = sum;
        }
        solveSubset(m, k);
        for (std::size_t a = 0; a < m; ++a)
        {
            force[subset_[a]] = work_[a];
        }
        // w = A f + q, and the rows that break complementarity (the first only, by the least-index rule).
        std::size_t flip = k;
        std::size_t violations = 0;
        for (std::size_t r = 0; r < k && (block || flip == k); ++r)
        {
            bool violated = false;
            if (active_[r] == Active)
            {
                violated = box(r) ? std::abs(force[r]) > friction[r] : force[r] < 0.0 && !equality(r);
            }
            else if (!box(r) || friction[r] > 0.0)
            {
                double w = rhs_[r];
                // PERF: only active rows carry force; the rest add exact zeros.
                for (std::size_t a = 0; a < m; ++a)
                {
                    w += delassus_[r * k + subset_[a]] * force[subset_[a]];
                }
                for (std::size_t b = 0; b < held; ++b)
                {
                    w += delassus_[r * k + bound_[b]] * force[bound_[b]];
                }
                slack_[r] = w;
                // Inactive or at the lower bound, w must not be negative; at the upper bound, not positive.
                violated = active_[r] == AtUpper ? w > tolerance : w < -tolerance;
            }
            if (violated)
            {
                flip = std::min(flip, r);
                violations += 1;
                flagged_[violations - 1U] = r;
            }
        }
        if (flip == k)
        {
            solution.activeRows = m + held;
            break;
        }
        if (++solution.pivots > maxPivots)
        {
            return invalidInput("soft contact pivoting did not terminate");
        }
        // A violating row's next state: a free friction-loss row to the bound it crossed, a held one free, others toggled.
        const auto next = [&](std::size_t r) -> std::uint8_t {
            if (box(r))
            {
                return active_[r] != Active ? Active : force[r] > 0.0 ? AtUpper : AtLower;
            }
            return active_[r] != Inactive ? Inactive : Active;
        };
        if (block)
        {
            if (violations < fewest)
            {
                fewest = violations;
                stalls = 3;
            }
            else if (--stalls < 0)
            {
                block = false; // Murty's rule from here terminates on its own
            }
        }
        if (block)
        {
            for (std::size_t v = 0; v < violations; ++v)
            {
                active_[flagged_[v]] = next(flagged_[v]);
            }
        }
        else
        {
            active_[flip] = next(flip);
        }
    }
    if (cones)
    {
        // Back from the edges' coordinates: f_n = g+ + g-, f_t = mu (g+ - g-).
        for (std::size_t r = 0; r < k; ++r)
        {
            if (isKind(r, SoftRowKind::ConeNormal))
            {
                const double plus = force[r];
                const double minus = force[r + 1];
                force[r] = plus + minus;
                force[r + 1] = friction[r + 1] * (plus - minus);
                ++r;
            }
        }
    }
    if (auto mapped = mapForces(model, jacobian, force, generalized, supports, k); !mapped)
    {
        return std::unexpected(mapped.error());
    }
    return solution;
}

std::expected<SoftConstraintSolution, Fault>
SoftConstraintSolver::solveBounded(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                                   std::span<const double> reference, std::span<const double> regularisation,
                                   std::span<const double> freeAcceleration, std::span<const double> bounds,
                                   std::span<double> force, std::span<double> generalized,
                                   std::span<const SoftRowSupport> supports, DelassusForm form,
                                   SoftRowSystem system) noexcept
{
    const std::size_t k = reference.size();
    if (auto checked = checkRows(model, workspace, jacobian, k, regularisation, freeAcceleration, force, generalized,
                                 supports, form, system);
        !checked)
    {
        return std::unexpected(checked.error());
    }
    if (bounds.size() != k)
    {
        return invalidInput("soft contact solve needs consistent sizes within the solver's capacity");
    }
    for (const double b : bounds)
    {
        if (!std::isfinite(b) || b < 0.0)
        {
            return invalidInput("bounded rows need a finite bound that is not negative");
        }
    }
    std::ranges::fill(generalized, 0.0);
    SoftConstraintSolution solution;
    if (k == 0)
    {
        return solution;
    }
    const auto formed = formRows(model, workspace, jacobian, reference, regularisation, freeAcceleration, supports, form,
                                 system, k);
    if (!formed)
    {
        return std::unexpected(formed.error());
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * *formed;
    constexpr std::uint8_t Free = 0U;
    constexpr std::uint8_t AtUpper = 1U;
    constexpr std::uint8_t AtLower = 2U;
    for (std::size_t r = 0; r < k; ++r)
    {
        active_[r] = bounds[r] > 0.0 ? Free : AtUpper; // a zero bound is pinned (and stays at zero force)
    }
    const std::uint32_t maxPivots = static_cast<std::uint32_t>(16U * k + 64U);
    std::size_t priorCount = 0;
    for (;;)
    {
        // Free rows: A_FF f_F = -(q_F + A_FB f_B), the others at their bounds.
        std::size_t m = 0;
        for (std::size_t r = 0; r < k; ++r)
        {
            if (active_[r] == Free)
            {
                subset_[m++] = r;
                force[r] = 0.0;
            }
            else
            {
                force[r] = active_[r] == AtUpper ? bounds[r] : -bounds[r];
            }
        }
        if (!factorSubset(m, k, priorCount))
        {
            return invalidInput("bounded row matrix is not positive definite");
        }
        for (std::size_t a = 0; a < m; ++a)
        {
            double sum = -rhs_[subset_[a]];
            for (std::size_t c = 0; c < k; ++c)
            {
                sum -= delassus_[subset_[a] * k + c] * force[c];
            }
            work_[a] = sum;
        }
        solveSubset(m, k);
        for (std::size_t a = 0; a < m; ++a)
        {
            force[subset_[a]] = work_[a];
        }
        // The first row that breaks the box's complementarity: a free row beyond a bound, or a row at a bound whose
        // gradient w = A f + q pushes it inside.
        std::size_t flip = k;
        std::uint8_t to = Free;
        for (std::size_t r = 0; r < k && flip == k; ++r)
        {
            if (active_[r] == Free)
            {
                if (force[r] > bounds[r])
                {
                    flip = r;
                    to = AtUpper;
                }
                else if (force[r] < -bounds[r])
                {
                    flip = r;
                    to = AtLower;
                }
                continue;
            }
            if (bounds[r] == 0.0)
            {
                continue;
            }
            double w = rhs_[r];
            for (std::size_t c = 0; c < k; ++c)
            {
                w += delassus_[r * k + c] * force[c];
            }
            if (active_[r] == AtUpper ? w > tolerance : w < -tolerance)
            {
                flip = r;
                to = Free;
            }
        }
        if (flip == k)
        {
            solution.activeRows = m;
            break;
        }
        if (++solution.pivots > maxPivots)
        {
            return invalidInput("bounded row pivoting did not terminate");
        }
        active_[flip] = to;
    }
    if (auto mapped = mapForces(model, jacobian, force, generalized, supports, k); !mapped)
    {
        return std::unexpected(mapped.error());
    }
    return solution;
}

} // namespace sub0physics::multibody
