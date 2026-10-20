// Constraints (docs/MULTIBODY.md section 6, M1c): loop closures and strict unilateral limits.
//
// A constrained step is the stepper's own free step plus impulses: with the step matrix A (M for semi-implicit
// Euler, M - h D - h^2 K~ for linearly implicit Euler) already factorised, the rows J (closures, and the limits
// engaged at the start) give the small dense Delassus matrix G = J A^-1 J^T; a mixed LCP (closures equalities,
// limits complementary) gives the impulses p, and qd+ = qd_free + A^-1 J^T p. The new positions are projected
// onto the closure manifold (Gauss-Newton in the metric A^-1 until the requested residual is met) and the
// velocities onto its tangent space. A limit crossed during the step splits it: the engagement is located by a
// root find on the step size, so every accepted state is inside its limits.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>

#include <sub0math/linear/mat.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/quat.hpp>
#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/linear/vec.hpp>

#include "multibody_joint.hpp"
#include "sub0physics/multibody/dynamics.hpp"
#include "sub0physics/multibody/kinematics.hpp"

namespace sub0physics::multibody
{

namespace
{

using sub0math::Mat3d;
using sub0math::Vec3d;

Vec3d plain(const sub0math::Vec3<sub0math::f64::Length>& v) noexcept
{
    return {v.x.si(), v.y.si(), v.z.si()};
}

// A site as a point on a body: its body (NoParent for the world) and world position. Needs worldToBody.
struct SitePoint
{
    std::size_t body;
    Vec3d position;
};

SitePoint sitePoint(const Model& model, const Workspace& ws, SiteId site) noexcept
{
    const Vec3d local = plain(model.siteFrame(site).translation);
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return {body, local};
    }
    const sub0math::PluckerTransformd& x = ws.worldToBody[body];
    return {body, x.translation + detail::bodyToWorld(x.rotation, local)};
}

// The site frame's axes in the world: columns are the site's x, y, z axes in world coordinates.
Mat3d siteRotation(const Model& model, const Workspace& ws, SiteId site) noexcept
{
    const Mat3d local = sub0math::toMat3(model.siteFrame(site).rotation);
    const std::size_t body = model.siteBodies()[site.value()];
    if (body == NoParent)
    {
        return local;
    }
    return sub0math::transpose(ws.worldToBody[body].rotation) * local;
}

// Adds `scale` times the translational Jacobian of the world point p, moved by `body`'s chain, to three rows.
void addPointJacobian(const Model& model, const Workspace& ws, std::size_t body, const Vec3d& p, double scale,
                      double* rx, double* ry, double* rz) noexcept
{
    const auto parents = model.parents();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const auto coordinates = model.coordinates();
    for (std::size_t j = body; j != NoParent; j = parents[j])
    {
        const sub0math::PluckerTransformd& x = ws.worldToBody[j];
        const Vec3d w = detail::bodyToWorld(x.rotation, axes[j]);
        const Vec3d column = types[j] == JointType::Revolute ? cross(w, p - x.translation) : w;
        const std::size_t c = coordinates[j];
        rx[c] += scale * column.x;
        ry[c] += scale * column.y;
        rz[c] += scale * column.z;
    }
}

// Adds `scale` times the angular-velocity Jacobian of `body` (world axes) to three rows.
void addAngularJacobian(const Model& model, const Workspace& ws, std::size_t body, double scale, double* rx,
                        double* ry, double* rz) noexcept
{
    const auto parents = model.parents();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const auto coordinates = model.coordinates();
    for (std::size_t j = body; j != NoParent; j = parents[j])
    {
        if (types[j] != JointType::Revolute)
        {
            continue;
        }
        const Vec3d w = detail::bodyToWorld(ws.worldToBody[j].rotation, axes[j]);
        const std::size_t c = coordinates[j];
        rx[c] += scale * w.x;
        ry[c] += scale * w.y;
        rz[c] += scale * w.z;
    }
}

// World velocity of the world point p moving with `body`. Needs bodyVelocity and worldToBody.
Vec3d pointVelocity(const Workspace& ws, std::size_t body, const Vec3d& p) noexcept
{
    if (body == NoParent)
    {
        return {};
    }
    const sub0math::PluckerTransformd& x = ws.worldToBody[body];
    const sub0math::SpatialMotiond& v = ws.bodyVelocity[body];
    const Vec3d r = x.rotation * (p - x.translation);
    return detail::bodyToWorld(x.rotation, v.linear + cross(v.angular, r));
}

// The body's spatial acceleration at zero qdd without gravity, in its own frame: what Jdot qd is made of. Needs the
// zero-qdd RNEA's bodyAcceleration (which carries -g) and worldToBody.
sub0math::SpatialMotiond biasAcceleration(const Model& model, const Workspace& ws, std::size_t body) noexcept
{
    const sub0math::PluckerTransformd& x = ws.worldToBody[body];
    const sub0math::SpatialMotiond& a = ws.bodyAcceleration[body];
    return {a.angular, a.linear - x.rotation * model.baseAcceleration().linear};
}

// World acceleration of the world point p moving with `body`, at zero qdd and without gravity.
Vec3d pointBiasAcceleration(const Model& model, const Workspace& ws, std::size_t body, const Vec3d& p) noexcept
{
    if (body == NoParent)
    {
        return {};
    }
    const sub0math::PluckerTransformd& x = ws.worldToBody[body];
    const sub0math::SpatialMotiond& v = ws.bodyVelocity[body];
    const sub0math::SpatialMotiond a = biasAcceleration(model, ws, body);
    const Vec3d r = x.rotation * (p - x.translation);
    const Vec3d local = a.linear + cross(a.angular, r) + cross(v.angular, v.linear + cross(v.angular, r));
    return detail::bodyToWorld(x.rotation, local);
}

Vec3d angularBiasAcceleration(const Model& model, const Workspace& ws, std::size_t body) noexcept
{
    if (body == NoParent)
    {
        return {};
    }
    return detail::bodyToWorld(ws.worldToBody[body].rotation, biasAcceleration(model, ws, body).angular);
}

// The orientation error of two frames as a vector: vee of the skew part of R1 R2^T, whose rate is w1 - w2 on the
// manifold (both axes sets in world coordinates).
Vec3d orientationError(const Mat3d& r1, const Mat3d& r2) noexcept
{
    const Mat3d e = r1 * sub0math::transpose(r2);
    return {0.5 * (e(2, 1) - e(1, 2)), 0.5 * (e(0, 2) - e(2, 0)), 0.5 * (e(1, 0) - e(0, 1))};
}

// The rows of a constraint set, assembled at the kinematics in the workspace.
struct Rows
{
    std::size_t count = 0;     // assembled rows
    std::size_t bilateral = 0; // of which closure rows (the first)
};

// Distance between a distance limit's sites: direction and length. False when they coincide.
struct Separation
{
    SitePoint first;
    SitePoint second;
    Vec3d direction;
    double length;
};

Separation separation(const Model& model, const Workspace& ws, const UnilateralLimit& limit) noexcept
{
    const SitePoint p1 = sitePoint(model, ws, limit.first);
    const SitePoint p2 = sitePoint(model, ws, limit.second);
    const Vec3d d = p1.position - p2.position;
    const double l = std::sqrt(dot(d, d));
    return {p1, p2, l > 0.0 ? d * (1.0 / l) : Vec3d{}, l};
}

// Gap of every limit at `q`; needs worldToBody at q when the model has distance limits. Faults when a distance
// limit's sites coincide.
std::expected<void, Fault> limitGaps(const Model& model, std::span<const double> q, Workspace& ws) noexcept
{
    const auto limits = model.limits();
    for (std::size_t i = 0; i < limits.size(); ++i)
    {
        const UnilateralLimit& limit = limits[i];
        double x = 0.0;
        if (limit.kind == LimitKind::Coordinate)
        {
            x = q[limit.coordinate];
        }
        else
        {
            const Separation s = separation(model, ws, limit);
            if (!(s.length > 0.0) || !std::isfinite(s.length))
            {
                return std::unexpected(Fault{FaultKind::OutsideDomain, "distance limit sites coincide"});
            }
            x = s.length;
        }
        ws.constraintScratch.gap[i] = limit.sign * (x - limit.bound);
    }
    return {};
}

// Assembles the closure rows and the rows of the limits flagged in `include` (none if null) into `jacobian`
// (rows x n), their residuals (closure phi, limit gap) into `residual` and their layout rows into `layoutRow`.
// Needs worldToBody (and limitGaps for the limits). With `bias`, also writes Jdot qd per row (needs bodyVelocity and
// the zero-qdd bodyAcceleration).
Rows assembleRows(const Model& model, Workspace& ws, const std::uint8_t* include, double* jacobian,
                  std::size_t* layoutRow, double* bias) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    const std::size_t n = model.dofCount();
    Rows rows;
    const auto newRow = [&](std::size_t layout) -> double* {
        double* row = jacobian + rows.count * n;
        std::fill_n(row, n, 0.0);
        layoutRow[rows.count] = layout;
        ++rows.count;
        return row;
    };
    for (const LoopClosure& closure : model.closures())
    {
        const SitePoint p1 = sitePoint(model, ws, closure.first);
        const SitePoint p2 = sitePoint(model, ws, closure.second);
        const std::size_t first = rows.count;
        double* rx = newRow(closure.firstRow);
        double* ry = newRow(closure.firstRow + 1);
        double* rz = newRow(closure.firstRow + 2);
        if (p1.body != NoParent)
        {
            addPointJacobian(model, ws, p1.body, p1.position, 1.0, rx, ry, rz);
        }
        if (p2.body != NoParent)
        {
            addPointJacobian(model, ws, p2.body, p2.position, -1.0, rx, ry, rz);
        }
        const Vec3d phi = p1.position - p2.position;
        sc.residual[first] = phi.x;
        sc.residual[first + 1] = phi.y;
        sc.residual[first + 2] = phi.z;
        if (bias != nullptr)
        {
            const Vec3d a = pointBiasAcceleration(model, ws, p1.body, p1.position) -
                            pointBiasAcceleration(model, ws, p2.body, p2.position);
            bias[first] = a.x;
            bias[first + 1] = a.y;
            bias[first + 2] = a.z;
        }
        if (closure.kind == ClosureKind::Frame)
        {
            double* ax = newRow(closure.firstRow + 3);
            double* ay = newRow(closure.firstRow + 4);
            double* az = newRow(closure.firstRow + 5);
            if (p1.body != NoParent)
            {
                addAngularJacobian(model, ws, p1.body, 1.0, ax, ay, az);
            }
            if (p2.body != NoParent)
            {
                addAngularJacobian(model, ws, p2.body, -1.0, ax, ay, az);
            }
            const Vec3d e =
                orientationError(siteRotation(model, ws, closure.first), siteRotation(model, ws, closure.second));
            sc.residual[first + 3] = e.x;
            sc.residual[first + 4] = e.y;
            sc.residual[first + 5] = e.z;
            if (bias != nullptr)
            {
                const Vec3d a =
                    angularBiasAcceleration(model, ws, p1.body) - angularBiasAcceleration(model, ws, p2.body);
                bias[first + 3] = a.x;
                bias[first + 4] = a.y;
                bias[first + 5] = a.z;
            }
        }
    }
    rows.bilateral = rows.count;
    if (include == nullptr)
    {
        return rows;
    }
    const auto limits = model.limits();
    for (std::size_t i = 0; i < limits.size(); ++i)
    {
        if (include[i] == 0U)
        {
            continue;
        }
        const UnilateralLimit& limit = limits[i];
        const std::size_t r = rows.count;
        double* row = newRow(limit.row);
        sc.residual[r] = sc.gap[i];
        if (limit.kind == LimitKind::Coordinate)
        {
            row[limit.coordinate] = limit.sign;
            if (bias != nullptr)
            {
                bias[r] = 0.0;
            }
            continue;
        }
        const Separation s = separation(model, ws, limit);
        // J_gap = sign u^T (J1 - J2): accumulate the three translational rows weighted by u.
        const auto addWeighted = [&](const SitePoint& p, double scale) {
            if (p.body == NoParent)
            {
                return;
            }
            const auto parents = model.parents();
            const auto types = model.jointTypes();
            const auto axes = model.jointAxes();
            for (std::size_t j = p.body; j != NoParent; j = parents[j])
            {
                const sub0math::PluckerTransformd& x = ws.worldToBody[j];
                const Vec3d w = detail::bodyToWorld(x.rotation, axes[j]);
                const Vec3d column = types[j] == JointType::Revolute ? cross(w, p.position - x.translation) : w;
                row[model.coordinates()[j]] += scale * dot(s.direction, column);
            }
        };
        addWeighted(s.first, limit.sign);
        addWeighted(s.second, -limit.sign);
        if (bias != nullptr)
        {
            const Vec3d dd = pointVelocity(ws, s.first.body, s.first.position) -
                             pointVelocity(ws, s.second.body, s.second.position);
            const double rate = dot(s.direction, dd);
            const Vec3d a = pointBiasAcceleration(model, ws, s.first.body, s.first.position) -
                            pointBiasAcceleration(model, ws, s.second.body, s.second.position);
            bias[r] = limit.sign * (dot(s.direction, a) + (dot(dd, dd) - rate * rate) / s.length);
        }
    }
    return rows;
}

// Solves G x = r over the rows `subset[0 .. m)` of the symmetric positive semidefinite G (stride `ld`), writing
// x[subset[k]] and leaving other entries of x alone. LDL^T with symmetric diagonal pivoting, the largest remaining
// pivot first (deterministic: the first of equal ones); pivots at or below `tolerance` times the largest diagonal
// are dependent rows, whose x is set to zero. Returns the rank.
std::size_t solvePivoted(std::size_t m, const double* g, std::size_t ld, const std::size_t* subset, const double* r,
                         double* x, ConstraintScratch& sc, double tolerance) noexcept
{
    double* work = sc.factor.data();
    std::size_t* perm = sc.permutation.data();
    double* y = sc.vector.data();
    double largest = 0.0;
    for (std::size_t i = 0; i < m; ++i)
    {
        perm[i] = subset[i];
        for (std::size_t j = 0; j < m; ++j)
        {
            work[i * m + j] = g[subset[i] * ld + subset[j]];
        }
        largest = std::max(largest, work[i * m + i]);
    }
    const double floor = tolerance * largest;
    std::size_t rank = 0;
    for (std::size_t k = 0; k < m; ++k)
    {
        std::size_t p = k;
        for (std::size_t i = k + 1; i < m; ++i)
        {
            if (work[i * m + i] > work[p * m + p])
            {
                p = i;
            }
        }
        if (!(work[p * m + p] > floor))
        {
            break;
        }
        if (p != k)
        {
            for (std::size_t j = 0; j < m; ++j)
            {
                std::swap(work[k * m + j], work[p * m + j]);
            }
            for (std::size_t i = 0; i < m; ++i)
            {
                std::swap(work[i * m + k], work[i * m + p]);
            }
            std::swap(perm[k], perm[p]);
        }
        const double d = work[k * m + k];
        for (std::size_t i = k + 1; i < m; ++i)
        {
            work[i * m + k] /= d;
        }
        for (std::size_t i = k + 1; i < m; ++i)
        {
            const double lik = work[i * m + k] * d;
            for (std::size_t j = k + 1; j <= i; ++j)
            {
                work[i * m + j] -= lik * work[j * m + k];
                work[j * m + i] = work[i * m + j];
            }
        }
        rank = k + 1;
    }
    for (std::size_t k = 0; k < rank; ++k)
    {
        double sum = r[perm[k]];
        for (std::size_t j = 0; j < k; ++j)
        {
            sum -= work[k * m + j] * y[j];
        }
        y[k] = sum;
    }
    for (std::size_t k = 0; k < rank; ++k)
    {
        y[k] /= work[k * m + k];
    }
    for (std::size_t k = rank; k-- > 0;)
    {
        double sum = y[k];
        for (std::size_t i = k + 1; i < rank; ++i)
        {
            sum -= work[i * m + k] * y[i];
        }
        y[k] = sum;
    }
    for (std::size_t k = 0; k < m; ++k)
    {
        x[perm[k]] = k < rank ? y[k] : 0.0;
    }
    return rank;
}

// Fills `solved` (A^-1 J^T by row) and `delassus` (J A^-1 J^T) for `count` rows of `jacobian`, with A applied by
// `solve(x)`.
template<typename Solve>
void delassusMatrix(std::size_t n, std::size_t count, const double* jacobian, ConstraintScratch& sc, Solve&& solve) noexcept
{
    for (std::size_t r = 0; r < count; ++r)
    {
        double* y = sc.solved.data() + r * n;
        std::copy_n(jacobian + r * n, n, y);
        solve(y);
    }
    for (std::size_t r = 0; r < count; ++r)
    {
        const double* jr = jacobian + r * n;
        for (std::size_t s = 0; s <= r; ++s)
        {
            const double* ys = sc.solved.data() + s * n;
            double sum = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                sum += jr[i] * ys[i];
            }
            sc.delassus[r * count + s] = sum;
            sc.delassus[s * count + r] = sum;
        }
    }
}

// Mixed LCP on `rows` with the Delassus matrix in sc.delassus (stride rows.count) and q in sc.rhs:
//   w = q + G p;  w_i = 0 on closure rows;  w_i >= 0, p_i >= 0, w_i p_i = 0 on limit rows.
// Principal pivoting by Murty's least-index rule: solve the active set's equalities, then flip the first limit row
// whose sign condition fails. Finite for a positive definite G; exact for a few rows. Writes p to sc.impulse and w to
// sc.slack; returns the number of pivots.
std::expected<std::uint32_t, Fault> solveMixedLcp(const Rows& rows, ConstraintScratch& sc, double tolerance) noexcept
{
    const std::size_t m = rows.count;
    double* q = sc.rhs.data();
    double* p = sc.impulse.data();
    double* w = sc.slack.data();
    double* negativeQ = sc.target.data(); // reused: the right-hand side -q
    double scaleQ = 0.0;
    for (std::size_t i = 0; i < m; ++i)
    {
        sc.active[i] = (i < rows.bilateral || q[i] < 0.0) ? 1U : 0U;
        negativeQ[i] = -q[i];
        scaleQ = std::max(scaleQ, std::abs(q[i]));
    }
    // Murty's rule terminates on a P-matrix in at most 2^k pivots; a handful of limits needs a handful.
    const std::uint32_t maxPivots = 8U * static_cast<std::uint32_t>(m) + 16U;
    for (std::uint32_t pivots = 0; pivots <= maxPivots; ++pivots)
    {
        std::size_t count = 0;
        for (std::size_t i = 0; i < m; ++i)
        {
            p[i] = 0.0;
            if (sc.active[i] != 0U)
            {
                sc.subset[count++] = i;
            }
        }
        solvePivoted(count, sc.delassus.data(), m, sc.subset.data(), negativeQ, p, sc, tolerance);
        double scaleP = 0.0;
        for (std::size_t i = 0; i < m; ++i)
        {
            double sum = q[i];
            for (std::size_t j = 0; j < m; ++j)
            {
                sum += sc.delassus[i * m + j] * p[j];
            }
            w[i] = sum;
            scaleP = std::max(scaleP, std::abs(p[i]));
        }
        const double tolerateW = tolerance * (scaleQ + 1e-300);
        const double tolerateP = tolerance * (scaleP + 1e-300);
        std::size_t flip = m;
        for (std::size_t i = rows.bilateral; i < m; ++i)
        {
            if ((sc.active[i] != 0U && p[i] < -tolerateP) || (sc.active[i] == 0U && w[i] < -tolerateW))
            {
                flip = i;
                break;
            }
        }
        if (flip == m)
        {
            for (std::size_t i = rows.bilateral; i < m; ++i)
            {
                p[i] = std::max(p[i], 0.0); // within tolerance of zero
            }
            return pivots;
        }
        sc.active[flip] = sc.active[flip] != 0U ? 0U : 1U;
    }
    return std::unexpected(Fault{FaultKind::OutsideDomain, "limit LCP did not converge"});
}

// How the projection measures distance: the step matrix A (its factor in implicitMatrix) or the mass matrix.
struct Metric
{
    bool mass;  // massFactor; else the step matrix
    bool dense; // the step matrix's factor is dense
};

void applyMetric(const Model& model, const Workspace& ws, const Metric& metric, double* x) noexcept
{
    if (metric.mass)
    {
        detail::solveMassFactor(model, ws, x);
    }
    else
    {
        detail::solveStepMatrix(model, ws, metric.dense, x);
    }
}

// Updates the transforms and world transforms (and the velocities with `velocities`) at a state.
void kinematicsAt(const Model& model, const State& state, Workspace& ws, bool velocities) noexcept
{
    if (velocities)
    {
        detail::updateTransforms<true>(model, state, ws);
    }
    else
    {
        detail::updateTransforms<false>(model, state, ws);
    }
    detail::composeWorldTransforms(model, ws);
}

bool hasDistanceLimits(const Model& model) noexcept
{
    return std::ranges::any_of(model.limits(), [](const UnilateralLimit& l) { return l.kind == LimitKind::Distance; });
}

// Energy of a state: kinetic, gravitational and the elements' potential.
std::expected<double, Fault> mechanicalEnergy(const Model& model, const State& state, Workspace& ws) noexcept
{
    const auto energy = computeEnergy(model, state, ws);
    if (!energy)
    {
        return std::unexpected(energy.error());
    }
    return energy->mechanical().si();
}

// Projects a state onto the constraints (docs/MULTIBODY.md section 6): Gauss-Newton on the closures and on the
// flagged limits that are violated, dq = -W J^T (J W J^T)^-1 phi with W = A^-1, until max |phi| meets the requested
// tolerance; coordinate limits still violated are then set to their bound (within the tolerance); finally the
// velocities are projected onto the closures' tangent space, qd -= W J^T (J W J^T)^-1 J qd. Writes the achieved
// residuals and the iteration count to the report.
std::expected<void, Fault> project(const Model& model, State& state, Workspace& ws, const Metric& metric,
                                   const std::uint8_t* includeLimits) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    ConstraintReport& report = ws.constraintReport;
    const ConstraintSettings& settings = model.constraintSettings();
    const std::size_t n = model.dofCount();
    const auto limits = model.limits();
    std::uint8_t* violated = sc.active.data(); // reused per iteration: the flagged limits with a negative gap
    const auto solve = [&](double* x) { applyMetric(model, ws, metric, x); };
    double residual = 0.0;
    Rows rows;
    for (std::uint32_t iteration = 0;; ++iteration)
    {
        kinematicsAt(model, state, ws, false);
        if (auto gaps = limitGaps(model, state.q(), ws); !gaps)
        {
            return gaps;
        }
        bool anyViolated = false;
        for (std::size_t i = 0; i < limits.size(); ++i)
        {
            // `active` has room for the limits: the layout has a row per limit.
            violated[i] = (includeLimits != nullptr && includeLimits[i] != 0U && sc.gap[i] < 0.0) ? 1U : 0U;
            anyViolated = anyViolated || violated[i] != 0U;
        }
        rows = assembleRows(model, ws, anyViolated ? violated : nullptr, sc.jacobian.data(),
                            sc.layoutRow.data(), nullptr);
        residual = 0.0;
        double closureResidual = 0.0;
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            residual = std::max(residual, std::abs(sc.residual[r]));
            if (r < rows.bilateral)
            {
                closureResidual = std::max(closureResidual, std::abs(sc.residual[r]));
            }
        }
        report.closureResidual = closureResidual;
        if (residual <= settings.closureTolerance || settings.maxProjectionIterations == 0U)
        {
            break;
        }
        if (iteration == settings.maxProjectionIterations)
        {
            return std::unexpected(Fault{FaultKind::OutsideDomain, "projection did not reach the closure tolerance"});
        }
        delassusMatrix(n, rows.count, sc.jacobian.data(), sc, solve);
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            sc.subset[r] = r;
            sc.rhs[r] = sc.residual[r];
        }
        // sc.vector is solvePivoted's scratch, so the multipliers go to sc.slack.
        solvePivoted(rows.count, sc.delassus.data(), rows.count, sc.subset.data(), sc.rhs.data(), sc.slack.data(), sc,
                     settings.solverTolerance);
        double* q = state.values.data();
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            const double* y = sc.solved.data() + r * n;
            for (std::size_t i = 0; i < n; ++i)
            {
                q[i] -= sc.slack[r] * y[i];
            }
        }
        ++report.projectionIterations;
    }
    if (settings.maxProjectionIterations == 0U)
    {
        return {};
    }
    // A coordinate limit is exact: set it to its bound (it is within the tolerance of it).
    for (std::size_t i = 0; i < limits.size(); ++i)
    {
        if (limits[i].kind == LimitKind::Coordinate && sc.gap[i] < 0.0 && includeLimits != nullptr &&
            includeLimits[i] != 0U)
        {
            state.q()[limits[i].coordinate] = limits[i].bound;
        }
    }
    if (rows.bilateral == 0)
    {
        return {};
    }
    // Velocities onto the closures' tangent space, with the rows at the projected positions.
    delassusMatrix(n, rows.bilateral, sc.jacobian.data(), sc, solve);
    double* qd = state.values.data() + n;
    for (std::size_t r = 0; r < rows.bilateral; ++r)
    {
        const double* jr = sc.jacobian.data() + r * n;
        double rate = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            rate += jr[i] * qd[i];
        }
        sc.rhs[r] = rate;
        sc.subset[r] = r;
    }
    solvePivoted(rows.bilateral, sc.delassus.data(), rows.bilateral, sc.subset.data(), sc.rhs.data(),
                 sc.slack.data(), sc, settings.solverTolerance);
    for (std::size_t r = 0; r < rows.bilateral; ++r)
    {
        const double* y = sc.solved.data() + r * n;
        for (std::size_t i = 0; i < n; ++i)
        {
            qd[i] -= sc.slack[r] * y[i];
        }
    }
    double velocityResidual = 0.0;
    for (std::size_t r = 0; r < rows.bilateral; ++r)
    {
        const double* jr = sc.jacobian.data() + r * n;
        double rate = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            rate += jr[i] * qd[i];
        }
        velocityResidual = std::max(velocityResidual, std::abs(rate));
    }
    report.closureVelocityResidual = velocityResidual;
    return {};
}

// Projection with the energy it changes measured (when the settings ask for it).
std::expected<double, Fault> projectMeasured(const Model& model, State& state, Workspace& ws, const Metric& metric,
                                             const std::uint8_t* includeLimits) noexcept
{
    const bool measure = model.constraintSettings().measureProjectionEnergy;
    double before = 0.0;
    if (measure)
    {
        const auto energy = mechanicalEnergy(model, state, ws);
        if (!energy)
        {
            return std::unexpected(energy.error());
        }
        before = *energy;
    }
    if (auto projected = project(model, state, ws, metric, includeLimits); !projected)
    {
        return std::unexpected(projected.error());
    }
    if (!measure)
    {
        return 0.0;
    }
    const auto after = mechanicalEnergy(model, state, ws);
    if (!after)
    {
        return std::unexpected(after.error());
    }
    return *after - before;
}

// One trial of a sub-step of `h` from sc.start into sc.trial. Needs prepareStepDynamics at sc.start, the start's
// element force in sc.elementForce and the start's rows in sc.stepJacobian.
struct Trial
{
    double closureWork = 0.0;
    double limitWork = 0.0;
    double projectionEnergy = 0.0;
    double minimumFreeGap = 0.0; // smallest gap among the limits not engaged at the start
    bool crossed = false;        // whether one of them is negative
};

std::expected<Trial, Fault> trialStep(const Model& model, const Forcing& forcing, double h, Workspace& ws,
                                      bool implicit, bool allowFallback) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    const ConstraintSettings& settings = model.constraintSettings();
    const std::size_t n = model.dofCount();
    const State& start = sc.start;
    const auto dense = detail::formStepMatrix(model, start, forcing.tau,
                                              model.elementCount() > 0 ? sc.elementForce.data() : nullptr, ws, h,
                                              implicit, allowFallback);
    if (!dense)
    {
        return std::unexpected(dense.error());
    }
    const Metric metric{false, *dense};
    double* v = sc.velocity.data();
    detail::solveStepMatrix(model, ws, *dense, ws.implicitRhs.data());
    const double* qd0 = start.values.data() + n;
    for (std::size_t i = 0; i < n; ++i)
    {
        v[i] = qd0[i] + ws.implicitRhs[i];
    }

    // Impulses: q_r = J_r v_free - target_r, the LCP, then v += A^-1 J^T p.
    Trial result;
    const Rows rows{sc.stepRows, sc.stepBilateral};
    if (rows.count > 0)
    {
        delassusMatrix(n, rows.count, sc.stepJacobian.data(), sc, [&](double* x) { applyMetric(model, ws, metric, x); });
        const auto limits = model.limits();
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            const double* jr = sc.stepJacobian.data() + r * n;
            double freeRate = 0.0;
            double startRate = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                freeRate += jr[i] * v[i];
                startRate += jr[i] * qd0[i];
            }
            double target = 0.0;
            if (r >= rows.bilateral)
            {
                const UnilateralLimit& limit = limits[sc.stepLayoutRow[r] - model.closureRowCount()];
                target = -limit.restitution * std::min(startRate, 0.0);
            }
            sc.residual[r] = freeRate; // kept for the impulse work
            sc.rhs[r] = freeRate - target;
        }
        const auto pivots = solveMixedLcp(rows, sc, settings.solverTolerance);
        if (!pivots)
        {
            return std::unexpected(pivots.error());
        }
        ws.constraintReport.lcpPivots += *pivots;
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            const double* y = sc.solved.data() + r * n;
            for (std::size_t i = 0; i < n; ++i)
            {
                v[i] += sc.impulse[r] * y[i];
            }
        }
        // Work of each impulse on the mean of the velocities before and after it: exact in the mass metric.
        for (std::size_t r = 0; r < rows.count; ++r)
        {
            const double* jr = sc.stepJacobian.data() + r * n;
            double rate = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                rate += jr[i] * v[i];
            }
            const double work = 0.5 * sc.impulse[r] * (sc.residual[r] + rate);
            (r < rows.bilateral ? result.closureWork : result.limitWork) += work;
        }
    }

    State& trial = sc.trial;
    const double* q0 = start.values.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(v[i]) || !std::isfinite(q0[i] + h * v[i]))
        {
            return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
        }
        trial.values[i] = q0[i] + h * v[i];
        trial.values[n + i] = v[i];
    }
    if (model.closures().size() > 0 || std::ranges::any_of(sc.engaged, [](std::uint8_t e) { return e != 0U; }))
    {
        const auto energy = projectMeasured(model, trial, ws, metric, sc.engaged.data());
        if (!energy)
        {
            return std::unexpected(energy.error());
        }
        result.projectionEnergy = *energy;
    }
    if (!model.limits().empty())
    {
        if (hasDistanceLimits(model))
        {
            kinematicsAt(model, trial, ws, false);
        }
        if (auto gaps = limitGaps(model, trial.q(), ws); !gaps)
        {
            return std::unexpected(gaps.error());
        }
        result.minimumFreeGap = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < model.limits().size(); ++i)
        {
            if (sc.engaged[i] == 0U)
            {
                result.minimumFreeGap = std::min(result.minimumFreeGap, sc.gap[i]);
            }
        }
        result.crossed = result.minimumFreeGap < 0.0;
    }
    return result;
}

// The start of a sub-step: dynamics at sc.start, the start's element force, the engaged limits and their rows.
std::expected<void, Fault> prepareSubstep(const Model& model, Workspace& ws, bool implicit) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    const ConstraintSettings& settings = model.constraintSettings();
    if (auto prepared = detail::prepareStepDynamics(model, sc.start, ws, implicit); !prepared)
    {
        return prepared;
    }
    std::copy(ws.elementForce.begin(), ws.elementForce.end(), sc.elementForce.begin());
    sc.startDissipationPower = model.elementCount() > 0 ? ws.elementDissipationPower : 0.0;
    if (auto gaps = limitGaps(model, sc.start.q(), ws); !gaps)
    {
        return gaps;
    }
    for (std::size_t i = 0; i < model.limits().size(); ++i)
    {
        sc.engaged[i] = sc.gap[i] <= settings.eventTolerance ? 1U : 0U;
    }
    const Rows rows = assembleRows(model, ws, sc.engaged.data(), sc.stepJacobian.data(),
                                   sc.stepLayoutRow.data(), nullptr);
    sc.stepRows = rows.count;
    sc.stepBilateral = rows.bilateral;
    return {};
}

} // namespace

namespace detail
{

std::expected<void, Fault> stepConstrained(const Model& model, const Forcing& forcing, double dt, State& state,
                                           Workspace& ws, bool implicit) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    ConstraintReport& report = ws.constraintReport;
    const ConstraintSettings& settings = model.constraintSettings();
    const std::size_t n = model.dofCount();
    report = {};
    std::ranges::fill(ws.constraintImpulse, 0.0);
    sc.start = state;

    // The invariant: an accepted state is inside its limits. A state given outside them is refused.
    if (!model.limits().empty())
    {
        kinematicsAt(model, state, ws, false);
        if (auto gaps = limitGaps(model, state.q(), ws); !gaps)
        {
            return gaps;
        }
        for (std::size_t i = 0; i < model.limits().size(); ++i)
        {
            if (sc.gap[i] < -settings.closureTolerance)
            {
                return invalidInput("state violates a limit");
            }
        }
    }

    // The tick is split at each engagement; `remaining` is the fraction of h still to go. Time stays tick * h: the
    // sub-steps only divide the integration, and nothing between them is exposed.
    double remaining = 1.0;
    double dissipated = 0.0;
    double closureWork = 0.0;
    double projectionEnergy = 0.0;
    for (;;)
    {
        if (auto prepared = prepareSubstep(model, ws, implicit); !prepared)
        {
            return prepared;
        }
        double startFreeGap = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < model.limits().size(); ++i)
        {
            if (sc.engaged[i] == 0U)
            {
                startFreeGap = std::min(startFreeGap, sc.gap[i]);
            }
        }
        const double h = remaining * dt;
        auto trial = trialStep(model, forcing, h, ws, implicit, true);
        if (!trial)
        {
            return std::unexpected(trial.error());
        }
        double fraction = 1.0;
        const bool split = trial->crossed;
        if (split)
        {
            if (report.substeps >= settings.maxEventsPerStep)
            {
                return std::unexpected(Fault{FaultKind::OutsideDomain, "too many limit engagements in one step"});
            }
            // Root find on the sub-step size (Illinois regula falsi): the smallest free gap is above eventTolerance
            // at 0 (the limits are not engaged) and negative at 1. Stop when it lands in [0, eventTolerance].
            double lo = 0.0;
            double hi = 1.0;
            double gapLo = startFreeGap;
            double gapHi = trial->minimumFreeGap;
            int side = 0;
            bool located = false;
            for (std::uint32_t iteration = 0; iteration < settings.maxEventIterations; ++iteration)
            {
                double theta = lo + (hi - lo) * gapLo / (gapLo - gapHi);
                if (!(theta > lo && theta < hi))
                {
                    theta = 0.5 * (lo + hi);
                }
                ++report.eventIterations;
                trial = trialStep(model, forcing, theta * h, ws, implicit, false);
                if (!trial)
                {
                    return std::unexpected(trial.error());
                }
                const double g = trial->minimumFreeGap;
                if (g >= 0.0 && g <= settings.eventTolerance)
                {
                    fraction = theta;
                    report.engagementGap = g;
                    located = true;
                    break;
                }
                if (g > settings.eventTolerance)
                {
                    lo = theta;
                    gapLo = g;
                    if (side == -1)
                    {
                        gapHi *= 0.5;
                    }
                    side = -1;
                }
                else
                {
                    hi = theta;
                    gapHi = g;
                    if (side == 1)
                    {
                        gapLo *= 0.5;
                    }
                    side = 1;
                }
            }
            if (!located)
            {
                return std::unexpected(Fault{FaultKind::OutsideDomain, "limit engagement not located to the event tolerance"});
            }
        }
        // Accept the trial: it is sc.trial, ending inside every limit.
        const double used = fraction * h;
        dissipated += used * sc.startDissipationPower - trial->limitWork;
        closureWork += trial->closureWork;
        projectionEnergy += trial->projectionEnergy;
        for (std::size_t r = 0; r < sc.stepRows; ++r)
        {
            ws.constraintImpulse[sc.stepLayoutRow[r]] += sc.impulse[r];
        }
        report.closureWork += trial->closureWork;
        report.limitWork += trial->limitWork;
        report.projectionEnergy += trial->projectionEnergy;
        ++report.substeps;
        std::copy(sc.trial.values.begin(), sc.trial.values.end(), sc.start.values.begin());
        if (!split)
        {
            break;
        }
        remaining -= fraction * remaining;
    }

    // Commit.
    double work = 0.0;
    const double* q0 = state.values.data();
    const double* q1 = sc.start.values.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        work += (forcing.tau.empty() ? 0.0 : forcing.tau[i]) * (q1[i] - q0[i]);
        ws.qdd[i] = (q1[n + i] - q0[n + i]) / dt;
    }
    report.minimumGap = std::numeric_limits<double>::infinity();
    if (!model.limits().empty())
    {
        if (auto gaps = limitGaps(model, sc.start.q(), ws); !gaps) // world transforms are at the end state
        {
            return gaps;
        }
        report.minimumGap = *std::ranges::min_element(sc.gap);
    }
    std::copy(sc.start.values.begin(), sc.start.values.end(), state.values.begin());
    state.dissipated += sub0math::f64::Energy::fromSi(dissipated);
    state.externalWork += sub0math::f64::Energy::fromSi(work);
    state.closureWork += sub0math::f64::Energy::fromSi(closureWork);
    state.projectionEnergy += sub0math::f64::Energy::fromSi(projectionEnergy);
    state.tick += 1U;
    return {};
}

std::expected<void, Fault> constrainAccelerations(const Model& model, const State& state, Workspace& ws) noexcept
{
    ConstraintScratch& sc = ws.constraintScratch;
    const ConstraintSettings& settings = model.constraintSettings();
    const std::size_t n = model.dofCount();
    detail::composeWorldTransforms(model, ws);
    std::ranges::fill(ws.constraintForce, 0.0);
    if (auto gaps = limitGaps(model, state.q(), ws); !gaps)
    {
        return gaps;
    }
    // A limit takes part when it is engaged (gap within eventTolerance) and not separating faster than
    // eventTolerance per second. An impact (closing speed) needs an impulse, which only a step can give.
    const auto limits = model.limits();
    for (std::size_t i = 0; i < limits.size(); ++i)
    {
        sc.engaged[i] = sc.gap[i] <= settings.eventTolerance ? 1U : 0U;
    }
    const Rows rows = assembleRows(model, ws, sc.engaged.data(), sc.jacobian.data(), sc.layoutRow.data(),
                                   sc.bias.data());
    // Drop the engaged limits that separate.
    std::size_t kept = rows.bilateral;
    const double* qd = state.values.data() + n;
    for (std::size_t r = rows.bilateral; r < rows.count; ++r)
    {
        const double* jr = sc.jacobian.data() + r * n;
        double rate = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            rate += jr[i] * qd[i];
        }
        if (rate > settings.eventTolerance)
        {
            continue;
        }
        if (kept != r)
        {
            std::copy_n(jr, n, sc.jacobian.data() + kept * n);
            sc.layoutRow[kept] = sc.layoutRow[r];
            sc.bias[kept] = sc.bias[r];
        }
        ++kept;
    }
    const Rows used{kept, rows.bilateral};
    if (used.count == 0)
    {
        return {};
    }
    delassusMatrix(n, used.count, sc.jacobian.data(), sc, [&](double* x) { detail::solveMassFactor(model, ws, x); });
    for (std::size_t r = 0; r < used.count; ++r)
    {
        const double* jr = sc.jacobian.data() + r * n;
        double rate = sc.bias[r];
        for (std::size_t i = 0; i < n; ++i)
        {
            rate += jr[i] * ws.qdd[i];
        }
        sc.rhs[r] = rate;
    }
    const auto pivots = solveMixedLcp(used, sc, settings.solverTolerance);
    if (!pivots)
    {
        return std::unexpected(pivots.error());
    }
    for (std::size_t r = 0; r < used.count; ++r)
    {
        const double* y = sc.solved.data() + r * n;
        for (std::size_t i = 0; i < n; ++i)
        {
            ws.qdd[i] += sc.impulse[r] * y[i];
        }
        ws.constraintForce[sc.layoutRow[r]] = sc.impulse[r];
    }
    return {};
}

} // namespace detail

std::expected<void, Fault> projectConstraints(const Model& model, State& state, Workspace& workspace) noexcept
{
    if (auto checked = detail::checkInputs(model, state, workspace); !checked)
    {
        return checked;
    }
    workspace.constraintReport = {};
    if (model.constraintCount() == 0)
    {
        return {};
    }
    if (auto mass = computeMassMatrix(model, state, workspace); !mass)
    {
        return mass;
    }
    if (auto factored = factorMassMatrix(model, workspace); !factored)
    {
        return factored;
    }
    ConstraintScratch& sc = workspace.constraintScratch;
    std::ranges::fill(sc.engaged, std::uint8_t{1});
    sc.trial = state;
    const auto energy = projectMeasured(model, sc.trial, workspace, Metric{true, false}, sc.engaged.data());
    if (!energy)
    {
        return std::unexpected(energy.error());
    }
    std::copy(sc.trial.values.begin(), sc.trial.values.end(), state.values.begin());
    state.projectionEnergy += sub0math::f64::Energy::fromSi(*energy);
    workspace.constraintReport.projectionEnergy = *energy;
    return {};
}

std::expected<void, Fault> computeReactions(const Model& model, const State& state, std::span<const double> tau,
                                            Workspace& workspace) noexcept
{
    if (auto accelerated = forwardDynamics(model, state, tau, workspace); !accelerated)
    {
        return accelerated;
    }
    Workspace& ws = workspace;
    ConstraintScratch& sc = ws.constraintScratch;
    const std::size_t n = model.dofCount();
    detail::composeWorldTransforms(model, ws);
    std::ranges::fill(sc.bodyExternal, sub0math::SpatialForced{});
    // A world wrench (force through world point p, plus a moment) on a body, in its frame about its origin.
    const auto addWrench = [&](std::size_t body, const Vec3d& p, const Vec3d& force, const Vec3d& moment) {
        if (body == NoParent)
        {
            return;
        }
        const sub0math::PluckerTransformd& x = ws.worldToBody[body];
        sc.bodyExternal[body] += sub0math::SpatialForced{x.rotation * (moment + cross(p - x.translation, force)),
                                                         x.rotation * force};
    };
    for (const ElementSiteWrench& w : ws.elementSiteWrenches)
    {
        const SitePoint p = sitePoint(model, ws, w.site);
        addWrench(p.body, p.position, w.force, w.moment);
    }
    for (const LoopClosure& closure : model.closures())
    {
        const double* lambda = ws.constraintForce.data() + closure.firstRow;
        const Vec3d force{lambda[0], lambda[1], lambda[2]};
        const Vec3d moment = closure.kind == ClosureKind::Frame ? Vec3d{lambda[3], lambda[4], lambda[5]} : Vec3d{};
        const SitePoint p1 = sitePoint(model, ws, closure.first);
        const SitePoint p2 = sitePoint(model, ws, closure.second);
        addWrench(p1.body, p1.position, force, moment);
        addWrench(p2.body, p2.position, force * -1.0, moment * -1.0);
        ws.constraintSiteWrenches[closure.wrenchSlot] = {closure.id, closure.first, force, moment};
        ws.constraintSiteWrenches[closure.wrenchSlot + 1] = {closure.id, closure.second, force * -1.0, moment * -1.0};
    }
    for (const UnilateralLimit& limit : model.limits())
    {
        if (limit.kind != LimitKind::Distance)
        {
            continue; // a coordinate limit acts inside its joint, like the applied force on it
        }
        const Separation s = separation(model, ws, limit);
        const Vec3d force = s.direction * (limit.sign * ws.constraintForce[limit.row]);
        addWrench(s.first.body, s.first.position, force, {});
        addWrench(s.second.body, s.second.position, force * -1.0, {});
        ws.constraintSiteWrenches[limit.wrenchSlot] = {limit.id, limit.first, force, {}};
        ws.constraintSiteWrenches[limit.wrenchSlot + 1] = {limit.id, limit.second, force * -1.0, {}};
    }
    detail::transmittedForces(model, state, ws.qdd.data(), sc.velocity.data(), sc.bodyExternal.data(), ws);
    for (std::size_t i = 0; i < n; ++i)
    {
        const Mat3d& e = ws.worldToBody[i].rotation;
        ws.jointReactions[i] = {detail::bodyToWorld(e, ws.bodyForce[i].linear),
                                detail::bodyToWorld(e, ws.bodyForce[i].angular)};
    }
    return {};
}

} // namespace sub0physics::multibody
