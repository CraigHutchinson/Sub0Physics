#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/workspace.hpp"

// The reference engine's soft-constraint law (docs/WHEEL.md section 3.2, "Soft-constraint emulation"; step M2c), written
// once for every row kind: a scalar row with a position (a contact's gap, an equality's violation), its rate and a
// diagonal approximation of `J M^-1 J^T` gives a reference acceleration and a regularisation, and one solve finds the
// forces of all rows at the start of a step. Transcribed from MuJoCo 3.16.0 (commit a892750); the derivation, the
// source lines and how it combines with the implicitfast step are in docs/WHEEL.md, "M2c status".

namespace sub0physics::multibody
{

/**
 * Parameters of the reference soft-constraint law (MuJoCo's `solref` in the standard format and `solimp`).
 */
struct SoftConstraintParameters
{
    double timeConstant = 0.02;  ///< solref[0] > 0 [s]: the reference's time constant (before the 2h safety bound).
    double dampingRatio = 1.0;   ///< solref[1] > 0.
    std::array<double, 5> impedance{0.9, 0.95, 0.001, 0.5, 2.0}; ///< solimp: d0, dwidth, width [m], midpoint, power.
};

/** One regularised row of the soft-constraint law at a state. */
struct SoftConstraintRow
{
    double referenceAcceleration = 0.0; ///< aref = -B v - K d(g) g [m/s^2].
    double regularisation = 0.0;        ///< R = (1 - d) / d * diagApprox [1/kg].
    double impedance = 0.0;             ///< d(g), in (0, 1).
    double stiffness = 0.0;             ///< K [1/s^2].
    double damping = 0.0;               ///< B [1/s].
};

/**
 * Evaluates the reference engine's soft-constraint law for one scalar row of any kind (MuJoCo 3.16.0, standard solref
 * format, a non-discrete integrator): a frictionless contact (`position` the signed gap), and equally a row of an
 * equality or a limit (`position` its violation), which differ only in their bounds (SoftRowKind). The time constant is first raised to at least `2 h`
 * (engine_core_constraint.c:1777-1780), the impedance `d(g)` follows solimp (getimpedance, :1851-1900), `K = 1 /
 * (dwidth^2 tc^2 dr^2)`, `B = 2 / (dwidth tc)` (mj_makeImpedance, :1947, :1956), `R = max(mjMINVAL, (1 - d) diagApprox
 * / d)` (:1930) and `aref = -B v - K d g` (mj_referenceConstraint, :3194-3196, the margin zero). Expression order
 * follows the source.
 *
 * @param parameters The law's parameters.
 * @param gap The row's position: the signed gap g of a contact [m], negative in overlap, or a row's violation.
 * @param rate Its rate v = J qd.
 * @param diagApprox The row's diagonal approximation of `J M^-1 J^T` [1/kg] (weighted node inverse weights).
 * @param step The integration step h [s].
 * @return The row.
 */
[[nodiscard]] SoftConstraintRow softConstraintRow(const SoftConstraintParameters& parameters, double gap, double rate,
                                                  double diagApprox, double step) noexcept;

/**
 * The tangential row of a contact on the reference's elliptic friction cone, from its normal row: the same law with no
 * position term. The impedance and damping are the normal row's, `K = 0` (mj_makeImpedance, :1940-1943), `R_t = R_n /
 * impratio` (:2012) and `aref = -B v` (mj_referenceConstraint, :3195, the position zero).
 *
 * @param normal The contact's normal row (softConstraintRow).
 * @param rate The tangential rate v = J_t qd [m/s].
 * @param impratio The reference's `opt.impratio` (its default is one).
 * @return The row.
 */
[[nodiscard]] SoftConstraintRow softConstraintTangentRow(const SoftConstraintRow& normal, double rate,
                                                         double impratio = 1.0) noexcept;

/**
 * A degree of freedom's friction-loss row (MuJoCo `dof_frictionloss`; the reference's brake, step M3c): the same law with
 * no position, from the dof's own `solref` and `solimp` (`solreffriction`, `solimpfriction`; engine_core_constraint.c
 * :1746-1749). Its position is zero, so the impedance is `solimp[0]` (getimpedance, :1866-1869), `K = 0` (mj_makeImpedance,
 * :1941-1943), `B = 2 / (dwidth tc)` with the time constant raised to at least `2 h`, `R = max(mjMINVAL, (1 - d)
 * diagApprox / d)` (:1930) and `aref = -B v` (mj_referenceConstraint, :3194-3196). The row's force is bounded by the
 * dof's capacity in the solve (SoftRowKind::FrictionLoss).
 *
 * @param parameters The dof's law (the reference's defaults are SoftConstraintParameters{}).
 * @param rate The dof's rate v = J qd (the joint's own rate) [rad/s or m/s].
 * @param diagApprox The dof's `dof_invweight0` (dofInverseWeight).
 * @param step The integration step h [s].
 * @return The row.
 */
[[nodiscard]] SoftConstraintRow softConstraintFrictionLossRow(const SoftConstraintParameters& parameters, double rate,
                                                              double diagApprox, double step) noexcept;

/**
 * A degree of freedom's inverse weight, MuJoCo's `dof_invweight0` for a hinge or slide joint (engine_setconst.c
 * :1261-1300): `(M^-1)_ii`, the diagonal of the inverse mass matrix, at the state the workspace's mass matrix was
 * factorised for (the reference evaluates it at the model's reference configuration).
 *
 * @param model The model.
 * @param workspace A workspace whose mass matrix is factorised (factorMassMatrix).
 * @param dof A coordinate.
 * @return The weight [1/kg or 1/(kg m^2)], or InvalidInput (another model, a coordinate out of range). Allocates scratch:
 *         a setup-time query.
 */
[[nodiscard]] std::expected<double, Fault> dofInverseWeight(const Model& model, const Workspace& workspace,
                                                            std::size_t dof);

/**
 * The exact limit of a tangent row (the improved model's Coulomb row, docs/WHEEL.md section 3.3): no regulariser and a
 * reference that brings the row's rate to zero at the end of the step, `J a = -v / h`, so a sticking row holds with zero
 * slip through the step when solved with the step's own matrix (SoftRowSystem::PreparedStep). `R` is a compliance of
 * 1e-14 of the row's diagonal, only so that dependent rows stay factorable: it leaves a slip of `h R f` (order 1e-14 m/s).
 *
 * @param rate The row's rate v = J qd [m/s].
 * @param diagApprox A scale of the row's diagonal of `J A^-1 J^T` [1/kg] (the node inverse weights, as the soft row).
 * @param step The integration step h [s].
 * @return The row (its stiffness and impedance unused).
 */
[[nodiscard]] SoftConstraintRow exactVelocityRow(double rate, double diagApprox, double step) noexcept;

/**
 * How a soft-constraint row is bounded (MuJoCo mj_constraintUpdate_impl, engine_core_constraint.c:3240-3408). Every kind is
 * a bound on one variable of the same problem, so all are solved together in the one solve.
 */
enum class SoftRowKind : std::uint8_t
{
    Unilateral,  ///< A contact or limit row: force `f >= 0`, zero while `J a - aref >= 0` (the row is satisfied).
    Equality,    ///< An equality row: always quadratic, force of either sign.
    ConeNormal,  ///< A contact's normal row on a friction cone; the next row is its ConeTangent.
    ConeTangent, ///< The tangential row after its ConeNormal, on the same support: `|f_t| <= mu f_n`.
    /**
     * A friction-loss row (a dof's dry friction, the reference's brake): `|f| <= b` with `b` given (`friction`), quadratic
     * inside the box and linear outside (:3257-3283): the force is `b` while `J a - aref <= -R b`, `-b` while `J a - aref
     * >= R b`, and `-(J a - aref) / R` between.
     */
    FrictionLoss
};

/**
 * How a soft-constraint solve forms its Delassus matrix `J M^-1 J^T`: the same matrix, two formulations that do
 * different work and map differently onto hardware. Both are kept behind the solve (docs/PERFORMANCE.md, "Delassus
 * A/B"; Sub0Math PERFORMANCE "Choosing between alternatives").
 */
enum class DelassusForm : std::uint8_t
{
    /**
     * Uniform: `Y = M^-1 J^T` by one full mass-matrix solve per row (through the ring condensation, O(n) each), then
     * `J Y`. Fixed shapes, any rows; the reference path's bits. With row supports the products `J Y`, `J a0` and
     * `J^T f` read only each row's support and all rows are solved in one batched sweep (bit for bit the per-row form).
     */
    MassSolvePerRow,
    /**
     * Structured: from the condensed Schur complement. A row on a ring edge is `j = (j_t, j_w)` (tree coordinates, the
     * edge's two nodes); with `A_i` a node's mass block, `B_i` its border and `S` the condensed tree matrix, `t = j_t -
     * J_h^T sum_i B_i^T A_i^-1 j_w,i` and `(J M^-1 J^T)_rc = t_r^T S^-1 t_c + sum_(shared nodes) j_r,i^T A_i^-1 j_c,i`:
     * one half sweep of the tree factor per row and O(k^2 tree) for the pairs, independent of the ring's node count.
     * Needs row supports (SoftRowSupport) and a model whose rows touch only the tree and those nodes. Equal to the
     * uniform form to rounding, not bit for bit.
     */
    Condensed
};

/**
 * Which matrix couples the rows of a solve.
 */
enum class SoftRowSystem : std::uint8_t
{
    MassMatrix,  ///< `J M^-1 J^T` at the start of the step (factorMassMatrix): the reference's forward pass.
    PreparedStep ///< `J A^-1 J^T` with the prepared step's own matrix (factorPreparedStepMatrix): rows hold exactly through the step.
};

/** The ring coordinates a row touches besides the tree's: one ring edge's two nodes (a ring contact's row). */
struct SoftRowSupport
{
    std::size_t ring = 0;  ///< Index in Model::rings().
    std::size_t node0 = 0; ///< First node (an edge's `e`).
    std::size_t node1 = 0; ///< Second node (`e + 1 mod N`).
};

/**
 * How a complementarity solve changes its active set: the same unique solution, reached differently (A/B; bit for bit the
 * same final solve when both end on one active set). docs/PERFORMANCE.md, "Friction rows A/B".
 */
enum class SoftPivotRule : std::uint8_t
{
    LeastIndex, ///< One row per iteration, the first that breaks complementarity (Murty): terminates on its own.
    Block       ///< Every such row per iteration (Judice and Pires), falling back to LeastIndex when the count stalls.
};

/** Outcome of a soft-constraint solve. */
struct SoftConstraintSolution
{
    std::uint32_t pivots = 0;    ///< Principal pivots taken.
    std::size_t activeRows = 0;  ///< Rows with a force; a cone counts its active edges (one sliding, two sticking), a friction-loss row counts inside its box or at a bound.
};

/**
 * Solves the reference engine's soft-constraint problem without friction at the start of a step (its forward pass,
 * with the mass matrix): minimise `1/2 (a - a0)^T M (a - a0) + sum_i s_i(J_i a - aref_i)`, where `s_i(x) = 1/2 x^2 /
 * R_i` for an equality row and, for a unilateral row, for `x < 0` and zero otherwise (mj_constraintUpdate_impl,
 * engine_core_constraint.c:3240-3300). Its unique minimiser solves the mixed linear complementarity problem with
 * matrix `J M^-1 J^T + R` and vector `J a0 - aref`: equality rows always active, unilateral rows with `f >= 0`,
 * `w >= 0`, `f w = 0`; row forces `f_i = -(J_i a - aref_i) / R_i` on the active rows. One law and one solve serve
 * every row kind (contacts now, equality closures in the reference path later). The reference reaches it with
 * Newton's method to a tolerance; here Murty's principal pivoting solves it exactly (the matrix is positive definite,
 * so the least-index rule terminates).
 *
 * **Friction loss.** A FrictionLoss row is a variable with two bounds, `-b <= f <= b`: free inside (`w = 0`), at its upper
 * bound while `w <= 0` and at its lower bound while `w >= 0`, where `w = A f + q`; that is the reference's Huber cost
 * (:3257-3283) as complementarity. The same principal pivoting moves it between its three states (Judice and Pires's box
 * form, as solveBounded), so brakes and contacts are one problem, solved together and exactly, as the reference's Newton
 * solves them together to a tolerance.
 *
 * **Friction cones.** In a plane a contact's cone `|f_t| <= mu f_n` is a wedge with two edges, so the reference's elliptic
 * cost (engine_core_constraint.c:3305-3408: a top zone with no force, a bottom zone quadratic in both rows, a middle zone
 * projecting onto the cone) is the same complementarity problem in the edges' coordinates: `f_n = g+ + g-`, `f_t = mu
 * (g+ - g-)`, `g+, g- >= 0`, with matrix `E^T (J M^-1 J^T + R) E` and vector `E^T (J a0 - aref)`. Both edges active is
 * the bottom zone, one the middle zone, none the top zone. The matrix stays positive definite, so the same pivoting
 * solves it exactly, with no iteration on the cone. (The reference's third contact row, out of the plane, has a zero
 * Jacobian here and carries no force. A cone with width needs its curved Hessian: M4.)
 *
 * Owns its scratch, sized once for a row capacity; solving never allocates.
 */
class SoftConstraintSolver
{
public:
    /**
     * Sizes the scratch.
     *
     * @param model The model it will serve.
     * @param capacity The largest number of rows a solve may take.
     */
    SoftConstraintSolver(const Model& model, std::size_t capacity);

    /**
     * Chooses the pivot rule of later solves, solve and solveBounded (the default is LeastIndex).
     *
     * @param rule The rule.
     */
    void setPivotRule(SoftPivotRule rule) noexcept
    {
        rule_ = rule;
    }

    /** @return The row capacity. */
    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return capacity_;
    }

    /**
     * Solves for the row forces and their generalized force `J^T f`.
     *
     * @param model The model.
     * @param workspace A workspace whose mass matrix is factorised at the step's start (factorMassMatrix).
     * @param jacobian `rows` rows of dofCount entries, row-major (ringContactJacobian).
     * @param reference aref per row.
     * @param regularisation R per row; positive.
     * @param freeAcceleration a0 = M^-1 (tau - c), dofCount entries.
     * @param kinds Each row's kind, or empty for all unilateral.
     * @param force Receives f per row (a unilateral row's never negative; a cone's in its normal and tangent rows) [N].
     * @param generalized Receives `J^T f`, dofCount entries.
     * @param supports Empty, or each row's ring support: then every row's entries outside the tree's coordinates and
     *        its two nodes must be zero (ringContactJacobian rows are), and the products read only the support.
     * @param form How the Delassus matrix is formed; Condensed needs supports.
     * @param friction Empty without cones or friction-loss rows, else one entry per row: mu on each ConeTangent row
     *        (finite and positive) and the bound `b` on each FrictionLoss row (finite, not negative; zero pins the force at
     *        zero); the other entries are not read.
     * @return The solution's counts, or InvalidInput (sizes, capacity, a regularisation that is not positive, a
     *         support outside the model, Condensed without supports, a cone without its pair or its mu, a friction-loss
     *         row without its bound) or NonFiniteState.
     * @param system The matrix that couples the rows; PreparedStep with the condensed form needs a block-diagonal ring
     *        block and a factor along the tree (DampingOnlyMassPattern).
     */
    [[nodiscard]] std::expected<SoftConstraintSolution, Fault>
    solve(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
          std::span<const double> reference, std::span<const double> regularisation,
          std::span<const double> freeAcceleration, std::span<const SoftRowKind> kinds, std::span<double> force,
          std::span<double> generalized, std::span<const SoftRowSupport> supports = {},
          DelassusForm form = DelassusForm::MassSolvePerRow, std::span<const double> friction = {},
          SoftRowSystem system = SoftRowSystem::MassMatrix) noexcept;

    /**
     * Solves bounded rows (docs/WHEEL.md section 10, item 2): `|f_i| <= b_i` with the bound given, not a function of the
     * solve, so a sliding row's force is exactly its bound and no normal row is lifted by it (a cone solved with hard
     * tangent rows would separate a sliding contact at `mu |v_t|`: the associated flow rule). The problem is the convex QP
     * minimise `1/2 f^T (J A^-1 J^T + R) f + (J a0 - aref)^T f` over the box, solved exactly by principal pivoting under
     * the chosen rule (setPivotRule; each row free, at its upper bound or at its lower bound; the matrix is positive
     * definite, so least index terminates; the solution is unique, so the rules agree to the pivoting's tolerance). A row
     * with a zero bound has zero force. Velocity-level only: no position term.
     *
     * @param model The model.
     * @param workspace A workspace whose chosen system is factorised.
     * @param jacobian `rows` rows of dofCount entries, row-major.
     * @param reference aref per row.
     * @param regularisation R per row; positive (exactVelocityRow).
     * @param freeAcceleration a0 under the chosen system, dofCount entries.
     * @param bounds b per row; finite and not negative [N].
     * @param force Receives f per row [N].
     * @param generalized Receives `J^T f`, dofCount entries.
     * @param supports As solve.
     * @param form As solve.
     * @param system As solve; the default, PreparedStep, holds exact stick through the step.
     * @return The solution's counts, or InvalidInput (as solve; a negative bound) or NonFiniteState.
     */
    [[nodiscard]] std::expected<SoftConstraintSolution, Fault>
    solveBounded(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                 std::span<const double> reference, std::span<const double> regularisation,
                 std::span<const double> freeAcceleration, std::span<const double> bounds, std::span<double> force,
                 std::span<double> generalized, std::span<const SoftRowSupport> supports = {},
                 DelassusForm form = DelassusForm::MassSolvePerRow,
                 SoftRowSystem system = SoftRowSystem::PreparedStep) noexcept;

private:
    /** Checks the sizes and supports every solve shares. */
    [[nodiscard]] std::expected<void, Fault> checkRows(const Model& model, const Workspace& workspace,
                                                       std::span<const double> jacobian, std::size_t k,
                                                       std::span<const double> regularisation,
                                                       std::span<const double> freeAcceleration,
                                                       std::span<const double> force, std::span<const double> generalized,
                                                       std::span<const SoftRowSupport> supports, DelassusForm form,
                                                       SoftRowSystem system) const noexcept;
    /** Forms `A = J S^-1 J^T + R` in `delassus_` and `q = J a0 - aref` in `rhs_`; returns the largest |q|. */
    [[nodiscard]] std::expected<double, Fault> formRows(const Model& model, const Workspace& workspace,
                                                        std::span<const double> jacobian,
                                                        std::span<const double> reference,
                                                        std::span<const double> regularisation,
                                                        std::span<const double> freeAcceleration,
                                                        std::span<const SoftRowSupport> supports, DelassusForm form,
                                                        SoftRowSystem system, std::size_t k) noexcept;
    /** The uniform form: Y = S^-1 J^T by full solves, then J Y (and q = J a0 - aref). */
    void formUniform(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                     std::span<const double> reference, std::span<const double> freeAcceleration,
                     std::span<const SoftRowSupport> supports, SoftRowSystem system, std::size_t k) noexcept;
    /** The structured form from the condensation (supports required). */
    void formCondensed(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                       std::span<const double> reference, std::span<const double> freeAcceleration,
                       std::span<const SoftRowSupport> supports, SoftRowSystem system, std::size_t k) noexcept;
    /** Cholesky of `delassus_` on `subset_[0..m)`, keeping the rows before the first that changed since the last call. */
    [[nodiscard]] bool factorSubset(std::size_t m, std::size_t k, std::size_t& priorCount) noexcept;
    /** Solves the factored subset in place: `work_[0..m)` is the right-hand side, then the solution. */
    void solveSubset(std::size_t m, std::size_t k) noexcept;
    /** Adds `J^T f` to `generalized`, over each row's support when given. */
    [[nodiscard]] std::expected<void, Fault> mapForces(const Model& model, std::span<const double> jacobian,
                                                       std::span<const double> force, std::span<double> generalized,
                                                       std::span<const SoftRowSupport> supports, std::size_t k) noexcept;

    std::size_t capacity_;
    std::size_t dof_;
    SoftPivotRule rule_ = SoftPivotRule::LeastIndex;
    std::vector<double> solved_;     // capacity x dof: M^-1 J^T, by row (uniform); per row the tree sweep and node terms (structured)
    std::vector<std::size_t> index_; // capacity x MaxSupport: each row's support coordinates, ascending
    std::vector<double> structured_; // capacity x (2 tree + 8): the structured form's per-row values
    std::vector<double> sweep_;      // 6 x capacity: the batched solve's host terms
    std::vector<std::uint64_t> key_; // 2 x capacity: the structured form's (ring, node) keys per row
    std::vector<std::uint8_t> touched_; // dofCount: the uniform form's ring nodes some row touches
    std::vector<double> delassus_;   // capacity^2
    std::vector<double> factor_;     // capacity^2
    std::vector<double> rhs_;        // capacity
    std::vector<double> slack_;      // capacity
    std::vector<double> work_;       // capacity
    std::vector<std::size_t> subset_; // capacity
    std::vector<std::size_t> prior_;  // capacity: the subset the factor rows were formed for
    std::vector<std::size_t> flagged_; // capacity: the rows the block rule flips
    std::vector<std::uint8_t> active_; // capacity: each row's state (inactive, free, at its upper or its lower bound)
    std::vector<std::size_t> bound_;   // capacity: the rows at a bound in this iteration (friction loss)
};

} // namespace sub0physics::multibody
