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

/** How a soft-constraint row is bounded (MuJoCo mj_constraintUpdate_impl, engine_core_constraint.c:3240-3300). */
enum class SoftRowKind : std::uint8_t
{
    Unilateral, ///< A contact or limit row: force `f >= 0`, zero while `J a - aref >= 0` (the row is satisfied).
    Equality    ///< An equality row: always quadratic, force of either sign.
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

/** The ring coordinates a row touches besides the tree's: one ring edge's two nodes (a ring contact's row). */
struct SoftRowSupport
{
    std::size_t ring = 0;  ///< Index in Model::rings().
    std::size_t node0 = 0; ///< First node (an edge's `e`).
    std::size_t node1 = 0; ///< Second node (`e + 1 mod N`).
};

/** Outcome of a soft-constraint solve. */
struct SoftConstraintSolution
{
    std::uint32_t pivots = 0;    ///< Principal pivots taken.
    std::size_t activeRows = 0;  ///< Rows with a positive force.
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
     * @param force Receives f per row (a unilateral row's never negative) [N].
     * @param generalized Receives `J^T f`, dofCount entries.
     * @param supports Empty, or each row's ring support: then every row's entries outside the tree's coordinates and
     *        its two nodes must be zero (ringContactJacobian rows are), and the products read only the support.
     * @param form How the Delassus matrix is formed; Condensed needs supports.
     * @return The solution's counts, or InvalidInput (sizes, capacity, a regularisation that is not positive, a
     *         support outside the model, Condensed without supports) or NonFiniteState.
     */
    [[nodiscard]] std::expected<SoftConstraintSolution, Fault>
    solve(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
          std::span<const double> reference, std::span<const double> regularisation,
          std::span<const double> freeAcceleration, std::span<const SoftRowKind> kinds, std::span<double> force,
          std::span<double> generalized, std::span<const SoftRowSupport> supports = {},
          DelassusForm form = DelassusForm::MassSolvePerRow) noexcept;

private:
    /** The uniform form: Y = M^-1 J^T by full solves, then J Y (and q = J a0 - aref). */
    void formUniform(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                     std::span<const double> reference, std::span<const double> freeAcceleration,
                     std::span<const SoftRowSupport> supports, std::size_t k) noexcept;
    /** The structured form from the condensation (supports required). */
    void formCondensed(const Model& model, const Workspace& workspace, std::span<const double> jacobian,
                       std::span<const double> reference, std::span<const double> freeAcceleration,
                       std::span<const SoftRowSupport> supports, std::size_t k) noexcept;

    std::size_t capacity_;
    std::size_t dof_;
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
    std::vector<std::uint8_t> active_; // capacity
};

} // namespace sub0physics::multibody
