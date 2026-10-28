#pragma once

#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/workspace.hpp"

// The tree part of the mass matrix's factor (factorMassMatrix), for solvers that work on the condensation directly
// (the soft-constraint solve's structured Delassus form). Internal.

namespace sub0physics::multibody::detail
{

// Solves the condensed tree matrix S x = b in place (treeDofCount entries) with the factor in Workspace::massFactor.
void solveMassTree(const Model& model, const Workspace& ws, double* x) noexcept;

// The first half of that solve, L^T y = b in place: then `b_r^T S^-1 b_c = sum_i y_r[i] y_c[i] / d_i` with d the
// factor's diagonal (Workspace::massFactor[i * n + i]).
void sweepMassTree(const Model& model, const Workspace& ws, double* x) noexcept;

// Solves M X = B for `columns` right-hand sides at once (X row-major dofCount x columns: coordinate i's values are
// contiguous), with the factor of factorMassMatrix: for each column the same operations in the same order as
// solveMassMatrix, so each column is bit for bit a single solve; the work is laid out so that the columns can share
// SIMD lanes. `scratch` holds 6 * columns values. `skipZeroRingRows`: the ring elimination skips nodes whose rows are
// zero in every column (their elimination adds exact zeros: still bit for bit). `touchedRows` (one flag per
// coordinate, a node's first) names the nodes that may be nonzero, replacing the scan; untouched nodes are not read.
void solveMassMatrixColumns(const Model& model, const Workspace& ws, double* x, std::size_t columns, double* scratch,
                            bool skipZeroRingRows, const std::uint8_t* touchedRows = nullptr) noexcept;

} // namespace sub0physics::multibody::detail
