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

} // namespace sub0physics::multibody::detail
