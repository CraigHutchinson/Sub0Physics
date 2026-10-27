// The soft-constraint solve's Delassus matrix, two forms (docs/PERFORMANCE.md, "Delassus A/B"): the uniform form (one
// full mass-matrix solve per row, then J Y) against the structured form (the condensed Schur complement, O(k^2)).
//
// Two shapes, because the forms may map differently onto hardware (Sub0Math PERFORMANCE, "Choosing between
// alternatives"):
// - **One wheel, the library's solve** (`soft solve 64-node ring, k rows, <form>`): the scalar consumer path.
// - **Many wheels at once, structure of arrays** (`Delassus W wheels x K rows SoA, <form>`): both forms written as
//   batched kernels with the wheel index innermost, so the compiler can map wheels onto SIMD lanes. The uniform form
//   is fixed-shape (dense vectors, no indexing); the structured form gathers each row's two nodes per lane and masks
//   the shared-node terms. Both kernels are checked against each other and against the library's matrix before
//   timing. Neither is a twin of the other: they are the two candidates, each written with the same care.

#include <algorithm>
#include <array>
#include <bench_runner.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "sub0physics/multibody.hpp"
#include "test_ring_models.hpp"

namespace
{

using namespace sub0physics;
using namespace sub0physics::multibody;

constexpr std::size_t Nodes = 64;

// Rows of a contact patch: k consecutive edges around node 0 (the bottom), a point on each edge, a normal near +z.
struct Rows
{
    std::vector<double> jacobian;
    std::vector<double> reference;
    std::vector<double> regularisation;
    std::vector<SoftRowSupport> supports;
};

std::array<double, 2> nodeWorld(const Model& model, const Workspace& ws, std::size_t node)
{
    const DeformableRing& ring = model.rings()[0];
    const RingWorkspace& rw = ws.rings[0];
    const PlanarTransform& host = ws.planar.worldToBody[ring.host];
    const double lx = rw.positionX[node];
    const double ly = rw.positionY[node];
    return {host.x + host.c * lx - host.s * ly, host.y + host.s * lx + host.c * ly};
}

sub0math::Vec3d toWorld(const Model& model, double i, double j)
{
    sub0math::Vec3d v{};
    v[model.planar().axisI] = i;
    v[model.planar().axisJ] = j;
    return v;
}

Rows patchRows(const Model& model, const Workspace& ws, std::size_t k)
{
    const std::size_t n = model.dofCount();
    Rows rows;
    rows.jacobian.assign(k * n, 0.0);
    rows.reference.assign(k, 0.0);
    rows.regularisation.assign(k, 0.0);
    rows.supports.resize(k);
    for (std::size_t r = 0; r < k; ++r)
    {
        const std::size_t edge = (Nodes - k / 2 + r) % Nodes;
        const auto a = nodeWorld(model, ws, edge);
        const auto b = nodeWorld(model, ws, (edge + 1) % Nodes);
        RingContactPoint p;
        p.edge = edge;
        p.weight0 = 0.3 + 0.4 * static_cast<double>(r % 3) / 2.0;
        p.weight1 = 1.0 - p.weight0;
        p.point = toWorld(model, p.weight0 * a[0] + p.weight1 * b[0], p.weight0 * a[1] + p.weight1 * b[1]);
        const double tilt = 0.05 * (static_cast<double>(r) - 0.5 * static_cast<double>(k));
        p.normal = toWorld(model, std::cos(tilt), std::sin(tilt)); // plane (i, j) = (z, x): mostly +z
        if (!ringContactJacobian(model, ws, p, std::span<double>(rows.jacobian.data() + r * n, n)))
        {
            std::fprintf(stderr, "Delassus bench: row failed\n");
            std::exit(1);
        }
        rows.supports[r] = {.ring = 0, .node0 = edge, .node1 = (edge + 1) % Nodes};
        rows.reference[r] = 2000.0 + 10.0 * static_cast<double>(r); // every row pushes (all active)
        rows.regularisation[r] = 0.3;
    }
    return rows;
}

// ---- Batched shape: W wheels, K rows each, structure of arrays (lane = wheel, innermost). ----

template<std::size_t W, std::size_t K, std::size_t T>
struct BatchedData
{
    static constexpr std::size_t N = Nodes;
    static constexpr std::size_t Dofs = T + 2 * N;
    // Per node: A^-1 (2 x 2 row-major), border B (2 x 3), Z = A^-1 B (2 x 3).
    std::vector<double> ainv = std::vector<double>(N * 4 * W);
    std::vector<double> border = std::vector<double>(N * 6 * W);
    std::vector<double> coupling = std::vector<double>(N * 6 * W);
    std::vector<double> hostMap = std::vector<double>(3 * T * W); // J_h: host motion (w, x, y) per tree coordinate
    std::vector<double> sinv = std::vector<double>(T * T * W);    // S^-1, dense
    // Rows: tree part, the four node entries and the two node indices.
    std::vector<double> rowTree = std::vector<double>(K * T * W);
    std::vector<double> rowNode = std::vector<double>(K * 4 * W);
    std::vector<std::uint32_t> rowIndex = std::vector<std::uint32_t>(K * 2 * W);
    // Scratch and results.
    std::vector<double> x = std::vector<double>(Dofs * W);
    std::vector<double> y = std::vector<double>(K * Dofs * W); // uniform: M^-1 J^T, all rows
    std::vector<double> jdense = std::vector<double>(K * Dofs * W);
    std::vector<double> tv = std::vector<double>(K * (2 * T + 4) * W); // structured: t, v, y_a, y_b
    std::vector<double> delassus = std::vector<double>(K * K * W);
};

// Uniform: Y = M^-1 J^T by the condensed full solve for every row, then D = J Y over every coordinate.
template<std::size_t W, std::size_t K, std::size_t T>
void uniformBatched(BatchedData<W, K, T>& d)
{
    constexpr std::size_t N = Nodes;
    constexpr std::size_t Dofs = T + 2 * N;
    for (std::size_t r = 0; r < K; ++r)
    {
        double* y = d.y.data() + r * Dofs * W;
        const double* j = d.jdense.data() + r * Dofs * W;
        std::array<double, 3 * W> g{};
        for (std::size_t i = 0; i < N; ++i)
        {
            const double* a = d.ainv.data() + i * 4 * W;
            const double* b = d.border.data() + i * 6 * W;
            const double* x0 = j + (T + 2 * i) * W;
            const double* x1 = x0 + W;
            double* y0 = y + (T + 2 * i) * W;
            double* y1 = y0 + W;
            for (std::size_t l = 0; l < W; ++l)
            {
                const double u = a[0 * W + l] * x0[l] + a[1 * W + l] * x1[l];
                const double v = a[2 * W + l] * x0[l] + a[3 * W + l] * x1[l];
                y0[l] = u;
                y1[l] = v;
                g[0 * W + l] += b[0 * W + l] * u + b[3 * W + l] * v;
                g[1 * W + l] += b[1 * W + l] * u + b[4 * W + l] * v;
                g[2 * W + l] += b[2 * W + l] * u + b[5 * W + l] * v;
            }
        }
        std::array<double, T * W> t{};
        for (std::size_t c = 0; c < T; ++c)
        {
            const double* h = d.hostMap.data() + 3 * c * W;
            for (std::size_t l = 0; l < W; ++l)
            {
                t[c * W + l] = j[c * W + l] - (h[0 * W + l] * g[0 * W + l] + h[1 * W + l] * g[1 * W + l] + h[2 * W + l] * g[2 * W + l]);
            }
        }
        for (std::size_t a = 0; a < T; ++a)
        {
            for (std::size_t l = 0; l < W; ++l)
            {
                double sum = 0.0;
                for (std::size_t c = 0; c < T; ++c)
                {
                    sum += d.sinv[(a * T + c) * W + l] * t[c * W + l];
                }
                y[a * W + l] = sum;
            }
        }
        std::array<double, 3 * W> hv{};
        for (std::size_t c = 0; c < T; ++c)
        {
            const double* h = d.hostMap.data() + 3 * c * W;
            for (std::size_t l = 0; l < W; ++l)
            {
                hv[0 * W + l] += h[0 * W + l] * y[c * W + l];
                hv[1 * W + l] += h[1 * W + l] * y[c * W + l];
                hv[2 * W + l] += h[2 * W + l] * y[c * W + l];
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            const double* z = d.coupling.data() + i * 6 * W;
            double* y0 = y + (T + 2 * i) * W;
            double* y1 = y0 + W;
            for (std::size_t l = 0; l < W; ++l)
            {
                y0[l] -= z[0 * W + l] * hv[0 * W + l] + z[1 * W + l] * hv[1 * W + l] + z[2 * W + l] * hv[2 * W + l];
                y1[l] -= z[3 * W + l] * hv[0 * W + l] + z[4 * W + l] * hv[1 * W + l] + z[5 * W + l] * hv[2 * W + l];
            }
        }
    }
    for (std::size_t r = 0; r < K; ++r)
    {
        for (std::size_t c = 0; c <= r; ++c)
        {
            const double* jr = d.jdense.data() + r * Dofs * W;
            const double* yc = d.y.data() + c * Dofs * W;
            std::array<double, W> sum{};
            for (std::size_t i = 0; i < Dofs; ++i)
            {
                for (std::size_t l = 0; l < W; ++l)
                {
                    sum[l] += jr[i * W + l] * yc[i * W + l];
                }
            }
            for (std::size_t l = 0; l < W; ++l)
            {
                d.delassus[(r * K + c) * W + l] = sum[l];
                d.delassus[(c * K + r) * W + l] = sum[l];
            }
        }
    }
}

// Structured: each row's two nodes gathered per lane, t = j_t - J_h^T B^T A^-1 j_w, v = S^-1 t, D = t . v + shared nodes.
template<std::size_t W, std::size_t K, std::size_t T>
void condensedBatched(BatchedData<W, K, T>& d)
{
    constexpr std::size_t Stride = 2 * T + 4;
    for (std::size_t r = 0; r < K; ++r)
    {
        double* row = d.tv.data() + r * Stride * W;
        const double* jn = d.rowNode.data() + r * 4 * W;
        const std::uint32_t* idx = d.rowIndex.data() + r * 2 * W;
        std::array<double, 3 * W> g{};
        for (std::size_t s = 0; s < 2; ++s)
        {
            for (std::size_t l = 0; l < W; ++l)
            {
                const std::size_t i = idx[s * W + l];
                const double* a = d.ainv.data() + i * 4 * W;
                const double* b = d.border.data() + i * 6 * W;
                const double x0 = jn[(2 * s) * W + l];
                const double x1 = jn[(2 * s + 1) * W + l];
                const double u = a[0 * W + l] * x0 + a[1 * W + l] * x1;
                const double v = a[2 * W + l] * x0 + a[3 * W + l] * x1;
                row[(2 * T + 2 * s) * W + l] = u;
                row[(2 * T + 2 * s + 1) * W + l] = v;
                g[0 * W + l] += b[0 * W + l] * u + b[3 * W + l] * v;
                g[1 * W + l] += b[1 * W + l] * u + b[4 * W + l] * v;
                g[2 * W + l] += b[2 * W + l] * u + b[5 * W + l] * v;
            }
        }
        for (std::size_t c = 0; c < T; ++c)
        {
            const double* h = d.hostMap.data() + 3 * c * W;
            const double* jt = d.rowTree.data() + (r * T + c) * W;
            for (std::size_t l = 0; l < W; ++l)
            {
                row[c * W + l] = jt[l] - (h[0 * W + l] * g[0 * W + l] + h[1 * W + l] * g[1 * W + l] + h[2 * W + l] * g[2 * W + l]);
            }
        }
        for (std::size_t a = 0; a < T; ++a)
        {
            for (std::size_t l = 0; l < W; ++l)
            {
                double sum = 0.0;
                for (std::size_t c = 0; c < T; ++c)
                {
                    sum += d.sinv[(a * T + c) * W + l] * row[c * W + l];
                }
                row[(T + a) * W + l] = sum;
            }
        }
    }
    for (std::size_t r = 0; r < K; ++r)
    {
        const double* rr = d.tv.data() + r * Stride * W;
        const double* jr = d.rowNode.data() + r * 4 * W;
        const std::uint32_t* ir = d.rowIndex.data() + r * 2 * W;
        for (std::size_t c = 0; c <= r; ++c)
        {
            const double* rc = d.tv.data() + c * Stride * W;
            const std::uint32_t* ic = d.rowIndex.data() + c * 2 * W;
            for (std::size_t l = 0; l < W; ++l)
            {
                double sum = 0.0;
                for (std::size_t i = 0; i < T; ++i)
                {
                    sum += rr[i * W + l] * rc[(T + i) * W + l];
                }
                for (std::size_t a = 0; a < 2; ++a)
                {
                    for (std::size_t b = 0; b < 2; ++b)
                    {
                        const double shared = jr[(2 * a) * W + l] * rc[(2 * T + 2 * b) * W + l] +
                                              jr[(2 * a + 1) * W + l] * rc[(2 * T + 2 * b + 1) * W + l];
                        sum += ir[a * W + l] == ic[b * W + l] ? shared : 0.0;
                    }
                }
                d.delassus[(r * K + c) * W + l] = sum;
                d.delassus[(c * K + r) * W + l] = sum;
            }
        }
    }
}

// Fills the batched data from the library at W different states of the ring rig, and the rows of a contact patch.
template<std::size_t W, std::size_t K, std::size_t T>
void fill(BatchedData<W, K, T>& d, const Model& model, Workspace& ws, std::vector<double>& libraryDelassus)
{
    constexpr std::size_t N = Nodes;
    constexpr std::size_t Dofs = T + 2 * N;
    const std::size_t n = model.dofCount();
    const DeformableRing& ring = model.rings()[0];
    libraryDelassus.assign(K * K * W, 0.0);
    for (std::size_t l = 0; l < W; ++l)
    {
        const State s = ringmodels::randomState(model, 100U + static_cast<unsigned>(l));
        if (!prepareLinearlyImplicitStep(model, s, ws) || !factorMassMatrix(model, ws))
        {
            std::fprintf(stderr, "Delassus bench: preparation failed\n");
            std::exit(1);
        }
        const RingWorkspace& rw = ws.rings[0];
        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t e = 0; e < 4; ++e)
            {
                d.ainv[(i * 4 + e) * W + l] = rw.mass.inversePivot[i][e];
            }
            for (std::size_t e = 0; e < 6; ++e)
            {
                d.border[(i * 6 + e) * W + l] = rw.mass.border[i][e];
                d.coupling[(i * 6 + e) * W + l] = rw.mass.coupling[i][e];
            }
        }
        for (std::size_t c = 0; c < T; ++c)
        {
            d.hostMap[(3 * c + 0) * W + l] = 0.0;
            d.hostMap[(3 * c + 1) * W + l] = 0.0;
            d.hostMap[(3 * c + 2) * W + l] = 0.0;
        }
        for (std::size_t c = 0; c < ring.hostChain.size(); ++c)
        {
            const std::size_t coordinate = model.coordinates()[ring.hostChain[c]];
            d.hostMap[(3 * coordinate + 0) * W + l] = rw.hostJacobian[c].w;
            d.hostMap[(3 * coordinate + 1) * W + l] = rw.hostJacobian[c].x;
            d.hostMap[(3 * coordinate + 2) * W + l] = rw.hostJacobian[c].y;
        }
        // S^-1 from tree-only right-hand sides (their ring rows eliminate to zero).
        for (std::size_t c = 0; c < T; ++c)
        {
            std::vector<double> e(n, 0.0);
            e[c] = 1.0;
            (void)solveMassMatrix(model, ws, e);
            for (std::size_t a = 0; a < T; ++a)
            {
                d.sinv[(a * T + c) * W + l] = e[a];
            }
        }
        const Rows rows = patchRows(model, ws, K);
        for (std::size_t r = 0; r < K; ++r)
        {
            const double* jr = rows.jacobian.data() + r * n;
            for (std::size_t i = 0; i < Dofs; ++i)
            {
                d.jdense[(r * Dofs + i) * W + l] = jr[i];
            }
            for (std::size_t c = 0; c < T; ++c)
            {
                d.rowTree[(r * T + c) * W + l] = jr[c];
            }
            const std::array<std::size_t, 2> nodes{rows.supports[r].node0, rows.supports[r].node1};
            for (std::size_t s2 = 0; s2 < 2; ++s2)
            {
                d.rowIndex[(r * 2 + s2) * W + l] = static_cast<std::uint32_t>(nodes[s2]);
                d.rowNode[(r * 4 + 2 * s2) * W + l] = jr[ring.firstCoordinate + 2 * nodes[s2]];
                d.rowNode[(r * 4 + 2 * s2 + 1) * W + l] = jr[ring.firstCoordinate + 2 * nodes[s2] + 1];
            }
        }
        // The library's matrix for this lane: J M^-1 J^T column by column.
        for (std::size_t c = 0; c < K; ++c)
        {
            std::vector<double> y(rows.jacobian.begin() + static_cast<std::ptrdiff_t>(c * n),
                                  rows.jacobian.begin() + static_cast<std::ptrdiff_t>((c + 1) * n));
            (void)solveMassMatrix(model, ws, y);
            for (std::size_t r = 0; r < K; ++r)
            {
                double sum = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                {
                    sum += rows.jacobian[r * n + i] * y[i];
                }
                libraryDelassus[(r * K + c) * W + l] = sum;
            }
        }
    }
}

template<std::size_t W, std::size_t K>
void benchBatched(sub0bench::BenchRunner& runner, const Model& model, Workspace& ws)
{
    constexpr std::size_t T = ringmodels::TreeDof;
    auto d = std::make_unique<BatchedData<W, K, T>>();
    std::vector<double> library;
    fill(*d, model, ws, library);
    uniformBatched(*d);
    const std::vector<double> uniform = d->delassus;
    condensedBatched(*d);
    for (std::size_t i = 0; i < library.size(); ++i)
    {
        const double scale = 1.0 + std::abs(library[i]);
        if (std::abs(uniform[i] - library[i]) > 1e-10 * scale || std::abs(d->delassus[i] - library[i]) > 1e-10 * scale)
        {
            std::fprintf(stderr, "Delassus bench: batched kernels disagree at %zu: library %.17g, uniform %.17g, condensed %.17g\n",
                         i, library[i], uniform[i], d->delassus[i]);
            std::exit(1);
        }
    }
    const std::string shape = "Delassus " + std::to_string(W) + " wheels x " + std::to_string(K) + " rows SoA, ";
    runner.bench().minEpochIterations(W >= 64 ? 50 : 500);
    runner.run((shape + "uniform (full solve per row)").c_str(), [&] {
        uniformBatched(*d);
        ankerl::nanobench::doNotOptimizeAway(d->delassus.data());
    });
    runner.run((shape + "structured (Schur, O(k^2))").c_str(), [&] {
        condensedBatched(*d);
        ankerl::nanobench::doNotOptimizeAway(d->delassus.data());
    });
}

} // namespace

void benchDelassus(sub0bench::BenchRunner& runner)
{
    const Model model = ringmodels::ringRig({.nodes = Nodes}, false, false);
    Workspace ws(model);
    const std::size_t n = model.dofCount();
    const State s = ringmodels::randomState(model, 3);
    if (!prepareLinearlyImplicitStep(model, s, ws) || !factorMassMatrix(model, ws))
    {
        std::fprintf(stderr, "Delassus bench: preparation failed\n");
        std::exit(1);
    }
    std::vector<double> a0(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        a0[i] = 0.01 * static_cast<double>(i % 7) - 0.03;
    }
    SoftConstraintSolver solver(model, 32);
    for (const std::size_t k : {std::size_t{4}, std::size_t{8}, std::size_t{16}})
    {
        const Rows rows = patchRows(model, ws, k);
        std::vector<double> force(k);
        std::vector<double> forceCondensed(k);
        std::vector<double> generalized(n);
        // The forms must agree before their times mean anything.
        if (!solver.solve(model, ws, rows.jacobian, rows.reference, rows.regularisation, a0, {}, force, generalized,
                          rows.supports) ||
            !solver.solve(model, ws, rows.jacobian, rows.reference, rows.regularisation, a0, {}, forceCondensed,
                          generalized, rows.supports, DelassusForm::Condensed))
        {
            std::fprintf(stderr, "Delassus bench: solve failed\n");
            std::exit(1);
        }
        for (std::size_t r = 0; r < k; ++r)
        {
            if (!(force[0] > 0.0) || std::abs(force[r] - forceCondensed[r]) > 1e-10 * (1.0 + std::abs(force[r])))
            {
                std::fprintf(stderr, "Delassus bench: forms disagree at row %zu: %.17g against %.17g\n", r, force[r],
                             forceCondensed[r]);
                std::exit(1);
            }
        }
        const std::string job = "soft solve 64-node ring, " + std::to_string(k) + " rows, ";
        runner.bench().minEpochIterations(2000);
        runner.run((job + "uniform (full solve per row)").c_str(), [&] {
            ankerl::nanobench::doNotOptimizeAway(solver.solve(model, ws, rows.jacobian, rows.reference,
                                                              rows.regularisation, a0, {}, force, generalized,
                                                              rows.supports));
        });
        runner.run((job + "structured (Schur, O(k^2))").c_str(), [&] {
            ankerl::nanobench::doNotOptimizeAway(solver.solve(model, ws, rows.jacobian, rows.reference,
                                                              rows.regularisation, a0, {}, force, generalized,
                                                              rows.supports, DelassusForm::Condensed));
        });
    }
    benchBatched<1, 8>(runner, model, ws);
    benchBatched<8, 8>(runner, model, ws);
    benchBatched<64, 8>(runner, model, ws);
}
