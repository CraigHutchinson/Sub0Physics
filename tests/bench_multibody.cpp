// Multibody dynamics on a fixed 10-DOF planar tree (testmodels::planarTree), each paired with the fastest
// hand-written code we know for the same job: planar spatial algebra (three components instead of six),
// the tree's parent array as a compile-time constant, and the joint types hard-coded. The library runs
// the general 3D formulation, by design (docs/MULTIBODY.md section 3: no separate 2D engine), so the gap
// is the price of generality and is tracked in docs/PERFORMANCE.md.

#include <algorithm>
#include <array>
#include <bench_runner.hpp>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"

namespace
{

using testmodels::PlanarTreeDof;
using testmodels::PlanarTreeParents;
constexpr std::size_t N = PlanarTreeDof;
constexpr std::size_t None = sub0physics::multibody::NoParent;

// Planar spatial vectors: motion (w, vx, vy), force (nz, fx, fy).
struct P3
{
    double w;
    double x;
    double y;
};

// Planar inertia about a frame origin: mass, first moment, rotational inertia about z.
struct PlanarInertia
{
    double m;
    double hx;
    double hy;
    double izz;
};

// Parent-to-body transform: rotation by -theta (c, s) after a shift by r.
struct PlanarX
{
    double c;
    double s;
    double rx;
    double ry;
};

inline P3 applyMotion(const PlanarX& x, const P3& m)
{
    const double ux = m.x - x.ry * m.w;
    const double uy = m.y + x.rx * m.w;
    return {m.w, x.c * ux + x.s * uy, -x.s * ux + x.c * uy};
}

inline P3 applyForceTranspose(const PlanarX& x, const P3& f) // child to parent
{
    const double fx = x.c * f.x - x.s * f.y;
    const double fy = x.s * f.x + x.c * f.y;
    return {f.w + x.rx * fy - x.ry * fx, fx, fy};
}

inline PlanarInertia applyInertiaTranspose(const PlanarX& x, const PlanarInertia& i)
{
    const double yx = x.c * i.hx - x.s * i.hy;
    const double yy = x.s * i.hx + x.c * i.hy;
    return {i.m, yx + i.m * x.rx, yy + i.m * x.ry,
            i.izz + 2.0 * (x.rx * yx + x.ry * yy) + i.m * (x.rx * x.rx + x.ry * x.ry)};
}

inline P3 times(const PlanarInertia& i, const P3& m)
{
    return {i.izz * m.w + i.hx * m.y - i.hy * m.x, i.m * m.x - i.hy * m.w, i.m * m.y + i.hx * m.w};
}

inline P3 crossMotion(const P3& v, const P3& m)
{
    return {0.0, -v.w * m.y + v.y * m.w, v.w * m.x - v.x * m.w};
}

inline P3 crossForce(const P3& v, const P3& f)
{
    return {v.x * f.y - v.y * f.x, -v.w * f.y, v.w * f.x};
}

// The tree of testmodels::planarTree(), hand-coded: body 0 slides along x, the rest turn about z.
struct HandPlanarTree
{
    std::array<PlanarInertia, N> inertia{};
    std::array<double, N> placementX{};
    std::array<double, N> placementY{};
    double gravity = 9.81;

    std::array<PlanarX, N> xup{};
    std::array<P3, N> v{};
    std::array<P3, N> a{};
    std::array<P3, N> f{};
    std::array<PlanarInertia, N> ic{};
    std::array<double, N * N> h{};
    std::array<double, N> qdd{};

    explicit HandPlanarTree(const sub0physics::multibody::Model& model)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            const auto& in = model.inertias()[i];
            inertia[i] = {in.mass, in.firstMoment.x, in.firstMoment.y, in.rotational(2, 2)};
            placementX[i] = model.treeTransforms()[i].translation.x;
            placementY[i] = model.treeTransforms()[i].translation.y;
        }
    }

    static double motionOf(std::size_t i, const P3& f) // S_i . f
    {
        return i == 0 ? f.x : f.w;
    }

    void transforms(const double* q)
    {
        xup[0] = {1.0, 0.0, q[0], 0.0};
        for (std::size_t i = 1; i < N; ++i)
        {
            xup[i] = {std::cos(q[i]), std::sin(q[i]), placementX[i], placementY[i]};
        }
    }

    void massMatrix()
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            ic[i] = inertia[i];
        }
        for (std::size_t i = N; i-- > 1;)
        {
            const PlanarInertia t = applyInertiaTranspose(xup[i], ic[i]);
            PlanarInertia& p = ic[PlanarTreeParents[i]];
            p = {p.m + t.m, p.hx + t.hx, p.hy + t.hy, p.izz + t.izz};
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            P3 force = times(ic[i], i == 0 ? P3{0.0, 1.0, 0.0} : P3{1.0, 0.0, 0.0});
            h[i * N + i] = motionOf(i, force);
            for (std::size_t j = i; PlanarTreeParents[j] != None;)
            {
                force = applyForceTranspose(xup[j], force);
                j = PlanarTreeParents[j];
                h[i * N + j] = h[j * N + i] = motionOf(j, force);
            }
        }
    }

    void factorAndSolve(double* x)
    {
        for (std::size_t k = N; k-- > 0;)
        {
            for (std::size_t i = PlanarTreeParents[k]; i != None; i = PlanarTreeParents[i])
            {
                const double s = h[k * N + i] / h[k * N + k];
                for (std::size_t j = i; j != None; j = PlanarTreeParents[j])
                {
                    h[i * N + j] -= s * h[k * N + j];
                }
                h[k * N + i] = s;
            }
        }
        for (std::size_t i = N; i-- > 0;)
        {
            for (std::size_t j = PlanarTreeParents[i]; j != None; j = PlanarTreeParents[j])
            {
                x[j] -= h[i * N + j] * x[i];
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            x[i] /= h[i * N + i];
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t j = PlanarTreeParents[i]; j != None; j = PlanarTreeParents[j])
            {
                x[i] -= h[i * N + j] * x[j];
            }
        }
    }

    void forwardDynamics(const double* q, const double* qd, const double* tau)
    {
        transforms(q);
        for (std::size_t i = 0; i < N; ++i)
        {
            const P3 vj = i == 0 ? P3{0.0, qd[0], 0.0} : P3{qd[i], 0.0, 0.0};
            const std::size_t p = PlanarTreeParents[i];
            const P3 vp = p == None ? P3{0.0, 0.0, 0.0} : v[p];
            const P3 ap = p == None ? P3{0.0, 0.0, gravity} : a[p];
            const P3 vt = applyMotion(xup[i], vp);
            v[i] = {vt.w + vj.w, vt.x + vj.x, vt.y + vj.y};
            const P3 at = applyMotion(xup[i], ap);
            const P3 ac = crossMotion(v[i], vj);
            a[i] = {at.w + ac.w, at.x + ac.x, at.y + ac.y};
            const P3 ia = times(inertia[i], a[i]);
            const P3 gyro = crossForce(v[i], times(inertia[i], v[i]));
            f[i] = {ia.w + gyro.w, ia.x + gyro.x, ia.y + gyro.y};
        }
        for (std::size_t i = N; i-- > 0;)
        {
            qdd[i] = tau[i] - motionOf(i, f[i]);
            if (const std::size_t p = PlanarTreeParents[i]; p != None)
            {
                const P3 t = applyForceTranspose(xup[i], f[i]);
                f[p] = {f[p].w + t.w, f[p].x + t.x, f[p].y + t.y};
            }
        }
        massMatrix();
        factorAndSolve(qdd.data());
    }
};

} // namespace

void benchMultibody(sub0bench::BenchRunner& runner)
{
    using namespace sub0physics::multibody;
    const Model model = testmodels::planarTree();
    Workspace ws(model);
    State state = State::atRest(model);
    std::array<double, N> tau{};
    for (std::size_t i = 0; i < N; ++i)
    {
        state.q[i] = 0.3 * std::sin(static_cast<double>(i) + 1.0);
        state.qd[i] = 0.5 * std::cos(static_cast<double>(i));
        tau[i] = 0.1 * static_cast<double>(i);
    }
    HandPlanarTree hand(model);

    // The reference must agree with the library before its time means anything.
    if (!forwardDynamics(model, state, tau, ws))
    {
        std::fprintf(stderr, "forwardDynamics failed\n");
        std::exit(1);
    }
    hand.forwardDynamics(state.q.data(), state.qd.data(), tau.data());
    for (std::size_t i = 0; i < N; ++i)
    {
        if (std::abs(hand.qdd[i] - ws.qdd[i]) > 1e-9 * (1.0 + std::abs(ws.qdd[i])))
        {
            std::fprintf(stderr, "hand-written planar tree disagrees at %zu: %.17g vs %.17g\n", i, hand.qdd[i], ws.qdd[i]);
            std::exit(1);
        }
    }

    runner.bench().minEpochIterations(20000);
    runner.run("planar tree 10-DOF mass matrix, hand-written planar", [&] {
        ankerl::nanobench::doNotOptimizeAway(state);
        hand.transforms(state.q.data());
        hand.massMatrix();
        ankerl::nanobench::doNotOptimizeAway(hand.h);
    });
    runner.run("planar tree 10-DOF mass matrix, computeMassMatrix", [&] {
        ankerl::nanobench::doNotOptimizeAway(state);
        ankerl::nanobench::doNotOptimizeAway(computeMassMatrix(model, state, ws).has_value());
        ankerl::nanobench::doNotOptimizeAway(ws.massMatrix.data());
    });

    std::array<double, N * N> massCopy{};
    std::copy(ws.massMatrix.begin(), ws.massMatrix.end(), massCopy.begin());
    std::array<double, N> x{};
    runner.run("planar tree 10-DOF LTDL factor and solve, hand-written", [&] {
        hand.h = massCopy;
        x = tau;
        ankerl::nanobench::doNotOptimizeAway(hand.h);
        hand.factorAndSolve(x.data());
        ankerl::nanobench::doNotOptimizeAway(x);
    });
    runner.run("planar tree 10-DOF LTDL factor and solve, factorMassMatrix + solveMassMatrix", [&] {
        x = tau;
        ankerl::nanobench::doNotOptimizeAway(ws.massMatrix.data());
        if (factorMassMatrix(model, ws))
        {
            ankerl::nanobench::doNotOptimizeAway(solveMassMatrix(model, ws, x).has_value());
        }
        ankerl::nanobench::doNotOptimizeAway(x);
    });

    runner.run("planar tree 10-DOF forward dynamics, hand-written planar", [&] {
        ankerl::nanobench::doNotOptimizeAway(state);
        hand.forwardDynamics(state.q.data(), state.qd.data(), tau.data());
        ankerl::nanobench::doNotOptimizeAway(hand.qdd);
    });
    runner.run("planar tree 10-DOF forward dynamics, forwardDynamics", [&] {
        ankerl::nanobench::doNotOptimizeAway(state);
        if (!forwardDynamics(model, state, tau, ws))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(ws.qdd.data());
    });
}
