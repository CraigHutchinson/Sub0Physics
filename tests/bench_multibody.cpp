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
#include <cstdint>
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


// The twin of stepLinearlyImplicitEuler on testmodels::planarTreeWithElements(): the same step on the same
// system, hand-written. Planar algebra and the compile-time tree from HandPlanarTree; the two elements
// written out (a cart spring-damper, and a preloaded strut between the tips of branches 0 and 2, in planar
// geometry with its exact Jacobians: the cart moves both tips equally, so it drops out of the strut); the
// implicit matrix factorised by a dense LDL^T of compile-time size. It produces the same state, dissipated
// energy and external work.
struct HandImplicitTree : HandPlanarTree
{
    testmodels::PlanarTreeElements e{};
    static constexpr std::array<std::size_t, 3> ChainA{1, 2, 3}; // branch 0, root to tip
    static constexpr std::array<std::size_t, 3> ChainB{7, 8, 9}; // branch 2
    static constexpr std::size_t M = 6;                          // strut support, cart excluded
    static constexpr std::array<std::size_t, M> Support{1, 2, 3, 7, 8, 9};

    std::array<double, N> c{};      // bias forces
    std::array<double, N> force{};  // element forces
    std::array<double, N * N> k{};  // dtau/dq
    std::array<double, N * N> d{};  // dtau/dqd
    std::array<double, N * N> am{}; // implicit matrix, then its factor
    std::array<double, N> b{};
    std::array<double, N> angle{};
    std::array<double, N> ox{};
    std::array<double, N> oy{};
    double dissipated = 0.0;
    double externalWork = 0.0;

    explicit HandImplicitTree(const sub0physics::multibody::Model& model) : HandPlanarTree(model)
    {
    }

    void biasForces(const double* qd)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            const P3 vj = i == 0 ? P3{0.0, qd[0], 0.0} : P3{qd[i], 0.0, 0.0};
            const std::size_t p = PlanarTreeParents[i];
            const P3 vp = p == None ? P3{0.0, 0.0, 0.0} : v[p];
            const P3 ap = p == None ? P3{0.0, 0.0, gravity} : a3[p];
            const P3 vt = applyMotion(xup[i], vp);
            v[i] = {vt.w + vj.w, vt.x + vj.x, vt.y + vj.y};
            const P3 at = applyMotion(xup[i], ap);
            const P3 ac = crossMotion(v[i], vj);
            a3[i] = {at.w + ac.w, at.x + ac.x, at.y + ac.y};
            const P3 ia = times(inertia[i], a3[i]);
            const P3 gyro = crossForce(v[i], times(inertia[i], v[i]));
            f[i] = {ia.w + gyro.w, ia.x + gyro.x, ia.y + gyro.y};
        }
        for (std::size_t i = N; i-- > 0;)
        {
            c[i] = motionOf(i, f[i]);
            if (const std::size_t p = PlanarTreeParents[i]; p != None)
            {
                const P3 t = applyForceTranspose(xup[i], f[i]);
                f[p] = {f[p].w + t.w, f[p].x + t.x, f[p].y + t.y};
            }
        }
    }

    std::array<P3, N> a3{}; // body accelerations at zero qdd

    // World angle and origin of each body.
    void worldPoses(const double* q)
    {
        angle[0] = 0.0;
        ox[0] = q[0];
        oy[0] = 0.0;
        for (std::size_t i = 1; i < N; ++i)
        {
            const std::size_t p = PlanarTreeParents[i];
            const double cp = std::cos(angle[p]);
            const double sp = std::sin(angle[p]);
            ox[i] = ox[p] + cp * placementX[i] - sp * placementY[i];
            oy[i] = oy[p] + sp * placementX[i] + cp * placementY[i];
            angle[i] = angle[p] + q[i];
        }
    }

    void elements(const double* q, const double* qd, double& power)
    {
        force = {};
        k = {};
        d = {};
        // Cart spring-damper on q0.
        const double x0 = q[0] - e.cartReference;
        force[0] = -e.cartStiffness * x0 - e.cartDamping * qd[0];
        k[0] = -e.cartStiffness;
        d[0] = -e.cartDamping;
        power = e.cartDamping * qd[0] * qd[0];

        // Strut between the tips (0, -L) of bodies 3 and 9.
        const double l3 = testmodels::PlanarTreeLink;
        const double p1x = ox[3] + l3 * std::sin(angle[3]);
        const double p1y = oy[3] - l3 * std::cos(angle[3]);
        const double p2x = ox[9] + l3 * std::sin(angle[9]);
        const double p2y = oy[9] - l3 * std::cos(angle[9]);
        const double dx = p1x - p2x;
        const double dy = p1y - p2y;
        const double l = std::sqrt(dx * dx + dy * dy);
        const double ux = dx / l;
        const double uy = dy / l;
        // Columns z x (p - o): (-(py - oy), px - ox); +J1 on branch 0, -J2 on branch 2.
        std::array<double, M> jx{};
        std::array<double, M> jy{};
        std::array<double, M> jl{};
        double rate = 0.0;
        for (std::size_t m = 0; m < 3; ++m)
        {
            const std::size_t i = ChainA[m];
            jx[m] = -(p1y - oy[i]);
            jy[m] = p1x - ox[i];
            const std::size_t j = ChainB[m];
            jx[m + 3] = p2y - oy[j];
            jy[m + 3] = -(p2x - ox[j]);
        }
        for (std::size_t m = 0; m < M; ++m)
        {
            jl[m] = ux * jx[m] + uy * jy[m];
            rate += jl[m] * qd[Support[m]];
        }
        const double s = e.strutReference - l; // compression
        const double sd = -rate;
        const double damping = e.strutDamping; // compression and rebound damping are equal in this model
        const double spread = e.strutPreload + e.strutStiffness * s + damping * sd; // spreading force
        // The law's force on compression is -f; the force conjugate to l is +f.
        for (std::size_t m = 0; m < M; ++m)
        {
            force[Support[m]] += jl[m] * spread;
        }
        power += damping * sd * sd;
        const double dfdl = -e.strutStiffness;
        const double dfdrate = -damping;
        // d2l/dq2: (Jd.Jd - jl jl) / l + u . (z x J_b) on one chain, with a the ancestor (all revolute here).
        std::array<double, M * M> hessian{};
        for (std::size_t m = 0; m < M; ++m)
        {
            for (std::size_t n = m; n < M; ++n)
            {
                double entry = (jx[m] * jx[n] + jy[m] * jy[n] - jl[m] * jl[n]) / l;
                if ((m < 3) == (n < 3))
                {
                    entry += ux * -jy[n] + uy * jx[n]; // u . (z x J_n), J_n already signed for its chain
                }
                hessian[m * M + n] = entry;
                hessian[n * M + m] = entry;
            }
        }
        std::array<double, M> g{};
        for (std::size_t m = 0; m < M; ++m)
        {
            double sum = 0.0;
            for (std::size_t n = 0; n < M; ++n)
            {
                sum += hessian[m * M + n] * qd[Support[n]];
            }
            g[m] = sum;
        }
        for (std::size_t m = 0; m < M; ++m)
        {
            for (std::size_t n = 0; n < M; ++n)
            {
                const std::size_t at = Support[m] * N + Support[n];
                k[at] += dfdl * jl[m] * jl[n] + spread * hessian[m * M + n] + dfdrate * jl[m] * g[n];
                d[at] += dfdrate * jl[m] * jl[n];
            }
        }
    }

    bool step(double* q, double* qd, const double* tau, double dt)
    {
        transforms(q);
        biasForces(qd);
        massMatrix();
        worldPoses(q);
        double power = 0.0;
        elements(q, qd, power);
        for (std::size_t i = 0; i < N; ++i)
        {
            double stiffnessRate = 0.0;
            for (std::size_t j = 0; j < N; ++j)
            {
                const double sym = 0.5 * (k[i * N + j] + k[j * N + i]);
                if (j <= i)
                {
                    am[i * N + j] = h[i * N + j] - dt * d[i * N + j] - dt * dt * sym;
                }
                stiffnessRate += sym * qd[j];
            }
            b[i] = dt * (tau[i] + force[i] - c[i] + dt * stiffnessRate);
        }
        // Dense LDL^T, then the three triangular sweeps.
        std::array<double, N> scratch{};
        for (std::size_t j = 0; j < N; ++j)
        {
            double pivot = am[j * N + j];
            for (std::size_t m = 0; m < j; ++m)
            {
                scratch[m] = am[j * N + m] * am[m * N + m];
                pivot -= am[j * N + m] * scratch[m];
            }
            if (!(pivot > 0.0))
            {
                return false;
            }
            am[j * N + j] = pivot;
            for (std::size_t i = j + 1; i < N; ++i)
            {
                double sum = am[i * N + j];
                for (std::size_t m = 0; m < j; ++m)
                {
                    sum -= am[i * N + m] * scratch[m];
                }
                am[i * N + j] = sum / pivot;
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            for (std::size_t m = 0; m < i; ++m)
            {
                b[i] -= am[i * N + m] * b[m];
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            b[i] /= am[i * N + i];
        }
        for (std::size_t i = N; i-- > 0;)
        {
            for (std::size_t m = i + 1; m < N; ++m)
            {
                b[i] -= am[m * N + i] * b[m];
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            if (!std::isfinite(qd[i] + b[i]) || !std::isfinite(q[i] + dt * (qd[i] + b[i])))
            {
                return false;
            }
        }
        double work = 0.0;
        for (std::size_t i = 0; i < N; ++i)
        {
            qd[i] += b[i];
            q[i] += dt * qd[i];
            work += tau[i] * dt * qd[i];
        }
        dissipated += dt * power;
        externalWork += work;
        return true;
    }
};

void benchImplicitStep(sub0bench::BenchRunner& runner)
{
    using namespace sub0physics::multibody;
    const Model model = testmodels::planarTreeWithElements();
    Workspace ws(model);
    State start = State::atRest(model);
    std::array<double, N> tau{};
    for (std::size_t i = 0; i < N; ++i)
    {
        start.q()[i] = 0.3 * std::sin(static_cast<double>(i) + 1.0);
        start.qd()[i] = 0.5 * std::cos(static_cast<double>(i));
        tau[i] = 0.1 * static_cast<double>(i);
    }
    const double h = 1e-3;
    HandImplicitTree hand(model);

    // The reference must agree with the library before its time means anything: 200 steps from the start.
    {
        State s = start;
        std::array<double, N> q{};
        std::array<double, N> qd{};
        std::copy_n(start.values.begin(), N, q.begin());
        std::copy_n(start.values.begin() + N, N, qd.begin());
        for (std::uint64_t tick = 0; tick < 200U; ++tick)
        {
            if (!stepLinearlyImplicitEuler(model, {tau}, tick, sub0math::seconds(h), s, ws) ||
                !hand.step(q.data(), qd.data(), tau.data(), h))
            {
                std::fprintf(stderr, "implicit step failed\n");
                std::exit(1);
            }
        }
        for (std::size_t i = 0; i < N; ++i)
        {
            if (std::abs(q[i] - s.q()[i]) > 1e-10 || std::abs(qd[i] - s.qd()[i]) > 1e-8)
            {
                std::fprintf(stderr, "hand-written implicit step disagrees at %zu: %.17g vs %.17g\n", i, q[i], s.q()[i]);
                std::exit(1);
            }
        }
        if (std::abs(hand.dissipated - s.dissipated.si()) > 1e-9 || std::abs(hand.externalWork - s.externalWork.si()) > 1e-9)
        {
            std::fprintf(stderr, "hand-written implicit step energy disagrees\n");
            std::exit(1);
        }
    }

    State s = start;
    std::array<double, N> q{};
    std::array<double, N> qd{};
    runner.bench().minEpochIterations(5000);
    runner.run("planar tree 10-DOF implicit Euler step, 2 elements, hand-written planar", [&] {
        std::copy_n(start.values.begin(), N, q.begin());
        std::copy_n(start.values.begin() + N, N, qd.begin());
        ankerl::nanobench::doNotOptimizeAway(q);
        if (!hand.step(q.data(), qd.data(), tau.data(), h))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(qd);
    });
    runner.run("planar tree 10-DOF implicit Euler step, 2 elements, stepLinearlyImplicitEuler", [&] {
        std::copy(start.values.begin(), start.values.end(), s.values.begin());
        s.tick = 0;
        ankerl::nanobench::doNotOptimizeAway(s.values.data());
        if (!stepLinearlyImplicitEuler(model, {tau}, 0U, sub0math::seconds(h), s, ws))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(s.values.data());
    });
    // The elements alone, with Jacobians accumulated into dense n x n matrices, from the coordinates (world
    // poses included on both sides).
    runner.run("planar tree 10-DOF elements with Jacobians, hand-written planar", [&] {
        ankerl::nanobench::doNotOptimizeAway(start.values.data());
        double power = 0.0;
        hand.worldPoses(start.values.data());
        hand.elements(start.values.data(), start.values.data() + N, power);
        ankerl::nanobench::doNotOptimizeAway(hand.k);
        ankerl::nanobench::doNotOptimizeAway(power);
    });
    runner.run("planar tree 10-DOF elements with Jacobians, evaluateElements", [&] {
        ankerl::nanobench::doNotOptimizeAway(start.values.data());
        if (!evaluateElements(model, start, ws))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(ws.elementStiffness.data());
    });
}

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
        state.q()[i] = 0.3 * std::sin(static_cast<double>(i) + 1.0);
        state.qd()[i] = 0.5 * std::cos(static_cast<double>(i));
        tau[i] = 0.1 * static_cast<double>(i);
    }
    HandPlanarTree hand(model);

    // The reference must agree with the library before its time means anything.
    if (!forwardDynamics(model, state, tau, ws))
    {
        std::fprintf(stderr, "forwardDynamics failed\n");
        std::exit(1);
    }
    hand.forwardDynamics(state.q().data(), state.qd().data(), tau.data());
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
        hand.transforms(state.q().data());
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
        hand.forwardDynamics(state.q().data(), state.qd().data(), tau.data());
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
    benchImplicitStep(runner);
}

