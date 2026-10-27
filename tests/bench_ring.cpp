// The deformable ring substructure (docs/WHEEL.md section 2.3): one linearly implicit Euler step of a small planar
// tree carrying two 64-node rings (ringmodels::ringRig), paired with the fastest hand-written code we know for the
// same step: the tree's structure and joint types hard-coded, the ring size a compile-time constant, planar
// three-component algebra, 2 x 2 blocks unrolled. Also the same description expanded into ordinary bodies and tendons
// (the oracle), which takes the dense path: the cost the substructure avoids.

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
#include "test_ring_models.hpp"

namespace
{

constexpr std::size_t Nodes = 64;
constexpr std::size_t Rings = 2;
constexpr std::size_t T = ringmodels::TreeDof; // travel, heave, pitch, rear spin, front spin
constexpr std::size_t Dof = T + Rings * 2 * Nodes;
constexpr std::array<std::size_t, T> Parent{static_cast<std::size_t>(-1), 0, 1, 2, 2};
constexpr std::size_t None = static_cast<std::size_t>(-1);

struct M3
{
    double w;
    double x;
    double y;
};

struct In
{
    double m;
    double hx;
    double hy;
    double i;
};

struct X
{
    double c;
    double s;
    double x;
    double y;
};

inline M3 apply(const X& t, const M3& m)
{
    const double ux = m.x - t.y * m.w;
    const double uy = m.y + t.x * m.w;
    return {m.w, t.c * ux + t.s * uy, t.c * uy - t.s * ux};
}

inline M3 applyT(const X& t, const M3& f)
{
    const double fx = t.c * f.x - t.s * f.y;
    const double fy = t.s * f.x + t.c * f.y;
    return {f.w + t.x * fy - t.y * fx, fx, fy};
}

inline In applyT(const X& t, const In& in)
{
    const double yx = t.c * in.hx - t.s * in.hy;
    const double yy = t.s * in.hx + t.c * in.hy;
    return {in.m, yx + in.m * t.x, yy + in.m * t.y, in.i + 2.0 * (t.x * yx + t.y * yy) + in.m * (t.x * t.x + t.y * t.y)};
}

inline M3 times(const In& in, const M3& m)
{
    return {in.i * m.w + in.hx * m.y - in.hy * m.x, in.m * m.x - in.hy * m.w, in.m * m.y + in.hx * m.w};
}

inline M3 crossM(const M3& v, const M3& m)
{
    return {0.0, v.y * m.w - v.w * m.y, v.w * m.x - v.x * m.w};
}

inline M3 crossF(const M3& v, const M3& f)
{
    return {v.x * f.y - v.y * f.x, -v.w * f.y, v.w * f.x};
}

inline M3 operator+(const M3& a, const M3& b)
{
    return {a.w + b.w, a.x + b.x, a.y + b.y};
}

// A 2 x 2 block, row-major.
struct B2
{
    double a;
    double b;
    double c;
    double d;
};

inline B2 inverse(const B2& m)
{
    const double off = 0.5 * (m.b + m.c);
    const double r = 1.0 / (m.a * m.d - off * off);
    return {m.d * r, -off * r, -off * r, m.a * r};
}

inline B2 mulT(const B2& x, const B2& y) // x^T y
{
    return {x.a * y.a + x.c * y.c, x.a * y.b + x.c * y.d, x.b * y.a + x.d * y.c, x.b * y.b + x.d * y.d};
}

inline B2 mul(const B2& x, const B2& y)
{
    return {x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d};
}

// The hand-written step.
struct HandRingStep
{
    // Configuration, read once from the model.
    std::array<In, T> inertia{};
    std::array<X, T> placement{};
    double gx = 0.0; // base acceleration (-g) in the plane
    double gy = 0.0;
    std::array<std::size_t, Rings> host{};
    std::array<double, 2> a0{};
    std::array<double, 2> a1{};
    std::array<std::array<double, Nodes>, Rings> restX{}, restY{}, rest0{}, rest1{};
    double nodeMass = 0.0;
    double nodeInertia = 0.0;
    double kf = 0.0;
    double ke = 0.0;
    double ce = 0.0;
    double l0 = 0.0;

    // Scratch.
    std::array<X, T> xup{};
    std::array<M3, T> v{}, acc{}, f{};
    std::array<In, T> ic{};
    std::array<double, T * T> h{};
    std::array<std::array<double, Nodes>, Rings> px{}, py{};
    std::array<std::array<M3, Nodes>, Rings> nv{};
    std::array<std::array<B2, Nodes>, Rings> kd{}, ku{}, dd{}, du{}, ad{}, au{}, pinv{}, ln{}, lw{};
    std::array<std::array<std::array<double, 6>, Nodes>, Rings> border{}, z{};
    std::array<std::array<M3, 4>, Rings> jac{};
    std::array<double, Dof> bias{}, force{}, b{};
    double dissipated = 0.0;

    explicit HandRingStep(const sub0physics::multibody::Model& model)
    {
        const auto& pm = model.planar();
        for (std::size_t i = 0; i < T; ++i)
        {
            const auto& body = pm.bodies[i];
            inertia[i] = {body.inertia.m, body.inertia.hx, body.inertia.hy, body.inertia.i};
            placement[i] = {body.placement.c, body.placement.s, body.placement.x, body.placement.y};
        }
        gx = -pm.gravityX;
        gy = -pm.gravityY;
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const auto& ring = model.rings()[r];
            host[r] = ring.host;
            a0 = ring.axis0;
            a1 = ring.axis1;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                restX[r][i] = ring.restX[i];
                restY[r][i] = ring.restY[i];
                rest0[r][i] = ring.rest0[i];
                rest1[r][i] = ring.rest1[i];
            }
            nodeMass = ring.mass[0];
            nodeInertia = ring.inertia[0];
            kf = ring.foundation[0].spec().stiffness.si();
            ke = ring.edgeLaws[0].spec().stiffness.si();
            ce = ring.edgeLaws[0].spec().damping.si();
            l0 = ring.edgeRestLength[0];
        }
    }

    // Linearly implicit Euler with every Jacobian (ImplicitTerms::Full), no applied forces.
    bool step(double* q, double* qd, double dt)
    {
        // Tree transforms and velocities: travel along j (x), heave along i (z), three hinges about +k.
        xup[0] = {1.0, 0.0, placement[0].x, placement[0].y + q[0]};
        xup[1] = {1.0, 0.0, placement[1].x + q[1], placement[1].y};
        for (std::size_t i = 2; i < T; ++i)
        {
            xup[i] = {std::cos(q[i]), std::sin(q[i]), placement[i].x, placement[i].y};
        }
        const auto joint = [&](std::size_t i, double rate) {
            return i == 0 ? M3{0.0, 0.0, rate} : i == 1 ? M3{0.0, rate, 0.0} : M3{rate, 0.0, 0.0};
        };
        const auto sdot = [](std::size_t i, const M3& m) { return i == 0 ? m.y : i == 1 ? m.x : m.w; };
        for (std::size_t i = 0; i < T; ++i)
        {
            const M3 vj = joint(i, qd[i]);
            const std::size_t p = Parent[i];
            v[i] = p == None ? vj : apply(xup[i], v[p]) + vj;
            acc[i] = apply(xup[i], p == None ? M3{0.0, gx, gy} : acc[p]) + crossM(v[i], vj);
            f[i] = times(inertia[i], acc[i]) + crossF(v[i], times(inertia[i], v[i]));
        }
        // Rings: node kinematics, Newton-Euler, composite inertia and border, elements.
        const double* w = q + T;
        const double* wd = qd + T;
        const In nodeIn{nodeMass, 0.0, 0.0, nodeInertia};
        for (std::size_t i = 0; i < T; ++i)
        {
            ic[i] = inertia[i];
        }
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const std::size_t hb = host[r];
            const double* wr = w + r * 2 * Nodes;
            const double* wdr = wd + r * 2 * Nodes;
            M3 total{0.0, 0.0, 0.0};
            In composite{0.0, 0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const double x = restX[r][i] + a0[0] * wr[2 * i] + a1[0] * wr[2 * i + 1];
                const double y = restY[r][i] + a0[1] * wr[2 * i] + a1[1] * wr[2 * i + 1];
                px[r][i] = x;
                py[r][i] = y;
                const X t{1.0, 0.0, x, y};
                const M3 rate{0.0, a0[0] * wdr[2 * i] + a1[0] * wdr[2 * i + 1], a0[1] * wdr[2 * i] + a1[1] * wdr[2 * i + 1]};
                const M3 vn = apply(t, v[hb]) + rate;
                nv[r][i] = vn;
                const M3 an = apply(t, acc[hb]) + crossM(vn, rate);
                const M3 fn = times(nodeIn, an) + crossF(vn, times(nodeIn, vn));
                bias[T + r * 2 * Nodes + 2 * i] = a0[0] * fn.x + a0[1] * fn.y;
                bias[T + r * 2 * Nodes + 2 * i + 1] = a1[0] * fn.x + a1[1] * fn.y;
                total = total + applyT(t, fn);
                const In ni = applyT(t, nodeIn);
                composite = {composite.m + ni.m, composite.hx + ni.hx, composite.hy + ni.hy, composite.i + ni.i};
                const double m = nodeMass;
                border[r][i] = {m * (a0[1] * x - a0[0] * y), m * a0[0], m * a0[1],
                                m * (a1[1] * x - a1[0] * y), m * a1[0], m * a1[1]};
            }
            f[hb] = f[hb] + total;
            ic[hb] = {ic[hb].m + composite.m, ic[hb].hx + composite.hx, ic[hb].hy + composite.hy, ic[hb].i + composite.i};
        }
        for (std::size_t i = T; i-- > 0;)
        {
            bias[i] = sdot(i, f[i]);
            if (Parent[i] != None)
            {
                f[Parent[i]] = f[Parent[i]] + applyT(xup[i], f[i]);
            }
        }
        for (std::size_t i = T; i-- > 0;)
        {
            if (Parent[i] != None)
            {
                const In t = applyT(xup[i], ic[i]);
                In& p = ic[Parent[i]];
                p = {p.m + t.m, p.hx + t.hx, p.hy + t.hy, p.i + t.i};
            }
        }
        for (std::size_t i = 0; i < T; ++i)
        {
            const In& c = ic[i];
            M3 fs = i == 0 ? M3{c.hx, 0.0, c.m} : i == 1 ? M3{-c.hy, c.m, 0.0} : M3{c.i, -c.hy, c.hx};
            h[i * T + i] = sdot(i, fs);
            for (std::size_t j = i; Parent[j] != None;)
            {
                fs = applyT(xup[j], fs);
                j = Parent[j];
                h[i * T + j] = sdot(j, fs);
            }
        }
        // Host Jacobians: hub, frame, heave, travel.
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const std::size_t chain[4] = {host[r], 2, 1, 0};
            for (std::size_t c = 0; c < 4; ++c)
            {
                M3 m = joint(chain[c], 1.0);
                for (std::size_t d = c; d-- > 0;)
                {
                    m = apply(xup[chain[d]], m);
                }
                jac[r][c] = m;
            }
        }
        // Ring elements with exact Jacobians, then the step's ring blocks.
        double power = 0.0;
        const double h2 = dt * dt;
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const double* wr = w + r * 2 * Nodes;
            const double* wdr = wd + r * 2 * Nodes;
            double* fr = force.data() + T + r * 2 * Nodes;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                fr[2 * i] = -kf * wr[2 * i];
                fr[2 * i + 1] = -kf * wr[2 * i + 1];
                kd[r][i] = {-kf, 0.0, 0.0, -kf};
                dd[r][i] = {0.0, 0.0, 0.0, 0.0};
            }
            for (std::size_t e = 0; e < Nodes; ++e)
            {
                const std::size_t n = e + 1 == Nodes ? 0 : e + 1;
                const double dx = (rest0[r][n] + wr[2 * n]) - (rest0[r][e] + wr[2 * e]);
                const double dy = (rest1[r][n] + wr[2 * n + 1]) - (rest1[r][e] + wr[2 * e + 1]);
                const double vx = wdr[2 * n] - wdr[2 * e];
                const double vy = wdr[2 * n + 1] - wdr[2 * e + 1];
                const double l = std::sqrt(dx * dx + dy * dy);
                const double il = 1.0 / l;
                const double ux = dx * il;
                const double uy = dy * il;
                const double rate = ux * vx + uy * vy;
                const double tension = -ke * (l - l0) - ce * rate;
                power += ce * rate * rate;
                fr[2 * n] += tension * ux;
                fr[2 * n + 1] += tension * uy;
                fr[2 * e] -= tension * ux;
                fr[2 * e + 1] -= tension * uy;
                const double p = tension * il;
                const double tx = (vx - rate * ux) * il;
                const double ty = (vy - rate * uy) * il;
                const B2 g{-ke * ux * ux + p * (1.0 - ux * ux) - ce * ux * tx, -ke * ux * uy - p * ux * uy - ce * ux * ty,
                           -ke * uy * ux - p * uy * ux - ce * uy * tx, -ke * uy * uy + p * (1.0 - uy * uy) - ce * uy * ty};
                const B2 d{-ce * ux * ux, -ce * ux * uy, -ce * uy * ux, -ce * uy * uy};
                kd[r][e] = {kd[r][e].a + g.a, kd[r][e].b + g.b, kd[r][e].c + g.c, kd[r][e].d + g.d};
                kd[r][n] = {kd[r][n].a + g.a, kd[r][n].b + g.b, kd[r][n].c + g.c, kd[r][n].d + g.d};
                ku[r][e] = {-g.a, -g.b, -g.c, -g.d};
                dd[r][e] = {dd[r][e].a + d.a, dd[r][e].b + d.b, dd[r][e].c + d.c, dd[r][e].d + d.d};
                dd[r][n] = {dd[r][n].a + d.a, dd[r][n].b + d.b, dd[r][n].c + d.c, dd[r][n].d + d.d};
                du[r][e] = {-d.a, -d.b, -d.c, -d.d};
            }
        }
        // Right-hand side h (tau_e - c + h K~ qd); the tree has no elements.
        for (std::size_t i = 0; i < T; ++i)
        {
            b[i] = dt * (0.0 - bias[i]);
        }
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const double* wdr = wd + r * 2 * Nodes;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const std::size_t pv = i == 0 ? Nodes - 1 : i - 1;
                const std::size_t nx = i + 1 == Nodes ? 0 : i + 1;
                const B2& k0 = kd[r][i];
                const B2& kn = ku[r][i];
                const B2& kp = ku[r][pv];
                const double o0 = 0.5 * (k0.b + k0.c);
                const double on = 0.5 * (kn.b + kn.c);
                const double op = 0.5 * (kp.b + kp.c);
                const double s0 = k0.a * wdr[2 * i] + o0 * wdr[2 * i + 1] + kn.a * wdr[2 * nx] + on * wdr[2 * nx + 1] +
                                  kp.a * wdr[2 * pv] + op * wdr[2 * pv + 1];
                const double s1 = o0 * wdr[2 * i] + k0.d * wdr[2 * i + 1] + on * wdr[2 * nx] + kn.d * wdr[2 * nx + 1] +
                                  op * wdr[2 * pv] + kp.d * wdr[2 * pv + 1];
                const std::size_t j = T + r * 2 * Nodes + 2 * i;
                b[j] = dt * (force[j] - bias[j] + dt * s0);
                b[j + 1] = dt * (force[j + 1] - bias[j + 1] + dt * s1);
            }
        }
        // Condense each ring: periodic block LDL^T, Z = A^-1 B, y = A^-1 b_w; tree -= J^T B^T Z J, b_t -= J^T B^T y.
        for (std::size_t r = 0; r < Rings; ++r)
        {
            auto& A = ad[r];
            auto& U = au[r];
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const B2& k0 = kd[r][i];
                const B2& kn = ku[r][i];
                const double o0 = 0.5 * (k0.b + k0.c);
                const double on = 0.5 * (kn.b + kn.c);
                A[i] = {nodeMass - dt * dd[r][i].a - h2 * k0.a, -dt * dd[r][i].b - h2 * o0, -dt * dd[r][i].c - h2 * o0,
                        nodeMass - dt * dd[r][i].d - h2 * k0.d};
                U[i] = {-(dt * du[r][i].a + h2 * kn.a), -(dt * du[r][i].b + h2 * on), -(dt * du[r][i].c + h2 * on),
                        -(dt * du[r][i].d + h2 * kn.d)};
            }
            constexpr std::size_t L = Nodes - 1;
            B2 piv = A[0];
            B2 last = A[L];
            B2 fill{U[L].a, U[L].c, U[L].b, U[L].d};
            for (std::size_t k = 0; k + 1 < L; ++k)
            {
                pinv[r][k] = inverse(piv);
                ln[r][k] = mulT(U[k], pinv[r][k]);
                lw[r][k] = mulT(fill, pinv[r][k]);
                const B2 nu = mul(ln[r][k], U[k]);
                const B2 lu = mul(lw[r][k], fill);
                const B2 fu = mul(ln[r][k], fill);
                piv = {A[k + 1].a - nu.a, A[k + 1].b - nu.b, A[k + 1].c - nu.c, A[k + 1].d - nu.d};
                last = {last.a - lu.a, last.b - lu.b, last.c - lu.c, last.d - lu.d};
                const B2 base = k + 2 == L ? U[k + 1] : B2{0.0, 0.0, 0.0, 0.0};
                fill = {base.a - fu.a, base.b - fu.b, base.c - fu.c, base.d - fu.d};
            }
            pinv[r][L - 1] = inverse(piv);
            lw[r][L - 1] = mulT(fill, pinv[r][L - 1]);
            const B2 lu = mul(lw[r][L - 1], fill);
            last = {last.a - lu.a, last.b - lu.b, last.c - lu.c, last.d - lu.d};
            pinv[r][L] = inverse(last);
            if (!std::isfinite(pinv[r][L].a))
            {
                return false;
            }
            // Four columns at once: B (three) and b_w.
            std::array<std::array<double, 8>, Nodes> x{};
            double* bw = b.data() + T + r * 2 * Nodes;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const auto& bb = border[r][i];
                x[i] = {bb[0], bb[1], bb[2], bw[2 * i], bb[3], bb[4], bb[5], bw[2 * i + 1]};
            }
            const auto sub = [](const B2& m, const std::array<double, 8>& from, std::array<double, 8>& to) {
                for (std::size_t c = 0; c < 4; ++c)
                {
                    to[c] -= m.a * from[c] + m.b * from[4 + c];
                    to[4 + c] -= m.c * from[c] + m.d * from[4 + c];
                }
            };
            const auto subT = [](const B2& m, const std::array<double, 8>& from, std::array<double, 8>& to) {
                for (std::size_t c = 0; c < 4; ++c)
                {
                    to[c] -= m.a * from[c] + m.c * from[4 + c];
                    to[4 + c] -= m.b * from[c] + m.d * from[4 + c];
                }
            };
            const auto app = [](const B2& m, const std::array<double, 8>& from) {
                std::array<double, 8> to{};
                for (std::size_t c = 0; c < 4; ++c)
                {
                    to[c] = m.a * from[c] + m.b * from[4 + c];
                    to[4 + c] = m.c * from[c] + m.d * from[4 + c];
                }
                return to;
            };
            for (std::size_t k = 0; k < L; ++k)
            {
                if (k + 1 < L)
                {
                    sub(ln[r][k], x[k], x[k + 1]);
                }
                sub(lw[r][k], x[k], x[L]);
            }
            x[L] = app(pinv[r][L], x[L]);
            for (std::size_t k = L; k-- > 0;)
            {
                std::array<double, 8> y = app(pinv[r][k], x[k]);
                if (k + 1 < L)
                {
                    subT(ln[r][k], x[k + 1], y);
                }
                subT(lw[r][k], x[L], y);
                x[k] = y;
            }
            std::array<double, 9> hm{};
            double g0 = 0.0;
            double g1 = 0.0;
            double g2 = 0.0;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const auto& bb = border[r][i];
                z[r][i] = {x[i][0], x[i][1], x[i][2], x[i][4], x[i][5], x[i][6]};
                for (std::size_t a = 0; a < 3; ++a)
                {
                    for (std::size_t c = a; c < 3; ++c)
                    {
                        hm[a * 3 + c] += bb[a] * x[i][c] + bb[3 + a] * x[i][4 + c];
                    }
                }
                bw[2 * i] = x[i][3];
                bw[2 * i + 1] = x[i][7];
                g0 += bb[0] * x[i][3] + bb[3] * x[i][7];
                g1 += bb[1] * x[i][3] + bb[4] * x[i][7];
                g2 += bb[2] * x[i][3] + bb[5] * x[i][7];
            }
            hm[3] = hm[1];
            hm[6] = hm[2];
            hm[7] = hm[5];
            const std::size_t chain[4] = {host[r], 2, 1, 0};
            for (std::size_t a = 0; a < 4; ++a)
            {
                const M3& ja = jac[r][a];
                const double q0 = hm[0] * ja.w + hm[1] * ja.x + hm[2] * ja.y;
                const double q1 = hm[3] * ja.w + hm[4] * ja.x + hm[5] * ja.y;
                const double q2 = hm[6] * ja.w + hm[7] * ja.x + hm[8] * ja.y;
                for (std::size_t c = a; c < 4; ++c)
                {
                    const M3& jc = jac[r][c];
                    h[chain[a] * T + chain[c]] -= jc.w * q0 + jc.x * q1 + jc.y * q2;
                }
                b[chain[a]] -= ja.w * g0 + ja.x * g1 + ja.y * g2;
            }
        }
        // Tree LTDL factor and solve (lower triangle).
        for (std::size_t k = T; k-- > 0;)
        {
            for (std::size_t i = Parent[k]; i != None; i = Parent[i])
            {
                const double s = h[k * T + i] / h[k * T + k];
                for (std::size_t j = i; j != None; j = Parent[j])
                {
                    h[i * T + j] -= s * h[k * T + j];
                }
                h[k * T + i] = s;
            }
        }
        for (std::size_t i = T; i-- > 0;)
        {
            for (std::size_t j = Parent[i]; j != None; j = Parent[j])
            {
                b[j] -= h[i * T + j] * b[i];
            }
        }
        for (std::size_t i = 0; i < T; ++i)
        {
            b[i] /= h[i * T + i];
        }
        for (std::size_t i = 0; i < T; ++i)
        {
            for (std::size_t j = Parent[i]; j != None; j = Parent[j])
            {
                b[i] -= h[i * T + j] * b[j];
            }
        }
        // Back-substitution: w = y - Z J dqd_tree.
        for (std::size_t r = 0; r < Rings; ++r)
        {
            const std::size_t chain[4] = {host[r], 2, 1, 0};
            double v0 = 0.0;
            double v1 = 0.0;
            double v2 = 0.0;
            for (std::size_t c = 0; c < 4; ++c)
            {
                v0 += jac[r][c].w * b[chain[c]];
                v1 += jac[r][c].x * b[chain[c]];
                v2 += jac[r][c].y * b[chain[c]];
            }
            double* bw = b.data() + T + r * 2 * Nodes;
            for (std::size_t i = 0; i < Nodes; ++i)
            {
                const auto& zz = z[r][i];
                bw[2 * i] -= zz[0] * v0 + zz[1] * v1 + zz[2] * v2;
                bw[2 * i + 1] -= zz[3] * v0 + zz[4] * v1 + zz[5] * v2;
            }
        }
        for (std::size_t i = 0; i < Dof; ++i)
        {
            if (!std::isfinite(b[i]))
            {
                return false;
            }
        }
        for (std::size_t i = 0; i < Dof; ++i)
        {
            qd[i] += b[i];
            q[i] += dt * qd[i];
        }
        dissipated += dt * power;
        return true;
    }
};

} // namespace

void benchRing(sub0bench::BenchRunner& runner)
{
    using namespace sub0physics::multibody;
    const ringmodels::RingParameters parameters{.nodes = Nodes};
    const Model model = ringmodels::ringRig(parameters, false);
    const Model expanded = ringmodels::ringRig(parameters, true);
    if (model.dofCount() != Dof)
    {
        std::fprintf(stderr, "ring rig has an unexpected size\n");
        std::exit(1);
    }
    Workspace ws(model);
    Workspace wx(expanded);
    const State start = ringmodels::randomState(model, 3);
    const double h = 31.25e-6;
    HandRingStep hand(model);

    // The twin and the oracle must agree with the library before their times mean anything: 200 steps.
    {
        State s = start;
        State e = start;
        std::array<double, Dof> q{};
        std::array<double, Dof> qd{};
        std::copy_n(start.values.begin(), Dof, q.begin());
        std::copy_n(start.values.begin() + Dof, Dof, qd.begin());
        for (std::uint64_t tick = 0; tick < 200U; ++tick)
        {
            if (!stepLinearlyImplicitEuler(model, {}, tick, sub0math::seconds(h), s, ws) ||
                !stepLinearlyImplicitEuler(expanded, {}, tick, sub0math::seconds(h), e, wx) ||
                !hand.step(q.data(), qd.data(), h))
            {
                std::fprintf(stderr, "ring step failed\n");
                std::exit(1);
            }
        }
        for (std::size_t i = 0; i < Dof; ++i)
        {
            const double scale = 1.0 + std::abs(s.qd()[i]);
            if (std::abs(q[i] - s.q()[i]) > 1e-12 * (1.0 + std::abs(s.q()[i])) || std::abs(qd[i] - s.qd()[i]) > 1e-10 * scale ||
                std::abs(e.q()[i] - s.q()[i]) > 1e-12 * (1.0 + std::abs(s.q()[i])) || std::abs(e.qd()[i] - s.qd()[i]) > 1e-10 * scale)
            {
                std::fprintf(stderr, "ring step disagrees at %zu: library %.17g, twin %.17g, oracle %.17g\n", i,
                             s.qd()[i], qd[i], e.qd()[i]);
                std::exit(1);
            }
        }
        if (std::abs(hand.dissipated - s.dissipated.si()) > 1e-12 * (1.0 + s.dissipated.si()))
        {
            std::fprintf(stderr, "ring twin dissipation disagrees\n");
            std::exit(1);
        }
    }

    State s = start;
    std::array<double, Dof> q{};
    std::array<double, Dof> qd{};
    runner.bench().minEpochIterations(5000);
    runner.run("two 64-node rings implicit Euler step, hand-written condensed", [&] {
        std::copy_n(start.values.begin(), Dof, q.begin());
        std::copy_n(start.values.begin() + Dof, Dof, qd.begin());
        ankerl::nanobench::doNotOptimizeAway(q);
        if (!hand.step(q.data(), qd.data(), h))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(qd);
    });
    runner.run("two 64-node rings implicit Euler step, stepLinearlyImplicitEuler", [&] {
        std::copy(start.values.begin(), start.values.end(), s.values.begin());
        s.tick = 0;
        ankerl::nanobench::doNotOptimizeAway(s.values.data());
        if (!stepLinearlyImplicitEuler(model, {}, 0U, sub0math::seconds(h), s, ws))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(s.values.data());
    });
    State e = start;
    runner.bench().minEpochIterations(30);
    runner.run("two 64-node rings implicit Euler step, expanded bodies (dense path)", [&] {
        std::copy(start.values.begin(), start.values.end(), e.values.begin());
        e.tick = 0;
        ankerl::nanobench::doNotOptimizeAway(e.values.data());
        if (!stepLinearlyImplicitEuler(expanded, {}, 0U, sub0math::seconds(h), e, wx))
        {
            std::abort();
        }
        ankerl::nanobench::doNotOptimizeAway(e.values.data());
    });
}
