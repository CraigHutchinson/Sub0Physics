#pragma once

// Terrain representation decision experiment (D0). Throwaway measurement code, not library API:
// it uses a minimal local vector type so the comparison measures representations, not vocabulary.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace d0
{

// ── Minimal vectors ─────────────────────────────────────────────────────────────────────────────

template<typename T>
struct V3
{
    T x{}, y{}, z{};
};

template<typename T> constexpr V3<T> operator+(V3<T> a, V3<T> b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
template<typename T> constexpr V3<T> operator-(V3<T> a, V3<T> b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
template<typename T> constexpr V3<T> operator*(V3<T> a, T s) { return {a.x * s, a.y * s, a.z * s}; }
template<typename T> constexpr T dot(V3<T> a, V3<T> b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template<typename T> constexpr V3<T> cross(V3<T> a, V3<T> b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template<typename T> T length(V3<T> a) { return std::sqrt(dot(a, a)); }
template<typename T> V3<T> normalized(V3<T> a)
{
    const T l = length(a);
    return l > T(0) ? a * (T(1) / l) : V3<T>{T(0), T(0), T(1)};
}

using V3d = V3<double>;
using V3f = V3<float>;

template<typename To, typename From> V3<To> cast(V3<From> v)
{
    return {static_cast<To>(v.x), static_cast<To>(v.y), static_cast<To>(v.z)};
}

struct V2d
{
    double x{}, y{};
};

constexpr double Pi = 3.14159265358979323846;

struct Box3
{
    V3d lo, hi;
};

/// Closest point on triangle abc to p (Ericson, Real-Time Collision Detection, 5.1.5).
template<typename T>
V3<T> closestOnTriangle(V3<T> p, V3<T> a, V3<T> b, V3<T> c)
{
    const V3<T> ab = b - a, ac = c - a, ap = p - a;
    const T d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= T(0) && d2 <= T(0)) return a;
    const V3<T> bp = p - b;
    const T d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= T(0) && d4 <= d3) return b;
    const T vc = d1 * d4 - d3 * d2;
    if (vc <= T(0) && d1 >= T(0) && d3 <= T(0)) return a + ab * (d1 / (d1 - d3));
    const V3<T> cp = p - c;
    const T d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= T(0) && d5 <= d6) return c;
    const T vb = d5 * d2 - d1 * d6;
    if (vb <= T(0) && d2 >= T(0) && d6 <= T(0)) return a + ac * (d2 / (d2 - d6));
    const T va = d3 * d6 - d5 * d4;
    if (va <= T(0) && (d4 - d3) >= T(0) && (d5 - d6) >= T(0))
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    const T denom = T(1) / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

/// Exact signed distance to a simple polygon; negative inside (Quilez).
inline double sdPolygon(const std::vector<V2d>& v, V2d p)
{
    const auto d2 = [](V2d a) { return a.x * a.x + a.y * a.y; };
    double d = d2({p.x - v[0].x, p.y - v[0].y});
    double s = 1.0;
    for (std::size_t i = 0, j = v.size() - 1; i < v.size(); j = i, ++i)
    {
        const V2d e{v[j].x - v[i].x, v[j].y - v[i].y};
        const V2d w{p.x - v[i].x, p.y - v[i].y};
        const double t = std::clamp((w.x * e.x + w.y * e.y) / (e.x * e.x + e.y * e.y), 0.0, 1.0);
        d = std::min(d, d2({w.x - e.x * t, w.y - e.y * t}));
        const bool c1 = p.y >= v[i].y, c2 = p.y < v[j].y, c3 = e.x * w.y > e.y * w.x;
        if ((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s = -s;
    }
    return s * std::sqrt(d);
}

// ── Exact features (E): oriented boxes and convex meshes ───────────────────────────────────────

template<typename T>
struct ObbT
{
    V3<T> c, half, ax, ay, az; ///< Centre, half extents, local axes in world.
};

template<typename T>
struct SdHit
{
    T d;         ///< Signed distance, positive outside.
    V3<T> n;     ///< Outward unit normal of the nearest boundary (direction from it to the point).
};

template<typename T>
SdHit<T> sdObb(const ObbT<T>& b, V3<T> p)
{
    const V3<T> d = p - b.c;
    const T l[3] = {dot(d, b.ax), dot(d, b.ay), dot(d, b.az)};
    const T h[3] = {b.half.x, b.half.y, b.half.z};
    const V3<T> axes[3] = {b.ax, b.ay, b.az};
    T q[3];
    bool outside = false;
    for (int k = 0; k < 3; ++k)
    {
        q[k] = std::abs(l[k]) - h[k];
        outside = outside || q[k] > T(0);
    }
    if (outside)
    {
        V3<T> closest = b.c;
        for (int k = 0; k < 3; ++k) closest = closest + axes[k] * std::clamp(l[k], -h[k], h[k]);
        const V3<T> v = p - closest;
        const T len = length(v);
        return {len, v * (T(1) / len)};
    }
    int k = 0;
    if (q[1] > q[k]) k = 1;
    if (q[2] > q[k]) k = 2;
    return {q[k], axes[k] * (l[k] < T(0) ? T(-1) : T(1))};
}

template<typename T>
struct ConvexMeshT
{
    std::vector<V3<T>> v;
    std::vector<std::array<std::uint32_t, 3>> tri;
    std::vector<V3<T>> n; ///< Outward face normals.
    std::vector<T> d;     ///< Plane offsets: dot(n, x) = d on the face.
};

template<typename T>
SdHit<T> sdConvexMesh(const ConvexMeshT<T>& m, V3<T> p)
{
    T maxPlane = -std::numeric_limits<T>::infinity();
    std::size_t face = 0;
    for (std::size_t f = 0; f < m.tri.size(); ++f)
    {
        const T s = dot(m.n[f], p) - m.d[f];
        if (s > maxPlane) { maxPlane = s; face = f; }
    }
    if (maxPlane <= T(0)) return {maxPlane, m.n[face]};
    T best = std::numeric_limits<T>::infinity();
    V3<T> bestC{};
    for (const auto& t : m.tri)
    {
        const V3<T> c = closestOnTriangle(p, m.v[t[0]], m.v[t[1]], m.v[t[2]]);
        const T dd = dot(p - c, p - c);
        if (dd < best) { best = dd; bestC = c; }
    }
    const T len = std::sqrt(best);
    return {len, (p - bestC) * (T(1) / len)};
}

template<typename To, typename From>
ConvexMeshT<To> castMesh(const ConvexMeshT<From>& m)
{
    ConvexMeshT<To> r;
    for (auto& x : m.v) r.v.push_back(cast<To>(x));
    r.tri = m.tri;
    for (auto& x : m.n) r.n.push_back(cast<To>(x));
    for (auto x : m.d) r.d.push_back(static_cast<To>(x));
    return r;
}

template<typename To, typename From>
ObbT<To> castObb(const ObbT<From>& b)
{
    return {cast<To>(b.c), cast<To>(b.half), cast<To>(b.ax), cast<To>(b.ay), cast<To>(b.az)};
}

/// Orients faces outward from an interior point and fills the planes.
inline void finishConvexMesh(ConvexMeshT<double>& m, V3d inside)
{
    m.n.clear();
    m.d.clear();
    for (auto& t : m.tri)
    {
        V3d n = normalized(cross(m.v[t[1]] - m.v[t[0]], m.v[t[2]] - m.v[t[0]]));
        if (dot(n, m.v[t[0]] - inside) < 0.0)
        {
            std::swap(t[1], t[2]);
            n = n * -1.0;
        }
        m.n.push_back(n);
        m.d.push_back(dot(n, m.v[t[0]]));
    }
}

/// A rock: an icosphere (one subdivision, 80 faces) scaled to an ellipsoid.
inline ConvexMeshT<double> makeEllipsoidRock(V3d centre, V3d radii)
{
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<V3d> v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                          {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    std::vector<std::array<std::uint32_t, 3>> f = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
        {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8},
        {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    std::unordered_map<std::uint64_t, std::uint32_t> mid;
    const auto midpoint = [&](std::uint32_t a, std::uint32_t b) {
        const std::uint64_t key = (std::uint64_t{std::min(a, b)} << 32) | std::max(a, b);
        if (auto it = mid.find(key); it != mid.end()) return it->second;
        v.push_back((v[a] + v[b]) * 0.5);
        const auto idx = static_cast<std::uint32_t>(v.size() - 1);
        mid.emplace(key, idx);
        return idx;
    };
    std::vector<std::array<std::uint32_t, 3>> g;
    for (auto& tri : f)
    {
        const auto a = midpoint(tri[0], tri[1]), b = midpoint(tri[1], tri[2]), c = midpoint(tri[2], tri[0]);
        g.push_back({tri[0], a, c});
        g.push_back({tri[1], b, a});
        g.push_back({tri[2], c, b});
        g.push_back({a, b, c});
    }
    ConvexMeshT<double> m;
    for (auto& p : v)
    {
        const V3d u = normalized(p);
        m.v.push_back({centre.x + u.x * radii.x, centre.y + u.y * radii.y, centre.z + u.z * radii.z});
    }
    m.tri = std::move(g);
    finishConvexMesh(m, centre);
    return m;
}

/// A convex polygon in (y, z) extruded along x.
inline ConvexMeshT<double> makePrismX(const std::vector<V2d>& yz, double x0, double x1)
{
    ConvexMeshT<double> m;
    const auto n = static_cast<std::uint32_t>(yz.size());
    for (auto& p : yz) m.v.push_back({x0, p.x, p.y});
    for (auto& p : yz) m.v.push_back({x1, p.x, p.y});
    for (std::uint32_t i = 1; i + 1 < n; ++i)
    {
        m.tri.push_back({0, i, i + 1});
        m.tri.push_back({n, n + i, n + i + 1});
    }
    for (std::uint32_t i = 0; i < n; ++i)
    {
        const std::uint32_t j = (i + 1) % n;
        m.tri.push_back({i, j, n + j});
        m.tri.push_back({i, n + j, n + i});
    }
    V3d c{};
    for (auto& p : m.v) c = c + p;
    finishConvexMesh(m, c * (1.0 / static_cast<double>(m.v.size())));
    return m;
}

inline ObbT<double> makeObb(V3d c, V3d half, double yaw, double pitch)
{
    const double cy = std::cos(yaw), sy = std::sin(yaw), cp = std::cos(pitch), sp = std::sin(pitch);
    return {c, half, {cy * cp, sy * cp, -sp}, {-sy, cy, 0.0}, {cy * sp, sy * sp, cp}};
}

// ── Scenes with exact references ────────────────────────────────────────────────────────────────

struct Scene
{
    std::string name;
    std::function<double(V3d)> sdf;          ///< Exact reference; positive outside the solid.
    std::function<double(double, double)> ground; ///< Single-valued ground top, without features.
    std::vector<ObbT<double>> boxes;         ///< Exact features (E).
    std::vector<ConvexMeshT<double>> meshes; ///< Exact features (E).
    Box3 domain;                             ///< Baked extent.
    Box3 queryRegion;                        ///< Where accuracy is sampled.
};

inline double featureSdf(const Scene& s, V3d p)
{
    double d = std::numeric_limits<double>::infinity();
    for (auto& b : s.boxes) d = std::min(d, sdObb(b, p).d);
    for (auto& m : s.meshes) d = std::min(d, sdConvexMesh(m, p).d);
    return d;
}

inline V3d referenceNormal(const std::function<double(V3d)>& sdf, V3d p)
{
    const double h = 1e-6;
    return normalized(V3d{sdf({p.x + h, p.y, p.z}) - sdf({p.x - h, p.y, p.z}),
                          sdf({p.x, p.y + h, p.z}) - sdf({p.x, p.y - h, p.z}),
                          sdf({p.x, p.y, p.z + h}) - sdf({p.x, p.y, p.z - h})});
}

/// Highest solid point on a vertical line, by sphere tracing down the exact SDF.
inline double topByRay(const std::function<double(V3d)>& sdf, double x, double y, double zTop, double zBottom)
{
    double z = zTop;
    for (int i = 0; i < 4096 && z > zBottom; ++i)
    {
        const double d = sdf({x, y, z});
        if (d < 1e-7) break;
        z -= std::max(d, 1e-7);
    }
    return z;
}

// ── B: heightfield page, 16-bit heights, fixed-diagonal triangulation, K deltas by chunk ─────────

class HeightField
{
public:
    static constexpr float Step = 1e-4f; ///< Quantisation: 0.1 mm.
    static constexpr int Chunk = 64;     ///< Dirty-delta chunk side in nodes.

    HeightField(double x0, double y0, double x1, double y1, float spacing)
        : x0_(x0), y0_(y0), h_(spacing),
          nx_(static_cast<int>(std::lround((x1 - x0) / spacing)) + 1),
          ny_(static_cast<int>(std::lround((y1 - y0) / spacing)) + 1),
          q_(static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_)),
          cx_((nx_ + Chunk - 1) / Chunk), cy_((ny_ + Chunk - 1) / Chunk),
          delta_(static_cast<std::size_t>(cx_ * cy_))
    {
    }

    template<typename F>
    void bake(F&& top)
    {
        std::vector<double> z(q_.size());
        double lo = std::numeric_limits<double>::infinity();
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
            {
                const double v = top(nodeX(i), nodeY(j));
                z[idx(i, j)] = v;
                lo = std::min(lo, v);
            }
        zOffset_ = static_cast<float>(lo);
        for (std::size_t k = 0; k < z.size(); ++k)
            q_[k] = static_cast<std::uint16_t>(std::clamp(std::lround((z[k] - lo) / Step), 0L, 65535L));
    }

    [[nodiscard]] double nodeX(int i) const { return x0_ + static_cast<double>(i) * h_; }
    [[nodiscard]] double nodeY(int j) const { return y0_ + static_cast<double>(j) * h_; }
    [[nodiscard]] int nx() const { return nx_; }
    [[nodiscard]] int ny() const { return ny_; }
    [[nodiscard]] float spacing() const { return h_; }

    [[nodiscard]] float z(int i, int j) const
    {
        float v = zOffset_ + static_cast<float>(q_[idx(i, j)]) * Step;
        if (const auto& c = delta_[chunkOf(i, j)]) v += static_cast<float>(c[local(i, j)]) * Step;
        return v;
    }

    /// Lowers node (i, j) to at most `target` through the sparse delta (K).
    void lowerTo(int i, int j, float target)
    {
        const float cur = z(i, j);
        if (target >= cur) return;
        auto& c = delta_[chunkOf(i, j)];
        if (!c)
        {
            c = std::make_unique<std::int16_t[]>(Chunk * Chunk);
            std::fill_n(c.get(), Chunk * Chunk, std::int16_t{0});
        }
        const long add = std::lround((target - cur) / Step);
        c[local(i, j)] = static_cast<std::int16_t>(std::clamp(static_cast<long>(c[local(i, j)]) + add, -32768L, 32767L));
    }

    [[nodiscard]] std::size_t baseBytes() const { return q_.size() * sizeof(std::uint16_t); }
    [[nodiscard]] std::size_t deltaBytes() const
    {
        std::size_t n = 0;
        for (auto& c : delta_) n += c ? Chunk * Chunk * sizeof(std::int16_t) : 0;
        return n + delta_.size() * sizeof(void*);
    }

    [[nodiscard]] float interp(float x, float y) const
    {
        const float lx = std::clamp((x - static_cast<float>(x0_)) / h_, 0.0f, static_cast<float>(nx_ - 1) - 1e-4f);
        const float ly = std::clamp((y - static_cast<float>(y0_)) / h_, 0.0f, static_cast<float>(ny_ - 1) - 1e-4f);
        const int i = static_cast<int>(lx), j = static_cast<int>(ly);
        const float fx = lx - static_cast<float>(i), fy = ly - static_cast<float>(j);
        const float za = z(i, j), zb = z(i + 1, j), zc = z(i, j + 1), zd = z(i + 1, j + 1);
        return fx >= fy ? za + (zb - za) * fx + (zd - zb) * fy : za + (zd - zc) * fx + (zc - za) * fy;
    }

    /// Exact signed distance to the triangulated surface (solid below) and its normal.
    [[nodiscard]] SdHit<float> query(V3f p) const
    {
        const float zs = interp(p.x, p.y);
        float best = std::abs(p.z - zs); // the vertical distance bounds the true distance
        V3f bestC{p.x, p.y, zs};
        const float fx0 = static_cast<float>(x0_), fy0 = static_cast<float>(y0_);
        const int i0 = std::max(0, static_cast<int>(std::floor((p.x - best - fx0) / h_)));
        const int i1 = std::min(nx_ - 2, static_cast<int>(std::floor((p.x + best - fx0) / h_)));
        const int j0 = std::max(0, static_cast<int>(std::floor((p.y - best - fy0) / h_)));
        const int j1 = std::min(ny_ - 2, static_cast<int>(std::floor((p.y + best - fy0) / h_)));
        for (int j = j0; j <= j1; ++j)
        {
            const float cy0 = fy0 + static_cast<float>(j) * h_;
            const float dy = std::max(0.0f, std::max(cy0 - p.y, p.y - (cy0 + h_)));
            if (dy * dy >= best * best) continue;
            for (int i = i0; i <= i1; ++i)
            {
                const float cx0 = fx0 + static_cast<float>(i) * h_;
                const float dx = std::max(0.0f, std::max(cx0 - p.x, p.x - (cx0 + h_)));
                const float dh = dx * dx + dy * dy;
                if (dh >= best * best) continue;
                const float za = z(i, j), zb = z(i + 1, j), zc = z(i, j + 1), zd = z(i + 1, j + 1);
                const float lo = std::min(std::min(za, zb), std::min(zc, zd));
                const float hi = std::max(std::max(za, zb), std::max(zc, zd));
                const float dz = std::max(0.0f, std::max(lo - p.z, p.z - hi));
                if (dh + dz * dz >= best * best) continue;
                const V3f a{cx0, cy0, za}, b{cx0 + h_, cy0, zb}, c{cx0, cy0 + h_, zc}, d{cx0 + h_, cy0 + h_, zd};
                for (const V3f& cp : {closestOnTriangle(p, a, b, d), closestOnTriangle(p, a, d, c)})
                {
                    const float dd = std::sqrt(dot(p - cp, p - cp));
                    if (dd < best) { best = dd; bestC = cp; }
                }
            }
        }
        const float sign = p.z >= zs ? 1.0f : -1.0f;
        const V3f n = best > 1e-7f ? (p - bestC) * (sign / best) : V3f{0, 0, 1};
        return {sign * best, n};
    }

private:
    [[nodiscard]] std::size_t idx(int i, int j) const
    {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(nx_) + static_cast<std::size_t>(i);
    }
    [[nodiscard]] std::size_t chunkOf(int i, int j) const
    {
        return static_cast<std::size_t>((j / Chunk) * cx_ + i / Chunk);
    }
    [[nodiscard]] static std::size_t local(int i, int j)
    {
        return static_cast<std::size_t>((j % Chunk) * Chunk + i % Chunk);
    }

    double x0_, y0_;
    float h_;
    int nx_, ny_;
    float zOffset_ = 0.0f;
    std::vector<std::uint16_t> q_;
    int cx_, cy_;
    std::vector<std::unique_ptr<std::int16_t[]>> delta_;
};

/// B+E: heightfield plus exact features, combined as a union.
struct HeightFieldWithFeatures
{
    const HeightField* field;
    std::vector<ObbT<float>> boxes;
    std::vector<ConvexMeshT<float>> meshes;

    [[nodiscard]] SdHit<float> query(V3f p) const
    {
        SdHit<float> best = field->query(p);
        for (auto& b : boxes)
            if (auto h = sdObb(b, p); h.d < best.d) best = h;
        for (auto& m : meshes)
            if (auto h = sdConvexMesh(m, p); h.d < best.d) best = h;
        return best;
    }
};

// ── H: sparse narrow-band SDF bricks (8^3 cells, 9^3 samples with a shared apron) ───────────────

class SdfBricks
{
public:
    static constexpr int Cells = 8;
    static constexpr int Samples = Cells + 1;
    static constexpr int BrickSamples = Samples * Samples * Samples;

    SdfBricks(float voxel, float band) : v_(voxel), band_(band), scale_(band / 32767.0f) {}

    template<typename F>
    void bake(F&& sdf, const Box3& region)
    {
        const double bs = static_cast<double>(Cells) * v_;
        const double halfDiag = std::sqrt(3.0) * bs * 0.5;
        for (int k = static_cast<int>(std::floor(region.lo.z / bs)); k <= static_cast<int>(std::floor(region.hi.z / bs)); ++k)
            for (int j = static_cast<int>(std::floor(region.lo.y / bs)); j <= static_cast<int>(std::floor(region.hi.y / bs)); ++j)
                for (int i = static_cast<int>(std::floor(region.lo.x / bs)); i <= static_cast<int>(std::floor(region.hi.x / bs)); ++i)
                {
                    const V3d c{(i + 0.5) * bs, (j + 0.5) * bs, (k + 0.5) * bs};
                    if (std::abs(sdf(c)) > band_ + halfDiag) continue;
                    bakeBrick(sdf, i, j, k, false);
                }
    }

    /// (Re)samples one brick from `sdf`; allocates it if missing (unless onlyIfNearSurface rejects it).
    template<typename F>
    bool bakeBrick(F&& sdf, int i, int j, int k, bool onlyIfNearSurface)
    {
        std::array<std::int16_t, BrickSamples> tmp{};
        bool near = false;
        for (int c = 0; c < Samples; ++c)
            for (int b = 0; b < Samples; ++b)
                for (int a = 0; a < Samples; ++a)
                {
                    const V3d p{(i * Cells + a) * static_cast<double>(v_), (j * Cells + b) * static_cast<double>(v_),
                                (k * Cells + c) * static_cast<double>(v_)};
                    const double d = sdf(p);
                    near = near || std::abs(d) < band_;
                    tmp[static_cast<std::size_t>((c * Samples + b) * Samples + a)] = quantise(static_cast<float>(d));
                }
        if (!near && (onlyIfNearSurface || !find(i, j, k))) return false;
        std::int16_t* dst = findOrAllocate(i, j, k, std::int16_t{0});
        std::copy(tmp.begin(), tmp.end(), dst);
        return true;
    }

    /// CSG subtraction in place: d = max(d, -cut(p)) over bricks overlapping `region`.
    /// Missing bricks under the old surface are allocated as "inside" (-band) when the cut reaches them.
    template<typename F, typename Inside>
    std::size_t subtract(F&& cut, const Box3& region, Inside&& isInside)
    {
        const double bs = static_cast<double>(Cells) * v_;
        std::size_t allocated = 0;
        for (int k = static_cast<int>(std::floor(region.lo.z / bs)); k <= static_cast<int>(std::floor(region.hi.z / bs)); ++k)
            for (int j = static_cast<int>(std::floor(region.lo.y / bs)); j <= static_cast<int>(std::floor(region.hi.y / bs)); ++j)
                for (int i = static_cast<int>(std::floor(region.lo.x / bs)); i <= static_cast<int>(std::floor(region.hi.x / bs)); ++i)
                {
                    const V3d c{(i + 0.5) * bs, (j + 0.5) * bs, (k + 0.5) * bs};
                    if (cut(c) > band_ + std::sqrt(3.0) * bs * 0.5) continue; // cut cannot reach this brick
                    std::int16_t* data = find(i, j, k);
                    if (!data)
                    {
                        if (!isInside(c)) continue; // empty air stays empty
                        data = findOrAllocate(i, j, k, quantise(-band_));
                        ++allocated;
                    }
                    for (int cz = 0; cz < Samples; ++cz)
                        for (int b = 0; b < Samples; ++b)
                            for (int a = 0; a < Samples; ++a)
                            {
                                const V3d p{(i * Cells + a) * static_cast<double>(v_), (j * Cells + b) * static_cast<double>(v_),
                                            (cz * 1.0 + k * Cells) * static_cast<double>(v_)};
                                auto& s = data[(cz * Samples + b) * Samples + a];
                                s = quantise(std::max(dequantise(s), static_cast<float>(-cut(p))));
                            }
                }
        return allocated;
    }

    struct Hit
    {
        float d;
        V3f n;
        bool known;
    };

    [[nodiscard]] Hit query(V3f p) const
    {
        const float bs = static_cast<float>(Cells) * v_;
        const int bi = static_cast<int>(std::floor(p.x / bs)), bj = static_cast<int>(std::floor(p.y / bs)),
                  bk = static_cast<int>(std::floor(p.z / bs));
        const std::int16_t* data = find(bi, bj, bk);
        if (!data) return {band_, {0, 0, 1}, false};
        const float lx = p.x / v_ - static_cast<float>(bi * Cells);
        const float ly = p.y / v_ - static_cast<float>(bj * Cells);
        const float lz = p.z / v_ - static_cast<float>(bk * Cells);
        const int a = std::clamp(static_cast<int>(lx), 0, Cells - 1), b = std::clamp(static_cast<int>(ly), 0, Cells - 1),
                  c = std::clamp(static_cast<int>(lz), 0, Cells - 1);
        const float fx = lx - static_cast<float>(a), fy = ly - static_cast<float>(b), fz = lz - static_cast<float>(c);
        const auto at = [&](int da, int db, int dc) {
            return dequantise(data[((c + dc) * Samples + (b + db)) * Samples + (a + da)]);
        };
        const float c000 = at(0, 0, 0), c100 = at(1, 0, 0), c010 = at(0, 1, 0), c110 = at(1, 1, 0);
        const float c001 = at(0, 0, 1), c101 = at(1, 0, 1), c011 = at(0, 1, 1), c111 = at(1, 1, 1);
        const float x00 = c000 + (c100 - c000) * fx, x10 = c010 + (c110 - c010) * fx;
        const float x01 = c001 + (c101 - c001) * fx, x11 = c011 + (c111 - c011) * fx;
        const float y0 = x00 + (x10 - x00) * fy, y1 = x01 + (x11 - x01) * fy;
        const float d = y0 + (y1 - y0) * fz;
        const float gx = ((c100 - c000) * (1 - fy) * (1 - fz) + (c110 - c010) * fy * (1 - fz) +
                          (c101 - c001) * (1 - fy) * fz + (c111 - c011) * fy * fz);
        const float gy = ((x10 - x00) * (1 - fz) + (x11 - x01) * fz);
        const float gz = (y1 - y0);
        return {d, normalized(V3f{gx, gy, gz}), true};
    }

    [[nodiscard]] std::size_t brickCount() const { return index_.size(); }
    [[nodiscard]] std::size_t bytes() const
    {
        return index_.size() * (BrickSamples * sizeof(std::int16_t) + 16); // 16 B: production page-table entry
    }
    [[nodiscard]] float voxel() const { return v_; }
    [[nodiscard]] float band() const { return band_; }

private:
    [[nodiscard]] std::int16_t quantise(float d) const
    {
        return static_cast<std::int16_t>(std::lround(std::clamp(d, -band_, band_) / scale_));
    }
    [[nodiscard]] float dequantise(std::int16_t s) const { return static_cast<float>(s) * scale_; }

    static std::uint64_t key(int i, int j, int k)
    {
        const auto u = [](int x) { return static_cast<std::uint64_t>(static_cast<std::uint32_t>(x + (1 << 20)) & 0x1FFFFF); };
        return (u(i) << 42) | (u(j) << 21) | u(k);
    }
    [[nodiscard]] const std::int16_t* find(int i, int j, int k) const
    {
        const auto it = index_.find(key(i, j, k));
        return it == index_.end() ? nullptr : &data_[static_cast<std::size_t>(it->second) * BrickSamples];
    }
    std::int16_t* find(int i, int j, int k)
    {
        const auto it = index_.find(key(i, j, k));
        return it == index_.end() ? nullptr : &data_[static_cast<std::size_t>(it->second) * BrickSamples];
    }
    std::int16_t* findOrAllocate(int i, int j, int k, std::int16_t fill)
    {
        if (auto* p = find(i, j, k)) return p;
        const auto idx = static_cast<std::uint32_t>(index_.size());
        index_.emplace(key(i, j, k), idx);
        data_.resize(data_.size() + BrickSamples, fill);
        return &data_[static_cast<std::size_t>(idx) * BrickSamples];
    }

    float v_, band_, scale_;
    std::unordered_map<std::uint64_t, std::uint32_t> index_;
    std::vector<std::int16_t> data_;
};

} // namespace d0
