#include "sub0physics/multibody/casing_models.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

#include <sub0math/number/det_math.hpp>

namespace sub0physics::multibody
{

namespace
{

// The moment about +y of a vector (fx, fz) at offset (dx, dz): z-component of r x f in the (z, x) plane convention.
double moment(double dx, double dz, double fx, double fz) noexcept
{
    return fx * dz - fz * dx;
}

} // namespace

TreadFaceResponse smoothedTreadFace(double stiffness, double damping, double band, double p0, double p1, double rate0,
                                    double rate1) noexcept
{
    TreadFaceResponse out;
    const double half = 0.5 * (p0 - p1);
    const double a = std::abs(half);
    const double sign = half >= 0.0 ? 1.0 : -1.0;
    double h = a;
    double slope = 1.0;
    if (a < band)
    {
        h = a * a * (2.0 * band - a) / (band * band);
        slope = a * (4.0 * band - 3.0 * a) / (band * band);
    }
    const double m = 0.5 * (p0 + p1) + h;
    if (!(m > 0.0))
    {
        return out;
    }
    const double w0 = 0.5 + 0.5 * slope * sign;
    const double w1 = 0.5 - 0.5 * slope * sign;
    const double rate = w0 * rate0 + w1 * rate1;
    const double elastic = stiffness * m;
    const double total = elastic + damping * rate;
    const double force = total > 0.0 ? total : 0.0;
    out.force0 = w0 * force;
    out.force1 = w1 * force;
    out.elastic0 = w0 * elastic;
    out.elastic1 = w1 * elastic;
    out.potential = 0.5 * elastic * m;
    out.dissipationPower = (force - elastic) * rate;
    out.penetration = m;
    return out;
}

TreadFaceResponse distributedTreadFace(double stiffness, double damping, double p0, double p1, double rate0,
                                       double rate1) noexcept
{
    TreadFaceResponse out;
    const bool flip = p1 > p0;
    const double pa = flip ? p1 : p0;
    const double pb = flip ? p0 : p1;
    if (!(pa > 0.0))
    {
        return out;
    }
    const double ra = flip ? rate1 : rate0;
    const double rb = flip ? rate0 : rate1;
    // The penetrating interval [0, a] measured from the deeper end; moments of (1 - s)^2, s (1 - s), s^2 over it.
    const double a = pb >= 0.0 ? 1.0 : pa / (pa - pb);
    const double b = 1.0 - a;
    const double maa = (1.0 - b * b * b) / 3.0;
    const double mab = 0.5 * a * a - a * a * a / 3.0;
    const double mbb = a * a * a / 3.0;
    const double ea = stiffness * (maa * pa + mab * pb);
    const double eb = stiffness * (mab * pa + mbb * pb);
    const double ta = ea + damping * (maa * ra + mab * rb);
    const double tb = eb + damping * (mab * ra + mbb * rb);
    const double fa = ta > 0.0 ? ta : 0.0;
    const double fb = tb > 0.0 ? tb : 0.0;
    out.potential = 0.5 * stiffness * (maa * pa * pa + 2.0 * mab * pa * pb + mbb * pb * pb);
    out.dissipationPower = (fa - ea) * ra + (fb - eb) * rb;
    out.penetration = pa;
    out.force0 = flip ? fb : fa;
    out.force1 = flip ? fa : fb;
    out.elastic0 = flip ? eb : ea;
    out.elastic1 = flip ? ea : eb;
    return out;
}

PlanarRigidState fitRigidRing(const NodalCasingView& casing, double ringInertia) noexcept
{
    const std::size_t n = casing.x.size();
    double mass = 0.0, cx = 0.0, cz = 0.0, vx = 0.0, vz = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double m = casing.mass[i];
        mass += m;
        cx += m * casing.x[i];
        cz += m * casing.z[i];
        vx += m * casing.vx[i];
        vz += m * casing.vz[i];
    }
    cx /= mass;
    cz /= mass;
    vx /= mass;
    vz /= mass;
    const auto sc = sub0math::det::sincos(casing.hostAngle);
    double sine = 0.0, cosine = 0.0, angular = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double m = casing.mass[i];
        const double ax = sc.cos * casing.restX[i] + sc.sin * casing.restZ[i];
        const double az = sc.cos * casing.restZ[i] - sc.sin * casing.restX[i];
        const double bx = casing.x[i] - cx, bz = casing.z[i] - cz;
        sine += m * moment(ax, az, bx, bz);
        cosine += m * (ax * bx + az * bz);
        angular += m * moment(bx, bz, casing.vx[i] - vx, casing.vz[i] - vz);
    }
    return {.x = cx,
            .z = cz,
            .angle = casing.hostAngle + sub0math::det::atan2(sine, cosine),
            .vx = vx,
            .vz = vz,
            .rate = (angular + casing.nodeSpin * casing.hostRate) / ringInertia};
}

double placeNodalCasing(const PlanarRigidState& ring, std::span<const double> restX, std::span<const double> restZ,
                        double hostInertia, double hostRate, double nodeSpin, std::span<double> x, std::span<double> z,
                        std::span<double> vx, std::span<double> vz) noexcept
{
    const auto sc = sub0math::det::sincos(ring.angle);
    for (std::size_t i = 0; i < restX.size(); ++i)
    {
        const double ox = sc.cos * restX[i] + sc.sin * restZ[i];
        const double oz = sc.cos * restZ[i] - sc.sin * restX[i];
        x[i] = ring.x + ox;
        z[i] = ring.z + oz;
        vx[i] = ring.vx + ring.rate * oz;
        vz[i] = ring.vz - ring.rate * ox;
    }
    return (hostInertia * hostRate + nodeSpin * ring.rate) / (hostInertia + nodeSpin);
}

} // namespace sub0physics::multibody
