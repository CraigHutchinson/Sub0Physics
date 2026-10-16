"""D1: vertical-height contact vs contact along the surface normal, in a 2D cross-section.

A tyre crown (a circle in the y-z plane, across the tread) is pressed into terrain by a fixed load.
Two contact models give the force per unit length on the tyre:

  vertical: penetration = terrain height under each crown sample minus the crown's bottom, force +z only
            (what a heightfield queried straight down gives).
  normal:   penetration along the line from the crown centre to each terrain point inside the crown,
            force along that line (pressure normal to the surfaces).

Both use the same linear penalty stiffness and are solved for vertical equilibrium. The ratio
|Fy| / Fz is the friction coefficient needed to stop the tyre sliding sideways: the rut's lateral push.
Standard library only.  python contact_projection.py
"""

import math

K = 1.0e8          # penalty stiffness per metre of contact length: a few mm of crown penetration at this load
LOAD = 400.0 / 0.05  # 400 N over a 50 mm long contact patch, as N per metre of length
RC = 0.029         # crown radius of a ~58 mm (2.3") tyre
MU = 0.6           # illustrative available friction, for the slide verdict only


def rut(y, width=0.08, depth=0.06):
    return -depth * 0.5 * (1 + math.cos(2 * math.pi * y / width)) if abs(y) < width / 2 else 0.0


def slope45(y):
    return max(0.0, y)  # flat, then a 45 degree berm face from y = 0


def forces_vertical(terrain, yc, zc, n=1000):
    fz = 0.0
    dy = 2 * RC / n
    for i in range(n):
        y = yc - RC + (i + 0.5) * dy
        bottom = zc - math.sqrt(max(0.0, RC * RC - (y - yc) ** 2))
        fz += K * max(0.0, terrain(y) - bottom) * dy
    return 0.0, fz


def forces_normal(terrain, yc, zc, n=1500, span=0.04):
    fy = fz = 0.0
    ys = [yc - span + 2 * span * i / n for i in range(n + 1)]
    for a, b in zip(ys, ys[1:]):
        za, zb = terrain(a), terrain(b)
        ym, zm = (a + b) / 2, (za + zb) / 2
        ds = math.hypot(b - a, zb - za)
        dx, dz = yc - ym, zc - zm
        dist = math.hypot(dx, dz)
        pen = RC - dist
        if pen > 0 and dist > 0:
            fy += K * pen * ds * dx / dist
            fz += K * pen * ds * dz / dist
    return fy, fz


def equilibrium(model, terrain, yc):
    # Walk down from clear of the ground to the first height that carries the load, then bisect.
    # (Bisecting the whole range fails: a fully buried crown has no terrain inside it.)
    hi = terrain(yc) + RC + 0.01
    lo = hi
    while model(terrain, yc, lo)[1] < LOAD:
        hi, lo = lo, lo - 0.0002
        if lo < terrain(yc) - 2 * RC:
            raise RuntimeError("load never carried; stiffness too low for this load")
    for _ in range(40):
        mid = (lo + hi) / 2
        _, fz = model(terrain, yc, mid)
        if fz > LOAD:
            lo = mid
        else:
            hi = mid
    zc = (lo + hi) / 2
    return zc, model(terrain, yc, zc)


def main():
    print("| Case | Crown offset from rut centre mm | Contact model | Lateral push / load (toward the rut centre, or into the corner) | Exceeds grip at mu = 0.6? |")
    print("|---|---|---|---|---|")
    for name, terrain, offsets in [("Rut 8 cm wide, 6 cm deep", rut, [0.0, 0.005, 0.01, 0.015, 0.02]),
                                   ("Berm face 45°", slope45, [0.05])]:
        for off in offsets:
            for label, model in [("vertical height", forces_vertical), ("surface normal", forces_normal)]:
                _, (fy, fz) = equilibrium(model, terrain, off)
                ratio = abs(fy) / fz if fz > 0 else 0.0
                print(f"| {name} | {off * 1000:.0f} | {label} | {ratio:.2f} | {'yes: the tyre slides' if ratio > MU else 'no'} |")


if __name__ == "__main__":
    main()
