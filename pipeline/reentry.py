"""Distal ejecta re-entry: where and when the fast ejecta fall back into the atmosphere, and how hot it gets.

Particles leave the impact site above the atmosphere at 2-11 km/s and 35-55 deg elevation, with the cratering
mass-velocity law M(>v) ~ v^-1.65 and a mild excess downrange (to the SW, for an impact from the NE). They fly
Kepler (two-body) orbits in an inertial frame; the Earth turns under them (the paleo frame's pole is the spin axis).
Re-entry is where they cross 100 km altitude. Particles faster than escape are lost.

Mass is scaled so that the distal layer (5,000-15,000 km from the crater) averages 3 mm (density 1,500 kg/m3),
the typical thickness of the distal K-Pg ejecta layer. Heating of the upper atmosphere = kinetic energy of the
re-entering mass per area and time. Ground irradiance depends on how much of that radiates down and on
self-shielding by the falling spherules (Melosh et al. 1990; Goldin & Melosh 2009); the viewer reports both.

Writes sim/data/reentry_f_*.png (360x180 frames every 240 s for 0-8 h, log10 W/m2 over 7 decades, three per pixel in RGB),
reentry_m.png (deposit thickness, log scale) and meta.reentry.

usage: python reentry.py <sim data dir>
"""
import json, math, os, sys, time
import numpy as np
from numba import njit, prange
from PIL import Image

OUT = sys.argv[1]
meta = json.load(open(os.path.join(OUT, "meta.json")))
LON0, LAT0 = meta["globe"]["impact"]["lon"], meta["globe"]["impact"]["lat"]
RE, GM, OMEGA = 6371e3, 3.986004e14, 7.2921e-5
R_IN = RE + 100e3                                  # re-entry interface
NP = 2_000_000
rng = np.random.default_rng(66)

# launch speeds from M(>v) ~ v^-1.65 between 2 and 11 km/s
k = 1.65; vmin, vmax = 2000.0, 11000.0
u = rng.random(NP)
v = (vmin ** -k - u * (vmin ** -k - vmax ** -k)) ** (-1 / k)
elev = np.radians(35 + 20 * rng.random(NP))
# azimuth (clockwise from north) with a downrange excess toward 225 deg (SW)
az = np.empty(NP); m = 0
while m < NP:
    a = rng.random(NP) * 2 * np.pi; keep = rng.random(NP) < (1 + 0.25 * np.cos(a - np.radians(225))) / 1.25
    a = a[keep][: NP - m]; az[m:m + len(a)] = a; m += len(a)

# launch site and local frame (x toward lon 0 on the equator, z = spin axis, y toward lon 90 E)
la, lo = math.radians(LAT0), math.radians(LON0)
up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
east = np.array([-math.sin(lo), math.cos(lo), 0.0])
north = np.cross(up, east)
r0 = up * R_IN
dirs = (np.cos(elev)[:, None] * (np.sin(az)[:, None] * east + np.cos(az)[:, None] * north) + np.sin(elev)[:, None] * up)
vel = v[:, None] * dirs + np.cross([0, 0, OMEGA], r0)[None, :]          # plus the Earth's own spin at the site


@njit(parallel=True, fastmath=True, cache=True)
def fly(vel, r0, dt, tmax, out):
    n = vel.shape[0]
    for i in prange(n):
        x, y, z = r0[0], r0[1], r0[2]; vx, vy, vz = vel[i, 0], vel[i, 1], vel[i, 2]
        t = 0.0; out[i, 0] = -1.0
        e = 0.5 * (vx * vx + vy * vy + vz * vz) - GM / math.sqrt(x * x + y * y + z * z)
        if e >= 0.0: continue                                              # escapes
        while t < tmax:
            r = math.sqrt(x * x + y * y + z * z); a = -GM / (r * r * r)
            vx += 0.5 * dt * a * x; vy += 0.5 * dt * a * y; vz += 0.5 * dt * a * z
            x += dt * vx; y += dt * vy; z += dt * vz
            r = math.sqrt(x * x + y * y + z * z); a = -GM / (r * r * r)
            vx += 0.5 * dt * a * x; vy += 0.5 * dt * a * y; vz += 0.5 * dt * a * z
            t += dt
            if t > 30.0 and r < R_IN:
                out[i, 0] = t; out[i, 1] = x; out[i, 2] = y; out[i, 3] = z
                out[i, 4] = math.sqrt(vx * vx + vy * vy + vz * vz); break

out = np.zeros((NP, 5))
tw = time.time(); fly(vel, r0, 2.0, 12 * 3600.0, out); print("flown in", round(time.time() - tw), "s")
ok = out[:, 0] > 0
t, X, Y, Z, ve = out[ok, 0], out[ok, 1], out[ok, 2], out[ok, 3], out[ok, 4]
lat = np.degrees(np.arcsin(Z / np.sqrt(X * X + Y * Y + Z * Z)))
lon = np.degrees(np.arctan2(Y, X) - OMEGA * t)                             # into the rotating frame
lon = (lon + 180) % 360 - 180
print("re-entered", ok.sum(), "of", NP, " median time", round(float(np.median(t)) / 60), "min, 95th", round(float(np.percentile(t, 95)) / 60), "min")

# grids (1 deg) and great-circle distance from the crater
NX, NY = 360, 180
gl = -180 + (np.arange(NX) + 0.5); gb = 90 - (np.arange(NY) + 0.5)
GLo, GLa = np.meshgrid(gl, gb)
area = (RE * math.radians(1)) ** 2 * np.cos(np.radians(GLa))
def gc(lo1, la1, lo2, la2):
    p1, p2 = np.radians(la1), np.radians(la2)
    return RE * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(np.radians(lo2 - lo1)), -1, 1))
DIST = gc(LON0, LAT0, GLo, GLa)
ix = np.clip(((lon + 180)).astype(int), 0, NX - 1); iy = np.clip(((90 - lat)).astype(int), 0, NY - 1)

# mass per particle: the distal layer averages 3 mm at 1,500 kg/m3
count = np.zeros((NY, NX)); np.add.at(count, (iy, ix), 1.0)
band = (DIST > 5e6) & (DIST < 15e6)
w_per = (0.003 * 1500 * area[band].sum()) / count[band].sum()             # kg per particle
thick_mm = count * w_per / area / 1500 * 1000
print("mass per particle", f"{w_per:.3g} kg", " total re-entering", f"{w_per * ok.sum():.3g} kg")

# heating of the upper atmosphere, W/m2, in 240 s bins over 8 h
DT, NT = 240.0, 120
tb = (t // DT).astype(int); sel = tb < NT
E = np.zeros((NT, NY, NX)); np.add.at(E, (tb[sel], iy[sel], ix[sel]), w_per * 0.5 * ve[sel] ** 2)
Q = E / area[None] / DT
# each particle stands for a swarm (one carries ~180 kW/m2 into a 1-deg cell for 4 min far out): smooth out the sampling
# noise over about 2 deg in space, wrapping in longitude; time is left sharp so arrival times keep
from scipy.ndimage import gaussian_filter
Q = np.maximum(gaussian_filter(Q, sigma=(0.0, 2.0, 2.0), mode=('nearest', 'nearest', 'wrap')), 0)
peak = Q.max(0); tpk = Q.argmax(0) * DT + DT / 2
for name, d in [("1,000 km", 1e6), ("3,000 km", 3e6), ("8,000 km", 8e6), ("15,000 km", 15e6), ("antipode", 19.5e6)]:
    kk = np.abs(DIST - d) < 0.4e6 if d < 19e6 else DIST > 19e6
    print(f"{name:10s} peak heating {np.median(peak[kk]) / 1e3:8.1f} kW/m2 at {np.median(tpk[kk]) / 60:5.0f} min, deposit {np.median(thick_mm[kk]):7.2f} mm")

def enc(q):                                       # 0 = nothing, else 1 + 254 * log10(q) / 7 (1 W/m2 .. 10 MW/m2)
    return np.where(q >= 1, np.clip(np.round(1 + 254 * np.log10(np.maximum(q, 1)) / 7), 1, 255), 0).astype(np.uint8)
frames = [enc(Q[i]) for i in range(NT)]
PR = 4; per = PR * PR * 3; files = []
for k0 in range(0, NT, per):
    img = np.zeros((PR * NY, PR * NX, 3), np.uint8)
    for kk, f in enumerate(frames[k0:k0 + per]):
        ch, cell = kk % 3, kk // 3; r_, c_ = divmod(cell, PR)
        img[r_ * NY:(r_ + 1) * NY, c_ * NX:(c_ + 1) * NX, ch] = f
    fn = f"reentry_f_{k0 // per}.png"; Image.fromarray(img).save(os.path.join(OUT, fn), optimize=True); files.append(fn)
tm = np.where(thick_mm > 0.001, np.clip(np.round(1 + 254 * (np.log10(np.maximum(thick_mm, 0.001)) + 3) / 6), 1, 255), 0).astype(np.uint8)   # 0.001 mm .. 1 m
Image.fromarray(np.dstack([tm, tm, tm])).save(os.path.join(OUT, "reentry_m.png"), optimize=True)
meta["reentry"] = {"nx": NX, "ny": NY, "dt": DT, "n": NT, "per_row": PR, "files": files,
                   "enc": "1 + 254 * log10(W/m2) / 7", "thick": "1 + 254 * (log10(mm) + 3) / 6",
                   "mass_kg": float(w_per * ok.sum()), "layer_mm": 3.0}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
