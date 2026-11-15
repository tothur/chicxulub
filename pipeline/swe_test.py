"""Checks of the shallow-water solver in swe.py against cases with known answers. Not part of the build.

1. A sea at rest over rough ground with islands stays at rest.
2. Dam break onto a dry bed: Ritter's solution (depth profile, front position).
3. Dam break onto shallow water: Stoker's solution (height and speed of the bore).
4. A hump of water collapsing over rough ground with shores, walled in: the volume of water does not change.

usage: python swe_test.py        (exits non-zero if a check fails)
"""
import math, os, sys
import numpy as np
from numba import njit, prange

src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "swe.py")).read()
exec(src[src.index("# ---- solver"):src.index("# ---- end of solver")].replace(", cache=True", ""))          # the solver section only (code run from a string cannot be cached)
g = 9.81
ok = True
def check(name, value, limit, unit=""):
    global ok
    good = value <= limit; ok &= good
    print(f"  {'ok  ' if good else 'FAIL'} {name}: {value:.3g}{unit} (limit {limit:g}{unit})")

def run(h, hu, hv, zb, dx, T):
    sol = Solver(h.shape[0]); t = 0.0
    while t < T - 1e-9: t += sol.step(h, hu, hv, zb, dx, T - t)
    return sol

n = 256; rng = np.random.default_rng(1)
def rough(amp):                                    # smooth random ground, m
    from scipy.ndimage import gaussian_filter
    a = gaussian_filter(rng.standard_normal((n, n)), 6); return amp * a / np.abs(a).max()

print("1. sea at rest over rough ground with islands")
zb = rough(1500.0) - 400.0                        # sea to ~1,900 m deep, islands to ~1,100 m
h = np.maximum(-zb, 0.0); hu = np.zeros((n, n)); hv = np.zeros((n, n))
wet0 = h > 0
run(h, hu, hv, zb, 2440.0, 3600.0)
check("largest change of the sea surface", float(np.abs((h + zb)[wet0]).max()), 1e-6, " m")
check("largest speed", float(np.abs(hu[wet0] / np.maximum(h[wet0], 1e-9)).max()), 1e-6, " m/s")
check("land that got wet", float(((h > 0) & ~wet0).sum()), 0, " cells")

print("2. dam break onto a dry bed (Ritter)")
dx = 10.0; x = (np.arange(n) + 0.5 - n / 2) * dx; h0 = 10.0; T = 30.0; c0 = math.sqrt(g * h0)
zb = np.zeros((n, n)); h = np.where(x < 0, h0, 0.0)[None, :].repeat(n, 0).copy(); hu = np.zeros((n, n)); hv = np.zeros((n, n))
V = h.sum(); sol = run(h, hu, hv, zb, dx, T)
exact = np.where(x < -c0 * T, h0, np.where(x > 2 * c0 * T, 0.0, (2 * c0 - x / T) ** 2 / (9 * g)))
check("mean depth error", float(np.abs(h[n // 2] - exact).mean()), 0.06, " m")
check("volume change", abs(float(h.sum() - V)) / V, 1e-12)
check("variation across the flow", float(np.abs(h - h[n // 2]).max()), 1e-9, " m")

print("3. dam break onto shallow water (Stoker)")
hl, hr = 10.0, 1.0
f = lambda hm: 2 * (math.sqrt(g * hl) - math.sqrt(g * hm)) - (hm - hr) * math.sqrt(0.5 * g * (hm + hr) / (hm * hr))
lo, hi = hr, hl
for _ in range(80):
    mid = 0.5 * (lo + hi)
    if f(mid) > 0: lo = mid
    else: hi = mid
hm = 0.5 * (lo + hi); um = 2 * (math.sqrt(g * hl) - math.sqrt(g * hm)); sbore = hm * um / (hm - hr)
h = np.where(x < 0, hl, hr)[None, :].repeat(n, 0).copy(); hu = np.zeros((n, n)); hv = np.zeros((n, n))
run(h, hu, hv, zb, dx, T)
row = h[n // 2]; front = x[np.argmin(np.abs(row - 0.5 * (hm + hr)) + 1e9 * (x < 0))]
plateau = row[(x > 0.3 * sbore * T) & (x < 0.8 * sbore * T)]
check("height behind the bore", float(np.abs(plateau - hm).max()), 0.02 * hm, " m")
check("bore position", abs(float(front) - sbore * T), 1.5 * dx, " m")

print("4. a hump of water collapsing over rough ground with shores, inside a wall of land")
zb = rough(300.0) - 60.0
zb[:4] = zb[-4:] = 2000.0; zb[:, :4] = zb[:, -4:] = 2000.0                 # nothing can leave
hump = np.exp(-((np.arange(n)[None, :] - n / 2) ** 2 + (np.arange(n)[:, None] - n / 2) ** 2) / 400.0)
h = np.maximum(200.0 * hump - zb, 0.0) * (hump > 0.05) + np.maximum(-zb, 0.0) * (hump <= 0.05)
hu = np.zeros((n, n)); hv = np.zeros((n, n)); V = h.sum(); dry0 = int((h == 0).sum())
sol = run(h, hu, hv, zb, 500.0, 900.0)
check("volume change", abs(float(h.sum() - V)) / V, 1e-12)
check("negative depths set to zero", float(sol.clip.sum()) / V, 1e-12)
moved = float(np.abs(hu).max())
print(f"       (the water did move: fastest flux {moved:.0f} m2/s, {dry0 - int((h == 0).sum())} dry cells got wet)")
ok &= moved > 1.0

print("all checks passed" if ok else "CHECKS FAILED")
sys.exit(0 if ok else 1)
