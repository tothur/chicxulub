"""Nonlinear shallow-water simulation of the resurge and the tsunami, 10 min to 10 h after impact.

Grid: the 66 Ma Gulf grid of gulf66.py (1024^2, 2.44 km), with the crater as it stands at 10 min (the viewer's
kinematic crater model, ported below) inside the block. Scheme: finite volumes on cell centres, second order
(limited slopes, Heun time stepping), HLL fluxes with the hydrostatic reconstruction of Audusse et al. 2004, so the
momentum equation is complete (bores form and lose energy as bores do), a sea at rest stays at rest against dry land,
and no water is created or thrown away: steep waves break by themselves, there is no limit on their height.
Manning friction, wetting and drying, open edges with a sponge. Heights near shore and inundation are limited by the
2.4 km grid and the 0.1-degree paleo-topography. pipeline/swe_test.py checks the solver against known solutions.

Initial state at 10 min (after Range et al. 2022): the sea is pushed back to 130 km from the impact point,
leaving the crater dry; a ring wave ~1.5 km high (less over the shallow shelf) at 220 km moves outward.

Writes sim/data/swe_far_*.png (512^2 frames every 240 s over the whole Gulf), swe_near_*.png (160^2 frames every
120 s around the crater, full resolution), rewrites gulf_w.png (first arrival, max height / flood depth) and
adds meta.gulf.swe.

usage: python swe.py <sim data dir>
"""
import json, math, os, sys, time
import numpy as np
from numba import njit, prange
from PIL import Image
from scipy.ndimage import map_coordinates

OUT = sys.argv[1]
meta = json.load(open(os.path.join(OUT, "meta.json")))
G = meta["gulf"]; N = G["N"]; HALF = G["half"]; dx = G["cell"]
g = 9.81

a = np.asarray(Image.open(os.path.join(OUT, "gulf_h.png"))).astype(np.int32)
z = ((a[..., 0] * 256 + a[..., 1]) * G["h"]["step"] + G["h"]["off"]).astype(np.float64)
u = -HALF + (np.arange(N) + 0.5) * dx
X, Zs = np.meshgrid(u, u)                     # x east, z south (row 0 = north edge)
R = np.hypot(X, Zs)

# ---- the crater at t = 600 s (mirrors impact(600), Wd, dexAt and deposits in the viewer) ----
def smooth(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1); return t * t * (3 - 2 * t)
I = dict(Hb=3500, Rb=62000, Cu=21000, Rc=22000, Lu=15000, U=300, Rr=95000, kC=1.0, dex=12000, Rex=46000,
         pkA=500, pkR=42000, melt=1.0, pool=1.0, brec=1.0)
def terrace(r):
    s = sum(1 - smooth(56000 + 10000 * k - 700, 56000 + 10000 * k + 700, r) for k in range(4))
    tilt = np.where((r > 56000) & (r < 96000), 0.25 * (((r - 56000) / 10000) % 1), 0)
    return (s + tilt) * 450 * smooth(44000, 54000, r)
def Wd(r, d):
    x = r / I["Rb"]; bowl = np.where(x < 1, np.power(np.clip(1 - x * x, 0, 1), 1.2), 0)
    y = (r - I["Rr"]) / (0.22 * I["Rr"] + 4000); c = r / I["Rc"]
    return (-I["Hb"] * bowl * np.exp(-d / 30000) + I["U"] * np.exp(-y * y) * np.exp(-d / 15000)
            + I["Cu"] * np.exp(-c * c) * np.exp(-d / I["Lu"]) - I["kC"] * terrace(r) * np.exp(-d / 14000))
dex = I["dex"] * np.maximum(0, 1 - (R / I["Rex"]) ** 2)
pk = I["pkA"] * np.exp(-((R - I["pkR"]) / 8000) ** 2)
bh = meta["half"]
inb = (np.abs(X) < bh) & (np.abs(Zs) < bh)
sf0 = np.where(inb, z, 0.0)                                             # 66 Ma seafloor inside the block
S = sf0 - dex + Wd(R, dex)
pool = np.clip(sf0 - 1300 - S, 0, 3000) * (1 - smooth(30000, 40000, R)) + 100 * (1 - smooth(50000, 58000, R))
inner = 250 * (1 - smooth(86000, 98000, R))
ej = 300 * (np.maximum(R, 90000) / 95000) ** -3 * smooth(82000, 96000, R) * (1 - smooth(110000, 165000, R))   # fades out inside the block
delta = -dex + Wd(R, dex) + pk + pool + inner + ej              # change of the surface from the 66 Ma seafloor
bh = meta["half"]
inb = (np.abs(X) < bh) & (np.abs(Zs) < bh)
z = np.where(inb, z + delta * np.clip(np.minimum(bh - np.abs(X), bh - np.abs(Zs)) / 15e3, 0, 1), z)

# ---- initial state at t0 = 600 s ----
T0, T1 = 600.0, 600.0 + float(os.environ.get("SWE_HOURS", 10)) * 3600
D0 = np.maximum(-z, 0)
eta = np.where((z < 0) & (R > 130e3), 0.0, z).astype(np.float64)        # pushed-back sea; dry crater
wet = eta > z + 0.01
fd = np.clip(np.sqrt((D0 + 30) / 1500), 0.15, 1.0)
A0, SIG = float(os.environ.get("SWE_A0", 900)), float(os.environ.get("SWE_SIG", 12e3))
ring = A0 * fd * np.exp(-((R - 220e3) / SIG) ** 2) * wet
eta = eta + ring
c0 = np.sqrt(g * np.maximum(eta - z, 0.01))
q = ring * c0                                                            # outward volume flux of a long wave
h = np.maximum(eta - z, 0.0)                                             # water depth; the state is (h, hu, hv) at cell centres
hu = q * X / np.maximum(R, 1); hv = q * Zs / np.maximum(R, 1)

SP = 48                                                                   # sponge width (cells)
k = np.minimum(np.minimum(np.arange(N), np.arange(N)[::-1])[None, :], np.minimum(np.arange(N), np.arange(N)[::-1])[:, None])
sponge = np.where(k < SP, 0.006 * ((SP - k) / SP) ** 2, 0.0)             # relaxation rate toward the sea at rest, 1/s

# ---- solver: finite volumes, second order (kept free of the grid above; swe_test.py runs it on test cases) ----
g = 9.81
HD = 0.05          # m: thinner water does not move (its velocity is not defined)
UMAX = 250.0       # m/s: cap on the speed of thin, fast sheets (keeps the time step finite; never reached by the waves)
TH = 1.5           # slope limiter: 1 = minmod (most diffusive), 2 = monotonised central

@njit(inline="always")
def lim(a, b):                                                            # limited slope from the one-sided differences
    if a > 0.0 and b > 0.0: return min(TH * a, 0.5 * (a + b), TH * b)
    if a < 0.0 and b < 0.0: return max(TH * a, 0.5 * (a + b), TH * b)
    return 0.0

@njit(inline="always")
def hll(hL, uL, vL, hR, uR, vR):                                          # HLL flux of (h, h un, h ut) through a face; un normal to it
    if hL <= 0.0 and hR <= 0.0: return 0.0, 0.0, 0.0
    cL = math.sqrt(g * hL); cR = math.sqrt(g * hR)
    if hL <= 0.0: sL = uR - 2.0 * cR; sR = uR + cR                        # dry bed on one side: the front runs at u -+ 2c
    elif hR <= 0.0: sL = uL - cL; sR = uL + 2.0 * cL
    else: sL = min(uL - cL, uR - cR); sR = max(uL + cL, uR + cR)
    a0 = hL * uL; a1 = a0 * uL + 0.5 * g * hL * hL; a2 = a0 * vL
    b0 = hR * uR; b1 = b0 * uR + 0.5 * g * hR * hR; b2 = b0 * vR
    if sL >= 0.0: return a0, a1, a2
    if sR <= 0.0: return b0, b1, b2
    w = 1.0 / (sR - sL)
    return ((sR * a0 - sL * b0 + sL * sR * (hR - hL)) * w, (sR * a1 - sL * b1 + sL * sR * (b0 - a0)) * w,
            (sR * a2 - sL * b2 + sL * sR * (hR * vR - hL * vL)) * w)

@njit(inline="always")
def speed(d, mu, mv):                                                     # velocity of a cell and its fastest wave
    if d > HD:
        a = mu / d; b = mv / d; sp = math.sqrt(a * a + b * b)
        if sp > UMAX: a *= UMAX / sp; b *= UMAX / sp; sp = UMAX
        return a, b, sp + math.sqrt(g * d)
    return 0.0, 0.0, 0.0

@njit(parallel=True, fastmath=True, cache=True)
def prim(h, hu, hv, u, v, rowmax):
    """Velocities and the fastest wave per row, for a state that advance() did not write."""
    n = h.shape[0]
    for j in prange(n):
        m = 0.0
        for i in range(n):
            u[j, i], v[j, i], c = speed(h[j, i], hu[j, i], hv[j, i])
            if c > m: m = c
        rowmax[j] = m

# The work arrays keep each cell's (or face's) numbers side by side: S[j, i, 0..3] the limited slopes of eta, h, u, v
# along x and S[j, i, 4..7] along z; FX[j, i, 0..2] the flux of mass, normal and transverse momentum through the face
# between cells i-1 and i, FX[j, i, 3..4] what the hydrostatic cut takes from the pressure on its low-index and
# high-index side (FZ[j, i] likewise for the face between rows j-1 and j). One pass then reads one stream of memory
# instead of eight; the run is limited by memory traffic, not arithmetic.
@njit(parallel=True, fastmath=True, cache=True)
def slopes(h, zb, u, v, S):
    n = h.shape[0]
    for j in prange(n):
        for i in range(n):
            if h[j, i] <= HD:                                             # dry cells stay flat: their bed is not tilted into the water
                for q in range(8): S[j, i, q] = 0.0
                continue
            e = h[j, i] + zb[j, i]
            if 0 < i < n - 1:
                S[j, i, 0] = lim(e - h[j, i - 1] - zb[j, i - 1], h[j, i + 1] + zb[j, i + 1] - e)
                S[j, i, 1] = lim(h[j, i] - h[j, i - 1], h[j, i + 1] - h[j, i])
                S[j, i, 2] = lim(u[j, i] - u[j, i - 1], u[j, i + 1] - u[j, i])
                S[j, i, 3] = lim(v[j, i] - v[j, i - 1], v[j, i + 1] - v[j, i])
            else:
                S[j, i, 0] = 0.0; S[j, i, 1] = 0.0; S[j, i, 2] = 0.0; S[j, i, 3] = 0.0
            if 0 < j < n - 1:
                S[j, i, 4] = lim(e - h[j - 1, i] - zb[j - 1, i], h[j + 1, i] + zb[j + 1, i] - e)
                S[j, i, 5] = lim(h[j, i] - h[j - 1, i], h[j + 1, i] - h[j, i])
                S[j, i, 6] = lim(u[j, i] - u[j - 1, i], u[j + 1, i] - u[j, i])
                S[j, i, 7] = lim(v[j, i] - v[j - 1, i], v[j + 1, i] - v[j, i])
            else:
                S[j, i, 4] = 0.0; S[j, i, 5] = 0.0; S[j, i, 6] = 0.0; S[j, i, 7] = 0.0

@njit(parallel=True, fastmath=True, cache=True)
def fluxes(h, zb, u, v, S, FX, FZ):
    """Face fluxes with hydrostatic reconstruction (Audusse et al. 2004): the water on each side is cut off at the higher
    of the two bed levels, so a sea at rest stays at rest, also against dry land."""
    n = h.shape[0]
    for j in prange(n):
        for i in range(n + 1):                                            # x-faces: face i lies between cells i-1 and i
            l = max(i - 1, 0); r = min(i, n - 1)                          # at the edges both sides are the edge cell (open boundary)
            hL = h[j, l] + 0.5 * S[j, l, 1]; eL = h[j, l] + zb[j, l] + 0.5 * S[j, l, 0]
            hR = h[j, r] - 0.5 * S[j, r, 1]; eR = h[j, r] + zb[j, r] - 0.5 * S[j, r, 0]
            zs = max(eL - hL, eR - hR)
            a = max(0.0, eL - zs); b = max(0.0, eR - zs)
            f0, f1, f2 = hll(a, u[j, l] + 0.5 * S[j, l, 2], v[j, l] + 0.5 * S[j, l, 3], b, u[j, r] - 0.5 * S[j, r, 2], v[j, r] - 0.5 * S[j, r, 3])
            FX[j, i, 0] = f0; FX[j, i, 1] = f1; FX[j, i, 2] = f2
            FX[j, i, 3] = 0.5 * g * (hL * hL - a * a); FX[j, i, 4] = 0.5 * g * (hR * hR - b * b)
    for j in prange(n + 1):
        jf = np.int64(j); jl = jf - 1 if jf > 0 else 0; jr = jf if jf < n else n - 1     # signed, and under their own names: the parallel loop's index is unsigned, and mixing the two makes floats
        for i in range(n):                                                # z-faces: face j lies between rows j-1 and j
            hL = h[jl, i] + 0.5 * S[jl, i, 5]; eL = h[jl, i] + zb[jl, i] + 0.5 * S[jl, i, 4]
            hR = h[jr, i] - 0.5 * S[jr, i, 5]; eR = h[jr, i] + zb[jr, i] - 0.5 * S[jr, i, 4]
            zs = max(eL - hL, eR - hR)
            a = max(0.0, eL - zs); b = max(0.0, eR - zs)
            f0, f1, f2 = hll(a, v[jl, i] + 0.5 * S[jl, i, 7], u[jl, i] + 0.5 * S[jl, i, 6], b, v[jr, i] - 0.5 * S[jr, i, 7], u[jr, i] - 0.5 * S[jr, i, 6])
            FZ[j, i, 0] = f0; FZ[j, i, 1] = f1; FZ[j, i, 2] = f2
            FZ[j, i, 3] = 0.5 * g * (hL * hL - a * a); FZ[j, i, 4] = 0.5 * g * (hR * hR - b * b)

@njit(parallel=True, fastmath=True, cache=True)
def advance(h, hu, hv, zb, S, FX, FZ, dt, dx, w, h0, hu0, hv0, ho, huo, hvo, u, v, rowmax, clip,
            fin, n2sea, n2land, sponge, eta0, t, arr, mx, lost):
    """out = w * U0 + (1 - w) * (U + dt * L(U)); clip[j] collects water added where a depth came out negative.
    With fin (the step's last stage, in a real run) the same pass also does what follows a step, cell by cell: Manning
    friction, the sponge at the open edges (lost[j]: water it takes out), first arrival and maximum.
    It leaves the new state's velocities in u, v and its fastest wave per row in rowmax, for the next stage."""
    n = h.shape[0]; k = dt / dx
    for j in prange(n):
        c = 0.0; s = 0.0; m = 0.0
        for i in range(n):
            e = h[j, i] + zb[j, i]
            hl = h[j, i] - 0.5 * S[j, i, 1]; hr = h[j, i] + 0.5 * S[j, i, 1]
            zl = e - 0.5 * S[j, i, 0] - hl; zr = e + 0.5 * S[j, i, 0] - hr
            hn = h[j, i] - 0.5 * S[j, i, 5]; hs = h[j, i] + 0.5 * S[j, i, 5]
            zn = e - 0.5 * S[j, i, 4] - hn; zz = e + 0.5 * S[j, i, 4] - hs
            d = h[j, i] - k * (FX[j, i + 1, 0] - FX[j, i, 0] + FZ[j + 1, i, 0] - FZ[j, i, 0])
            mu = hu[j, i] - k * (FX[j, i + 1, 1] + FX[j, i + 1, 3] - FX[j, i, 1] - FX[j, i, 4] - 0.5 * g * (hl + hr) * (zl - zr)
                                 + FZ[j + 1, i, 2] - FZ[j, i, 2])
            mv = hv[j, i] - k * (FZ[j + 1, i, 1] + FZ[j + 1, i, 3] - FZ[j, i, 1] - FZ[j, i, 4] - 0.5 * g * (hn + hs) * (zn - zz)
                                 + FX[j, i + 1, 2] - FX[j, i, 2])
            d = w * h0[j, i] + (1.0 - w) * d; mu = w * hu0[j, i] + (1.0 - w) * mu; mv = w * hv0[j, i] + (1.0 - w) * mv
            if d < 0.0: c -= d; d = 0.0
            if d <= HD: mu = 0.0; mv = 0.0
            if fin:
                zc = zb[j, i]
                if d > HD:
                    sp = math.sqrt(mu ** 2 + mv ** 2) / d
                    f = 1.0 / (1.0 + dt * g * (n2land if zc > 0.0 else n2sea) * sp / d ** 1.3333)
                    mu *= f; mv *= f
                if sponge[j, i] > 0.0 and zc < 0.0:
                    f = math.exp(-dt * sponge[j, i])
                    dn = max((d + zc) * f - zc, 0.0)
                    s += d - dn; d = dn; mu *= f; mv *= f
                if d > HD:
                    e = d + zc
                    if arr[j, i] > 1e9 and (abs(e - eta0[j, i]) > 1.0 or (zc > 0 and d > 0.3)):
                        arr[j, i] = t
                    val = d if zc > 0 else e       # land: flood depth; sea: surface height
                    if val > mx[j, i]: mx[j, i] = val
            ho[j, i] = d; huo[j, i] = mu; hvo[j, i] = mv
            u[j, i], v[j, i], cc = speed(d, mu, mv)
            if cc > m: m = cc
        clip[j] += c; rowmax[j] = m
        if fin: lost[j] += s

class Solver:
    """One step = Heun's second-order Runge-Kutta over the reconstruction above. CFL 0.3 on the fastest wave.
    after: (n2sea, n2land, sponge, eta0, arr, mx, lost) switches on friction, the edge sponge and the records of first
    arrival and maximum. The state must only be changed through step() once it has been called."""
    def __init__(self, n, after=None):
        self.u = np.zeros((n, n)); self.v = np.zeros((n, n)); self.S = np.zeros((n, n, 8)); self.rowmax = np.zeros(n)
        self.FX = np.zeros((n, n + 1, 5)); self.FZ = np.zeros((n + 1, n, 5)); self.clip = np.zeros(n)
        self.h1 = np.zeros((n, n)); self.hu1 = np.zeros((n, n)); self.hv1 = np.zeros((n, n)); self.u1 = np.zeros((n, n)); self.v1 = np.zeros((n, n))
        self.started = False; self.rm1 = np.zeros(n); self.cfl = 0.3
        e = np.zeros((1, 1)); self.fin = after is not None
        self.after = after if self.fin else (0.0, 0.0, e, e, e, e, np.zeros(1))
    def step(self, h, hu, hv, zb, dx, dt_max, t=0.0):
        if not self.started: prim(h, hu, hv, self.u, self.v, self.rowmax); self.started = True
        dt = min(self.cfl * dx / max(float(self.rowmax.max()), 1e-6), dt_max)
        a = self.after
        slopes(h, zb, self.u, self.v, self.S)
        fluxes(h, zb, self.u, self.v, self.S, self.FX, self.FZ)
        advance(h, hu, hv, zb, self.S, self.FX, self.FZ, dt, dx, 0.0, h, hu, hv, self.h1, self.hu1, self.hv1, self.u1, self.v1, self.rm1, self.clip,
                False, a[0], a[1], a[2], a[3], t + dt, a[4], a[5], a[6])
        slopes(self.h1, zb, self.u1, self.v1, self.S)
        fluxes(self.h1, zb, self.u1, self.v1, self.S, self.FX, self.FZ)
        advance(self.h1, self.hu1, self.hv1, zb, self.S, self.FX, self.FZ, dt, dx, 0.5, h, hu, hv, h, hu, hv, self.u, self.v, self.rowmax, self.clip,
                self.fin, a[0], a[1], a[2], a[3], t + dt, a[4], a[5], a[6])
        return dt
# ---- end of solver ----


ETA_MAX = 1600.0
def encode(e, w):                      # 0 = dry, else 128 + 127 sign sqrt(|eta| / ETA_MAX)
    v = 128 + 127 * np.sign(e) * np.sqrt(np.minimum(np.abs(e), ETA_MAX) / ETA_MAX)
    return np.where(w, np.clip(np.round(v), 1, 255), 0).astype(np.uint8)

NB = 160; c0i = N // 2 - NB // 2                                          # near window: 160 cells around the crater
far, near = [], []
def snapshot(t, do_far, do_near):
    w = eta > z + 0.05
    if do_far:
        ws = w.reshape(N // 2, 2, N // 2, 2); es = np.where(w, eta, 0).reshape(N // 2, 2, N // 2, 2)
        cnt = ws.sum(axis=(1, 3)); mean = es.sum(axis=(1, 3)) / np.maximum(cnt, 1)
        far.append(encode(mean, cnt > 0))
    if do_near:
        sl = (slice(c0i, c0i + NB), slice(c0i, c0i + NB))
        near.append(encode(eta[sl], w[sl]))

arr = np.full((N, N), 2e9); mx = np.full((N, N), -1e9)
eta0 = np.where(wet, 0.0, eta)
lost = np.zeros(N); sol = Solver(N, (0.025 ** 2, 0.045 ** 2, sponge, eta0, arr, mx, lost))      # Manning n: sea, land
sol.cfl = float(os.environ.get("SWE_CFL", sol.cfl))
V0 = float(h.sum()) * dx * dx; ring_km3 = float(ring.sum()) * dx * dx / 1e9
t = T0; step = 0; DF, DN = 240.0, 120.0
nf, nn = 0, 0
snapshot(t, True, True); nf, nn = 1, 1
t_wall = time.time()
while t < T1 - 1e-6:
    nxt = min(T0 + nf * DF, T0 + nn * DN)
    dt = sol.step(h, hu, hv, z, dx, nxt - t, t)
    t += dt; step += 1
    if abs(t - (T0 + nn * DN)) < 1e-3 or t >= T1 - 1e-6:
        eta = h + z
        df = abs(t - (T0 + nf * DF)) < 1e-3
        snapshot(t, df, True); nn += 1
        if df: nf += 1
    if step % int(os.environ.get("SWE_LOG", 1000)) == 0:
        print(f"t={t/3600:5.2f} h  step {step}  dt {dt:4.2f}s  wall {time.time()-t_wall:5.0f}s  frames {len(far)}/{len(near)}", flush=True)
eta = h + z
print("done", step, "steps", round(time.time() - t_wall), "s", len(far), "far", len(near), "near")
# water balance: nothing is created or thrown away inside the grid; the open edges and their sponge let it out (or in)
dV = (float(h.sum()) * dx * dx - V0) / 1e9
print(f"water: ring wave {ring_km3:.0f} km3; change over the run {dV:+.0f} km3, of which the edge sponge {-lost.sum() * dx * dx / 1e9:+.0f} km3, "
      f"negative depths set to zero {sol.clip.sum() * dx * dx / 1e9:+.3f} km3, the rest through the open edges")

def atlas(frames, tile, per_row, name):     # three frames per pixel in R, G, B; alpha stays opaque (browsers premultiply it)
    per_img = per_row * per_row * 3; files = []
    for k0 in range(0, len(frames), per_img):
        chunk = frames[k0:k0 + per_img]
        img = np.zeros((per_row * tile, per_row * tile, 3), np.uint8)
        for k, f in enumerate(chunk):
            ch, cell = k % 3, k // 3; r, c = divmod(cell, per_row)
            img[r * tile:(r + 1) * tile, c * tile:(c + 1) * tile, ch] = f
        fn = f"{name}_{k0 // per_img}.png"; Image.fromarray(img).save(os.path.join(OUT, fn), optimize=True); files.append(fn)
    return files

farf = atlas(far, N // 2, 4, "swe_far")
nearf = atlas(near, NB, 8, "swe_near")

# arrival and maximum, for the readout and the labels
land = z > 0
A = np.where(land, np.maximum(mx, 0) / 2, np.maximum(mx, 0) / 8)            # land: flood depth / 2 m; sea: max height / 8 m
qt = np.where(arr < 1e9, np.clip(np.round(arr / 2), 0, 65534), 65535).astype(np.uint32)
Image.fromarray(np.dstack([qt >> 8, qt & 255, np.clip(np.round(A), 0, 255)]).astype(np.uint8)).save(os.path.join(OUT, "gulf_w.png"), optimize=True)

E, Nn = X, -Zs
for L in G["labels"]:
    if L["kind"] != "town": continue
    kk = np.hypot(X - L["x"], Zs - L["z"]) < 80e3
    sea_k = kk & ~land & (arr < 1e9); land_k = kk & land & (arr < 1e9)
    L["arrive_s"] = round(float(arr[sea_k | land_k].min())) if (sea_k | land_k).any() else None
    L["height_m"] = round(float(mx[sea_k & (z > -200)].max())) if (sea_k & (z > -200)).any() else (round(float(mx[sea_k].max())) if sea_k.any() else None)
    if L["height_m"] is None and land_k.any(): L["height_m"] = round(float(mx[land_k].max()))      # no sea within 80 km: the deepest flood on land
    L["flood_km"] = round(float(np.sqrt(land_k.sum()) * dx / 1000)) if land_k.any() else 0
    print(L["text"], L["arrive_s"] and f'{L["arrive_s"]/3600:.2f} h', L["height_m"], "m")
deep = (z < -2000) & (R > 600e3) & (R < 1000e3) & (arr < 1e9)
print("open Gulf 600-1000 km, deep water: max height median", round(float(np.median(mx[deep]))), "m, 90th", round(float(np.percentile(mx[deep], 90))), "m")
flooded = land & (mx > 0.5)
print("flooded land:", round(flooded.sum() * dx * dx / 1e6), "km2")
G["swe"] = {"t0": T0, "t1": T1, "far": {"files": farf, "n": len(far), "dt": DF, "size": N // 2, "per_row": 4},
            "near": {"files": nearf, "n": len(near), "dt": DN, "size": NB, "per_row": 8, "half": NB / 2 * dx},
            "eta_max": ETA_MAX, "w_land_step": 2, "w_sea_step": 8}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
