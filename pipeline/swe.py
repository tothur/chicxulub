"""Nonlinear shallow-water simulation of the resurge and the tsunami, 10 min to 10 h after impact.

Grid: the 66 Ma Gulf grid of gulf66.py (1024^2, 2.44 km), with the crater as it stands at 10 min (the viewer's
kinematic crater model, ported below) inside the block. Scheme: staggered C-grid, explicit, flux form
(as in TUNAMI-N2 / COMCOT), Manning friction, wetting and drying, sponge layers at the edges. Advection terms are
left out, so bores are smeared; heights near shore and inundation are limited by the 0.2-degree paleo-topography.

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
A0, SIG = float(os.environ.get("SWE_A0", 465)), float(os.environ.get("SWE_SIG", 12e3))
ring = A0 * fd * np.exp(-((R - 220e3) / SIG) ** 2) * wet
eta = eta + ring
c0 = np.sqrt(g * np.maximum(eta - z, 0.01))
q = ring * c0                                                            # outward volume flux of a long wave
Mx = np.zeros((N, N + 1)); Nz = np.zeros((N + 1, N))
Mx[:, 1:N] = 0.5 * ((q * X / np.maximum(R, 1))[:, :-1] + (q * X / np.maximum(R, 1))[:, 1:])
Nz[1:N, :] = 0.5 * ((q * Zs / np.maximum(R, 1))[:-1, :] + (q * Zs / np.maximum(R, 1))[1:, :])

SP = 48                                                                   # sponge width (cells)
k = np.minimum(np.minimum(np.arange(N), np.arange(N)[::-1])[None, :], np.minimum(np.arange(N), np.arange(N)[::-1])[:, None])
damp = np.where(k < SP, 1 - 0.03 * ((SP - k) / SP) ** 2, 1.0)
dampM = np.concatenate([damp, damp[:, -1:]], axis=1)                      # on x-faces
dampN = np.concatenate([damp, damp[-1:, :]], axis=0)                      # on z-faces


@njit(parallel=True, fastmath=True, cache=True)
def momentum(eta, z, Mx, Nz, dt, dx, g, n2sea, n2land, vmax):
    n = eta.shape[0]
    for j in prange(n):
        for i in range(1, n):
            el = eta[j, i - 1]; er = eta[j, i]; zl = z[j, i - 1]; zr = z[j, i]
            D = max(el, er) - max(zl, zr)
            if D < 0.05:
                Mx[j, i] = 0.0; continue
            m = Mx[j, i] - g * D * dt / dx * (er - el)
            n2 = n2land if max(zl, zr) > 0.0 else n2sea
            m = m / (1.0 + dt * g * n2 * abs(m) / D ** 2.3333)
            cap = vmax * D
            m = min(max(m, -cap), cap)
            if m > 0 and el - zl < 0.05: m = 0.0
            if m < 0 and er - zr < 0.05: m = 0.0
            Mx[j, i] = m
    for j in prange(1, n):
        for i in range(n):
            el = eta[j - 1, i]; er = eta[j, i]; zl = z[j - 1, i]; zr = z[j, i]
            D = max(el, er) - max(zl, zr)
            if D < 0.05:
                Nz[j, i] = 0.0; continue
            m = Nz[j, i] - g * D * dt / dx * (er - el)
            n2 = n2land if max(zl, zr) > 0.0 else n2sea
            m = m / (1.0 + dt * g * n2 * abs(m) / D ** 2.3333)
            cap = vmax * D
            m = min(max(m, -cap), cap)
            if m > 0 and el - zl < 0.05: m = 0.0
            if m < 0 and er - zr < 0.05: m = 0.0
            Nz[j, i] = m


@njit(parallel=True, fastmath=True, cache=True)
def mass(eta, z, Mx, Nz, dt, dx, damp, dampM, dampN, eta0, t, arr, mx):
    n = eta.shape[0]
    for j in prange(n):
        for i in range(n):
            e = eta[j, i] - dt / dx * (Mx[j, i + 1] - Mx[j, i] + Nz[j + 1, i] - Nz[j, i])
            if e < z[j, i]: e = z[j, i]
            if z[j, i] < 0.0 and e > 0.8 * (-z[j, i]) + 30.0: e = 0.8 * (-z[j, i]) + 30.0     # breaking: a crest cannot stand much higher than the water is deep
            d = damp[j, i]
            if d < 1.0: e = z[j, i] + (e - z[j, i]) if z[j, i] > 0 else e * d
            eta[j, i] = e
            depth = e - z[j, i]
            if depth > 0.05:
                if arr[j, i] > 1e9 and (abs(e - eta0[j, i]) > 1.0 or (z[j, i] > 0 and depth > 0.3)):
                    arr[j, i] = t
                v = depth if z[j, i] > 0 else e          # land: flood depth; sea: surface height
                if v > mx[j, i]: mx[j, i] = v
    for j in prange(n):
        for i in range(n + 1):
            Mx[j, i] *= dampM[j, i]
    for j in prange(n + 1):
        for i in range(n):
            Nz[j, i] *= dampN[j, i]


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
t = T0; step = 0; DF, DN = 240.0, 120.0
nf, nn = 0, 0
snapshot(t, True, True); nf, nn = 1, 1
t_wall = time.time()
while t < T1 - 1e-6:
    Dmax = float(np.max(eta - z))
    dt = min(0.45 * dx / math.sqrt(g * max(Dmax, 1.0)), 10.0)
    nxt = min(T0 + nf * DF, T0 + nn * DN)
    dt = min(dt, nxt - t)
    momentum(eta, z, Mx, Nz, dt, dx, g, 0.025 ** 2, 0.045 ** 2, 60.0)
    mass(eta, z, Mx, Nz, dt, dx, damp, dampM, dampN, eta0, t + dt, arr, mx)
    t += dt; step += 1
    if abs(t - (T0 + nn * DN)) < 1e-3 or t >= T1 - 1e-6:
        df = abs(t - (T0 + nf * DF)) < 1e-3
        snapshot(t, df, True); nn += 1
        if df: nf += 1
    if step % 500 == 0:
        print(f"t={t/3600:5.2f} h  step {step}  dt {dt:4.1f}s  wall {time.time()-t_wall:5.0f}s  frames {len(far)}/{len(near)}", flush=True)
print("done", step, "steps", round(time.time() - t_wall), "s", len(far), "far", len(near), "near")

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
