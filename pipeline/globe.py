"""Globe overview: the 66 Ma Earth and a global tsunami, for the viewer's globe mode.

Topography: PALEOMAP 65 Ma (Scotese & Wright 2018, the 0.1 deg PaleoDEM), in its own (66 Ma) coordinates: averaged
to 0.25 deg for the tsunami run (topo.py writes the full 0.1 deg for drawing, globe_t.png). As in the Gulf grid,
low land within 450 km of the impact site is the Yucatán shelf (40 m of water).

Tsunami: linear shallow-water equations on the sphere (lat-lon C-grid, Manning friction, coasts as walls, a breaking
cap near shore), from 30 min to 48 h. It starts from an outgoing ring 400 km from the impact site whose height is
calibrated so the open Gulf sees the same crests as the regional simulation (swe.py). No inundation at this scale.

Writes sim/data/globe_h.png (elevation at
0.25 deg, same encoding; B = land), globe_w.png (first arrival (R*256+G)*4 s,
B = max height in a sqrt scale), globe_f_*.png (720x360 frames, three per pixel in RGB) and meta.globe.

usage: python globe.py <paleodem 65 Ma .nc> <sim data dir> [calibrate]
"""
import json, math, os, sys, time
import numpy as np
import netCDF4
from numba import njit, prange
from PIL import Image
from scipy.ndimage import map_coordinates

NC, OUT = sys.argv[1], sys.argv[2]
CAL = len(sys.argv) > 3
meta = json.load(open(os.path.join(OUT, "meta.json")))
fr = meta["gulf"]["frame"]
LON0, LAT0 = fr["paleo_lon"], fr["paleo_lat"]                     # the impact site at 66 Ma
RE, g = 6371e3, 9.81

NX, NY = 1440, 720
dl = 2 * math.pi / NX
lon = -180 + (np.arange(NX) + 0.5) * 360 / NX
lat = 90 - (np.arange(NY) + 0.5) * 180 / NY                       # row 0 = north
LON, LAT = np.meshgrid(lon, lat)

ds = netCDF4.Dataset(NC)
var = lambda pre: next(k for k in ds.variables if k.lower().startswith(pre))      # lon/lat or longitude/latitude
L, B, Z = ds[var("lon")][:], ds[var("lat")][:], np.asarray(ds["z"][:], dtype=np.float64) - 20.0  # PALEOMAP heights come in 40 m steps, its 0 m cells being sea: sea level is put midway between 0 and 40 m
from scipy.ndimage import gaussian_filter
dsrc = abs(float(L[1] - L[0]))
Zs = gaussian_filter(Z, 0.4 * 0.25 / dsrc, mode=("nearest", "wrap")) if dsrc < 0.2 else Z    # area-average a finer source to 0.25 deg
h = map_coordinates(Zs, [(LAT - B[0]) / (B[1] - B[0]), (LON - L[0]) / (L[1] - L[0])], order=1, mode="nearest")

def gc(lon1, lat1, lon2, lat2):                                    # great-circle distance (m)
    p1, p2, dlon = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    return RE * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dlon), -1, 1))
R = gc(LON0, LAT0, LON, LAT)
h[(R < 450e3) & (h > 0) & (h < 400)] = -40.0                       # the Yucatán shelf
land = h > 0
wall = land | (np.abs(LAT) > 78)
D = np.where(wall, 0.0, -h)

cosc = np.cos(np.radians(lat))                                     # at cell centres (per row)
cosn = np.cos(np.radians(lat + 90 / NY))                           # at the face north of each row
dphi = math.pi / NY

# ---- initial state at T0: an outgoing ring ----
T0, T1 = 1800.0, 48 * 3600.0
A0 = float(os.environ.get("GLOBE_A0", 205)); R0, SIG = 400e3, 60e3
fd = np.clip(np.sqrt((D + 30) / 1500), 0.15, 1.0)
eta = np.where(wall, 0.0, A0 * fd * np.exp(-((R - R0) / SIG) ** 2))
q = eta * np.sqrt(g * np.maximum(D, 1.0))
# outward direction from the distance field (east and south components)
dRx = (np.roll(R, -1, 1) - np.roll(R, 1, 1)) / (2 * RE * cosc[:, None] * dl)
dRy = (np.vstack([R[1:], R[-1:]]) - np.vstack([R[:1], R[:-1]])) / (2 * RE * dphi)   # toward increasing row = south
nrm = np.maximum(np.hypot(dRx, dRy), 1e-9)
qe, qs = q * dRx / nrm, q * dRy / nrm
Mx = 0.5 * (qe + np.roll(qe, -1, 1))                               # face east of each cell
Ny = np.zeros((NY + 1, NX)); Ny[1:NY] = 0.5 * (qs[:-1] + qs[1:])    # face north of row j (j = 1..NY-1), + = south


@njit(parallel=True, fastmath=True, cache=True)
def momentum(eta, D, wall, Mx, Ny, dt, cosc, n2):
    ny, nx = eta.shape
    for j in prange(ny):
        ax = g * dt / (RE * cosc[j] * dl)
        for i in range(nx):
            ii = (i + 1) % nx
            if wall[j, i] or wall[j, ii]:
                Mx[j, i] = 0.0; continue
            d = 0.5 * (D[j, i] + D[j, ii])
            m = Mx[j, i] - ax * d * (eta[j, ii] - eta[j, i])
            Mx[j, i] = m / (1.0 + dt * g * n2 * abs(m) / max(d, 1.0) ** 2.3333)
    for j in prange(1, ny):
        ay = g * dt / (RE * dphi)
        for i in range(nx):
            if wall[j, i] or wall[j - 1, i]:
                Ny[j, i] = 0.0; continue
            d = 0.5 * (D[j, i] + D[j - 1, i])
            m = Ny[j, i] - ay * d * (eta[j, i] - eta[j - 1, i])
            Ny[j, i] = m / (1.0 + dt * g * n2 * abs(m) / max(d, 1.0) ** 2.3333)


@njit(parallel=True, fastmath=True, cache=True)
def mass(eta, D, wall, Mx, Ny, dt, cosc, cosn, t, arr, mx):
    ny, nx = eta.shape
    for j in prange(ny):
        cn = cosn[j]; cs = cosn[j + 1] if j + 1 < ny else 0.0
        bx = dt / (RE * cosc[j] * dl); by = dt / (RE * cosc[j] * dphi)
        for i in range(nx):
            if wall[j, i]: continue
            e = eta[j, i] - bx * (Mx[j, i] - Mx[j, (i - 1) % nx]) - by * (Ny[j + 1, i] * cs - Ny[j, i] * cn)
            cap = 0.8 * D[j, i] + 30.0                               # breaking near shore
            if e > cap: e = cap
            if e < -D[j, i]: e = -D[j, i]
            eta[j, i] = e
            if arr[j, i] > 1e9 and abs(e) > 0.2: arr[j, i] = t
            if e > mx[j, i]: mx[j, i] = e

# cosn[j] = cos of the latitude of the face north of row j (index 0..NY)
cosn = np.cos(np.radians(90 - np.arange(NY + 1) * 180 / NY))

ETA_MAX = 400.0
def encode(e, w):
    v = 128 + 127 * np.sign(e) * np.sqrt(np.minimum(np.abs(e), ETA_MAX) / ETA_MAX)
    return np.where(w, np.clip(np.round(v), 1, 255), 0).astype(np.uint8)

frames = []; DF = 1200.0
def snap():
    w = ~wall
    ws = w.reshape(NY // 2, 2, NX // 2, 2); es = np.where(w, eta, 0).reshape(NY // 2, 2, NX // 2, 2)
    cnt = ws.sum(axis=(1, 3)); frames.append(encode(es.sum(axis=(1, 3)) / np.maximum(cnt, 1), cnt > 0))

arr = np.full((NY, NX), 2e9); mx = np.zeros((NY, NX))
cmax = math.sqrt(g * float(D.max()))
dxmin = RE * math.cos(math.radians(78)) * dl
dt0 = 0.45 * dxmin / cmax
t = T0; nf = 1; snap()
tw = time.time(); step = 0
Tend = T0 + 3 * 3600 if CAL else T1
while t < Tend - 1e-6:
    dt = min(dt0, T0 + nf * DF - t)
    momentum(eta, D, wall, Mx, Ny, dt, cosc, 0.025 ** 2)
    mass(eta, D, wall, Mx, Ny, dt, cosc, cosn, t + dt, arr, mx)
    t += dt; step += 1
    if abs(t - (T0 + nf * DF)) < 1e-3:
        snap(); nf += 1
        if nf % 12 == 0: print(f"t={t/3600:5.1f} h  steps {step}  wall {time.time()-tw:4.0f}s", flush=True)
print("done", step, "steps", round(time.time() - tw), "s,", len(frames), "frames")

gulf = (R > 600e3) & (R < 1000e3) & (D > 2000) & (np.abs(LON - LON0) < 15) & (np.abs(LAT - LAT0) < 12)
print("open Gulf 600-1000 km: max height median", round(float(np.median(mx[gulf]))), "m")
if CAL: sys.exit()

def pack16(x, step, off):
    q = np.round((x - off) / step); assert q.min() >= 0 and q.max() < 65536; return q.astype(np.uint32)
qh = pack16(np.clip(h, -11000, 9000), 0.5, -11000)
Image.fromarray(np.dstack([qh >> 8, qh & 255, land * 255]).astype(np.uint8)).save(os.path.join(OUT, "globe_h.png"), optimize=True)
# (the full-resolution topography for drawing, globe_t.png, is written by topo.py)
qt = np.where(arr < 1e9, np.clip(np.round(arr / 4), 0, 65534), 65535).astype(np.uint32)
qa = np.clip(np.round(255 * np.sqrt(np.clip(mx, 0, ETA_MAX) / ETA_MAX)), 0, 255).astype(np.uint32)
Image.fromarray(np.dstack([qt >> 8, qt & 255, qa]).astype(np.uint8)).save(os.path.join(OUT, "globe_w.png"), optimize=True)

TW, TH, PR = NX // 2, NY // 2, 4                                   # 720x360 tiles, 4x4 per image, 3 frames per pixel
files = []
per = PR * PR * 3
for k0 in range(0, len(frames), per):
    img = np.zeros((PR * TH, PR * TW, 3), np.uint8)
    for k, f in enumerate(frames[k0:k0 + per]):
        ch, cell = k % 3, k // 3; r, c = divmod(cell, PR)
        img[r * TH:(r + 1) * TH, c * TW:(c + 1) * TW, ch] = f
    fn = f"globe_f_{k0 // per}.png"; Image.fromarray(img).save(os.path.join(OUT, fn), optimize=True); files.append(fn)

# ocean labels at 66 Ma: the deepest water near each guess
def ocean(name, lo, la):
    k = (np.abs(LON - lo) < 12) & (np.abs(LAT - la) < 10) & (D > 3000)
    if not k.any(): return None
    i = np.argmax(np.where(k, -np.hypot(LON - lo, LAT - la), -1e9)); return {"text": name, "lon": round(float(LON.flat[i]), 2), "lat": round(float(LAT.flat[i]), 2)}
labels = [x for x in [ocean("Gulf of Mexico", LON0 - 1, LAT0 + 5), ocean("North Atlantic", -42, 30), ocean("South Atlantic", -15, -28),
                      ocean("Pacific Ocean", -130, 5), ocean("Tethys Ocean", 45, 15), ocean("Indian Ocean", 75, -25), ocean("Southern Pacific", -120, -45)] if x]
# today's coastlines, moved to 66 Ma (cached GPlates output from gulf66.py), simplified to ~0.15 deg
from skimage.measure import approximate_polygon
today = []
cj = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "raw", "gws", "coastlines.json")     # the GPlates cache in the project
for f in json.load(open(cj)).get("features", []):
    gm = f["geometry"]; cs = gm["coordinates"]
    rings = [cs] if gm["type"] == "LineString" else cs if gm["type"] in ("Polygon", "MultiLineString") else [r for p_ in cs for r in p_]
    for ring in rings:
        a_ = np.array(ring, float)
        if len(a_) < 6: continue
        a_ = approximate_polygon(a_, tolerance=0.15)
        if len(a_) >= 3: today.append([[round(x, 2), round(y, 2)] for x, y in a_])
meta["globe"] = {"nx": NX, "ny": NY, "h": {"step": 0.5, "off": -11000}, "w": {"tstep": 4, "eta_max": ETA_MAX},
                 "frames": {"files": files, "n": len(frames), "dt": DF, "t0": T0, "w": TW, "h": TH, "per_row": PR},
                 "impact": {"lon": LON0, "lat": LAT0}, "labels": labels, "a0": A0, "today": today}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
fin = arr < 1e9
print("reached", round(100 * fin[~wall].mean()), "% of ocean cells;", "median arrival", round(float(np.median(arr[fin])) / 3600, 1), "h;", labels)
