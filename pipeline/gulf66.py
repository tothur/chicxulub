"""Gulf-scale sea and tsunami for the 3D viewer, on a 66 Ma paleogeography.

Geography: PALEOMAP PaleoDEM for 65 Ma (Scotese & Wright 2018, 0.1 deg), which is in 66 Ma-ish paleo-coordinates.
The grid is built in that paleo frame, centred on the impact site's paleo position and turned so that the block's
north (today's north at the crater) lines up: rotations come from the GPlates Web Service, PALEOMAP model, 66 Ma.
Inside the block's 172 km disc the detailed 66 Ma target from paleo.py replaces the PaleoDEM, blended out to 235 km.

Tsunami (kinematic, after Range et al. 2022): rim wave at 220 km, 10 min, ~1.5 km high; crest arrival by fast
marching at sqrt(g h); height ~1/r with Green's law shoaling, capped near shore (same as gulf.py).

Also writes today's coastlines and a few of today's cities moved to where their rocks were at 66 Ma.

usage: python gulf66.py <paleodem 65Ma .nc> <cache dir> <sim data dir>
"""
import json, math, os, sys, urllib.request, urllib.parse
import numpy as np
import netCDF4, skfmm
from PIL import Image
from pyproj import Proj, Geod
from scipy.ndimage import map_coordinates, gaussian_filter
from skimage.measure import find_contours, approximate_polygon

NC, CACHE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(CACHE, exist_ok=True)
meta = json.load(open(os.path.join(OUT, "meta.json")))
GWS = "https://gws.gplates.org/reconstruct/"
AGE, MODEL = 66, "PALEOMAP"


def gws(path, params, cache):
    f = os.path.join(CACHE, cache)
    if not os.path.exists(f):
        url = GWS + path + "?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=180) as r: open(f, "wb").write(r.read())
    return json.load(open(f))


def rotate_points(lonlat, tag):
    out = []
    for i in range(0, len(lonlat), 100):
        chunk = lonlat[i:i + 100]
        pts = ",".join(f"{lo:.4f},{la:.4f}" for lo, la in chunk)
        out += gws("reconstruct_points/", {"points": pts, "time": AGE, "model": MODEL}, f"pts_{tag}_{i}.json")["coordinates"]
    return out

# ---- the paleo frame: impact site and the local turn of today's north ----
clat, clon = meta["crater"]["lat"], meta["crater"]["lon"]
(pc, pn) = rotate_points([(clon, clat), (clon, clat + 1.0)], "frame")
geod = Geod(ellps="WGS84")
theta = math.radians(geod.inv(pc[0], pc[1], pn[0], pn[1])[0])          # bearing of today's north in the paleo frame
pproj = Proj(proj="aeqd", lat_0=pc[1], lon_0=pc[0], ellps="WGS84")
ct, st = math.cos(theta), math.sin(theta)
to_paleo_xy = lambda x, n: (x * ct + n * st, -x * st + n * ct)            # block frame (east, north) -> paleo AEQD
to_block_xy = lambda X, Y: (X * ct - Y * st, X * st + Y * ct)
print("impact site at 66 Ma:", pc, "today's north turned by", round(math.degrees(theta), 1), "deg")

# ---- grid ----
N, HALF = 1024, 1250e3
cell = 2 * HALF / N
u = -HALF + (np.arange(N) + 0.5) * cell
v = HALF - (np.arange(N) + 0.5) * cell
E, Nn = np.meshgrid(u, v)
X, Y = to_paleo_xy(E, Nn)
plon, plat = pproj(X, Y, inverse=True)

ds = netCDF4.Dataset(NC)
names = list(ds.variables)
zname = [k for k in names if ds.variables[k].ndim == 2][0]
lons, lats = ds.variables[[k for k in names if k.lower().startswith("lon")][0]][:], ds.variables[[k for k in names if k.lower().startswith("lat")][0]][:]
Z = np.asarray(ds.variables[zname][:], dtype=np.float32) - 20.0  # PALEOMAP heights come in 40 m steps, its 0 m cells being sea: sea level is put midway between 0 and 40 m
ci = (plon - lons[0]) / (lons[1] - lons[0]); ri = (plat - lats[0]) / (lats[1] - lats[0])
h = map_coordinates(Z, [ri, ci], order=1, mode="nearest")

R = np.hypot(E, Nn)
# The Yucatán platform around the block continues the block's ramp (paleo.py): ~40 m of water to the SSW, ~2 km to the
# NNE, depth contours curving around a point 450 km SSW; where PALEOMAP's sea is deeper, it is kept. PALEOMAP draws
# Yucatán as low land at this resolution; the platform drowns it. The platform fades into PALEOMAP between ~430 and
# 650 km, along an edge that wanders by ~100 km so it is not a circle.
from scipy.ndimage import gaussian_filter
def sstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1); return t * t * (3 - 2 * t)
_rw = np.random.default_rng(66)
_w = gaussian_filter(_rw.standard_normal((N, N)), 40e3 / cell); _w /= np.abs(_w).max()
B, RP = math.radians(22.5), 450e3
s_r = (np.hypot(E + RP * math.sin(B), Nn + RP * math.cos(B)) - RP) / 1e3 + 18 * _w     # km down the ramp, + toward NNE
ramp = 40 + 1960 * sstep(-60, 140, s_r)
plat = -np.maximum(ramp, np.where(h < 0, -h, 0.0))
wp = 1 - sstep(430e3, 650e3, R + 110e3 * _w)
h = plat * wp + h * (1 - wp)

# inside the block: the detailed 66 Ma target
a = np.asarray(Image.open(os.path.join(OUT, "paleo.png"))).astype(np.int32)
d66 = (a[..., 0] * 256 + a[..., 1]) * 0.1
bn, bh = meta["N"], meta["half"]
inb = (np.abs(E) < bh) & (np.abs(Nn) < bh)
col = (E + bh) / (2 * bh) * bn - 0.5; row = (bh - Nn) / (2 * bh) * bn - 0.5
blk = -map_coordinates(d66, [row, col], order=1, mode="nearest")             # beyond the block this repeats its edge
# the viewer draws the block's surface as a disc of 172 km (R_BLK): exactly the block inside 165 km, blended out to 235 km
w = np.clip(1 - (R - 165e3) / 70e3, 0, 1) ** 1.0
w = w * w * (3 - 2 * w)
h = blk * w + h * (1 - w)

land = h > 0
D = np.where(land, 0.0, -h)

R0, T0 = 220e3, 600.0
phi = np.ma.MaskedArray(R - R0, mask=land)
tt = skfmm.travel_time(phi, np.sqrt(9.81 * np.maximum(D, 10.0)), dx=cell)
T = np.where(R <= R0, 180.0 + np.maximum(R - 100e3, 0) / 290.0, T0 + np.ma.filled(tt, np.inf))
T = np.where(land | ~np.isfinite(T), np.inf, T)
A = 1500.0 * np.minimum(1.0, R0 / np.maximum(R, 1.0)) * np.clip((1000.0 / np.maximum(D, 1.0)) ** 0.25, 0.6, 3.0)
A = np.where(land, 0, np.minimum(A, 0.8 * D + 30.0))

def pack16(x, step, off):
    q = np.round((x - off) / step); assert q.min() >= 0 and q.max() < 65536, (q.min(), q.max()); return q.astype(np.uint32)

qh = pack16(np.clip(h, -8000, 8000), 0.25, -8000)
Image.fromarray(np.dstack([qh >> 8, qh & 255, land * 255]).astype(np.uint8)).save(os.path.join(OUT, "gulf_h.png"), optimize=True)
qt = np.where(np.isfinite(T), np.clip(np.round(T / 2), 0, 65534), 65535).astype(np.uint32)
qa = np.clip(np.round(A / 8), 0, 255).astype(np.uint32)
Image.fromarray(np.dstack([qt >> 8, qt & 255, qa]).astype(np.uint8)).save(os.path.join(OUT, "gulf_w.png"), optimize=True)

def grid_lines(field, level, tol):
    out = []
    for c in find_contours(field, level):
        if len(c) < 20: continue
        c = approximate_polygon(c, tolerance=tol)
        out.append([[round(-HALF + (q[1] + 0.5) * cell), round(-(HALF - (q[0] + 0.5) * cell))] for q in c])
    return out

coast66 = grid_lines(gaussian_filter(h, 0.8), 0.0, 0.8)

# today's coastlines, moved to where their rocks were at 66 Ma
gj = gws("coastlines/", {"time": AGE, "model": MODEL}, "coastlines.json")
today = []
def add_ring(ring):
    pts = []
    for lo, la in ring:
        X1, Y1 = pproj(lo, la); x, n = to_block_xy(X1, Y1)
        if abs(x) < HALF * 1.02 and abs(n) < HALF * 1.02: pts.append([round(x), round(-n)])
        elif len(pts) > 1: today.append(pts); pts = []
        else: pts = []
    if len(pts) > 1: today.append(pts)
for f in gj.get("features", []):
    g = f["geometry"]; cs = g["coordinates"]
    rings = [cs] if g["type"] == "LineString" else cs if g["type"] in ("Polygon", "MultiLineString") else [r for p in cs for r in p]
    for ring in rings: add_ring(ring)

def block_of(lo, la):
    X1, Y1 = pproj(lo, la); x, n = to_block_xy(X1, Y1); return x, -n

def arrival(x, z, rad=120e3):
    k = (np.hypot(E - x, Nn + z) < rad) & np.isfinite(T)
    if not k.any(): return None, None
    i = np.argmin(np.where(k, T, np.inf)); sh = k & (D < 150)
    return float(T.flat[i]), float(np.max(A[sh])) if sh.any() else float(A.flat[i])

towns = [("Houston", 29.76, -95.37), ("New Orleans", 29.95, -90.07), ("Tampa", 27.95, -82.46), ("Havana", 23.11, -82.37),
         ("Veracruz", 19.17, -96.13), ("Monterrey", 25.67, -100.31), ("Dallas", 32.78, -96.80), ("Memphis", 35.15, -90.05)]
rot = rotate_points([(lo, la) for _, la, lo in towns], "places")
labels = []
for (name, la, lo), (plo, pla) in zip(towns, rot):
    x, z = block_of(plo, pla)
    if not (np.isfinite(x) and np.isfinite(z)) or max(abs(x), abs(z)) > HALF * 0.97: continue
    t, amp = arrival(x, z)
    labels.append({"text": name, "x": round(x), "z": round(z), "kind": "town", "arrive_s": None if t is None else round(t), "height_m": None if amp is None else round(amp)})
# the Gulf's oceanic crust is not on a plate that GWS rotates; label the middle of the deep basin instead
deep = (h < -2500) & (Nn > 150e3) & (R < 1000e3)
labels.append({"text": "Gulf of Mexico", "x": round(float(E[deep].mean())), "z": round(float(-Nn[deep].mean())), "kind": "water"})

meta["gulf"] = {"N": N, "half": HALF, "cell": cell, "h": {"step": 0.25, "off": -8000}, "w": {"tstep": 2, "astep": 8},
                "R0": R0, "T0": T0, "coast": coast66, "today": today, "labels": labels,
                "frame": {"paleo_lon": pc[0], "paleo_lat": pc[1], "turn_deg": round(math.degrees(theta), 2), "age": AGE, "model": MODEL}}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
for L in labels:
    if L["kind"] == "town": print(L["text"], None if L["arrive_s"] is None else f'{L["arrive_s"]/3600:.2f} h', L["height_m"], "m")
print("today's coast lines:", len(today), " 66 Ma coast lines:", len(coast66))
