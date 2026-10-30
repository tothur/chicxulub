"""Gulf-scale sea and tsunami for the 3D viewer.

Grid: 2500 km square, azimuthal equidistant around the crater (same frame as the block), 1024^2 cells.
Depth: GMRT (today) as a stand-in for the 66 Ma Gulf; inside the block the 66 Ma model from paleo.py;
today's Yucatán peninsula, which was a shallow shelf at 66 Ma, is turned into 40 m of water.

Tsunami (kinematic, after Range et al. 2022): a rim wave pushed out by the ejecta curtain reaches 220 km at
10 min about 1.5 km high; from there the crest travels at the shallow-water speed sqrt(g h) (arrival times by
fast marching), its height falling ~1/r and rising over shallow water (Green's law, h^-1/4), capped near shore.

Writes sim/data/gulf_h.png (elevation: (R*256+G)*0.25 - 8000 m, B = land), sim/data/gulf_w.png
(arrival: (R*256+G)*2 s, 65535 = never; B = crest height / 8 m) and adds "gulf" to sim/data/meta.json.

usage: python gulf.py <gmrt_gulf.tif> <sim data dir>
"""
import json, math, os, sys
import numpy as np
import rasterio, skfmm
from PIL import Image
from pyproj import Proj
from scipy.ndimage import map_coordinates, gaussian_filter, binary_dilation
from skimage.measure import find_contours, approximate_polygon

TIF, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
m = meta["map"]
proj = Proj(proj="aeqd", lat_0=m["lat0"], lon_0=m["lon0"], ellps="WGS84")
CX, CY = m["cx_km"] * 1000, m["cy_km"] * 1000              # crater centre in map metres

N, HALF = 1024, 1250e3
cell = 2 * HALF / N
u = -HALF + (np.arange(N) + 0.5) * cell                      # east (world x)
v = HALF - (np.arange(N) + 0.5) * cell                       # north; row 0 = north (world z = -north)
E, Nn = np.meshgrid(u, v)
lon, lat = proj(E + CX, Nn + CY, inverse=True)

src = rasterio.open(TIF)
g = src.read(1).astype(np.float32); g[np.isnan(g)] = 0
c, r = ~src.transform * (lon, lat)
today = map_coordinates(g, [r - 0.5, c - 0.5], order=1, mode="nearest")
h = today.copy()

# Yucatán peninsula: shallow sea at 66 Ma
pen = (today > 0) & (lon > -92.0) & (lon < -86.4) & (lat > 17.6) & (lat < 22.0)
h[pen] = -40.0
# inside the block: the 66 Ma seafloor of paleo.py
a = np.asarray(Image.open(os.path.join(OUT, "paleo.png"))).astype(np.int32)
d66 = (a[..., 0] * 256 + a[..., 1]) * 0.1
bn, bh = meta["N"], meta["half"]
inb = (np.abs(E) < bh) & (np.abs(Nn) < bh)
col = (E + bh) / (2 * bh) * bn - 0.5; row = (bh - Nn) / (2 * bh) * bn - 0.5
h[inb] = -map_coordinates(d66, [row[inb], col[inb]], order=1, mode="nearest")

land = h > 0
D = np.where(land, 0.0, -h)
R = np.hypot(E, Nn)

# arrival time of the crest
R0, T0 = 220e3, 600.0
phi = np.ma.MaskedArray(R - R0, mask=land)
speed = np.sqrt(9.81 * np.maximum(D, 10.0))
tt = skfmm.travel_time(phi, speed, dx=cell)
T = np.where(R <= R0, 180.0 + np.maximum(R - 100e3, 0) / 290.0, T0 + np.ma.filled(tt, np.inf))
T = np.where(land | ~np.isfinite(T), np.inf, T)

# crest height
A = 1500.0 * np.minimum(1.0, R0 / np.maximum(R, 1.0)) * np.clip((1000.0 / np.maximum(D, 1.0)) ** 0.25, 0.6, 3.0)
A = np.minimum(A, 0.8 * D + 30.0)
A = np.where(land, 0, A)

def pack16(x, step, off):
    q = np.round((x - off) / step); assert q.min() >= 0 and q.max() < 65536, (q.min(), q.max()); return q.astype(np.uint32)

qh = pack16(np.clip(h, -8000, 8000), 0.25, -8000)
Image.fromarray(np.dstack([qh >> 8, qh & 255, land * 255]).astype(np.uint8)).save(os.path.join(OUT, "gulf_h.png"), optimize=True)
qt = np.where(np.isfinite(T), np.clip(np.round(T / 2), 0, 65534), 65535).astype(np.uint32)
qa = np.clip(np.round(A / 8), 0, 255).astype(np.uint32)
Image.fromarray(np.dstack([qt >> 8, qt & 255, qa]).astype(np.uint8)).save(os.path.join(OUT, "gulf_w.png"), optimize=True)

# today's coast (for orientation), world metres (x east, z south)
coast = []
for cc in find_contours(gaussian_filter(today, 0.7), 0.0):
    if len(cc) < 20: continue
    cc = approximate_polygon(cc, tolerance=0.8)
    coast.append([[round(-HALF + (q[1] + 0.5) * cell), round(-(HALF - (q[0] + 0.5) * cell))] for q in cc])

def world(la, lo):
    x, y = proj(lo, la); return [x - CX, -(y - CY)]

def arrival(la, lo, rad=60e3):        # first crest at sea within rad of a coastal place, and its height
    x, z = world(la, lo); k = (np.hypot(E - x, Nn + z) < rad) & np.isfinite(T)
    if not k.any(): return None, None
    i = np.argmin(np.where(k, T, np.inf)); sh = k & (D < 150)              # height where the water is shallow, near shore
    return float(T.flat[i]), float(np.max(A[sh])) if sh.any() else float(A.flat[i])

places = [("Houston", 29.76, -95.37), ("New Orleans", 29.95, -90.07), ("Tampa", 27.95, -82.46), ("Havana", 23.11, -82.37),
          ("Veracruz", 19.17, -96.13), ("Tampico", 22.25, -97.86), ("Belize City", 17.50, -88.20)]
labels = []
for name, la, lo in places:
    x, z = world(la, lo); t, amp = arrival(la, lo)
    labels.append({"text": name, "x": round(x), "z": round(z), "kind": "town", "arrive_s": None if t is None else round(t), "height_m": None if amp is None else round(amp)})
for name, la, lo in [("Gulf of Mexico", 25.6, -90.6), ("Caribbean Sea", 17.2, -82.5), ("Florida Straits", 24.2, -81.4), ("Bay of Campeche", 19.9, -94.6)]:
    x, z = world(la, lo); labels.append({"text": name, "x": round(x), "z": round(z), "kind": "water"})

meta["gulf"] = {"N": N, "half": HALF, "cell": cell, "h": {"step": 0.25, "off": -8000}, "w": {"tstep": 2, "astep": 8},
                "R0": R0, "T0": T0, "coast": coast, "labels": labels}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
for L in labels:
    if "arrive_s" in L: print(L["text"], None if L["arrive_s"] is None else f'{L["arrive_s"]/3600:.2f} h', L["height_m"], "m")
fin = np.isfinite(T); print("max T in sea (h):", np.percentile(T[fin], [50, 90, 99]) / 3600)
