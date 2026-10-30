"""66 Ma target for the 3D Chicxulub viewer: paleo seafloor depth and Mesozoic sediment thickness.

Model (all on the regional map's azimuthal-equidistant grid, re-centred on the crater centre):
- A carbonate ramp deepening from tens of metres in the SSW to ~2 km in the NNE
  (Gulick et al. 2008; Collins et al. 2008), along bearing 22.5 deg; its depth contours curve gently around the platform.
- Where today's seafloor is deep (off the platform) the deeper of ramp and today's depth is kept.
- Mesozoic sediments ~3 km thick in the SW, ~4 km in the NE (Collins et al. 2008 end-members),
  on crystalline crust; Moho ~33 km (Collins et al. 2020).
Small fractal noise is added so the surfaces are not machine-smooth; it carries no data.

Writes sim/data/paleo.png (depth, 0.1 m), sim/data/sed.png (thickness, 0.1 m), sim/data/meta.json.

usage: python paleo.py <map data dir> <features.json> <out data dir>
"""
import json, math, os, sys
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates, gaussian_filter, zoom
from skimage.measure import find_contours, approximate_polygon

MAP, FEAT, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(66)

mm = json.load(open(os.path.join(MAP, "meta.json")))
feat = json.load(open(FEAT))
a = np.asarray(Image.open(os.path.join(MAP, "dem.png"))).astype(np.int32)
H = (a[..., 0] * 256 + a[..., 1]) * mm["encode"]["step_m"] + mm["encode"]["offset_m"]
X0, X1, Y0, Y1 = mm["extent_km"]; RES = mm["res_km"]

CX, CY = feat["centre"]["xy"]          # crater centre in map km
HALF = 180.0                           # block half-size, km
N = 1024
cell = 2 * HALF / N
u = -HALF + (np.arange(N) + 0.5) * cell            # east, km from centre
v = HALF - (np.arange(N) + 0.5) * cell             # north, km from centre (row 0 = north)
E, Nn = np.meshgrid(u, v)
col = (E + CX - X0) / RES - 0.5
row = (Y1 - (Nn + CY)) / RES - 0.5
today = map_coordinates(H, [row, col], order=1, mode="nearest")
D = np.maximum(0, -today)


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1)
    return t * t * (3 - 2 * t)


def fractal(octaves):
    out = np.zeros((N, N))
    for g, amp in octaves:
        r = rng.standard_normal((g, g))
        out += amp * zoom(gaussian_filter(r, 0.8), N / g, order=3)[:N, :N]
    return out / max(1e-9, np.abs(out).max())


BEAR = math.radians(22.5)
# distance down the ramp, + toward NNE. The depth contours curve around the platform (arcs about a point 450 km to the
# SSW) and wander by ~15 km, so the shelf margin is not ruler-straight; at the crater s is the plain NNE distance.
RP = 450.0
s = np.hypot(E + RP * math.sin(BEAR), Nn + RP * math.cos(BEAR)) - RP
_rw = np.random.default_rng(67)
_w = sum(a * zoom(gaussian_filter(_rw.standard_normal((g, g)), 0.8), N / g, order=3)[:N, :N] for g, a in [(4, 1.0), (9, 0.45)])
s = s + 15.0 * _w / np.abs(_w).max() * smoothstep(20, 80, np.hypot(E, Nn))
k = smoothstep(-60, 140, s)
ramp = 40 + 1960 * k
nz = fractal([(6, 1.0), (16, 0.5), (48, 0.25), (128, 0.12)])
ramp = ramp + nz * (12 + 0.06 * ramp)                  # a few metres on the inner platform, tens on the slope
off = smoothstep(200, 800, D)
d66 = np.maximum(8, ramp * (1 - off) + np.maximum(ramp, D) * off)

T = 3000 + 1000 * k + fractal([(5, 1.0), (14, 0.4)]) * 90

def save16(arr, name, step):
    q = np.round(arr / step)
    assert q.min() >= 0 and q.max() < 65536, (name, arr.min(), arr.max())
    q = q.astype(np.uint32)
    Image.fromarray(np.dstack([q >> 8, q & 255, np.zeros_like(q)]).astype(np.uint8)).save(os.path.join(OUT, name), optimize=True)

save16(d66, "paleo.png", 0.1)
save16(T, "sed.png", 0.1)

# today's coastline inside the block, as world-metre polylines (x east, z south)
coast = []
for c in find_contours(today, 0.0):
    if len(c) < 25: continue
    c = approximate_polygon(c, tolerance=0.6)
    xs = -HALF + (c[:, 1] + 0.5) * cell
    ys = HALF - (c[:, 0] + 0.5) * cell
    coast.append([[round(x * 1000), round(-y * 1000)] for x, y in zip(xs, ys)])

def world(xy):
    return [round((xy[0] - CX) * 1000), round(-(xy[1] - CY) * 1000)]

labels = []
SKIP = {"Chicxulub Puerto", "Sisal", "Dzilam de Bravo"}     # too close to the impact-point label
for p in feat["places"]:
    if p["name"] in SKIP: continue
    w = world(p["xy"])
    if max(abs(w[0]), abs(w[1])) < HALF * 1000 - 8000:
        labels.append({"text": p["name"], "x": w[0], "z": w[1], "kind": "town"})
for st in feat["sites"]:
    w = world(st["xy"]); labels.append({"text": st["name"], "x": w[0], "z": w[1], "kind": "site"})

# layers, youngest first: fraction of the Mesozoic pile, age span (Ma) over which each is laid down
LAYERS = [
    {"id": "marl", "name": "Marly limestone", "age": "Campanian–Maastrichtian", "frac": 0.10, "ma": [72, 66.04]},
    {"id": "lime", "name": "Limestone and dolomite", "age": "Late Cretaceous", "frac": 0.24, "ma": [100, 74]},
    {"id": "evap", "name": "Anhydrite and dolomite", "age": "Early Cretaceous", "frac": 0.38, "ma": [145, 102]},
    {"id": "red", "name": "Red beds and salt", "age": "Jurassic", "frac": 0.28, "ma": [195, 150]},
]
ic = N // 2
centre = {"water_m": round(float(d66[ic, ic])), "sed_m": round(float(T[ic, ic])),
          "layers_m": [round(float(T[ic, ic]) * L["frac"]) for L in LAYERS]}
meta = {
    "N": N, "half": HALF * 1000, "cell": cell * 1000,
    "map": {"lat0": mm["proj"]["lat0"], "lon0": mm["proj"]["lon0"], "cx_km": CX, "cy_km": CY},
    "crater": {"lat": feat["centre"]["lat"], "lon": feat["centre"]["lon"], "rim": [90000, 100000], "peak": [40000, 45000]},
    "ramp": {"bearing": 22.5, "s0_km": -60, "s1_km": 140, "shallow_m": 40, "deep_m": 2000},
    "sed": {"sw_m": 3000, "ne_m": 4000},
    "crust": {"conrad_m": 20000, "moho_m": 34000, "bottom_m": 40000},
    "layers": LAYERS, "centre": centre,
    "range": {"water": [float(d66.min()), float(d66.max())], "sed": [float(T.min()), float(T.max())]},
    "coast": coast, "labels": labels,
}
old = json.load(open(os.path.join(OUT, "meta.json"))) if os.path.exists(os.path.join(OUT, "meta.json")) else {}
old.update(meta)                          # keep what later steps added (gulf, globe, reentry, climate)
json.dump(old, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
print({k: meta[k] for k in ("centre", "range")}, len(coast), "coast lines", sum(len(c) for c in coast), "pts")
