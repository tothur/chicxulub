"""The 66 Ma topography for drawing the globe: the PALEOMAP grid at its full 0.1 deg, with its 40 m steps smoothed out.

PALEOMAP stores every height (land and sea floor) in steps of 40 m, so a gentle lowland is a staircase of flat
terraces, and relief shading lights only the step edges: from low orbit the land looks covered in contour lines. The
true height of a cell lies somewhere within its step (+-20 m of the stored value, after globe.py's -20 m shift), so
the steps are removed by smoothing that never leaves that range: repeated Gaussian smoothing at falling scales, each
pass clipped back into every cell's own range (a projection onto the heights the data allow). Terraces become slopes,
nothing moves by more than 20 m, and land stays land and sea stays sea (cells keep their sign). Mountains, where
neighbours differ by far more than a step, are unchanged to within 20 m.

The tsunami runs (globe.py, swe.py) keep their own inputs; this is for drawing only (and for landcover.py's colours).

Writes sim/data/globe_t.png (3600x1800 elevation for drawing, (R*256+G)*0.5 - 11000 m; cell centres at 0.1 deg, the
source having nodes on whole tenths) and meta.globe.topo. Run after globe.py, which rewrites meta.globe.

usage: python topo.py <paleodem 65 Ma .nc> <sim data dir>
"""
import json, math, os, sys
import numpy as np
import netCDF4
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates

NC, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
RE, STEP = 6371e3, 40.0

ds = netCDF4.Dataset(NC)
var = lambda pre: next(k for k in ds.variables if k.lower().startswith(pre))      # lon/lat or longitude/latitude
L, B, Z = ds[var("lon")][:], ds[var("lat")][:], np.asarray(ds["z"][:], dtype=np.float64) - 20.0  # as in globe.py: sea level midway between 0 and 40 m
# (its -180 and 180 columns are not always equal, so the grid is not wrapped: smoothing and sampling stop at its edges,
# and the GPU's wrap draws the seam as before)


def dequant(z):
    """Smooth the 40 m steps away, keeping every cell within its own step and on its own side of sea level."""
    lo, hi = z - STEP / 2, z + STEP / 2
    lo = np.where(z > 0, np.maximum(lo, 0.5), lo); hi = np.where(z > 0, hi, np.minimum(hi, -0.5))
    x = z.copy()
    for s in (16, 8, 4, 2, 1):                                          # source cells (0.1 deg)
        for _ in range(6):
            x = np.clip(gaussian_filter(x, s, mode="nearest"), lo, hi)
    return x


Zd = dequant(Z)
print("dequantised: max change", round(float(np.abs(Zd - Z).max()), 1), "m; land cells", int((Z > 0).sum()), "->", int((Zd > 0).sum()))
TX, TY = 3600, 1800
tlon = -180 + (np.arange(TX) + 0.5) * 360 / TX; tlat = 90 - (np.arange(TY) + 0.5) * 180 / TY
ht = np.empty((TY, TX), np.float32)
for r0 in range(0, TY, 300):
    TLo, TLa = np.meshgrid(tlon, tlat[r0:r0 + 300])
    ht[r0:r0 + 300] = map_coordinates(Zd, [(TLa - B[0]) / (B[1] - B[0]), (TLo - L[0]) / (L[1] - L[0])], order=1, mode="nearest")
fr = meta["gulf"]["frame"]; LON0, LAT0 = fr["paleo_lon"], fr["paleo_lat"]
TLo, TLa = np.meshgrid(tlon, tlat)
p1, p2 = math.radians(LAT0), np.radians(TLa)
R = RE * np.arccos(np.clip(math.sin(p1) * np.sin(p2) + math.cos(p1) * np.cos(p2) * np.cos(np.radians(TLo - LON0)), -1, 1))
ht[(R < 450e3) & (ht > 0) & (ht < 400)] = -40.0                         # the Yucatán shelf, as in globe.py
q = np.round((np.clip(ht, -11000, 9000) + 11000) / 0.5).astype(np.uint32)
Image.fromarray(np.dstack([q >> 8, q & 255, np.zeros_like(q)]).astype(np.uint8)).save(os.path.join(OUT, "globe_t.png"), optimize=True)
meta["globe"]["topo"] = {"file": "globe_t.png", "nx": TX, "ny": TY, "dequantised": "40 m steps smoothed within +-20 m"}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
