"""Fine coastlines for the viewer: where the shore runs, on grids finer than the ones it draws.

The Gulf grid (2.4 km cells) and the globe (0.25 deg) decide land and sea per cell, so close up their coasts are
polygons and staircases. This writes the signed distance to the shore on finer grids; the viewer uses its sign for
land and sea (the relief and the colours by height stay on the coarse grids):

  coast_gulf.png    4096 x 4096 over the Gulf grid (610 m): its heights (gulf_h.png), resampled smoothly
  coast_globe.png   5760 x 2880 (0.0625 deg): the 65 Ma PALEOMAP grid (0.1 deg) resampled smoothly, with the Yucatán
                    shelf as in globe.py

Bilinear resampling and a light Gaussian turn the cell steps into smooth curves (a cubic spline rings, and lifts the
many shallow seas a few metres deep next to land above sea level). The 66 Ma map knows nothing finer than its 0.1 deg
(about 11 km), so the small wiggles of the coast below that (bays, spits, lagoons) are illustrative: a few metres of fractal
relief (seeded, up to about 5 m in the Gulf, 7 m on the globe), scaled down with the ground's slope so the shoreline
moves by no more than ~1.5 km in the Gulf and ~6 km on the globe: less than a cell of the grids it refines.
A distance rather than a height, because on a steep coast the height jumps from deep sea to high land within one
cell, and the GPU's blend of such a near-binary field draws the shore as scallops along the cells. Near the shore
the distance is height / slope (sub-cell accurate), further out the distance to the nearest cell across the shore.
Encoding: grey value = 128 + 126 d / (6 cells), d > 0 on land, clipped (2..254).

usage: python coast.py <paleodem 65 Ma .nc> <sim data dir>
"""
import json, math, os, sys
import numpy as np
import netCDF4
from PIL import Image
from scipy.ndimage import zoom, gaussian_filter, map_coordinates, distance_transform_edt

NC, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
RE = 6371e3
rng = np.random.default_rng(66)


def fractal(shape, sigmas, amps):
    """Seeded fractal relief: smoothed noise at a few scales (in cells), each scaled to its amplitude (m, std)."""
    out = np.zeros(shape, np.float32)
    for s, a in zip(sigmas, amps):
        n = gaussian_filter(rng.standard_normal(shape).astype(np.float32), s, mode="wrap")
        out += a * n / n.std()
    return out


def wiggle(h, cell, dmax, sigmas, amps):
    """Fractal relief whose amplitude is capped so it shifts the zero contour by at most ~dmax (m) on this slope."""
    gy, gx = np.gradient(h, cell)
    cap = dmax * np.hypot(gx, gy)
    return fractal(h.shape, sigmas, amps) / max(1e-9, float(np.sqrt(np.sum(np.square(amps))))) * np.minimum(float(np.sqrt(np.sum(np.square(amps)))), cap)


def enc(h):
    """Signed distance to the zero contour of h, in cells (land positive), encoded in a byte."""
    gy, gx = np.gradient(h); g = np.hypot(gx, gy)
    near = h / np.maximum(g, 1e-6)
    land = h > 0
    far = np.where(land, distance_transform_edt(land) - 0.5, -(distance_transform_edt(~land) - 0.5))
    d = np.where(np.abs(near) < 1.5, near, far)
    return np.clip(np.round(128 + 126 * np.clip(d / 6.0, -1, 1)), 2, 254).astype(np.uint8)


# ---- the Gulf grid: 1024 -> 4096, the block's disc left alone ----
G = meta["gulf"]; N, HALF = G["N"], G["half"]
a = np.asarray(Image.open(os.path.join(OUT, "gulf_h.png")).convert("RGB")).astype(np.float64)
h = (a[..., 0] * 256 + a[..., 1]) * G["h"]["step"] + G["h"]["off"]
K = 4; NF = N * K; cell = 2 * HALF / NF
hf = gaussian_filter(zoom(h, K, order=1, grid_mode=True, mode="nearest"), 0.6 * K / 2, mode="nearest").astype(np.float32)
u = -HALF + (np.arange(NF) + 0.5) * cell
R = np.hypot(u[None, :], u[:, None])
keep = np.clip((R - 240e3) / 60e3, 0, 1)                                 # inside ~240 km the block's own surface is drawn
hf += wiggle(hf, cell, 1500.0, [3, 6, 12, 24], [1.5, 2.5, 3.0, 3.5]) * keep
Image.fromarray(enc(hf)).save(os.path.join(OUT, "coast_gulf.png"), optimize=True)
print("gulf coast:", NF, "cells of", round(cell), "m; land share", round(float((hf > 0).mean()), 3), "(coarse", round(float((h > 0).mean()), 3), ")")

# ---- the globe: 0.2 deg PALEOMAP -> 0.0625 deg ----
ds = netCDF4.Dataset(NC)
var = lambda pre: next(k for k in ds.variables if k.lower().startswith(pre))      # lon/lat or longitude/latitude
L, B, Z = ds[var("lon")][:], ds[var("lat")][:], np.asarray(ds["z"][:], dtype=np.float64) - 20.0  # PALEOMAP heights come in 40 m steps, its 0 m cells being sea: sea level is put midway between 0 and 40 m
Z = Z[:, :-1]; L = L[:-1]                                               # the last column repeats the first (-180 = 180)
P = 4; Zp = np.concatenate([Z[:, -P:], Z, Z[:, :P]], axis=1)            # wrap in longitude
NX, NY = 5760, 2880
lon = -180 + (np.arange(NX) + 0.5) * 360 / NX
lat = 90 - (np.arange(NY) + 0.5) * 180 / NY
dl, db = L[1] - L[0], B[1] - B[0]
ci = (lon - L[0]) / dl + P; ri = (lat - B[0]) / db
hg = np.empty((NY, NX), np.float32)
for r0 in range(0, NY, 360):                                            # in strips, to keep memory modest
    rr = ri[r0:r0 + 360]
    RI, CI = np.meshgrid(rr, ci, indexing="ij")
    hg[r0:r0 + 360] = map_coordinates(Zp, [RI, CI], order=1, mode="nearest")
hg = gaussian_filter(hg, (0.6 * abs(db) / (180 / NY),) * 2, mode=("nearest", "wrap"))         # about 0.6 source cells
fr = G["frame"]; LON0, LAT0 = fr["paleo_lon"], fr["paleo_lat"]
LONm, LATm = np.meshgrid(lon, lat)
p1, p2 = math.radians(LAT0), np.radians(LATm)
Rg = RE * np.arccos(np.clip(math.sin(p1) * np.sin(p2) + math.cos(p1) * np.cos(p2) * np.cos(np.radians(LONm - LON0)), -1, 1))
hg[(Rg < 450e3) & (hg > 0) & (hg < 400)] = -40.0                        # the Yucatán shelf, as in globe.py
hg += wiggle(hg, RE * math.radians(360 / NX), 6000.0, [1.5, 3, 6], [3.0, 4.0, 5.0])        # (cells near the poles are narrower; the cap is looser there)
Image.fromarray(enc(hg)).save(os.path.join(OUT, "coast_globe.png"), optimize=True)
print("globe coast:", NX, "x", NY, "; land share", round(float((hg > 0).mean()), 3))

meta["coast_fine"] = {"gulf": {"file": "coast_gulf.png", "n": NF}, "globe": {"file": "coast_globe.png", "nx": NX, "ny": NY},
                 "enc": "128 + 126 * signed distance to the shore / 6 cells (land > 128)"}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
