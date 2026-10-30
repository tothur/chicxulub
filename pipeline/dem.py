"""Regional elevation grid for the Chicxulub map.

Merges Copernicus GLO-90 (land) with GMRT (seafloor), resamples both onto an
azimuthal-equidistant grid centred on the crater, so rings are true circles and
distances from the centre are true. Writes data/dem.png (16-bit elevation packed
into R/G, land flag in B) and data/meta.json.

usage: python dem.py <raw dir> <out data dir>
  raw dir holds glo90/*.tif and gmrt/gmrt.tif
"""
import glob, json, os, sys
import numpy as np
import rasterio
from rasterio.merge import merge
from scipy.ndimage import map_coordinates, uniform_filter
from pyproj import Proj
from PIL import Image

RAW, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

LAT0, LON0 = 21.40, -89.55          # projection centre (near the crater centre)
X0, X1, Y0, Y1 = -420.0, 380.0, -360.0, 380.0   # km
RES = 0.4                           # km per pixel
STEP = 0.1                          # m per encoded unit

proj = Proj(proj="aeqd", lat_0=LAT0, lon_0=LON0, ellps="WGS84")
nx, ny = int(round((X1 - X0) / RES)), int(round((Y1 - Y0) / RES))
xs = (X0 + (np.arange(nx) + 0.5) * RES) * 1000
ys = (Y1 - (np.arange(ny) + 0.5) * RES) * 1000
X, Y = np.meshgrid(xs, ys)
lon, lat = proj(X, Y, inverse=True)


def sample(arr, transform, lon, lat):
    inv = ~transform
    c, r = inv * (lon, lat)
    return map_coordinates(arr, [r - 0.5, c - 0.5], order=1, mode="nearest")


# land: GLO-90 mosaic, box-averaged to ~the output pixel size before sampling
srcs = [rasterio.open(p) for p in sorted(glob.glob(os.path.join(RAW, "glo90", "*.tif")))]
glo, gt = merge(srcs, nodata=-9999)
glo = glo[0].astype(np.float32)
have = (glo > -9000).astype(np.float32)
glo[glo < -9000] = 0
k = max(1, int(round(RES * 1000 / 90)))
glo_s = uniform_filter(glo, k)
have_s = sample(uniform_filter(have, k), gt, lon, lat)
land_h = sample(glo_s, gt, lon, lat)

# sea: GMRT
g = rasterio.open(os.path.join(RAW, "gmrt", "gmrt.tif"))
gm = g.read(1).astype(np.float32)
gm[np.isnan(gm)] = np.nanmean(gm)
sea_h = sample(gm, g.transform, lon, lat)

covered = have_s > 0.99
land = np.where(covered, land_h > 0.25, sea_h > 0)
h = np.where(land, np.where(covered, land_h, sea_h), np.minimum(sea_h, -0.5))

off = float(np.floor(h.min()))
v = np.round((h - off) / STEP)
assert v.max() < 65536, f"range too large: {h.min()}..{h.max()}"
v = v.astype(np.uint32)
rgb = np.dstack([(v >> 8) & 255, v & 255, land * 255]).astype(np.uint8)
Image.fromarray(rgb, "RGB").save(os.path.join(OUT, "dem.png"), optimize=True)

meta = {
    "proj": {"type": "aeqd", "lat0": LAT0, "lon0": LON0, "ellps": "WGS84"},
    "extent_km": [X0, X1, Y0, Y1], "res_km": RES, "size": [nx, ny],
    "encode": {"offset_m": off, "step_m": STEP},
    "range_m": [float(h.min()), float(h.max())],
    "sources": ["Copernicus GLO-90 DEM (land)", "GMRT v4 (Ryan et al. 2009), seafloor"],
}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
print(meta)
