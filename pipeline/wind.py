"""Prevailing winds on the day of the impact, for the viewer's fires: which way they spread and where their smoke goes.

The winds of a climate model run on the same 66 Ma geography: HadCM3L, the Bristol BRIDGE group's Phanerozoic series
on Scotese's PALEOMAP maps (Valdes, Scotese & Lunt 2021, Clim. Past; run texpn1, 66 Ma, Foster et al. 2017 CO2;
www.paleo.bristol.ac.uk/ummodel/data/texpn1/climate/). Its land and sea match globe_h.png in 98% of its 3.75 x 2.5 deg
cells, with no shift. The model's April and May climate means (the last decades of a multi-millennial run) are blended
for 5 May, the date the viewer assumes (solar longitude 45 deg; During et al. 2022): 1/3 April, 2/3 May.

  surface   the 10 m wind, which drives the spread of a surface fire
  smoke     the mean of the 925, 850 and 700 hPa winds (about 0.8 to 3 km up), where most wildfire smoke travels

These are monthly means, the prevailing winds: steady in the trade-wind belts, where the mean is close to any day's
wind, but in middle latitudes the weather of the day (fronts, cyclones) blows from almost any side, and the mean is the
average of it. The impact itself (the blast, the heating of the upper air) is not in them.

Writes sim/data/wind.png (360 x 360 RGB, the 1 deg grid of the fire codes: rows 0-179 the surface wind, rows 180-359
the smoke layer; R east, G north, 128 + 127 * speed / scale, scales 15 and 25 m/s) and meta.wind.

usage: python wind.py <hadcm3 climate dir> <sim data dir>
"""
import json, os, sys
import numpy as np
import netCDF4
from PIL import Image
from scipy.ndimage import map_coordinates

RAW, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
MONTHS = {"apr": 1 / 3, "may": 2 / 3}
S10, SSM = 15.0, 25.0


def field(kind, name, lev=None):
    out = 0
    for m, w in MONTHS.items():
        d = netCDF4.Dataset(os.path.join(RAW, f"texpn1a.{kind}cl{m}.nc"))
        a = np.ma.filled(d[name][0], np.nan).astype(np.float64)
        out = out + w * (a[0] if lev is None else a[lev].mean(0))
        lat, lon = d["latitude_1"][:].astype(np.float64), d["longitude_1"][:].astype(np.float64)   # the wind (u) grid
    bad = np.isnan(out)                                             # a stray missing cell or two: the mean of its neighbours
    if bad.any():
        k = np.nan_to_num(out); n = (~bad).astype(float)
        sm = lambda a: sum(np.roll(np.roll(a, i, 0), j, 1) for i in (-1, 0, 1) for j in (-1, 0, 1))
        out = np.where(bad, sm(k) / np.maximum(sm(n), 1), out)
    print(name, "missing cells filled:", int(bad.sum()))
    return out, lat, lon


def regrid(a, lat, lon):
    """Cubic resampling to 1 deg cell centres (89.5 N first, -179.5 E first), wrapping in longitude."""
    P = 3; ap = np.concatenate([a[:, -P:], a, a[:, :P]], axis=1)
    la = 89.5 - np.arange(180); lo = -179.5 + np.arange(360)
    ri = (la - lat[0]) / (lat[1] - lat[0])                          # lat runs north to south
    ci = ((lo - lon[0]) % 360) / (lon[1] - lon[0]) + P
    RI, CI = np.meshgrid(ri, ci, indexing="ij")
    return map_coordinates(ap, [RI, CI], order=3, mode="nearest")


d = netCDF4.Dataset(os.path.join(RAW, "texpn1a.pcclmay.nc")); P = list(d["p"][:].astype(int))
LEV = [P.index(925), P.index(850), P.index(700)]
u10, lat, lon = field("pd", "u_mm_10m"); v10, _, _ = field("pd", "v_mm_10m")
us, _, _ = field("pc", "u_mm_p", LEV); vs, _, _ = field("pc", "v_mm_p", LEV)
enc = lambda a, s: np.clip(np.round(128 + 127 * a / s), 1, 255).astype(np.uint8)
img = np.zeros((360, 360, 3), np.uint8)
for r0, (u, v, s) in enumerate([(u10, v10, S10), (us, vs, SSM)]):
    U, V = regrid(u, lat, lon), regrid(v, lat, lon)
    img[r0 * 180:(r0 + 1) * 180, :, 0] = enc(U, s); img[r0 * 180:(r0 + 1) * 180, :, 1] = enc(V, s)
    sp = np.hypot(U, V); band = lambda a, b: round(float(sp[int(90 - b):int(90 - a)].mean()), 1)
    print(["surface", "smoke  "][r0], "max", round(float(sp.max()), 1), "m/s; mean speed 0-30 deg N", band(0, 30), " 30-60 N", band(30, 60),
          " 30-0 S", band(-30, 0), " 60-30 S", band(-60, -30), "; zonal mean u at 15 N", round(float(U[74].mean()), 1), "45 N", round(float(U[44].mean()), 1))
Image.fromarray(img).save(os.path.join(OUT, "wind.png"), optimize=True)
meta["wind"] = {"file": "wind.png", "nx": 360, "ny": 180, "scale": [S10, SSM], "date": "5 May (1/3 April, 2/3 May means)",
                "source": "HadCM3L texpn1 (66 Ma, Valdes et al. 2021)", "rows": "0-179 10 m wind, 180-359 925-700 hPa mean; R east, G north"}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
