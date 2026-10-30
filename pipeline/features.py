"""Map features for the Chicxulub regional map, projected to the DEM's grid (km).

Fits a circle to the mapped cenotes (OpenStreetMap natural=sinkhole) to locate the
cenote ring from data, and projects places, drill sites and crater geometry.

usage: python features.py <osm sinkholes.json> <data dir>
"""
import json, math, os, sys
import numpy as np
from pyproj import Proj
from scipy.optimize import least_squares

OSM, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
pp = meta["proj"]
proj = Proj(proj="aeqd", lat_0=pp["lat0"], lon_0=pp["lon0"], ellps=pp["ellps"])


def km(lat, lon):
    x, y = proj(lon, lat)
    return [round(x / 1000, 2), round(y / 1000, 2)]


# cenotes
pts = []
for e in json.load(open(OSM))["elements"]:
    la, lo = (e["lat"], e["lon"]) if "lat" in e else (e["center"]["lat"], e["center"]["lon"])
    pts.append(km(la, lo))
P = np.array(pts)

# robust circle fit, tightening the band around the current circle each pass
c, r = np.zeros(2), 85.0
for it in range(8):
    d = np.hypot(*(P - c).T)
    sel = np.abs(d - r) < (15 if it < 3 else 6)
    q = least_squares(lambda q: np.hypot(*(P[sel] - q[:2]).T) - q[2], [*c, r],
                      loss="soft_l1", f_scale=2).x
    c, r = q[:2], q[2]
resid = np.hypot(*(P[sel] - c).T) - r
ang = np.degrees(np.arctan2(P[sel, 1] - c[1], P[sel, 0] - c[0]))
lo, la = proj(c[0] * 1000, c[1] * 1000, inverse=True)

# structural centre: gravity-defined, as commonly cited (21°24'N 89°31'W)
centre = km(21.40, -89.5167)

places = [  # name, lat, lon, kind
    ("Mérida", 20.967, -89.623, "city"), ("Progreso", 21.283, -89.663, "town"),
    ("Chicxulub Puerto", 21.293, -89.607, "town"), ("Celestún", 20.859, -90.400, "town"),
    ("Sisal", 21.165, -90.030, "town"), ("Dzilam de Bravo", 21.393, -88.892, "town"),
    ("Campeche", 19.845, -90.524, "city"), ("Cancún", 21.161, -86.851, "city"),
    ("Valladolid", 20.690, -88.202, "town"), ("Chetumal", 18.500, -88.300, "city"),
    ("Ciudad del Carmen", 18.650, -91.810, "town"), ("Cozumel", 20.420, -86.920, "town"),
    ("Tulum", 20.211, -87.465, "town"), ("Chichén Itzá", 20.684, -88.568, "site"),
]
sites = [
    {"id": "M0077A", "name": "M0077A", "lat": 21.45, "lon": -89.95,
     "what": "IODP–ICDP Expedition 364 (2016). Cored the peak ring from 505.7 to 1334.7 m below the seafloor; shocked granite from the deep crust sits under impact breccia and melt."},
    {"id": "Yax-1", "name": "Yaxcopoil-1", "lat": 20.7440, "lon": -89.7185,
     "what": "ICDP borehole (2001–02), 1511 m deep, between the peak ring and the crater rim. About 100 m of impactites under ~800 m of post-impact limestone."},
]
for s in sites:
    s["xy"] = km(s["lat"], s["lon"])
    s["r_km"] = round(float(np.hypot(*(np.array(s["xy"]) - centre))), 1)

feat = {
    "centre": {"xy": centre, "lat": 21.40, "lon": -89.5167},
    "cenote_ring": {"xy": [round(float(c[0]), 2), round(float(c[1]), 2)], "r_km": round(float(r), 1),
                    "lat": round(la, 4), "lon": round(lo, 4), "n": int(sel.sum()),
                    "sd_km": round(float(resid.std()), 2),
                    "arc_deg": [round(float(ang.min()), 1), round(float(ang.max()), 1)],
                    "offset_km": round(float(np.hypot(*(c - centre))), 1)},
    "cenotes": [[round(x, 1), round(y, 1)] for x, y in pts],
    "places": [{"name": n, "xy": km(a, b), "kind": k} for n, a, b, k in places],
    "sites": sites,
}
json.dump(feat, open(os.path.join(OUT, "features.json"), "w"), separators=(",", ":"))
print({k: v for k, v in feat.items() if k not in ("cenotes", "places")})
