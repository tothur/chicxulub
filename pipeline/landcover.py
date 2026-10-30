"""Land cover for the 66 Ma Earth, in true colour: what a satellite would see on a cloud-free day.

Not a vegetation model: a rule from three things the map knows (latitude, distance from the sea, height), following
the reconstructions of latest Cretaceous vegetation (Upchurch, Otto-Bliesner & Scotese 1999; Chumakov et al. 1995):
a greenhouse world, mostly forested, with no ice sheets. Evergreen rainforest in the wet tropics; broad-leaved and
conifer forests at middle latitudes; deciduous forests reaching the poles; open woodland, scrub and deserts in the
subtropical belts, more so far inland. There were no grasslands yet (grass-dominated plains came some 40 Myr later),
so open country is drawn as scrub and bare ground, not grass. High ground goes to rock; snow only where it is cold
enough (rare at 66 Ma).

Moisture index P (0..1): latitude belts (wet tropics, dry subtropics, wet storm tracks, moderately wet polar
regions), raised near the sea (within ~1,000 km) and on high ground. Temperature: the climate model's Maastrichtian
annual mean (29 - 24 sin^2 lat degC) less 6.5 degC/km. Colours are sRGB, picked from true-colour satellite mosaics of
today's equivalent landscapes; seeded noise at a few scales breaks up the belts. The boundaries are illustrative:
their positions follow the rule, not mapped fossil sites.

Writes sim/data/globe_c.png (3600x1800, the grid of globe_t.png; sea cells carry the nearest land's colour so that
blending at the shore stays clean) and meta.globe.cover.

usage: python landcover.py <sim data dir>
"""
import json, math, os, sys
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, distance_transform_edt
from scipy.spatial import cKDTree

OUT = sys.argv[1]
meta = json.load(open(os.path.join(OUT, "meta.json")))
G = meta["globe"]; T = G["topo"]
a = np.asarray(Image.open(os.path.join(OUT, T["file"])).convert("RGB")).astype(np.float64)
h = (a[..., 0] * 256 + a[..., 1]) * G["h"]["step"] + G["h"]["off"]
NY, NX = h.shape
lat = 90 - (np.arange(NY) + 0.5) * 180 / NY; lon = -180 + (np.arange(NX) + 0.5) * 360 / NX
LO, LA = np.meshgrid(np.radians(lon), np.radians(lat))
land = h > 0
rng = np.random.default_rng(66)


def unit(lo, la):
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def noise(sig_deg, amp):
    """Seeded smooth noise, sigma in degrees of latitude (wraps in longitude), std = amp."""
    s = sig_deg * NY / 180
    n = gaussian_filter(rng.standard_normal((NY, NX)), (s, s), mode=("nearest", "wrap"))
    return amp * n / n.std()


# distance from the sea (km), on the sphere: nearest shore cell
shore = ~land & (gaussian_filter(land.astype(float), 1.0, mode=("nearest", "wrap")) > 0.01)
tree = cKDTree(unit(LO[shore], LA[shore]))
dist = np.zeros((NY, NX))
chord, _ = tree.query(unit(LO[land], LA[land]), workers=-1)
dist[land] = 2 * np.arcsin(np.minimum(chord / 2, 1)) * 6371.0
mar = np.exp(-dist / 1000.0)

# moisture: latitude belts, the sea's reach, high ground
phi = np.abs(lat)[:, None] * np.ones((1, NX))
belt = 0.9 * np.exp(-(phi / 14) ** 2) + 0.75 * np.exp(-((phi - 50) / 16) ** 2) + 0.45 * np.exp(-((phi - 80) / 14) ** 2)
P = 0.15 + belt * (0.4 + 0.6 * mar) + 0.3 * mar + 0.1 * np.clip((h - 500) / 1500, 0, 1)
P = np.clip(P + noise(4, 0.07) + noise(1, 0.04), 0, 1)
Tm = 29 - 24 * np.sin(np.radians(lat))[:, None] ** 2 - 6.5 * np.maximum(h, 0) / 1000 + noise(3, 1.0)


def ss(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1); return t * t * (3 - 2 * t)


C = lambda r, g, b: np.array([r, g, b])[None, None, :]
rain, warm, cool, polar = C(.09, .19, .07), C(.13, .22, .09), C(.11, .18, .10), C(.16, .19, .11)
wood, scrub, desert, red = C(.30, .30, .16), C(.50, .44, .29), C(.76, .64, .45), C(.71, .50, .32)
rock, snow, wet = C(.42, .39, .35), C(.92, .93, .95), C(.08, .15, .07)
t_ = Tm[..., None]
forest = np.where(t_ > 22, rain, np.where(t_ > 15, rain + (warm - rain) * ss(22, 15, t_), warm + (cool - warm) * ss(15, 5, t_)))
forest = forest + (polar - forest) * ss(4, -2, t_)
p_ = P[..., None]
dry = desert + (red - desert) * np.clip(0.5 + noise(6, 0.35), 0, 1)[..., None]       # some deserts redder
col = dry + (scrub - dry) * ss(0.18, 0.3, p_)
col = col + (wood - col) * ss(0.3, 0.45, p_)
col = col + (forest - col) * ss(0.45, 0.62, p_)
col = col + (wet - col) * (ss(40, 5, h) * ss(0.55, 0.75, P))[..., None]            # swampy coastal lowlands
col = col + (rock - col) * ss(1800, 3000, h)[..., None]
col = col + (snow - col) * (ss(-1, -6, Tm) * ss(1500, 2500, h))[..., None]
col = col * (1 + noise(2, 0.07) + noise(0.4, 0.05))[..., None]                       # patchiness
col = np.clip(col, 0, 1)

# sea cells take the nearest land colour (clean blending along the shore)
idx = distance_transform_edt(~land, return_distances=False, return_indices=True)
col = col[idx[0], idx[1]]
Image.fromarray((col * 255 + 0.5).astype(np.uint8)).save(os.path.join(OUT, "globe_c.png"), optimize=True)

w = np.cos(np.radians(lat))[:, None] * land
share = lambda m: round(float((m * w).sum() / w.sum()), 3)
print("land: forest", share(P > 0.55), " woodland/scrub", share((P > 0.25) & (P <= 0.55)), " desert", share(P <= 0.25),
      " rock", share(h > 2400), " snow", share((Tm < -3) & (h > 1500)))
G["cover"] = {"file": "globe_c.png", "nx": NX, "ny": NY}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
