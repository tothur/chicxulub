"""The night sky behind the globe: today's real stars and Milky Way.

The Milky Way is NASA's Deep Star Maps 2020 background layer (NASA/Goddard Space Flight Center Scientific
Visualization Studio; Gaia DR2: ESA/Gaia/DPAC; public domain), from Gaia DR2 star counts with the bright stars taken
out: plate carree in J2000 right ascension and declination, RA increasing to the left with 0h at the centre
(checked against the galactic centre and the Large Magellanic Cloud). The stars are the Yale Bright Star Catalogue,
5th edition (Hoffleit and Warren 1991, CDS V/50): the 9,000-odd stars to magnitude 6.5, with J2000 positions,
V magnitude and B-V colour, drawn as points in the viewer.

This is today's sky. At 66 Ma none of it stood where it does now: over 66 million years every naked-eye star has
moved hundreds of light-years, and many of today's brightest (Orion's, for one) had not yet formed. Only the band
of the Milky Way would have looked broadly the same.

Writes sim/data/sky_mw.webp (4096x2048 sRGB, lossy: a smooth background), sky_stars.json and meta.sky.

usage: python sky.py <raw/sky dir> <sim data dir>
"""
import gzip, json, os, sys
import numpy as np
import OpenEXR
from PIL import Image

RAW, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))

# ---- the Milky Way: linear light (the map peaks at 1) -> sRGB bytes, a little brighter so the faint glow keeps steps
p = np.asarray(OpenEXR.File(os.path.join(RAW, "milkyway_2020_4k.exr")).channels()["RGB"].pixels, np.float32)
GAIN = 2.0
v = np.clip(p * GAIN, 0, 1)
srgb = np.where(v <= 0.0031308, 12.92 * v, 1.055 * np.power(v, 1 / 2.4) - 0.055)
Image.fromarray((srgb * 255 + 0.5).astype(np.uint8)).save(os.path.join(OUT, "sky_mw.webp"), "WEBP", quality=88, method=6)
print("milky way:", p.shape[1], "x", p.shape[0], round(os.path.getsize(os.path.join(OUT, "sky_mw.webp")) / 1e6, 2), "MB")

# ---- the bright stars (fixed-width records; ReadMe in raw/sky)
stars = []
for line in gzip.open(os.path.join(RAW, "bsc5_catalog.gz"), "rt", encoding="latin-1"):
    try:
        if not line[75:77].strip() or not line[102:107].strip(): continue          # no position (novae, removed entries)
        ra = 15 * (int(line[75:77]) + int(line[77:79]) / 60 + float(line[79:83]) / 3600)
        de = (int(line[84:86]) + int(line[86:88]) / 60 + int(line[88:90]) / 3600) * (-1 if line[83] == "-" else 1)
        vm = float(line[102:107]); bv = float(line[109:114]) if line[109:114].strip() else 0.6
    except ValueError:
        continue
    if vm <= 6.5: stars.append([round(ra, 3), round(de, 3), round(vm, 2), round(bv, 2)])
stars.sort(key=lambda s: s[2])
json.dump({"cols": ["ra", "dec", "vmag", "bv"], "stars": stars}, open(os.path.join(OUT, "sky_stars.json"), "w"), separators=(",", ":"))
print("stars:", len(stars), "brightest", stars[0], " to V 6.5")
meta["sky"] = {"mw": {"file": "sky_mw.webp", "gain": GAIN, "uv": "u = 0.5 - RA/360, v = 0.5 - Dec/180"}, "stars": "sky_stars.json",
               "credit": "NASA/GSFC Scientific Visualization Studio (Deep Star Maps 2020; Gaia DR2: ESA/Gaia/DPAC); Yale Bright Star Catalogue"}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
