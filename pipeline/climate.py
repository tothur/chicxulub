"""Impact winter and the long warm-up: darkness, cold and warmth over days to 100,000 years, for the viewer's globe.

Not a climate model of record: a small energy-balance model run with the aerosol loads of three published studies,
which disagree about which agent mattered most. Each scenario's aerosol optical depth over time is tuned so the
model gives that study's headline darkness and cooling (printed below, next to the published numbers):

  dust     Senel et al. 2023 (Nat. Geosci.): fine silicate dust with sulfur and soot; photosynthesis shut down for
           ~1.7-2 yr, up to ~15 degC global cooling, pre-impact temperatures back after ~15-20 yr.
  soot     Bardeen et al. 2017 (PNAS), Tabor et al. 2020 (GRL): soot from global fires; light below 1% of normal
           for ~2 yr, land up to -28 degC, ocean -11 degC, global -16.5 degC.
  sulfate  Brugger et al. 2017 (GRL): sulfate aerosols; global annual mean at least -26 degC, 3-16 yr below
           freezing, recovery longer than 30 yr.

Spreading: each agent starts where the impact put it (dust and sulfur with the distal ejecta, from reentry.py; soot
over land where re-entry heating was strongest) and is carried by stratospheric winds and eddy mixing on a 5 deg
grid for the first months, by which time it covers the globe.

Climate (2.5 deg): anomalies of an ocean mixed layer (150 m), a deep ocean (1,000 m, per latitude band), and land
(a thin surface plus boundary layer), with Budyko radiation (B = 1.8 W/m2/K), heat exchange between land, sea and air
in each latitude band, continentality (inland cells are less tied to the sea) and meridional diffusion. Light at the
surface = exp(-absorption) / (1 + 0.225 * scattering optical depth). Baseline Maastrichtian climate (annual means):
sea surface 29 degC at the equator to 5 degC at the poles, land 6.5 degC/km colder with height.

Seasons: the sunlight the aerosols take away is the day's own (daily-mean insolation for the latitude and date,
obliquity 23.4 deg), starting in mid boreal spring (solar longitude 45 deg, about 5 May; During et al. 2022). So the
darkness costs the north its first summer's sunshine, while the south goes into winter, when there is little sunlight
to lose. The model's anomalies are relative to the normal for that time of year. Frames are shown as they are for the
first ~2 years and blend into running annual means by 4 years (later frames are too far apart to sample the seasons).

After the aerosols: about 5 degC of warming over ~100,000 yr, found in fish-tooth oxygen isotopes at El Kef
(MacLeod et al. 2018). Its cause is debated (impact CO2 alone is too little); it is added as observed, not modelled.

Writes sim/data/climate_<id>.png (frames at 144x72, 20 per row: R light 1 + 254*(log10 + 6)/6, G sea and B land
anomaly 128 + 3*degC) and meta.climate (times, scenarios, global series for the chart).

usage: python climate.py <sim data dir>
"""
import json, math, os, sys
import numpy as np
from numba import njit
from PIL import Image
from scipy.ndimage import gaussian_filter, zoom

OUT = sys.argv[1]
meta = json.load(open(os.path.join(OUT, "meta.json")))
DAY, YR = 86400.0, 365.25 * 86400.0
RE = 6371e3

# ---------- inputs: 66 Ma land and height (0.25 deg), the distal ejecta layer and re-entry heating (1 deg) ----------
def px(name):
    return np.array(Image.open(os.path.join(OUT, name)).convert("RGB")).astype(np.int64)
gh = px("globe_h.png")
H = (gh[..., 0] * 256 + gh[..., 1]) * meta["globe"]["h"]["step"] + meta["globe"]["h"]["off"]
LAND = gh[..., 2] > 127
NX, NY = 144, 72                                       # 2.5 deg
fL = LAND.reshape(NY, 10, NX, 10).mean((1, 3))
hL = np.where(LAND, np.maximum(H, 0), 0).reshape(NY, 10, NX, 10).sum((1, 3)) / np.maximum(LAND.reshape(NY, 10, NX, 10).sum((1, 3)), 1)
lat = 90 - 2.5 * (np.arange(NY) + 0.5); phi = np.radians(lat)
wgt = np.cos(phi)[:, None] * np.ones((NY, NX)); wgt /= wgt.sum()
# continentality: how much sea lies within ~1,000 km
seaNear = 1 - gaussian_filter(fL, sigma=(4, 4), mode=("nearest", "wrap"))

rm = px("reentry_m.png")[..., 0].astype(float)
thick = np.where(rm > 0, 10 ** ((rm - 1) / 254 * 6 - 3), 0.0)          # mm, 1 deg
RX = meta["reentry"]; peak = np.zeros((180, 360))
for k, f in enumerate(RX["files"]):
    im = px(f); per = RX["per_row"]
    for c in range(per * per):
        for ch in range(3):
            n = k * per * per * 3 + c * 3 + ch
            if n >= RX["n"]: continue
            r0, c0 = divmod(c, per)
            q = im[r0 * 180:(r0 + 1) * 180, c0 * 360:(c0 + 1) * 360, ch].astype(float)
            peak = np.maximum(peak, np.where(q > 0, 10 ** ((q - 1) / 254 * 7), 0))

# ---------- spreading on a 5 deg grid ----------
TX, TY = 72, 36
def to5(a):  # 1 deg -> 5 deg block mean
    return a.reshape(TY, 5, TX, 5).mean((1, 3))
land1 = LAND.reshape(180, 4, 360, 4).mean((1, 3)) > 0.5
src = {
    "dust": to5(np.minimum(thick, 10.0)),                               # the distal layer; proximal debris never reaches the stratosphere
    "sulfate": to5(np.minimum(thick, 10.0)),                            # sulfur vapour travelled with the plume
    "soot": to5(land1 * np.clip((np.log10(np.maximum(peak, 1)) - 4.0) / 1.0, 0, 1) + 0.05 * land1),   # fires: land under strong heating
}
lat5 = np.radians(90 - 5 * (np.arange(TY) + 0.5))
u5 = -12 * np.cos(lat5) ** 2 + 28 * np.sin(2 * lat5) ** 2                # m/s: tropical easterlies, mid-latitude westerlies
KX, KY = 2.0e7, 5.0e6                                                    # eddy mixing, m2/s

@njit(cache=True)
def spread(c, u, lat, dt, nsteps, KX, KY):
    ny, nx = c.shape
    dy = math.pi / ny * 6371e3; dl = 2 * math.pi / nx
    out = c.copy(); nw = c.copy()
    for s in range(nsteps):
        for j in range(ny):
            cj = math.cos(lat[j]); dx = max(dl * 6371e3 * cj, 2.0e5)          # polar rows are mixed below; keeps diffusion stable
            for i in range(nx):
                ip = (i + 1) % nx; im = (i - 1) % nx
                a = out[j, i]
                adv = -u[j] * ((a - out[j, im]) if u[j] > 0 else (out[j, ip] - a)) / dx
                dfx = KX * (out[j, ip] - 2 * a + out[j, im]) / (dx * dx)
                jn = max(j - 1, 0); js = min(j + 1, ny - 1)
                cn = math.cos(0.5 * (lat[j] + lat[jn])) if j > 0 else 0.0
                cs = math.cos(0.5 * (lat[j] + lat[js])) if j < ny - 1 else 0.0
                dfy = KY * (cn * (out[jn, i] - a) + cs * (out[js, i] - a)) / (cj * dy * dy)
                nw[j, i] = a + dt * (adv + dfx + dfy)
        out[:, :] = nw
        # the polar rows mix around the pole
        out[0, :] = out[0, :].mean(); out[ny - 1, :] = out[ny - 1, :].mean()
    return out

# ---------- frame times: 6 h to 300 yr from the model, then the long warm-up to 100 kyr ----------
T_MODEL = np.geomspace(6 * 3600, 300 * YR, 170)
T_LONG = np.geomspace(330 * YR, 1e5 * YR, 30)
TIMES = np.concatenate([T_MODEL, T_LONG])

def patterns():
    """Normalised column pattern of each agent at each model frame time (stops changing after ~5 months)."""
    w5 = np.cos(lat5)[:, None] * np.ones((TY, TX)); w5 /= w5.sum()
    out = {}
    for k, c0 in src.items():
        c = c0 / (c0 * w5).sum(); t = 0.0; frames = []
        for T in T_MODEL:
            T = min(T, 150 * DAY)
            n = int(max(0, (T - t)) // 300)
            if n: c = spread(c, u5, lat5, 300.0, n, KX, KY); t += n * 300.0
            frames.append(zoom(c / (c * w5).sum(), 2, order=1, mode="grid-wrap"))
        out[k] = np.array(frames)
    return out

# ---------- the energy-balance model ----------
S0, OBL, LS0 = 1360.0, math.radians(23.4), math.radians(45.0)     # solar constant; obliquity; solar longitude at impact
A_SW = 0.7                                                          # share of the sunlight that the surface and air absorb

@njit(cache=True)
def insol(phi, t):
    """Daily-mean sunlight at the top of the atmosphere, W/m2, t seconds after an impact at solar longitude LS0."""
    ls = LS0 + 2 * math.pi * t / (365.25 * 86400.0)
    dec = math.asin(math.sin(OBL) * math.sin(ls))
    c = -math.tan(phi) * math.tan(dec); h0 = math.pi if c < -1 else (0.0 if c > 1 else math.acos(c))
    return S0 / math.pi * (h0 * math.sin(phi) * math.sin(dec) + math.cos(phi) * math.cos(dec) * math.sin(h0))
T0_SEA = 29 - 24 * np.sin(phi) ** 2
T0 = {"sea": T0_SEA[:, None] * np.ones((NY, NX)), "land": T0_SEA[:, None] - 6.5 * hL / 1000.0}
B, C_O, C_D, C_L = 1.8, 4.2e6 * 150, 4.2e6 * 2000, 2.0e7
GAMMA, GCONV, KA, KL, DM = 0.7, 6.0, 2.0, 7.0, 0.6     # sea-deep exchange (stable / convecting), air, land, meridional

@njit(cache=True)
def ebm(times, F, fL, seaNear, A_SW, B, C_O, C_D, C_L, GAMMA, GCONV, KA, KL, DM, outO, outL, outOy, outLy):
    """F[k] = lost sunlight fraction at times[k] (interpolated in between); the day's own sunlight. Anomalies, degC:
    as they are (outO, outL) and as running means over the past year (outOy, outLy)."""
    ny, nx = fL.shape
    To = np.zeros((ny, nx)); Tl = np.zeros((ny, nx)); Td = np.zeros(ny)
    dphi = math.pi / ny
    S = np.zeros(ny)
    NB = 365; ringO = np.zeros((NB, ny, nx)); ringL = np.zeros((NB, ny, nx)); sumO = np.zeros((ny, nx)); sumL = np.zeros((ny, nx))
    nday = 0; tday = 86400.0
    t = 0.0; k = 0; kf = 0
    while kf < len(times):
        for j in range(ny):
            S[j] = A_SW * insol(math.pi / 2 - (j + 0.5) * dphi, t)
        dt = 3600.0 if t < 10 * 86400.0 else (21600.0 if t < 2 * 365.25 * 86400.0 else 86400.0)
        while k + 1 < len(times) and times[k + 1] <= t: k += 1
        w = 0.0 if t <= times[0] else min(1.0, (t - times[k]) / (times[k + 1] - times[k])) if k + 1 < len(times) else 1.0
        # band means of sea, land and air
        mo = np.zeros(ny); ma = np.zeros(ny)
        for j in range(ny):
            so = 0.0; sw = 0.0; sa = 0.0
            for i in range(nx):
                fo = 1 - fL[j, i]; so += fo * To[j, i]; sw += fo; sa += fo * To[j, i] + fL[j, i] * Tl[j, i]
            mo[j] = so / sw if sw > 0 else 0.0; ma[j] = sa / nx
        for j in range(ny):
            pj = math.pi / 2 - (j + 0.5) * dphi; cj = math.cos(pj)
            jn = max(j - 1, 0); js = min(j + 1, ny - 1)
            cn = math.cos(pj + dphi / 2) if j > 0 else 0.0; cs = math.cos(pj - dphi / 2) if j < ny - 1 else 0.0
            dif = DM * (cn * (mo[jn] - mo[j]) + cs * (mo[js] - mo[j])) / (cj * dphi * dphi)
            for i in range(nx):
                lost = S[j] * ((1 - w) * F[k, j, i] + w * F[min(k + 1, len(times) - 1), j, i])
                if fL[j, i] < 0.999:
                    g = GCONV if To[j, i] < Td[j] else GAMMA            # a surface colder than the water below overturns
                    To[j, i] += dt / C_O * (-lost - B * To[j, i] + g * (Td[j] - To[j, i]) + KA * (ma[j] - To[j, i]) + dif)
                air = (1 - seaNear[j, i]) * ma[j] + seaNear[j, i] * (To[j, i] if fL[j, i] < 0.999 else mo[j])
                kl = KL * (0.35 + 0.65 * seaNear[j, i])
                Tl[j, i] += dt / C_L * (-lost - B * Tl[j, i] + kl * (air - Tl[j, i]))
            Td[j] += dt / C_D * (GCONV if mo[j] < Td[j] else GAMMA) * (mo[j] - Td[j])
        t += dt
        while t >= tday:                                                   # one sample a day for the running annual mean
            b = nday % NB
            sumO += To - ringO[b]; sumL += Tl - ringL[b]; ringO[b] = To; ringL[b] = Tl
            nday += 1; tday += 86400.0
        while kf < len(times) and times[kf] <= t:
            m = max(min(nday, NB), 1)
            outO[kf] = To; outL[kf] = Tl; outOy[kf] = sumO / m; outLy[kf] = sumL / m; kf += 1

# ---------- scenarios: optical depth of each agent over time (global mean), tuned to each study ----------
def plateau(t, tau0, tc, n):            # stays near tau0, then falls away from about tc
    return tau0 / (1 + (t / tc) ** n)
def rise(t, tr):
    return 1 - np.exp(-t / tr)
DUST = lambda t: (plateau(t, 1000, 0.45 * YR, 3) + plateau(t, 1.2, 4 * YR, 2)) * rise(t, 0.3 * DAY)    # coarse and fine dust
SCEN = [
    {"id": "dust", "name": "Fine dust", "ref": "Senel et al. 2023",
     "abs": lambda t: 0.14 * DUST(t),                                   # silicate dust absorbs about a seventh of what it intercepts
     "sca": lambda t: 0.86 * DUST(t) + plateau(t, 1.5, 2.5 * YR, 2) * rise(t, 20 * DAY),
     "mix": {"dust": 0.75, "sulfate": 0.24, "soot": 0.01},
     "paper": "photosynthesis off ~1.7-2 yr; global cooling up to ~15 degC; recovery ~15-20 yr"},
    {"id": "soot", "name": "Soot from global fires", "ref": "Bardeen et al. 2017; Tabor et al. 2020",
     "abs": lambda t: plateau(t, 160, 0.6 * YR, 3) * rise(t, 2 * DAY),
     "sca": lambda t: plateau(t, 3, 1.0 * YR, 2) * rise(t, 20 * DAY),
     "mix": {"soot": 0.9, "sulfate": 0.1, "dust": 0.0},
     "paper": "light <1% for ~2 yr; land up to -28 degC, ocean -11 degC; global -16.5 degC (Tabor)"},
    {"id": "sulfate", "name": "Sulfate haze", "ref": "Brugger et al. 2017",
     "abs": lambda t: 0 * t,
     "sca": lambda t: plateau(t, 20, 3.2 * YR, 2.0) * rise(t, 25 * DAY),
     "mix": {"sulfate": 1.0, "dust": 0.0, "soot": 0.0},
     "paper": "global annual mean at least -26 degC; 3-16 yr below freezing; recovery >30 yr"},
]
EL_KEF_LAT = 25.0
W_LONG = lambda t: 5.0 * np.clip(np.log10(np.maximum(t, 1) / (100 * YR)) / 3.0, 0, 1)   # +5 degC by 100 kyr, schematic

def run(sc, P):
    # local optical depth = global value x the spread pattern of the agents in the mix
    pat = sum(sc["mix"][k] * P[k] for k in P if sc["mix"].get(k, 0) > 0) / sum(v for v in sc["mix"].values())
    ta = sc["abs"](T_MODEL)[:, None, None] * pat; ts = sc["sca"](T_MODEL)[:, None, None] * pat
    light = np.exp(-ta) / (1 + 0.225 * ts)
    # of the light that does not get through, a share is reflected to space and the rest absorbed in the layer, whose
    # heat is re-radiated, about half of it down to the surface. Reflected share: single scattering for a thin layer,
    # the two-stream reflectance of a semi-infinite one for a thick layer (asymmetry g = 0.7)
    om = ts / np.maximum(ta + ts, 1e-9); om1 = om * 0.3 / (1 - 0.7 * om)
    rinf = (1 - np.sqrt(1 - om1)) / (1 + np.sqrt(1 - om1))
    r = om1 + (rinf - om1) * (1 - np.exp(-(ta + 0.3 * ts)))
    F = ((1 - light) * (r + 0.5 * (1 - r))).astype(np.float64)
    oO = np.zeros((len(T_MODEL), NY, NX)); oL = np.zeros_like(oO); oOy = np.zeros_like(oO); oLy = np.zeros_like(oO)
    ebm(T_MODEL, F, fL, seaNear, A_SW, B, C_O, C_D, C_L, GAMMA, GCONV, KA, KL, DM, oO, oL, oOy, oLy)
    # as they are for the first ~2 years, then running annual means (frames further apart than the seasons)
    wy = np.clip((T_MODEL / YR - 1.5) / 2.5, 0, 1); wy = (wy * wy * (3 - 2 * wy))[:, None, None]
    oO = oO * (1 - wy) + oOy * wy; oL = oL * (1 - wy) + oLy * wy
    # the long warm-up: what remains of the model anomaly fades; the observed warming, stronger toward the poles
    # stronger toward the poles, scaled so the El Kef sea (about 25 deg N at 66 Ma) sees the recorded 5 degC; land 1.2x
    amp = ((1 + 0.6 * np.sin(phi) ** 2) / (1 + 0.6 * math.sin(math.radians(EL_KEF_LAT)) ** 2))[:, None]
    fade = np.exp(-(T_LONG[:, None, None] - T_MODEL[-1]) / (300 * YR))
    fadeO, fadeL = oO[-1] * fade, oL[-1] * fade
    lo = W_LONG(T_LONG)[:, None, None] * amp
    sea = np.concatenate([oO + W_LONG(T_MODEL)[:, None, None] * amp, fadeO + lo])
    land = np.concatenate([oL + W_LONG(T_MODEL)[:, None, None] * amp * 1.2, fadeL + lo * 1.2])
    light = np.concatenate([light, np.ones((len(T_LONG), NY, NX))])
    return light, sea, land

def metrics(light, sea, land):
    land_w = wgt * fL; sea_w = wgt * (1 - fL)
    g = (sea * sea_w + land * land_w).sum((1, 2)); gl = (land * land_w).sum((1, 2)) / land_w.sum(); go = (sea * sea_w).sum((1, 2)) / sea_w.sum()
    lm = (light * wgt).sum((1, 2))
    dark = (((light < 0.01) * wgt).sum((1, 2)))                       # share of the globe below 1% light
    Tabs = ((T0["sea"] + sea) * sea_w + (T0["land"] + land) * land_w).sum((1, 2))
    return g, gl, go, lm, dark, Tabs

def span(t, cond):
    idx = np.where(cond)[0]
    return (t[idx[-1]] - t[idx[0]]) / YR if len(idx) else 0.0

if __name__ == "__main__":
    P = patterns()
    out = {"times": [float(x) for x in TIMES], "nx": NX, "ny": NY, "per_row": 20,
           "enc": {"light": "1 + 254 * (log10 + 6) / 6", "dT": "128 + 3 * degC (G sea, B land)"},
           "base": {"sea": "29 - 24 sin^2(lat)", "lapse": 6.5}, "scenarios": []}
    for sc in SCEN:
        light, sea, land = run(sc, P)
        g, gl, go, lm, dark, Tabs = metrics(light, sea, land)
        t = TIMES
        rec = t[np.where((t > t[np.argmin(g)]) & (g > -1.0))[0]]
        print(f"== {sc['id']}: {sc['paper']}")
        print(f"   model: global min {g.min():.1f} at {t[np.argmin(g)] / YR:.2f} yr; land min {gl.min():.1f}, sea min {go.min():.1f}; "
              f"half the globe <1% light for {span(t, dark > 0.5):.2f} yr, 90% for {span(t, dark > 0.9):.2f} yr; "
              f"global mean below 0 degC for {span(t, Tabs < 0):.1f} yr; back within 1 degC at {rec[0] / YR if len(rec) else -1:.0f} yr; "
              f"baseline global {Tabs[0]:.1f} degC")
        # frames -> atlas
        enc_l = np.where(light > 1e-6, np.clip(np.round(1 + 254 * (np.log10(np.maximum(light, 1e-6)) + 6) / 6), 1, 255), 0).astype(np.uint8)
        enc_o = np.clip(np.round(128 + 3 * sea), 0, 255).astype(np.uint8); enc_g = np.clip(np.round(128 + 3 * land), 0, 255).astype(np.uint8)
        PR = out["per_row"]; rows = math.ceil(len(t) / PR)
        img = np.zeros((rows * NY, PR * NX, 3), np.uint8)
        for k in range(len(t)):
            r, c = divmod(k, PR)
            img[r * NY:(r + 1) * NY, c * NX:(c + 1) * NX] = np.dstack([enc_l[k], enc_o[k], enc_g[k]])
        fn = f"climate_{sc['id']}.png"; Image.fromarray(img).save(os.path.join(OUT, fn), optimize=True)
        r2 = lambda a, n=2: [round(float(x), n) for x in a]
        out["scenarios"].append({"id": sc["id"], "name": sc["name"], "ref": sc["ref"], "paper": sc["paper"], "file": fn,
                                 "global": r2(g), "land": r2(gl), "sea": r2(go), "light": [float(f"{x:.3g}") for x in lm],
                                 "dark": r2(dark, 3), "abs": r2(Tabs)})
    meta["climate"] = out
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
