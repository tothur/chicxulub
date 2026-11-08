# Chicxulub regional map: data pipeline

Raw inputs (re-download into a scratch `raw/` dir; `osm_sinkholes.json` is kept in `chicxulub/raw/`):

- `glo90/*.tif`: Copernicus GLO-90 tiles N17–N23, W086–W094, from `https://copernicus-dem-90m.s3.amazonaws.com/Copernicus_DSM_COG_30_N{lat}_00_W{lon}_00_DEM/...tif`
- `gmrt/gmrt.tif`: GMRT GridServer, `minlongitude=-95&maxlongitude=-84&minlatitude=16&maxlatitude=27&format=geotiff&resolution=high&layer=topo`
- `osm_sinkholes.json`: Overpass, `natural=sinkhole` (nodes + ways, `out center`) and cave entrances named "cenote", bbox 19.3,-91.2,22.0,-86.7

Build (venv with numpy scipy rasterio pillow pyproj):

```
python -I pipeline/dem.py <raw dir> data
python -I pipeline/features.py raw/osm_sinkholes.json data
```

The grid is azimuthal equidistant on 21.40 N, 89.55 W at 400 m. `dem.png` packs elevation as (R·256+G)·0.1 m + offset; B = land flag.

## 3D viewer (`../sim/`)

```
python -I pipeline/paleo.py data data/features.json sim/data
```

Builds the 66 Ma target on a 360 km block (1024², 352 m) centred on the crater: `paleo.png` = water depth at 66 Ma,
`sed.png` = Mesozoic pile thickness (both (R·256+G)·0.1 m), `meta.json` = layer split, crust depths, today's coast,
labels. Ramp: 40 m → 2000 m along bearing 22.5° (smoothstep −60…140 km); pile 3 → 4 km along the same axis; off-platform
areas keep today's depth if deeper. `sim/index.html` was assembled once from the Mátra viewer's CSS plus new code
(block shaders, time column of the target's build-up); edit it directly. Impact phases are meant to follow the 'eve' band.

Impact phases (2026-10-24) live entirely in `sim/index.html`, no extra data: the time column's upper half (p 0.5–1)
runs −6 s → 10 min. `impact(t)` gives the crater shape parameters; `Wd(r, d)` (JS and GLSL copies, keep in sync) moves
rock that started d m below the seafloor; rock shallower than `dexAt(r)` is excavated; peak-ring sheet, melt and
breccia/ejecta lie on top (`deposits`). Face/wall fragments invert z(d) by bisection to find the starting layer.
Ejecta are stateless ballistic GPU points launched when the excavation front reaches them. Timings follow
Morgan et al. 2016 (T = 1, 3, 4, 5, 10 min snapshots); amplitudes are tuned by eye, not simulated.

Sea return (2026-10-24): `python -I pipeline/gulf.py <gmrt_gulf_wide.tif> sim/data` (GMRT med, lon −104…−74, lat 7…35).
2500 km Gulf grid (1024², 2.44 km) around the crater: today's bathymetry as a stand-in, 66 Ma block inside, Yucatán
peninsula → 40 m shelf. Tsunami crest: rim wave at 220 km at 10 min, 1.5 km high (Range et al. 2022), then fast
marching (scikit-fmm) at sqrt(g h); height ∝ 1/r × Green's law, capped at 0.8 h + 30 m. Viewer: time column now
build-up 0–0.4, impact 0.4–0.8, sea 0.8–1 (10 min → 1 day); resurge from the NNE, bore, seiches, steam, and
resurge/seiche suevite (+10 m by 1 h, +80 m by 12 h; Gulick et al. 2019).

66 Ma Gulf (replaces gulf.py output): `python -I pipeline/rzip.py <EarthByte paleoDEM_0.2Deg_grids.zip URL>
paleoDEM_0.2Deg_grids/Map16_PALEOMAP_1deg_KT_Boundary_65Ma.nc <dir>` fetches just the 65 Ma PALEOMAP grid via range
requests (Zenodo has no range support; its 0.1° archive is 207 MB), then
`python -I pipeline/gulf66.py <Map16 ... 65Ma.nc> raw/gws sim/data`. The grid is in the paleo frame: AEQD on the
impact site's 66 Ma position (−71.07, 25.49) turned 13° so block north matches; rotations, today's coastlines and city
positions come from gws.gplates.org (PALEOMAP, 66 Ma), cached in raw/gws. Low PALEOMAP land within 450 km of the
impact is turned into 40 m shelf (Yucatán platform, as in the block). gulf.py (today's geography) is kept for comparison.

v5 (2026-10-25): time column −3 h → +1 day (build-up removed). `python -I pipeline/swe.py sim/data` (after gulf66.py)
runs a nonlinear shallow-water model (numba; C-grid, Manning 0.025 sea / 0.045 land, wet/dry, breaking cap
η ≤ 0.8·depth + 30 m, sponge edges) from 10 min to 10 h on the 1024² Gulf grid, starting from the 10-min crater (ported
from the viewer: keep Hb 3500, Cu 21000, Rc 22000, U 300, pkA 500, melt level sf0−1300, breccia 250 in sync) and a ring
wave A0 = 500 m (SWE_A0), σ = 12 km at 220 km, calibrated to ~110 m median crests in the open Gulf (Range et al. 2022).
Frames: swe_far_*.png (512², 240 s), swe_near_*.png (160² around the crater, 120 s), 3 frames per pixel in RGB —
never in alpha (canvas premultiplication wipes it). Results: centre of the crater floods at ~60 min; ~114,000 km² of
coastal plain flooded; New Orleans 1.5 h, Houston 2.7 h.

v6 (2026-10-25): cutaway is a half-disc of 150 km radius (R_CUT, GLSL + JS) with a curved strata wall; ejecta blanket
300 m·(r/95 km)^-3, faded out between 110 and 165 km (viewer and swe.py in sync) so the crater merges into the map;
fine/coarse simulation frames blend over 25 km; land the tsunami reached is tinted (toggle "Land the tsunami reached");
"Fly along the flooded north coast" follows the 66 Ma north shore east→west from 1.2 h to 3.3 h.

v7 (2026-10-25), globe overview: `python -I pipeline/globe.py <Map16 ... 65Ma.nc> sim/data [calibrate]` builds the
65 Ma PALEOMAP Earth at 0.25° (paleo frame; Yucatán shelf fix) and runs a linear shallow-water model on the sphere
(lat-lon C-grid, walls at land and |lat| > 78°, breaking cap) from 30 min to 48 h, starting from a ring 400 km out
(GLOBE_A0 = 165 m, σ 60 km), calibrated to the same ~110 m open-Gulf crests as swe.py. Outputs globe_h.png, globe_w.png
(first arrival, max crest), globe_f_*.png (720×360, 20 min, RGB-packed) and meta.globe (incl. today's coastlines at
66 Ma from raw/gws/coastlines.json). Viewer: "Globe" dock button / G key / "The whole Earth"; globe data load lazily;
double-click the globe to dive into the crater. Time column now ends at 2 days.

v8 (2026-10-25), re-entry pulse, air blast, seismic waves: `python -I pipeline/reentry.py sim/data` (after globe.py)
flies 2 million ejecta particles (M(>v) ~ v^-1.65 at 2–11 km/s, 35–55° elevation, mild SW downrange excess) on Kepler
orbits from 100 km altitude over the turning Earth until they re-enter at 100 km. Particle mass is scaled so the distal
layer (5,000–15,000 km) averages 3 mm at 1,500 kg/m³. Heating of the upper atmosphere = kinetic energy per area per 240 s,
smoothed over ~2° (each particle stands for a swarm; one carries ~180 kW/m² into a 1° cell far out), not in time.
Outputs reentry_f_*.png (360×180, 120 frames to 8 h, 1 + 254·log10(W/m²)/7, RGB-packed), reentry_m.png (spherule layer,
1 + 254·(log10 mm + 3)/6) and meta.reentry. Results: peaks ~13 MW/m² at 1,000 km (10 min), ~1.3 MW/m² at 3,000 km,
~90 kW/m² at 8,000 km (38 min), ~18 kW/m² at the antipode (2.4 h); 15–65× sunlight far out, below Melosh et al. 1990's 50–150×. The ground
share is not modelled; the viewer quotes half of it unshielded and Goldin & Melosh 2009's ~6 kW/m² with self-shielding.
Air blast and seismic shaking are computed in the viewer from the Earth Impact Effects Program (Collins et al. 2005):
E = 4.76e23 J, M ≈ 10, overpressure/wind from the surface-burst fit, blast arrival by integrating the shock speed
(Lamb wave at 318 m/s far out, as measured for Hunga Tonga 2022, where the overpressure fit is flagged unreliable beyond 3,000 km); approximate
IASP91 P/S travel times to 100° (shadow 103–142°, core-refracted PKP ~19–20 min beyond); Rayleigh surface waves at 4 km/s (EIEP uses 5, which would
outrun S). Viewer: Effects switches, fronts as rings on the globe and the Gulf, blast wall in the regional view, sky
glow overlay from the re-entry flux above the camera target, and a probe (click the globe) summarising every effect.

v9 (2026-10-25), impact winter and the long warm-up: `python -I pipeline/climate.py sim/data` (after reentry.py).
Time column extended past 2 days on a log scale: weeks, months, years, 30 yr → 100 kyr (first two days = 78%).
Three scenarios, because the studies disagree about what darkened the sky; each scenario's aerosol optical depth over
time is tuned so a small energy-balance model gives the study's headline numbers:
dust (Senel et al. 2023) → model 1.6 yr under 1% light, −19 °C global (published ~1.7–2 yr, up to −15 °C);
soot (Bardeen et al. 2017, Tabor et al. 2020) → 1.9 yr, land −29, ocean −11, global −16 °C (published ~2 yr, −28/−11/−16.5);
sulfate (Brugger et al. 2017) → −26 °C, global mean below 0 °C for 6.7 yr, >30 yr recovery (published ≥ −26, 3–16 yr, >30 yr).
Spreading: 5° advection–diffusion (stratospheric winds, eddy mixing) of each agent from where the impact put it (distal
ejecta layer for dust and sulfur, strongly heated land for soot). Climate at 2.5°: 150 m ocean layer that overturns
into a 2,000 m deep ocean when chilled, land with continentality, Budyko B = 1.8, meridional diffusion; absorbing layers
re-radiate half of what they absorb downward; reflected share by two-stream reflectance (g = 0.7). Baseline
Maastrichtian sea surface 29 → 5 °C equator → pole, land 6.5 °C/km lapse. After ~300 yr: +5 °C warming over 100 kyr as
recorded at El Kef (MacLeod et al. 2018), 1.6× toward the poles, 1.2× over land — drawn as observed, not modelled.
Outputs climate_<id>.png (200 frames at 144×72, 20 per row; R light 1 + 254·(log10 + 6)/6, G sea / B land anomaly
128 + 3·°C) and meta.climate (times, global/land/sea/light/dark series). Viewer: globe veil (dust brown, soot black,
sulfate white), snow and sea ice below freezing, optional temperature-change map; regional view dims, hazes and
frosts; climate card with scenario switch and chart; probe gives darkness, coldest temperature and freeze duration.

v11 (2026-10-25), science audit. Checked against the literature: impactor energy (12 km, 2,630 kg/m³, 20 km/s →
4.76e23 J ≈ 1.1e8 Mt); 60° from the NE (Collins 2020); cavity, uplift and peak-ring timing (Morgan 2016); resurge within
the hour, ~10 m then ~80 m suevite, reflected tsunami within a day (Gulick 2019); open-Gulf crests ~100 m (Range 2022);
distal layer 2–3 mm, ~65 cm near 1,000 km; EIEP magnitude, Meff (Δ in radians, the only unit continuous at 700 km),
Mercalli table and overpressure; IASP91 travel times; impact-winter scenarios as in v9. Corrected:
far-field blast speed 330 → 318 m/s (Hunga Tonga Lamb wave); probe now gives core-refracted PKP (~19–20 min) beyond
142° instead of "core shadow"; re-entry heating described as below Melosh et al. 1990 (15–65× vs 50–150× sunlight)
rather than "the same order"; El Kef warming scaled to the site's ~25° N paleolatitude (was the equator; +13% too warm);
climate.py now fades the model's land anomaly after 300 yr instead of dropping it.

v12 (2026-10-25), one zoomable globe: the globe and the regional model share one scene, camera and uniforms. The
regional map is AEQD about the crater, so every regional vertex shader wraps it exactly onto a sphere of 6,371 km
(GLSL bend(), JS bendV/unbendV); the Earth's centre is 6,371 km below the crater and the globe mesh is rotated so the
66 Ma impact site is at the top with paleo north turned by the frame angle (GM_ROT). The globe skips the Gulf square
and uses the regional colour, light, fog and sea-surface model; the regional meshes use the globe's climate veil,
snow, re-entry glow and (far out) log-scale tsunami tint, so they meet without a seam. Navigation: OrbitControls in
metres; beyond ~3,000 km altitude the pivot sinks to the Earth's centre (fully at 12,000 km) so the globe turns, and
coming back in the camera swings round to the impact region. Globe relief exaggeration rises from the regional value
near the ground to 18× far out. Views and the coast tour are defined on the flat map and wrapped (flyF, bendV).

v13 (2026-10-25), round, blended boundaries. Viewer: the block's surface and sea are a disc of 172 km (R_BLK) inside
its square; the Gulf grid ends at a square of half-size 1,240 km with corners rounded by 550 km (gulfSD), and across its
outer 220 km the Gulf meshes hand over to the globe's data (height, land, global tsunami, relief exaggeration), so the
globe, which skips exactly that shape, meets an identical surface; the block's sea now gets the far-view tsunami tint.
Data: paleo.py bends the ramp's depth contours into arcs about a point 450 km SSW, with ~15 km of low-frequency wander
(same depths at the crater), and merges into meta.json instead of rewriting it; gulf66.py blends the block over a ring
165–235 km (was a 40 km square ring), and builds the Yucatán shelf from a smoothed, noise-wiggled mask of PALEOMAP land
(all land within 380 km; low land to ~520 km), sloping into deep water instead of 1° steps. swe.py re-calibrated:
SWE_A0 = 290 m (open-Gulf median crest 120 m; New Orleans 1.55 h, Houston 2.7 h; 121,000 km² flooded).

v16 (2026-10-25), the asteroid's arrival and wildfires (viewer only, no new data files).
Arrival: the impactor is a dark, lumpy carbonaceous body (Fischer-Gödde et al. 2024) lit by the Sun until it enters the
atmosphere (100 km, ~5.8 s before contact); then a shock cap on its leading face, an air tunnel behind it back to the
entry point, and its light on the sea, all scaled by air density exp(-h/8 km) (the visible glow also ramps from 100 km,
as meteors do). Playback runs at 0.22× through the last 8 s.
Wildfires: fireCodes() turns the re-entry heating frames into ignition times per 1° cell: ground heat ≥ 12 kW/m² for
4 min (dry litter) or 20 min (live forest, after Belcher et al. 2015), the ground getting 0.1 (Patchy: spherule
self-shielding, Goldin & Melosh 2009) or 0.5 (Firestorm: Melosh et al. 1990) of the upper-air heating. Patchy: litter
on ~32% of vegetated land, almost no forest; Firestorm: litter ~88%, forest ~47%. Near the crater the pulse is intense
but ~12 min long, so forest does not ignite there in either case. Drawn per fragment on the globe and the Gulf: patchy
fuel, a front spreading over hours, smouldering for days, smoke tint, burnt ground greening over 2–25 years. Options →
Wildfires (Off / Patchy / Firestorm; the soot climate scenario switches Patchy to Firestorm); probe section.
Pitfall: tNoise tiles every 1/divisor; divisors of a few km show as a regular lattice from afar — use long, unmatched
periods (e.g. 260 km, 61 km, 17.3 km).

v17 (2026-10-25), fixes. gulf66.py: the Yucatán platform around the block now continues the block's ramp (curved
contours about a point 450 km SSW, ±18 km wander; 40 m SSW → 2 km NNE, PALEOMAP kept where deeper), fading into
PALEOMAP between ~430 and 650 km along a wandering edge. The v13 rule "all land within 380 km → 40 m shelf" had made a
shallow ridge beyond the deep NNE ramp that glowed cyan. swe.py re-calibrated to SWE_A0 = 260 m: open-Gulf median
crest 127 m, 90th percentile 353 m (the median is not monotonic in A0 here: breaking dominates); New Orleans 1.55 h,
Houston 2.7 h; 117,000 km² flooded. Viewer: the far-out point of light hands over to the asteroid body at STAR_D
(650 km, 32.5 s) instead of both showing between 650 and 800 km; zooming in from the globe no longer swings the view
back to the crater: far out, the ground pivot follows the point below the camera, so zooming goes straight down onto
what is in view (double-click or G flies back to the crater); the camera stays at least 3 km above the globe.

v18 (2026-10-25): no fixed night side. Every place is lit as the crater region is (sun 38° up in the SE): the regional
meshes use uSunDir in their own flat-map frame, and the globe carries the same local sun round in that frame
(fromFlat). Season: boreal spring / austral autumn (During et al. 2022, Tanis fish), stated in the story and About;
the climate model stays annual-mean.

v19–20 (2026-10-25), navigation: drag grabs the map. The orbit controls no longer handle the left button or one
finger; a grab handler turns camera, pivot, ground pivot and camera.up together about the Earth's centre so the point
under the cursor stays there (a slide near the ground, a turn of the globe far out). Right-drag / Shift-drag still
orbits; two fingers turn and zoom. Ground-following corrections run only near the ground (they had pulled the camera
when the pivot was the Earth's centre), and the view re-levels gently toward the local vertical when not grabbing,
not far out and not looking straight down. Tested: grabbed point stays within 0 km of the cursor at 46 km to 15,000 km.

v21 (2026-10-25): the far asteroid is no longer a glowing stand-in held at 650 km. It is drawn at its true position all
the way in: a sharp point of reflected sunlight (about magnitude +1 three hours out) while its 12 km span fewer than
~5 pixels, then the dark body. Its label uses the true world position (not the map wrap). It comes in 60° up in the
NE, so the opening view no longer shows it; the story text says so.

v22 (2026-10-26)
Reproducible build: the raw inputs now live in raw/ (GLO-90 tiles, GMRT, Gulf GMRT subsets, the 65 Ma PALEOMAP grid,
OSM sinkholes, the GPlates cache; ~150 MB) and the environment is pinned in pipeline/requirements.txt.
  python3 -m venv .venv && .venv/bin/pip install -r pipeline/requirements.txt
  pipeline/build.sh                      # dem → features → paleo → gulf66 → swe → globe → reentry → climate → webp
Verified: a full rebuild into a scratch copy reproduced every data file bit for bit. globe.py now finds the GPlates
cache next to the script (it looked relative to the output folder).
WebP: pipeline/webp.py writes lossless, bit-checked .webp copies (27.0 → 21.8 MB); the viewer loads those.
Streaming: the first view needs meta.json, paleo, sed, gulf_h/w and globe_h/w (~2.9 MB); the regional tsunami, the
global tsunami and re-entry frames, and the climate frames stream in afterwards (each effect switches on when its
data arrive; a note shows if time runs ahead of them).
Quality: Options → Quality (Auto adapts resolution, phones may drop to 0.75×; High; Low = 1×, no bloom).
Tsunami check vs Range et al. 2022 (globe run): Gulf open-water median 131 m (>100), Pacific via the seaway 4.7 h
(~4), Indian Ocean 21.9 h (~24), North Atlantic coasts ≥10 m 59% (“many”), S Pacific 21% (“some”), sheltered Tethys /
Indian / S Atlantic medians 1.0 / 1.4 / 2.2 m. No recalibration needed.
Viewer: turning/tilting (right-drag, Shift-drag, two-finger twist, Q/E/R/F) is about the local vertical at the pivot;
overlapping labels give way by rank; a sky dome (blue by day, darkened and tinted by the impact-winter haze, red during
re-entry, black with stars above ~160 km, dark underground) also sets the distance-fog colour.

v23 (2026-10-26): fixed the Gulf drawn as a flat pale wash with a visible rounded-square edge. The v22 streaming patch
had put a comment in front of the uSw2 uniform (near-frame count, grid sizes, wave-height scale), so the Gulf shaders
read zeros and sampled a single texel. Also: the Gulf sea now calms after the simulated 10 h as the block's sea does
(it had kept the last frame's heights, a pale ring round the block); tsunami mud on flooded land is grey-brown and
overgrown within ~12 years; burn scars fade between 1.5 and 10 years (the fern spike), not 2 to 25.

v24 (2026-10-26)
Seasons in the impact winter (climate.py): the sunlight the aerosols take away is each day's own (daily-mean
insolation by latitude and date, obliquity 23.4°), starting at solar longitude 45° (mid boreal spring, ~5 May; During
et al. 2022). Anomalies are relative to the normal for the time of year; frames are instantaneous for ~2 years and
blend into running annual means by 4 years (later frames are too far apart to sample the seasons). Effect, dust
scenario: northern land 40 °C below its normal summer 3 months in, while the south (into winter) is 12 °C below;
the hemispheres then alternate. Global numbers barely move (dust −19.6 °C at 2.2 yr, was −19.1 at 2.6 yr; soot land
min −32, was −29; sulfate unchanged at −26.2).
Debris cloud above the atmosphere (viewer): 40,000 points (16,000 on phones) drawn from reentry.py's launch law (2–11
km/s, M(>v) ~ v^-1.65, 35–55° up, SW excess), flying exact two-body orbits in the vertex shader (Kepler's equation,
checked against leapfrog integration to <0.25 km) in the turning Earth's frame, from the crater (through the air column
the plume clears) until they are back at 100 km; the slowest land near the crater. Parcels fade within a few of their
own widths of the camera, where the crater's curtain and fireball show the debris in finer detail. Hot for
the first minute (the condensing vapour plume), then a dusty veil that spreads at ~400 m/s per parcel, glowing as it
falls back in. Opacity floor set from the distal layer's optical depth (~10 when spread globally; Goldin & Melosh
2009).
Fine coastlines (new pipeline/coast.py, after globe.py): signed distance to the shore at 610 m over the Gulf
(coast_gulf.png, 4096²) and 0.0625° for the globe (coast_globe.png), from cubic fits to the Gulf grid and the 0.2°
PALEOMAP grid, plus seeded sub-grid wiggles capped by slope (shore moves ≤1.5 km / ≤6 km; illustrative). Distances,
not heights: near-binary heights on steep coasts made texel-scale scallops, and whole-metre heights made staircases
on flat ones. The viewer takes land/sea from the fine sign and heights from the coarse grids; where the 2.4 km
tsunami grid had land but the fine coast has sea, the water stays at rest. meta.coast_fine (meta.coast is paleo.py's).

v26 (2026-10-26): higher-resolution continents. The input was PALEOMAP's 1° map run through a 400 km Gaussian
(grdfilter -Fg400) and resampled to 0.2° (raw/paleodem02, EarthByte's "paleoDEM_0.2Deg_grids"): ~24 m rms of real
sub-degree relief. Now the 0.1° (6 arc-minute) PALEOMAP PaleoDEM for 66 Ma (Scotese & Wright 2018, Zenodo 5460860,
CC BY 4.0; Map16_PALEOMAP_6min_KT_Boundary_65Ma.nc in raw/paleodem6min, fetched with rzip.py): ~156 m rms below 1°,
mountains to 3,600 m (was 2,700). Its heights come in 40 m steps with 0 m cells being sea; all scripts subtract 20 m
so sea level falls midway (land > 0 had grown every coast by up to half a cell). globe.py now area-averages the source
to 0.25° for the tsunami and writes globe_t.png (3600x1800, the source at full resolution) for drawing; the viewer
shades the globe from it per pixel (normals from the 0.1° heights at the drawn exaggeration, the step widening to the
pixel footprint and handing back to the mesh normal far out). coast.py: bilinear + light Gaussian instead of a cubic
spline (the spline's ringing lifted shallow seas above sea level).
Re-calibration on the new geography: SWE_A0 = 465 m (open-Gulf median crest 128 m, 90th 342 m; New Orleans 1.8 h,
Houston 3.0 h, Tampa 2.1 h, Veracruz 1.2 h; 74,000 km² flooded); GLOBE_A0 = 205 m (open-Gulf median 135 m). Range 2022
checks (same script on old and new data): Pacific via the seaway 5.1 h (was 4.7), Indian Ocean 24.2 h (24.2), North
Atlantic coasts ≥10 m 46% (53%), S American Pacific 18% (29%), Tethys / Indian / S Atlantic medians 1.2 / 1.6 / 2.5 m.
Climate unchanged within 0.1 °C.

v27 (2026-10-26): true-colour view. New pipeline/landcover.py (after coast.py) writes globe_c.png (3600x1800 sRGB
land cover on globe_t's grid): moisture from latitude belts (wet tropics, dry subtropics, wet storm tracks, moderately
wet poles), raised within ~1,000 km of the sea (great-circle distance to the nearest shore, KD-tree) and on high
ground; Maastrichtian annual-mean temperature less 6.5 °C/km; biome colours from today's true-colour mosaics (no
grassland: scrub instead); seeded patchiness; sea cells carry the nearest land colour. Result: forest ~80% of land,
scrub/woodland ~19%, true desert <1%, rock and snow only on the highest ground. Viewer: land on the globe and in the
Gulf from it (sRGB texture; the old height palette until it loads). Water: absorption of clear seawater, two-way
(0.22/0.045/0.018 per m, was ~10x too weak in red, which showed pale seabeds through tens of metres), and wet-sediment
seabed at half the rock's albedo: shelves turquoise, deep sea dark blue. Sun glint off from far out (with each place
lit by the crater's own sun it traced an arc). Fixed the atmosphere shell: its rim term mixed a normal in the globe's
own frame with a world-space view direction, which drew a pale band across the disc.

v28 (2026-10-26): a real sky. New pipeline/sky.py (inputs in raw/sky: NASA SVS Deep Star Maps 2020 milkyway_2020_4k.exr,
public domain; the Yale Bright Star Catalogue BSC5, CDS V/50) writes sky_mw.webp (4096x2048 lossy sRGB, gain 2) and
sky_stars.json (8,404 stars to V 6.5: J2000 RA/Dec, V, B-V). The NASA map is plate carree with RA increasing to the
left and 0h at the centre (checked: galactic centre at u 0.75, v 0.66; the LMC stands out only under this reading).
Viewer: the sky dome samples the Milky Way in space (replacing the random dots); the stars are points far out (sized
and brightened by flux, coloured by B-V; the globe hides them by depth), shown above ~25-160 km. Celestial pole = the
66 Ma spin axis, RA 0h over the paleo 0 deg meridian (time of day unknown; the sky is held still). It is today's sky:
at 66 Ma the stars stood elsewhere and many had not formed; the About text says so. Needs OpenEXR (requirements.txt).

v29 (2026-10-26): wildfires redrawn (viewer only). The old look (noise patches with a hard threshold, each patch
glowing orange at once, flat grey smoke) is replaced by spreading fires: ignition points jittered in cells (14 km in
Firestorm, 26 km in Patchy) light within the first hour after the heat makes the ground ignitable (fire codes as
before); the front arrives at distance / rate of spread (~0.42 m/s ≈ 1.5 km/h, varying), on domain-warped coordinates
so fronts are ragged and merge; fuel is a warped fractal with soft edges (same burnt shares as before). Flames are a
~25 min (~600 m) band behind the front, widened and dimmed where the band is narrower than a pixel (fwidth of the
arrival time, passed in from the fragment shaders: GL_COMMON also goes into vertex shaders); smouldering patches for
~12 h; charred ground with ash drifts, fading as before. Smoke: a smooth regional burning share gathered from 0-66 km
upwind (trade winds west in the tropics, westerlies beyond), the plume widening downwind, drawn out into streaks
along the wind; only land counts as a source.

v30 (2026-10-26): the wildfires' continental pattern. The noise texture tiles at each layer's period, so the wide
layers (fuel at 181/41 km, warps, the smoke's cover and streaks) repeated across a continent as a lattice of tufts and
a diagonal hatching. Layers of 40 km and wider now use a period-free hash value noise (lnz: two octaves, the texture's
statistics); the wind meanders ±35° over hundreds of km; streaks are built from 181/61/19 km layers on warped
coordinates, the finest fading once a pixel spans more than a few km (pixel size passed in from the fragment
shaders); a ~900 km field leaves some regions under heavy pall and others clear; the burnt region's edge is jittered.
Dropped the smoke term shaped like each burnt patch (it drew pale stars).

v31 (2026-10-27): prevailing winds for the fires. New step wind.py: the April and May climate means of the BRIDGE
HadCM3L 66 Ma run (texpn1; Valdes, Scotese & Lunt 2021, on Scotese's PALEOMAP maps: its land and sea match
globe_h.png in 98% of its cells with no shift), downloaded to raw/hadcm3_texpn1/ from
www.paleo.bristol.ac.uk/ummodel/data/texpn1/climate/ (texpn1a.pdcl*.nc, texpn1a.pccl*.nc), blended 1/3 : 2/3 for
5 May, resampled cubically to the fire codes' 1 deg grid: the 10 m wind and the 925-850-700 hPa mean (wind.png,
meta.wind). One missing 10 m cell is filled from its neighbours. The viewer turns east/north into the regional
map's frame per pixel (a 2x2 Jacobian of uvOfFlat). Spread: each ignition grows as an ellipse with its ignition near
the back (Van Wagner 1969), length-to-breadth from the 10 m wind (Alexander 1985), head rate 0.2 m/s x
exp(0.05039 U km/h) (the Canadian FBP wind function), back rate from the ellipse; closed-form arrival time. The
ignition search is 4x4 cells, two upwind, so fast heads do not stop at cell borders. Flames scale with the front's
rate^0.46 (Byram). Smoke: drifts with the smoke-layer wind, reaching further in a stronger wind, wandering about it
more where the mean is weak (middle latitudes, where the monthly mean averages changeable weather).

v32 (2026-10-27): glitches and an Atmosphere switch.
- New step topo.py (after globe.py, which no longer writes it) writes globe_t.png with PALEOMAP's 40 m height steps
  smoothed out: repeated Gaussian smoothing at falling scales, each pass clipped back into the cell's own step (+-20 m,
  land and sea keeping their sign). The source's -180 and 180 columns differ, so the grid is not wrapped. landcover.py
  rerun on it. Drawing only; the tsunami inputs are unchanged.
- The "wood grain" contour rings were in the smoke, not the ground: the streak pattern rotated map coordinates thousands
  of km long by the (meandering) wind direction, shearing it into rings, and its warp changed faster than its own
  size (lnz has cells an eighth of its nominal scale). Streaks are now one pattern laid along three fixed directions,
  blended by closeness to the wind; warp and meander use vnz at their nominal scales.
- Far from the crater the fires' map is the azimuthal equidistant map about the antipode, folded onto the crater's at
  90 deg (squeeze <= pi/2 instead of 3x at 135 deg). Wind directions on the map are computed from the sphere (the
  texture-difference Jacobian had ~0.1% noise), and the globe passes the map coordinate from its vertices (per-pixel
  atan put steps of a few hundred metres into it far out); both had turned the fine streaks into bricks.
- Dark dots along the Gulf/globe seam: the Gulf's sea surface and its seabed (and the globe's seabed) each tested
  "inside the Gulf" at their own depth; seabeds now test where the line of sight crosses sea level (seaCrossXZ).
- Re-entry glow shimmer: smooth triplanar noise instead of a random value per 0.125 deg texel (drew squares).
- Smoke sources follow the fine coast (soft over ~2 km) instead of the 0.25 deg land grid (stair steps).
- Options > Atmosphere (uAir): off hides the impact winter's veil and darkness, fire smoke, distance fog, the sky's
  haze and the limb glow.

v33 (2026-10-27): fossil sites and K-Pg boundary sections.
- New step sites.py: 17 sites (Tanis, Hell Creek, Raton Basin, Brazos River, El Mimbral, Beloc, Blake Nose ODP 1049,
  Demerara Rise ODP 1259, Gubbio, Stevns Klint, Caravaca, Zumaia, El Kef, Maastricht, Seymour Island, Woodside Creek,
  the Lameta beds), today's positions rounded to ~0.1 deg, moved to 66 Ma with the PALEOMAP plate model (rotations and
  static polygons from the GPlates plate model repository, Zenodo 10.5281/zenodo.10596610, kept in raw/paleomap_model)
  by pygplates (added to requirements.txt). The GPlates Web Service was down (503); the local reconstruction reproduces
  its cached results for the crater and the towns to 0.0001 deg. Writes meta.sites with each site's record in a few
  sentences, its references, its setting then (land/sea) and the PaleoDEM height under its 66 Ma position.
- The map and the rocks disagree at nine of the twelve marine sites: the PaleoDEM has land under them (Brazos River,
  Zumaia, El Kef, Maastricht, Woodside Creek; Gubbio, Stevns Klint, El Mimbral and Beloc with sea within 30 km).
  PALEOMAP's coastlines are approximate. Markers stay at the reconstructed positions; the viewer's card says so.
- Viewer: site labels (diamonds = boundary sections, circles = fossil sites) on the globe and close up; click one for
  its record above the probe's model results for the same place. Options > Fossil sites and boundary rocks; Views >
  Fossil sites and boundary rocks.

v35 (2026-11-08): navigation bug fixes.
- Zooming toward the cursor no longer uses the orbit controls' cursor zoom, which moved the camera sideways and re-aimed
  it, so the view swung round and could spin. The camera now moves along the line of sight and the Earth is carried, as
  in a grab, so the point under the cursor stays there with tilt and roll unchanged.
- The view's up can no longer lie along the line of sight, where it stops fixing the roll. Tilting from straight down
  keeps the heading, and the compass no longer jitters looking straight down.
- Keyboard moves work anywhere on the Earth (they had been undone by the pivot each frame). North up keeps the same place
  in view on the globe. The relief switch leaves the camera alone far from the crater. Views arrive level. Leaving Fly
  mode keeps the view, and choosing Fly on the globe goes to the crater first. Grabbing works from below the ground.

v36 (2026-11-08): Hungarian.
- sim/hu.js maps the viewer's English text to Hungarian (long texts by their first 60 characters); site records and names
  are in it too, so sites.py and meta.json stay English. Values in generated text go through tr() with named placeholders;
  the {x|val} and {x|ig} filters add Hungarian case endings to times. The 2D map keeps its smaller dictionary inline.

v37 (2026-11-08): guided tour.
- A captioned run through the whole event in about 3.5 minutes, started only by the Guided tour button or Views > Watch
  (GUIDE in sim/index.html: each scene has a camera shot, a time span and a caption). Navigating by hand pauses it; Play
  returns to the scene's view. One scene steps back in time, to the first hour at the crater, and its caption says so.
