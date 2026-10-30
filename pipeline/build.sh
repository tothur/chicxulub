#!/bin/sh
# Rebuild every data file of the 2D map (data/) and the 3D viewer (sim/data/) from the raw inputs in raw/.
# Order matters: gulf66 rewrites meta.gulf, swe adds meta.gulf.swe; globe, topo, coast, landcover, reentry, climate, wind and sites add their own keys.
#   usage:  pipeline/build.sh                (from anywhere; uses .venv/bin/python unless PY is set)
#           MAP=/tmp/m SIM=/tmp/s pipeline/build.sh   to build into other folders (they need a copy of data/ and sim/data/)
set -e
cd "$(dirname "$0")/.."
PY=${PY:-.venv/bin/python}
MAP=${MAP:-data}; SIM=${SIM:-sim/data}
PALEO=raw/paleodem6min/Map16_PALEOMAP_6min_KT_Boundary_65Ma.nc    # Scotese & Wright 2018, 0.1 deg (Zenodo 5460860, CC BY 4.0)
step() { echo "== $*"; "$PY" -I "$@"; }
step pipeline/dem.py raw "$MAP"                                  # today's terrain for the 2D map (GLO-90 + GMRT)
step pipeline/features.py raw/osm_sinkholes.json "$MAP"          # cenote ring, crater centre, places
step pipeline/paleo.py "$MAP" "$MAP/features.json" "$SIM"        # the 66 Ma target block
step pipeline/gulf66.py "$PALEO" raw/gws "$SIM"                  # the 66 Ma Gulf grid (GPlates rotations cached in raw/gws)
step pipeline/swe.py "$SIM"                                      # the regional tsunami
step pipeline/globe.py "$PALEO" "$SIM"                           # the 66 Ma Earth and the global tsunami
step pipeline/topo.py "$PALEO" "$SIM"                            # the 0.1 deg topography for drawing, its 40 m steps smoothed out
step pipeline/coast.py "$PALEO" "$SIM"                           # fine coastlines (signed distance) for the Gulf and the globe
step pipeline/landcover.py "$SIM"                                # land cover in true colour (latitude belts, distance from the sea, height)
step pipeline/sky.py raw/sky "$SIM"                              # today's stars (Yale BSC) and Milky Way (NASA SVS Deep Star Maps 2020)
step pipeline/reentry.py "$SIM"                                  # ejecta re-entry heating
step pipeline/climate.py "$SIM"                                  # impact winter scenarios
step pipeline/wind.py raw/hadcm3_texpn1 "$SIM"                    # prevailing winds of 5 May for the fires (HadCM3L 66 Ma, Valdes et al. 2021)
step pipeline/sites.py raw/paleomap_model "$SIM"                  # fossil sites and K-Pg boundary sections at their 66 Ma positions (PALEOMAP, pygplates)
step pipeline/webp.py "$SIM"                                     # lossless WebP copies the viewer loads
echo "done: $MAP $SIM"
