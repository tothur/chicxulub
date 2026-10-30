# Chicxulub

Interactive reconstructions of the Chicxulub impact, 66 million years ago, in the browser (three.js, no build step).

- **`sim/`**: the 3D viewer. One zoomable scene from the crater out to the whole Earth as it was then: the target rocks,
  the crater forming, the sea's return and tsunami, ejecta re-entry, air blast, seismic waves, wildfires driven by the
  prevailing winds, the impact winter and the long warm-up, plus fossil sites and K–Pg boundary sections at their
  66 Ma positions. Open `sim/index.html` from a local web server (e.g. `python3 -m http.server` in `sim/`).
- **`index.html`** + **`data/`**: the 2D regional map of today's crater (terrain, cenote ring, drill sites).
- **`pipeline/`**: the Python scripts that build every data file from the raw inputs; `pipeline/build.sh` runs them in
  order, and `pipeline/README.md` describes the sources, methods and each version.

The viewer's About panel lists the models, their limits and the sources.

## Rebuilding the data

```sh
python3 -m venv .venv && .venv/bin/pip install -r pipeline/requirements.txt
pipeline/build.sh
```

Most raw inputs are not in the repository (download locations are in `pipeline/README.md`); `raw/` keeps only the
OpenStreetMap cenotes, the cached GPlates reconstructions and the PALEOMAP plate model.

## Data credits

PALEOMAP PaleoDEM and plate model (Scotese & Wright 2018, CC BY 4.0), Copernicus GLO-90, GMRT, OpenStreetMap
contributors (ODbL), HadCM3L run texpn1 (Valdes et al. 2021), NASA SVS Deep Star Maps 2020, Yale Bright Star Catalogue.

## License

The code is under the MIT License (see `LICENSE`). Data derived from the sources above keeps their terms
(PALEOMAP CC BY 4.0, OpenStreetMap ODbL, and so on).
