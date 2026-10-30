"""Famous fossil sites and K-Pg boundary sections, at the places their rocks stood 66 Ma ago.

Today's positions (rounded to about 0.1 deg; some excavations publish no exact spot) are moved to 66 Ma with the
PALEOMAP plate model (Scotese; rotations and plate polygons from the GPlates plate model repository, in
raw/paleomap_model), reconstructed here with pygplates: the same model and age the GPlates Web Service used for the
crater and the Gulf's towns (gulf66.py), which this reproduces to 0.0001 deg, so they sit on the same 66 Ma map. Each
site's setting then (land or sea) is checked against the PaleoDEM under its 66 Ma position, as a test of both.

kind  'boundary'  a section through the boundary layer itself (iridium, spherules, shocked quartz, the extinction)
      'fossil'    a fossil site known for the life of the last Cretaceous days

Writes meta.sites (name, place today, today's and 66 Ma positions, kind, setting, distance from the impact, the
record in a few sentences, references).

usage: python sites.py <PALEOMAP plate model dir> <sim data dir>
"""
import json, math, os, sys
import numpy as np
import pygplates
from PIL import Image

MODEL_DIR, OUT = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(OUT, "meta.json")))
AGE, MODEL = 66, "PALEOMAP"

# name, place today, lat, lon, kind, setting then, record, references
SITES = [
    ("Tanis", "Hell Creek Formation, North Dakota", 46.1, -103.4, "fossil", "land",
     "A river bank buried within the first hour or so: sturgeon and paddlefish piled on a sandbar under about 2 m of "
     "sediment, impact spherules up to 5 mm in their gill rakers, spherule pits in the mud and spherules caught in amber. "
     "Its excavators read it as surges of river water set off by the seismic shaking (a seiche), not a tsunami, which "
     "could not have come this far in time. Bone growth in its fish points to spring in the north. The site has been "
     "studied mostly by one team, and parts of its interpretation are disputed.",
     "DePalma et al. 2019 (PNAS); During et al. 2022 (Nature)"),
    ("Hell Creek", "Garfield County, Montana", 47.5, -106.9, "fossil", "land",
     "The best-known record of the last dinosaurs: Tyrannosaurus, Triceratops, Edmontosaurus and Ankylosaurus in river "
     "and floodplain sands and muds, with mammals, lizards, turtles and birds. At the top lies the boundary clay, "
     "centimetres thick, with iridium and shocked quartz under a coal seam; no dinosaur bones are found above it, and more "
     "than half the kinds of leaves disappear.",
     "Johnson 1992; Clemens & Hartman 2014 (GSA Special Paper 503)"),
    ("Raton Basin", "near Trinidad, Colorado", 37.1, -104.5, "boundary", "land",
     "Where the iridium layer was first found in rocks laid down on land: a claystone about a centimetre thick in the coal "
     "swamps' beds, with iridium, shocked quartz and spherules. Above it the pollen of flowering plants gives way to fern "
     "spores (the 'fern spike'): ferns recolonised the bare ground first.",
     "Orth et al. 1981 (Science); Izett 1990"),
    ("Brazos River", "Falls County, Texas", 31.1, -96.8, "boundary", "sea",
     "A muddy shelf of the Gulf coast then, tens of metres to a hundred metres deep. At the boundary lies a sandstone bed "
     "with torn-up mud clasts, shell debris and impact spherules, read as laid down by the tsunami and the currents the "
     "impact set off. Whether it lies exactly at the boundary has been argued.",
     "Bourgeois et al. 1988 (Science); Keller et al. 2007"),
    ("El Mimbral", "Tamaulipas, Mexico", 23.2, -98.8, "boundary", "sea",
     "A layered clastic unit up to 3 m thick between the marls of the deep Gulf floor: impact spherules with glass cores at "
     "its base, then sandstone, then rippled sands, under an iridium-rich clay. Read as the work of the impact's slumps, "
     "turbidity currents and tsunami backwash; Keller and Stinnesbeck argued it built up over a longer time.",
     "Smit et al. 1992 (Geology); Stinnesbeck et al. 1993"),
    ("Beloc", "Haiti", 18.4, -72.6, "boundary", "sea",
     "Deep Caribbean sea floor then: a bed of impact spherules up to tens of centimetres thick, some still with cores of "
     "black and yellow glass whose make-up matches melted Yucatán rock and its evaporites. With shocked quartz and "
     "iridium it helped point to Chicxulub as the crater.",
     "Sigurdsson et al. 1991 (Nature); Izett 1991; Maurrasse & Sen 1991"),
    ("Blake Nose, ODP 1049", "off Florida, western Atlantic", 30.1, -76.1, "boundary", "sea",
     "A drill core from the Atlantic slope: a graded bed of green spherules some 10 to 17 cm thick, with shocked quartz and "
     "an iridium layer on top; the plankton below it are Cretaceous, those above it Paleocene, with almost nothing "
     "surviving across.",
     "Norris et al. 1999 (Geology); Klaus et al. 2000"),
    ("Demerara Rise, ODP 1259", "off Suriname, western Atlantic", 9.3, -54.2, "boundary", "sea",
     "Drill cores from a submarine plateau: a spherule bed about a centimetre or two thick, the sudden loss of most "
     "plankton above it, and the slow return of life to the open sea in the cores that follow.",
     "MacLeod et al. 2007 (GSA Bulletin)"),
    ("Gubbio", "Bottaccione Gorge, Italy", 43.4, 12.6, "boundary", "sea",
     "Pink limestone of the deep Tethys sea floor, broken by a clay about a centimetre thick. Walter and Luis Alvarez's "
     "team measured iridium in the clay at some 30 times the level in the limestone and proposed in 1980 that an asteroid "
     "had struck and caused the extinction. The foraminifera of the limestone change abruptly across it.",
     "Alvarez et al. 1980 (Science)"),
    ("Stevns Klint", "Zealand, Denmark", 55.3, 12.4, "boundary", "sea",
     "A 15 km sea cliff in the chalk of a warm northern sea: white Cretaceous chalk, the dark Fish Clay (Fiskeler) a few "
     "centimetres thick with fish remains, iridium (about 160 times background, one of the three 1980 sites) and soot "
     "from the fires, then Paleocene limestone built by bryozoans. A UNESCO World Heritage Site since 2014.",
     "Alvarez et al. 1980 (Science); Wolbach et al. 1985 (Nature); Surlyk et al. 2013"),
    ("Caravaca", "Barranco del Gredero, Spain", 38.1, -1.9, "boundary", "sea",
     "One of the most complete marine sections: a dark boundary clay a few centimetres thick with a reddish basal layer "
     "rich in iridium and nickel-rich spinel crystals, where most planktic foraminifera vanish. Its iridium was found the "
     "same year as Gubbio's, independently.",
     "Smit & Hertogen 1980 (Nature)"),
    ("Zumaia", "Basque coast, Spain", 43.3, -2.3, "boundary", "sea",
     "Flysch cliffs, deep sea floor then, with millions of years in tilted beds: the boundary is a thin dark clay with "
     "iridium over the last Cretaceous marls. The same cliffs hold the reference sections of two later Paleocene stages.",
     "Smit & ten Kate 1982; Schmitz et al. 2011 (Episodes)"),
    ("El Kef", "Tunisia", 36.2, 8.6, "boundary", "sea",
     "The official reference section for the start of the Paleogene: the boundary is set at the base of a dark clay with a "
     "thin reddish layer rich in iridium, where most planktic foraminifera die out. Oxygen isotopes of fish teeth from "
     "here record the warming that followed (about 5 °C), which this viewer's long warm-up follows.",
     "Molina et al. 2006 (Episodes); MacLeod et al. 2018 (Science)"),
    ("Maastricht", "Sint-Pietersberg, Netherlands", 50.8, 5.7, "fossil", "sea",
     "The town the Maastrichtian, the last age of the Cretaceous, is named after. Its chalk quarries gave up the first "
     "mosasaur skull in the 1760s and 1770s, later named Mosasaurus hoffmanni: the sea's top predators, marine reptiles "
     "up to some 15 m long, which died out at the boundary.",
     "Dumont 1849; Mulder et al. 2004"),
    ("Seymour Island", "Antarctic Peninsula", -64.3, -56.8, "boundary", "sea",
     "A shelf sea far south then: the López de Bertodano Formation, rich in ammonites, snails, clams and marine reptiles up "
     "to the boundary, which carries iridium and a bed of fish remains just above it. Even this far from the crater the "
     "extinction was sudden.",
     "Elliot et al. 1994 (Geology); Witts et al. 2016 (Nat. Commun.)"),
    ("Woodside Creek", "Marlborough, New Zealand", -41.9, 174.2, "boundary", "sea",
     "Limestone of the open southern ocean, among the farthest boundary sites from the crater: the third of the 1980 "
     "iridium sites, with a thin dark clay that also holds soot.",
     "Alvarez et al. 1980 (Science); Wolbach et al. 1985 (Nature)"),
    ("Lameta beds", "Jabalpur, India", 23.2, 79.9, "fossil", "land",
     "Nesting grounds of the last dinosaurs of India: titanosaur eggs and bones in the Lameta beds, beneath and between the "
     "first lavas of the Deccan Traps. The Deccan eruptions were flooding India with lava before and after the impact; "
     "how much they added to the extinction is still argued.",
     "Mohabey 1996; Schoene et al. 2019 (Science); Sprain et al. 2019 (Science)"),
]


def reconstruct(lonlat):
    """Today's points to AGE: each takes the plate of the static polygon it falls in (None where there is none)."""
    rot = pygplates.RotationModel(os.path.join(MODEL_DIR, "Rotations", "PALEOMAP_PlateModel.rot"))
    polys = pygplates.FeatureCollection(os.path.join(MODEL_DIR, "StaticPolygons", "PALEOMAP_PlatePolygons.gpml"))
    out = []
    for lo, la in lonlat:
        f = pygplates.Feature(); f.set_geometry(pygplates.PointOnSphere(la, lo))
        part = pygplates.partition_into_plates(polys, rot, [f], properties_to_copy=[pygplates.PartitionProperty.reconstruction_plate_id])
        if not part: out.append(None); continue
        g = part[0].get_geometry(); pla, plo = (rot.get_rotation(AGE, part[0].get_reconstruction_plate_id()) * g).to_lat_lon()
        out.append((plo, pla))
    return out


rot = reconstruct([(lo, la) for _, _, la, lo, *_ in SITES])
G = meta["globe"]; T = G["topo"]
a = np.asarray(Image.open(os.path.join(OUT, T["file"])).convert("RGB")).astype(np.float64)
H = (a[..., 0] * 256 + a[..., 1]) * G["h"]["step"] + G["h"]["off"]
imp = G["impact"]
out = []
for (name, place, la, lo, kind, env, rec, refs), pr in zip(SITES, rot):
    if pr is None: sys.exit(f"{name}: no plate polygon at {la}, {lo}")
    plo, pla = pr
    p1, p2, dl = math.radians(imp["lat"]), math.radians(pla), math.radians(plo - imp["lon"])
    d = 6371.0 * math.acos(min(1, math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)))
    i = min(int((plo + 180) / 360 * T["nx"]), T["nx"] - 1); j = min(int((90 - pla) / 180 * T["ny"]), T["ny"] - 1)
    h = float(H[j, i]); hm = H[max(j - 3, 0):j + 4, max(i - 3, 0):i + 4]
    match = (h > 0) == (env == "land")
    near = bool(((hm > 0) == (env == "land")).any())                  # within ~30 km
    print(f"{name:24s} 66 Ma {plo:8.2f} {pla:7.2f}  {d:7.0f} km  map {h:6.0f} m  {env}: {'ok' if match else 'within 30 km' if near else 'MISMATCH'}")
    out.append({"name": name, "place": place, "lat": la, "lon": lo, "plat": round(pla, 3), "plon": round(plo, 3), "kind": kind, "env": env,
                "km": round(d), "map_h": round(h), "record": rec, "refs": refs})
meta["sites"] = {"age": AGE, "model": MODEL, "list": out}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
