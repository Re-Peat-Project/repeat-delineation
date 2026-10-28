"""Test sites. dtm may depend on the resolution; boundary None = stand-in (valid-data footprint of the DTM)."""
import os

from _common import DATA

# folder of repeat-dem runs (repeat-dem/runs); only the sites other than "padang" need it
RUNS = os.environ.get("REPEAT_DEM_RUNS", "")

SITES = {
    "padang":    {"dtm": lambda res: str(DATA / f"dtm_{res}m.tif"), "boundary": None, "res": (10, 30)},
    "ks":        {"dtm": lambda res: RUNS + r"\ks_full_30m\dtm_smooth.tif",   "boundary": RUNS + r"\aoi_ks.gpkg", "res": (30,)},
    "mendahara": {"dtm": lambda res: RUNS + r"\mendahara_30m\dtm_smooth.tif", "boundary": RUNS + r"\aoi_mendahara.gpkg", "res": (30,)},
    "bengkalis": {"dtm": lambda res: RUNS + r"\bengkalis_30m\dtm_smooth.tif", "boundary": RUNS + r"\aoi_bengkalis.gpkg", "res": (30,)},
}
