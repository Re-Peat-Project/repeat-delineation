"""Run delineate.delineate on the Pulau Padang DTM.
usage: python harness.py <site> <res_m> ['{"outlet_order": 3}'] [tag]
"""
import sys, logging, time, json
from _common import DATA, OUT
from sites import SITES
import numpy as np, rasterio, geopandas as gpd
from rasterio.features import shapes
from scipy import ndimage
from shapely.geometry import shape
from repeat_delineation import delineate as D
from repeat_common.grid import make_grid

logging.basicConfig(level=logging.INFO, format="%(message)s")
site = sys.argv[1]; res = int(sys.argv[2])
extra = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
tag = sys.argv[4] if len(sys.argv) > 4 else ""
S = SITES[site]; dtm = S["dtm"](res)

if S["boundary"]:
    geom_ll = gpd.read_file(S["boundary"]).to_crs(4326).union_all()
else:   # stand-in KHG = valid-data footprint of the DTM
    with rasterio.open(DATA / "dtm_30m.tif") as s:
        m = np.isfinite(s.read(1)); t = s.transform; crs = s.crs
    m = ndimage.binary_fill_holes(ndimage.binary_closing(m, iterations=2))
    lab, n = ndimage.label(m); sz = ndimage.sum(m, lab, range(1, n + 1))
    keep = np.isin(lab, 1 + np.flatnonzero(sz > 0.05 * sz.max()))
    g = gpd.GeoSeries([shape(x) for x, _ in shapes(keep.astype("uint8"), mask=keep, transform=t)],
                      crs=crs).union_all().simplify(30)
    geom_ll = gpd.GeoSeries([g], crs=crs).to_crs(4326).iloc[0]
grid = make_grid(geom_ll, res, 1000)
cache = OUT / f"{site}_{res}{tag}"; (cache / "static").mkdir(parents=True, exist_ok=True)
cfg = {"units": {"subkhg": extra}, "rasters": {"dtm": dtm}, "vectors": {}}

orig = D._units_from_outlets          # report cells that reach no outlet (filled by nearest neighbour later)
def wrap(fdir, outlet, khg, min_px, split_px=None):
    lab, n = orig(fdir, outlet, khg, min_px, split_px)
    print(f"  unassigned before fill: {((lab == 0) & khg).sum() / khg.sum():.1%} (n={n})")
    return lab, n
D._units_from_outlets = wrap

t0 = time.time()
gdf, info = D.delineate(cfg, site, geom_ll, grid, cache)
print("time", round(time.time() - t0), "s;", info)
gp = gdf.to_crs(grid.crs); print((gp.area / 1e6).describe().round(1).to_string())
