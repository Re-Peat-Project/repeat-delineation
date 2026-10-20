"""Quick-look map: hillshade + sub-KHG boundaries + streams. usage: python plot.py <site_res[tag]> [more...]  e.g. ks_30 padang_30"""
import sys
from _common import DATA, OUT
from sites import SITES
import numpy as np, rasterio, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

names = sys.argv[1:]
fig, axs = plt.subplots(1, len(names), figsize=(8 * len(names), 11), squeeze=False)
for ax, name in zip(axs[0], names):
    site, res = name.split("_")[:2]
    with rasterio.open(SITES[site]["dtm"](int(res))) as s:
        z = s.read(1, out_shape=(min(s.height, 1800), min(s.width, 1800))); b = s.bounds
        z = np.where(z == s.nodata, np.nan, z) if s.nodata is not None else z
    zz = np.nan_to_num(z, nan=np.nanmin(z))
    ext = (b.left, b.right, b.bottom, b.top)
    ax.imshow(LightSource(315, 45).hillshade(zz, vert_exag=20), cmap="gray", extent=ext)
    so = rasterio.open(OUT / name / "static" / "stream_order.tif")
    o = so.read(1); bb = so.bounds
    ax.imshow(np.where(o >= 3, 1.0, np.nan), cmap="winter", extent=(bb.left, bb.right, bb.bottom, bb.top), alpha=.8)
    g = gpd.read_file(OUT / name / "subkhg.gpkg").to_crs(so.crs)
    g.boundary.plot(ax=ax, color="red", lw=1.2)
    ax.set_xlim(bb.left, bb.right); ax.set_ylim(bb.bottom, bb.top)
    ax.set_title(f"{name}: {len(g)} sub-KHG"); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(OUT / ("quicklook_" + "_".join(names) + ".png"), dpi=90)
