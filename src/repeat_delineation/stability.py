#!/usr/bin/env python3
"""How stable are the sub-KHG units? Delineate one KHG many times with only the terrain or the settings changed.

  python -m repeat_delineation.stability --config ../re-peat/sites/16.02.01_sugihan_lumpur/delineation_stability.yaml

The config names a base delineation config, a set of DTMs (`flat` = a constant surface, so only the mapped channels
shape the units), optional noise members (the reference DTM plus spatially correlated gaussian noise, made on the
delineation grid) and a set of settings (overrides of units.subkhg). Every DTM and noise member runs with every setting.
Writes to <out>/:
  units/<setting>__<member>.gpkg      the units of each run
  members.csv                         per run: units, area range, best-match IoU and adjusted Rand index against the
                                      reference DTM with the same setting
  summary.csv                         per setting: agreement among the noise members and between the DTMs
  boundary_frequency_<setting>.tif    share of runs (DTMs and noise members) in which a cell lies on a unit boundary
  boundary_frequency.png              the same, side by side
The adjusted Rand index compares two partitions cell by cell (1 = identical, 0 = no better than chance), so it is not
thrown off by different unit counts. Best-match IoU: for each reference unit, the IoU of the run's unit that matches it
best, averaged.
"""
from __future__ import annotations

import argparse
import copy
import itertools
import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml
from rasterio.features import rasterize
from scipy import ndimage

from repeat_common import grid as G
from repeat_common.config import deep_update, expand_env

from .config import load_config
from .delineate import delineate

log = logging.getLogger("repeat_delineation.stability")


def label_raster(units, grid):
    u = units.to_crs(grid.crs).reset_index(drop=True)
    return rasterize([(g, i + 1) for i, g in enumerate(u.geometry)], out_shape=grid.shape, transform=grid.transform,
                     fill=0, dtype="int32")


def boundaries(lab):
    b = np.zeros(lab.shape, bool)
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    return b & (lab > 0)


def ari(a, b, mask):
    """Adjusted Rand index of two label rasters over mask."""
    x, y = a[mask].astype(np.int64), b[mask].astype(np.int64)
    _, xi = np.unique(x, return_inverse=True)
    _, yi = np.unique(y, return_inverse=True)
    pair = xi * (yi.max() + 1) + yi
    nij = np.bincount(pair).astype(np.float64)
    ai = np.bincount(xi).astype(np.float64)
    bj = np.bincount(yi).astype(np.float64)
    c2 = lambda v: (v * (v - 1) / 2).sum()
    s_ij, s_a, s_b, n = c2(nij), c2(ai), c2(bj), float(x.size)
    exp = s_a * s_b / (n * (n - 1) / 2)
    mx = (s_a + s_b) / 2
    return float((s_ij - exp) / (mx - exp)) if mx != exp else 1.0


def best_iou(ref, other):
    """Mean over reference units of the IoU of the best-matching unit (label rasters, 0 = outside)."""
    out = []
    for r in np.unique(ref[ref > 0]):
        m = ref == r
        cand, cnt = np.unique(other[m & (other > 0)], return_counts=True)
        best = 0.0
        for c, i in zip(cand, cnt):
            best = max(best, i / float((m | (other == c)).sum()))
        out.append(best)
    return float(np.mean(out))


def noise_field(shape, res, sigma, corr, rng):
    f = ndimage.gaussian_filter(rng.standard_normal(shape), corr / res / 2.0)
    return (f / f.std() * sigma).astype("float32")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    cpath = Path(a.config).resolve()
    sc = expand_env(yaml.safe_load(open(cpath, encoding="utf-8")))
    base = load_config(str(cpath.parent / sc["base"]))
    out = Path(sc["out"])
    (out / "units").mkdir(parents=True, exist_ok=True)
    khg = gpd.read_file(base["khg"]).to_crs(4326)
    kid = str(khg[base["khg_id_field"]].iloc[0]) if base["khg_id_field"] in khg else "KHG01"
    geom = khg.union_all()
    grid = G.make_grid(geom, float(base["grid_res_m"]))
    inside = rasterize([(gpd.GeoSeries([geom], crs=4326).to_crs(grid.crs).iloc[0], 1)], out_shape=grid.shape,
                       transform=grid.transform, fill=0, dtype="uint8").astype(bool)

    # members: DTM files, the flat control, noise members of the reference (written once on the delineation grid)
    members = {}
    ref_name = sc["reference"]
    for name, path in sc["dtms"].items():
        if path == "flat":
            z = G.warp_to_grid(sc["dtms"][ref_name], grid)
            path = out / "dtm" / "flat.tif"
            G.write(path, np.where(np.isfinite(z), 5.0, np.nan).astype("float32"), grid)
        members[name] = str(path)
    nz = sc.get("noise") or {}
    if nz.get("members"):
        rng = np.random.default_rng(int(nz.get("seed", 1)))
        z = G.warp_to_grid(sc["dtms"][ref_name], grid)
        for i in range(int(nz["members"])):
            f = out / "dtm" / f"{ref_name}_noise{i + 1}.tif"
            if not f.exists():
                G.write(f, z + noise_field(z.shape, grid.res, float(nz["sigma_m"]), float(nz["corr_length_m"]), rng), grid)
            members[f"{ref_name}_noise{i + 1}"] = str(f)
    noise_names = [m for m in members if "_noise" in m]

    rows, summary, freq = [], [], {}
    for sname, over in sc["settings"].items():
        labs = {}
        for mname, path in members.items():
            gp = out / "units" / f"{sname}__{mname}.gpkg"
            if gp.exists():
                units = gpd.read_file(gp)
            else:
                cfg = copy.deepcopy(base)
                cfg["rasters"]["dtm"] = path
                deep_update(cfg["units"]["subkhg"], copy.deepcopy(over))
                cache = out / "cache" / f"{sname}__{mname}"
                cache.mkdir(parents=True, exist_ok=True)
                units, info = delineate(cfg, kid, geom, grid, cache)
                units.to_file(gp, driver="GPKG")
            labs[mname] = label_raster(units, grid)
            log.info("%s / %s: %d units", sname, mname, len(units))
        ref = labs[ref_name]
        mask = inside & (ref > 0)
        for mname, lab in labs.items():
            area = np.bincount(lab[lab > 0]) * grid.res ** 2 / 1e6
            area = area[area > 0]
            rows.append({"setting": sname, "member": mname, "units": int(len(area)),
                         "area_min_km2": round(float(area.min()), 1), "area_median_km2": round(float(np.median(area)), 1),
                         "area_max_km2": round(float(area.max()), 1),
                         "ari_vs_reference": round(ari(ref, lab, mask & (lab > 0)), 3),
                         "best_iou_vs_reference": round(best_iou(ref, lab), 3)})
        dtm_names = [m for m in members if m not in noise_names and m != ref_name]
        pair = lambda names: [ari(labs[x], labs[y], mask & (labs[x] > 0) & (labs[y] > 0)) for x, y in itertools.combinations(names, 2)]
        nz_pairs = pair(noise_names + [ref_name]) if noise_names else []
        summary.append({"setting": sname,
                        "ari_noise_members_mean": round(float(np.mean(nz_pairs)), 3) if nz_pairs else None,
                        "ari_dtms_vs_reference_mean": round(float(np.mean([r["ari_vs_reference"] for r in rows
                                                                          if r["setting"] == sname and r["member"] in dtm_names])), 3),
                        "ari_flat_vs_reference": next((r["ari_vs_reference"] for r in rows if r["setting"] == sname and r["member"] == "flat"), None),
                        "units_range": f"{min(r['units'] for r in rows if r['setting'] == sname)}-{max(r['units'] for r in rows if r['setting'] == sname)}"})
        f = np.mean([boundaries(l) for l in labs.values()], axis=0).astype("float32")
        f[~inside] = np.nan
        freq[sname] = f
        G.write(out / f"boundary_frequency_{sname}.tif", f, grid)
    pd.DataFrame(rows).to_csv(out / "members.csv", index=False)
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False)
    plot(freq, grid, base, sc, out)
    print(pd.DataFrame(summary).to_string(index=False))
    print(pd.DataFrame(rows).to_string(index=False))


def plot(freq, grid, base, sc, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    canals = gpd.read_file(base["vectors"]["canals"]).to_crs(grid.crs)
    field = base["units"]["subkhg"].get("canal_class_field")
    ext = (grid.transform.c, grid.transform.c + grid.width * grid.res, grid.transform.f - grid.height * grid.res, grid.transform.f)
    fig, axs = plt.subplots(1, len(freq), figsize=(7.5 * len(freq), 8.5), dpi=150)
    axs = np.atleast_1d(axs)
    for ax, (sname, f) in zip(axs, freq.items()):
        outl = sc["settings"][sname].get("outlet_canal_classes") or []
        canals.plot(ax=ax, color="#d9d9d9", linewidth=0.25)                  # all canals, context
        if field:
            canals[canals[field].isin(outl)].plot(ax=ax, color="#9e9e9e", linewidth=0.5)   # outlet canals
        # boundaries are one cell wide: thickened for display only (max over 5 x 5 cells); the GeoTIFF keeps the cell values
        show = ndimage.maximum_filter(np.nan_to_num(f, nan=0), 5)
        show = np.where(np.isfinite(f) & (show >= 0.15), show, np.nan)
        im = ax.imshow(show, extent=ext, cmap="magma_r", vmin=0, vmax=1, interpolation="nearest", zorder=3)
        ax.set_title(f"outlets: {', '.join(outl)}", fontsize=11)
        ax.set_axis_off()
    fig.colorbar(im, ax=axs, shrink=0.5, label="share of runs with a unit boundary here (DTMs and noise members)")
    fig.savefig(out / "boundary_frequency.png", bbox_inches="tight")


if __name__ == "__main__":
    main()
