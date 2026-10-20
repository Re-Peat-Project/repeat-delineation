#!/usr/bin/env python3
"""Sub-KHG delineation only: no cloud fetching, no indicators. Needs rasters.dtm; rivers and canals optional.

  python -m repeat_delineation --config delineate.yaml
  python -m repeat_delineation --config delineate.yaml --khg-id KHG01
  python -m repeat_delineation --config delineate.yaml --sweep 5,20,35,50,70,100

--sweep also counts the sub-KHG each minimum area gives, with the adjacent merge and the earlier nearest-unit merge.

Writes <out>/subkhg_units.gpkg (khg_id, unit_id, geometry), <out>/delineation_report.json (the settings used, the
inputs, pieces before and units after merging, unit areas; read by the prioritisation dashboard) and, per KHG, cache/<khg>/subkhg.gpkg with
static/stream_order.tif and static/flow_accumulation.tif. Check the polygons, then score them with
units: {mode: vector, path: <out>/subkhg_units.gpkg} in a sub-KHG run of repeat-prioritisation.
"""
import argparse
import json
import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from repeat_common import grid as G

from .config import load_config
from .delineate import SUBKHG_DEFAULTS, delineate, sweep


def _report(cfg, kid, units, info) -> dict:
    """The settings and counts of one delineation, for the report and the explanatory page."""
    from repeat_common.grid import auto_utm
    p = {**SUBKHG_DEFAULTS, **(cfg["units"].get("subkhg") or {})}
    c = units.geometry.union_all().centroid
    area = units.to_crs(auto_utm(c.x, c.y)).area / 1e6
    inputs = {}
    if cfg["vectors"].get("canals"):
        cg = gpd.read_file(cfg["vectors"]["canals"])
        inputs["canal_segments"] = int(len(cg))
        f = p["canal_class_field"]
        if f and f in cg.columns:
            inputs["canal_segments_by_class"] = {str(k): int(v) for k, v in cg[f].astype(str).value_counts().items()}
            inputs["outlet_canal_segments"] = int(cg[f].astype(str).str.lower().isin([x.lower() for x in p["outlet_canal_classes"]]).sum())
    for key in ("rivers", "lakes"):
        if cfg["vectors"].get(key):
            inputs[f"{key}_segments"] = int(len(gpd.read_file(cfg["vectors"][key])))
    return {"khg_id": kid, "params": {k: v for k, v in p.items()}, **info, "inputs": inputs,
            "area_km2": {"min": round(float(area.min()), 1), "median": round(float(area.median()), 1), "max": round(float(area.max()), 1),
                         "total": round(float(area.sum()), 1)}}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config")
    ap.add_argument("--khg", help="KHG boundary file (overrides config)")
    ap.add_argument("--khg-id", action="append", help="only this KHG id (repeatable)")
    ap.add_argument("--out", help="output folder (overrides config)")
    ap.add_argument("--sweep", help="comma-separated minimum areas in km2, e.g. 5,20,50,100: count the sub-KHG each gives")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    cfg = load_config(a.config, a.khg, a.out)
    if not cfg["khg"]:
        raise SystemExit("set 'khg' (KHG boundary file) in the config or with --khg")
    khg = gpd.read_file(cfg["khg"]).to_crs(4326)
    idf = cfg["khg_id_field"]
    if idf not in khg.columns:
        khg[idf] = [f"KHG{i + 1:02d}" for i in range(len(khg))]
    khg[idf] = khg[idf].astype(str)
    ids = a.khg_id or cfg["khg_filter"]
    if ids:
        khg = khg[khg[idf].isin([str(x) for x in ids])]

    frames, reports, sweeps = [], [], []
    for _, r in khg.iterrows():
        kid = r[idf]
        cache = Path(cfg["cache_dir"]) / kid
        cache.mkdir(parents=True, exist_ok=True)
        grid = G.make_grid(r.geometry, float(cfg["grid_res_m"]))
        units, info = delineate(cfg, kid, r.geometry, grid, cache)
        logging.getLogger("repeat_delineation").info("KHG %s: %s", kid, info)
        frames.append(units)
        reports.append(_report(cfg, kid, units, info))
        if a.sweep:
            sweeps.append(sweep(cfg, kid, r.geometry, grid, cache, [float(x) for x in a.sweep.split(",")]))
    if not frames:
        raise SystemExit("no KHG processed")
    out = Path(cfg["out"])
    out.mkdir(parents=True, exist_ok=True)
    allu = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=4326)
    allu.to_file(out / "subkhg_units.gpkg", driver="GPKG")
    (out / "delineation_report.json").write_text(json.dumps(reports, indent=2), encoding="utf-8")
    if sweeps:
        (out / "delineation_sweep.json").write_text(json.dumps(sweeps[0] if len(sweeps) == 1 else sweeps, indent=2), encoding="utf-8")
        for s in sweeps:
            print("minimum area km2:", s["min_area_km2"], "| pieces:", s["pieces"], "| adjacent merge:", s["units_adjacent"], "| nearest merge:", s["units_nearest"])
    print(f"{len(allu)} sub-KHG in {len(frames)} KHG -> {out / 'subkhg_units.gpkg'}")


if __name__ == "__main__":
    main()
