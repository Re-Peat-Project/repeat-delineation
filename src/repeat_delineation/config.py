"""Config of a delineation run: defaults, merge with the user's YAML, ${NAME} expansion."""
from __future__ import annotations

import copy

from repeat_common.config import deep_update, expand_env

DEFAULTS = {
    "khg": None,                      # vector file with KHG boundaries
    "khg_id_field": "khg_id",         # created from row number if missing
    "khg_filter": [],                 # run only these KHG ids (empty = all)
    "grid_res_m": 30,
    "cache_dir": "cache",             # per KHG: <cache_dir>/<khg_id>/subkhg.gpkg, static/stream_order.tif, static/flow_accumulation.tif
    "out": "out_delineation",
    "rasters": {"dtm": None},         # dtm_smooth.tif from repeat-dem (or lidar); a DSM gives wrong catchments
    "vectors": {"canals": None, "rivers": None, "lakes": None},
    # method settings: see SUBKHG_DEFAULTS in delineate.py (kept under units.subkhg so existing configs keep working)
    "units": {"subkhg": {}},
}


def load_config(path=None, khg=None, out=None) -> dict:
    cfg = copy.deepcopy(DEFAULTS)
    if path:
        import yaml
        with open(path, encoding="utf-8") as f:
            deep_update(cfg, expand_env(yaml.safe_load(f) or {}))
    if khg:
        cfg["khg"] = khg
    if out:
        cfg["out"] = out
    return cfg
