# repeat-delineation

Sub-KHG delineation for RE-PEAT (GitHub org Re-Peat-Project). Takes a ground DTM, rivers and canals, and writes sub-KHG polygons that the neraca and prioritisation modules read. Follows the Kahayan-Sebangau study (2021), automated.

```
dtm_smooth.tif (repeat-dem) ─┐
rivers, canals, KHG polygon ─┴─► repeat-delineation ─► subkhg_units.gpkg (khg_id, unit_id, geometry)
                                                      delineation_report.json
```

File formats: `re-peat/docs/11-data-contracts.md`, section C3.

## Install and run

```bash
pip install -e ../repeat-common
pip install -e .
python -m repeat_delineation --config configs/pulau_padang.yaml          # smoke test, about 1 minute
python -m repeat_delineation --config my.yaml --khg-id KHG01
python -m repeat_delineation --config my.yaml --sweep 5,20,50,100        # count the sub-KHG each minimum area gives
python -m tests.test_delineate_merge ; python -m tests.test_config
```

Copy [configs/example.yaml](configs/example.yaml) for a real KHG. Output: `<out>/subkhg_units.gpkg`, `<out>/delineation_report.json` (settings, inputs, pieces before and units after merging, areas), and per KHG `<cache_dir>/<khg>/subkhg.gpkg` with `static/stream_order.tif` and `static/flow_accumulation.tif` to check the result. Units are named `<khg_id>_S01`, `_S02`, north to south. Check the polygons, then give them to `repeat-prioritisation` as `units: {mode: vector, path: <out>/subkhg_units.gpkg}`; edit them first if needed, keeping `khg_id` and `unit_id`.

## Inputs

| Input | Config | Role |
|---|---|---|
| KHG polygon | `khg` | Clips the result |
| DTM | `rasters.dtm` (required) | Flow routing |
| Rivers, lakes | `vectors.rivers`, `vectors.lakes` | Always boundaries |
| Canals | `vectors.canals` | Burned into the DTM; primary canals are boundaries |
| Canal class column | `units.subkhg.canal_class_field` | Which canals are outlet channels (default: all) |

## Method

1. Rivers and canals are burned 3 m into the DTM (canals are narrower than a 30 m cell).
2. Pits, depressions and flats are resolved. D8 flow direction, accumulation and Strahler order follow (`pysheds`).
3. Outlet channels = mapped rivers + outlet canals + DEM streams of order >= k. Streams start at 2 km2 of contributing area. k is picked automatically so the KHG gets 6-30 sub-KHG, or set `outlet_order`.
4. Each cell is traced downstream to the first outlet cell. A sub-KHG is the catchment of one channel reach on one bank, so divides and outlet channels are the boundaries.
5. Coast-draining catchments over 50 km2 are cut into compact pieces.
6. Pieces under `min_area_km2` are merged, smallest first, into the neighbour they share the longest divide border with (`merge: adjacent`; an outlet-channel border only if there is no divide neighbour). `min_area_km2: auto` raises the minimum 10% at a time from 5 km2 until at most `target_count[1]` units remain. `merge: nearest` is the earlier behaviour (small pieces dropped, cells go to the nearest piece). Channel cells go to the nearest unit.

Without flow directions (`use_dem_streams: false`, or the hydrology fails), units are land enclosed by mapped channels. The count target is then ignored.

Parameters (`units.subkhg.*`, defaults in `SUBKHG_DEFAULTS`): `burn_depth_m`, `stream_threshold_km2`, `outlet_order`, `target_count`, `min_area_km2` (number or `auto`), `merge`, `coastal_split_km2`, `canal_class_field`, `outlet_canal_classes`, `use_dem_streams`.

Field checks in the 2021 study found sub-KHG boundaries are the wide outlet canals, so mapped rivers and primary canals are always boundaries; DEM-derived streams only add boundaries where no mapped channel exists.

## Examples

[examples/pulau_padang/](examples/pulau_padang/README.md): scripts that run the delineation on the Pulau Padang LiDAR DTM at 10 and 30 m and plot it. Only the 30 m DTM (5 MB) is included; `aggregate.py` rebuilds the 10 m one from the 2 m source (set `PULAU_PADANG_DTM_2M`; sites other than `padang` read `REPEAT_DEM_RUNS`).(set `PULAU_PADANG_DTM_2M`; sites other than `padang` read `REPEAT_DEM_RUNS`).

## Notes for the next steps

- Config key `units.subkhg` is kept from the old combined program so existing site configs work; it can become `delineation:` when the site configs move to the hub.
- Terrain limits the result: see the accuracy notes in `repeat-dem`. Delineation on a smoothed DEMNAS gives regional divides, not decimetre-scale ones.
