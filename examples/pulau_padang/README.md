# Pulau Padang sub-KHG delineation test

Run from this folder with `repeat-delineation` and `repeat-common` installed (`pip install -e .`).

| file | purpose |
|---|---|
| `inspect_dtm.py` | CRS / nodata / value check of the 2 m source DTM |
| `aggregate.py` | average-resample the 2 m DTM to `data/dtm_10m.tif`, `data/dtm_30m.tif` |
| `harness.py` | run `delineate.delineate` (stand-in KHG = DTM footprint); writes `output/cache_<res><tag>/` |
| `plot.py` | `output/quicklook.png` from one or more cache folders |
| `diagnose_unassigned.py` | classifies cells left unassigned before the nearest-neighbour fill; writes `output/unassigned_30.png` |

    python harness.py 30                          # automatic outlet order
    python harness.py 30 "{\"outlet_order\":3}" _o3
    python harness.py 10
    python plot.py cache_30_o3 cache_30

`output/cache_*/subkhg.gpkg` = sub-KHG polygons; `static/` = stream order and flow accumulation rasters.
