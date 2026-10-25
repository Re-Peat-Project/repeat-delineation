"""Shared paths for the Pulau Padang test scripts."""
import os
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "data"
OUT = HERE / "output"
# the 2 m LiDAR DTM, only needed by aggregate.py and inspect_dtm.py (not in the repo): set PULAU_PADANG_DTM_2M,
# else $REPEAT_DATA (default: the workspace _data) / lidar_pulau_padang_2016
_DATA_ROOT = Path(os.environ.get("REPEAT_DATA") or HERE.parents[2] / "_data")
SRC_DTM = Path(os.environ.get("PULAU_PADANG_DTM_2M", _DATA_ROOT / "lidar_pulau_padang_2016" / "DTM_2mFloat.tif"))
