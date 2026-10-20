"""Shared paths for the Pulau Padang test scripts."""
import os
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "data"
OUT = HERE / "output"
# the 2 m LiDAR DTM, only needed by aggregate.py and inspect_dtm.py (not in the repo): set PULAU_PADANG_DTM_2M
SRC_DTM = Path(os.environ.get("PULAU_PADANG_DTM_2M", "DTM_2mFloat.tif"))
