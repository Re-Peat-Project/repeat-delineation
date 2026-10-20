"""Config loading: defaults are filled in, the user's values win, ${NAME} is expanded.

    python -m tests.test_config
"""
import os
import tempfile
from pathlib import Path

from repeat_delineation.config import DEFAULTS, load_config


def test_defaults_and_overrides():
    os.environ["REPEAT_TEST_ROOT"] = "D:/data"
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "c.yaml"
        p.write_text("khg: k.gpkg\nrasters: {dtm: '${REPEAT_TEST_ROOT}/dtm.tif'}\nunits: {subkhg: {merge: nearest}}\n", encoding="utf-8")
        cfg = load_config(p, out="o")
    assert cfg["rasters"]["dtm"] == "D:/data/dtm.tif"
    assert cfg["units"]["subkhg"] == {"merge": "nearest"}
    assert cfg["out"] == "o" and cfg["grid_res_m"] == 30 and cfg["vectors"]["canals"] is None
    assert DEFAULTS["rasters"]["dtm"] is None, "defaults must not be modified"
    print("ok test_defaults_and_overrides")


if __name__ == "__main__":
    test_defaults_and_overrides()
