"""Draw the sub-KHG delineation as a draw.io diagram: inputs, processing steps, outputs, the modules that read them.

    python docs/delineation_diagram.py [KHG folder name]     # default 16.02.01_sugihan_lumpur; writes docs/delineation_flow.drawio

Steps follow src/repeat_delineation/delineate.py; parameter defaults are read from SUBKHG_DEFAULTS and DEFAULTS, the
example numbers from the KHG's last run (_khg/<khg>/02_delineation/delineation_report.json) and stability test
(02_delineation/stability/summary.csv, python -m repeat_delineation.stability). Rerun after changing the method.
Layout: one row grid; inputs enter their step from the left on the step's own row, outputs leave to the right, the
steps run top to bottom, and every arrow has fixed exit and entry points so draw.io does not re-route them.
"""
import json
import os
import sys
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
WS = REPO.parent
sys.path[:0] = [str(REPO / "src"), str(WS / "repeat-common" / "src")]
from repeat_delineation.config import DEFAULTS  # noqa: E402
from repeat_delineation.delineate import SUBKHG_DEFAULTS as P  # noqa: E402

KHG = sys.argv[1] if len(sys.argv) > 1 else "16.02.01_sugihan_lumpur"
KDIR = Path(os.environ.get("REPEAT_KHG") or WS / "_khg") / KHG
OUT = HERE / "delineation_flow.drawio"

STY = {
    "lane": "rounded=1;whiteSpace=wrap;html=1;fillColor={f};strokeColor={s};verticalAlign=top;fontStyle=1;fontSize=14;arcSize=2;",
    "in": "shape=cylinder3;whiteSpace=wrap;html=1;boundedLbl=1;size=7;fillColor=#dae8fc;strokeColor=#6c8ebf;fontSize=11;align=left;spacingLeft=6;",
    "step": "rounded=1;whiteSpace=wrap;html=1;fillColor=#e1d5e7;strokeColor=#9673a6;fontSize=11;align=left;spacingLeft=8;verticalAlign=middle;",
    "alt": "rounded=1;whiteSpace=wrap;html=1;fillColor=#f5f5f5;strokeColor=#999999;dashed=1;fontSize=11;align=left;spacingLeft=8;",
    "out": "shape=cylinder3;whiteSpace=wrap;html=1;boundedLbl=1;size=7;fillColor=#ffe6cc;strokeColor=#d79b00;fontSize=11;align=left;spacingLeft=6;",
    "use": "rounded=1;whiteSpace=wrap;html=1;fillColor=#d5e8d4;strokeColor=#82b366;fontSize=11;align=left;spacingLeft=6;",
    "note": "shape=note;whiteSpace=wrap;html=1;size=12;fillColor=#ffffcc;strokeColor=#b3b300;fontSize=10;align=left;spacingLeft=6;verticalAlign=top;spacingTop=4;",
    "title": "text;html=1;fontSize=20;fontStyle=1;align=left;verticalAlign=middle;",
    "sub": "text;html=1;fontSize=11;align=left;verticalAlign=top;whiteSpace=wrap;fontColor=#52514e;",
}
EDGE = "edgeStyle=orthogonalEdgeStyle;rounded=1;html=1;endArrow=block;endFill=1;strokeColor=#52514e;fontSize=10;labelBackgroundColor=#ffffff;"
BR = "<br>"
X = {"A": 20, "B": 330, "C": 900, "D": 1220}
W = {"A": 270, "B": 520, "C": 280, "D": 260}
Y0, RH, BH = 130, 96, 74
cells, boxes = [], {}


def ry(r):
    return Y0 + r * RH


def box(name, kind, col, r0, label, r1=None, x=None, w=None, h=None, y=None):
    r1 = r0 if r1 is None else r1
    x = X[col] if x is None else x
    w = W[col] if w is None else w
    y = ry(r0) if y is None else y
    h = (ry(r1) + BH - y) if h is None else h
    boxes[name] = (x, y, w, h)
    cells.append(f'<mxCell id="{name}" value="{escape(label)}" style="{STY[kind]}" vertex="1" parent="1">'
                 f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')


def edge(a, b, ya=None, yb=None, side=("R", "L"), label="", dashed=False):
    """Arrow from box a to box b; ya/yb are absolute y of the exit/entry points (default: box middle)."""
    def pt(n, yy, s):
        x, y, w, h = boxes[n]
        if s in "LR":
            return (0 if s == "L" else 1), ((yy if yy is not None else y + h / 2) - y) / h
        return 0.5, (0 if s == "T" else 1)
    ex, ey = pt(a, ya, side[0])
    nx, ny = pt(b, yb, side[1])
    st = EDGE + f"exitX={ex};exitY={ey:.4f};exitDx=0;exitDy=0;entryX={nx};entryY={ny:.4f};entryDx=0;entryDy=0;" + ("dashed=1;" if dashed else "")
    cells.append(f'<mxCell id="e{len(cells)}" value="{escape(label)}" style="{st}" edge="1" parent="1" source="{a}" target="{b}">'
                 f'<mxGeometry relative="1" as="geometry"/></mxCell>')


def mid(r):
    return ry(r) + BH / 2


def main():
    rep = json.loads((KDIR / "02_delineation" / "delineation_report.json").read_text(encoding="utf-8"))[0]
    st_f = KDIR / "02_delineation" / "stability" / "summary.csv"
    st = pd.read_csv(st_f) if st_f.exists() else None
    res, ar = DEFAULTS["grid_res_m"], rep["area_km2"]

    cells.append(f'<mxCell id="t" value="Sub-KHG delineation (repeat-delineation)" style="{STY["title"]}" vertex="1" parent="1">'
                 f'<mxGeometry x="20" y="10" width="900" height="34" as="geometry"/></mxCell>')
    sub = ("python -m repeat_delineation --config delineate.yaml. Steps as in src/repeat_delineation/delineate.py; values in "
           "brackets are the defaults (units.subkhg in the config overrides them). Blue: inputs. Purple: steps, top to bottom. "
           "Orange: files written. Green: modules that read them.")
    cells.append(f'<mxCell id="s" value="{escape(sub)}" style="{STY["sub"]}" vertex="1" parent="1">'
                 f'<mxGeometry x="20" y="46" width="1460" height="40" as="geometry"/></mxCell>')
    for col, title, f, s in (("A", "Inputs", "#f4f8fd", "#6c8ebf"), ("B", "Processing", "#faf6fc", "#9673a6"),
                             ("C", "Outputs", "#fffaf3", "#d79b00"), ("D", "Read by", "#f6fbf5", "#82b366")):
        cells.append(f'<mxCell id="lane{col}" value="{title}" style="{STY["lane"].format(f=f, s=s)}" vertex="1" parent="1">'
                     f'<mxGeometry x="{X[col] - 10}" y="{Y0 - 34}" width="{W[col] + 20}" height="{ry(12) + BH + 20 - Y0 + 34}" as="geometry"/></mxCell>')

    # inputs
    box("khg", "in", "A", 0, f"<b>KHG polygon</b>{BR}khg ({KHG}/00_inputs/khg.gpkg)")
    box("dtm", "in", "A", 1, f"<b>Ground DTM</b>{BR}rasters.dtm: dtm_smooth.tif from repeat-dem, or LiDAR. A DSM gives wrong catchments")
    box("riv", "in", "A", 2, f"<b>Rivers, lakes</b>{BR}vectors.rivers, vectors.lakes: always boundaries")
    box("can", "in", "A", 3, f"<b>Canals</b>{BR}vectors.canals, class column canal_class_field (here KET_KANAL: Primer, Sekunder, Tersier)")

    # steps
    S = [
        ("g", 0, 0, f"<b>1. Grid</b>{BR}UTM grid at grid_res_m ({res} m) over the KHG plus 1 km"),
        ("w", 1, 1, f"<b>2. DTM onto the grid</b>{BR}bilinear warp; KHG mask; cells without DTM leave the KHG (warning if more than 1%)"),
        ("ch", 2, 3, f"<b>3. Channel masks</b>{BR}rivers + lakes; all canals; <b>outlet canals</b> = canals whose class is in "
                     f"outlet_canal_classes ({', '.join(P['outlet_canal_classes'])}; site config: Kanal Primer). "
                     f"Outlet set = rivers + outlet canals"),
        ("b", 4, 4, f"<b>4. Burn</b>{BR}DTM lowered by burn_depth_m ({P['burn_depth_m']:g} m) under rivers and all canals "
                    f"(canals are narrower than a cell)"),
        ("d8", 5, 5, f"<b>5. Condition and route</b> (pysheds){BR}fill pits, fill depressions, resolve flats; D8 flow direction; flow accumulation"),
        ("so", 6, 6, f"<b>6. DEM streams</b>{BR}accumulation ≥ stream_threshold_km2 ({P['stream_threshold_km2']:g} km²); Strahler order"),
        ("os", 7, 7, f"<b>7. Choose the outlet set</b>{BR}first the mapped outlets alone; if the unit count is outside target_count "
                     f"({P['target_count'][0]}-{P['target_count'][1]}), add DEM streams of order ≥ k, k from the highest order down "
                     f"(or outlet_order fixed)"),
        ("pc", 8, 8, f"<b>8. Pieces</b>{BR}outlet network cut into reaches at junctions; every cell traced down D8 to its first outlet cell "
                     f"= the catchment of that reach; split at divides into connected banks. Flow leaving the KHG without an outlet = "
                     f"coastal catchment, cut into pieces of about coastal_split_km2 ({P['coastal_split_km2']:g} km²)"),
        ("mg", 9, 10, f"<b>9. Merge small pieces</b> (merge: {P['merge']}){BR}pieces under min_area_km2 join, smallest first, the "
                      f"neighbour they share the longest <i>divide</i> border with (a border across an outlet channel only if there is "
                      f"no divide neighbour). min_area_km2: auto = start at 5 km², raise 10% at a time until at most "
                      f"{P['target_count'][1]} units remain"),
        ("fn", 11, 12, f"<b>10. Units</b>{BR}channel cells go to the nearest unit; clip to the KHG; polygonise; named "
                       f"&lt;khg_id&gt;_S01, _S02 ... north to south"),
    ]
    for n, r0, r1, lab in S:
        box(n, "step", "B", r0, lab, r1)
    box("fb", "alt", "A", 6, f"<b>Fallback</b>{BR}no flow directions (use_dem_streams: false, or step 5 fails): units = land "
                             f"enclosed by the mapped outlets; target_count ignored", r1=7)

    # outputs and readers
    box("acc", "out", "C", 5, f"<b>flow_accumulation.tif</b>{BR}cache/&lt;khg&gt;/static/")
    box("ord", "out", "C", 6, f"<b>stream_order.tif</b>{BR}cache/&lt;khg&gt;/static/")
    box("un", "out", "C", 11, f"<b>subkhg_units.gpkg</b>{BR}&lt;out&gt;/ (khg_id, unit_id, geometry); also cache/&lt;khg&gt;/subkhg.gpkg")
    box("rp", "out", "C", 12, f"<b>delineation_report.json</b>{BR}settings, inputs, pieces, units, areas")
    box("nr", "use", "D", 11, f"<b>repeat-neraca</b>{BR}water balance per unit (units: in neraca.yaml)")
    box("pr", "use", "D", 12, f"<b>repeat-prioritisation subkhg</b>{BR}units: {{mode: vector, path: subkhg_units.gpkg}}; dashboard")

    # edges: inputs on their rows, the step chain, outputs
    edge("khg", "g")
    edge("dtm", "w")
    edge("riv", "ch", yb=mid(2))
    edge("can", "ch", yb=mid(3))
    chain = ["g", "w", "ch", "b", "d8", "so", "os", "pc", "mg", "fn"]
    for a, b in zip(chain, chain[1:]):
        edge(a, b, side=("B", "T"))
    edge("d8", "acc")
    edge("so", "ord")
    edge("d8", "fb", side=("L", "R"), ya=mid(5) + 18, yb=mid(6), label="fails", dashed=True)
    edge("fn", "un", ya=mid(11), yb=mid(11))
    edge("fn", "rp", ya=mid(12), yb=mid(12))
    edge("un", "nr")
    edge("un", "pr", ya=mid(11) + 20, yb=mid(12), side=("R", "L"))

    # notes below the lanes
    ny = ry(13) + 50
    pieces, units = rep["pieces_before_merge"], rep["n_subkhg"]
    n1 = (f"<b>{KHG}, last run</b>{BR}{rep['inputs'].get('canal_segments', '?')} canal segments "
          f"({', '.join(f'{k} {v}' for k, v in rep['inputs'].get('canal_segments_by_class', {}).items())}); outlets: "
          f"{', '.join(rep['params']['outlet_canal_classes'])}{BR}{pieces} pieces before merging → {units} units; minimum area "
          f"{rep['min_area_km2']} km² (auto); unit area {ar['min']} to {ar['max']} km² (median {ar['median']})")
    box("n1", "note", "A", 0, n1, x=20, w=560, y=ny, h=72)
    if st is not None:
        lines = [f"{r.setting}: noise members agree ARI {r.ari_noise_members_mean}; DTMs vs LiDAR {r.ari_dtms_vs_reference_mean}; "
                 f"flat surface vs LiDAR {r.ari_flat_vs_reference}; {r.units_range} units" for r in st.itertuples()]
        n2 = (f"<b>Stability</b> (02_delineation/stability/, python -m repeat_delineation.stability){BR}" + BR.join(lines) +
              f"{BR}ARI: adjusted Rand index of two partitions, 1 = identical, 0 = chance.")
        box("n2", "note", "B", 0, n2, x=600, w=560, y=ny, h=72)
    n3 = (f"<b>Where the boundaries come from</b>{BR}mapped outlets (rivers, outlet canals) are fixed boundaries; inside them "
          f"the DTM only decides the divides (steps 4-8), and step 9 decides which pieces join. Without outlets inside the KHG, "
          f"small DTM changes move the divides and the merges")
    box("n3", "note", "C", 0, n3, x=1180, w=330, y=ny, h=72)

    xml = ('<mxfile host="app.diagrams.net"><diagram name="sub-KHG delineation"><mxGraphModel dx="1600" dy="1200" grid="1" '
           'gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1600" '
           'pageHeight="1700" math="0" shadow="0"><root><mxCell id="0"/><mxCell id="1" parent="0"/>'
           + "".join(cells) + "</root></mxGraphModel></diagram></mxfile>")
    OUT.write_text(xml, encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
