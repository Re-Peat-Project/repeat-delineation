"""
Automatic sub-KHG delineation, following the approach of the Kahayan-Sebangau study (2021):

  1. burn rivers and canals into the DTM (canals are narrower than a cell)
  2. D8 flow direction and flow accumulation (pits, depressions and flats resolved)
  3. stream network from accumulation, Strahler order
  4. outlet channels = rivers + selected canal classes + streams of order >= outlet_order
     (outlet_order chosen automatically so the number of sub-KHG falls in a target range)
  5. sub-KHG = catchment of each outlet-channel reach (cells traced downstream along D8 to the
     first outlet cell), split into the two banks of the channel; divides and channels are
     therefore the boundaries. Pieces smaller than min_area_km2 are merged, smallest first, into the
     adjacent piece they share the longest border with (merge: "adjacent"). A border across a divide
     is preferred to a border across an outlet channel, so outlet canals stay boundaries where they can.
     merge: "nearest" is the earlier behaviour: small pieces are dropped and their cells go to the
     nearest remaining piece. Channel cells are assigned to the nearest unit.

Field checks in that study found sub-KHG boundaries are the wide outlet canals, so mapped
rivers and primary canals are always treated as boundaries; DEM-derived streams only add
boundaries where no mapped channel exists.
"""
from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import numpy as np
from rasterio.features import rasterize, shapes
from scipy import ndimage

from repeat_common.grid import Grid, vector_to_grid_mask, warp_to_grid, write

log = logging.getLogger("repeat_delineation")

SUBKHG_DEFAULTS = {
    "burn_depth_m": 3.0,              # how deep mapped channels are cut into the DTM
    "stream_threshold_km2": 2.0,      # contributing area that starts a DEM-derived stream
    "outlet_order": None,             # None = automatic; or an integer Strahler order
    "target_count": [6, 30],          # wanted number of sub-KHG per KHG (automatic mode)
    "min_area_km2": 5.0,              # smaller pieces are merged into a neighbour; "auto" = the smallest value that gives at most target_count[1] units
    "merge": "adjacent",              # "adjacent": into the neighbour with the longest shared border; "flow": prefer the divide
                                      # neighbour that drains to the same channel line; "nearest": earlier behaviour
    "coastal_split_km2": 50.0,        # coast-draining catchments (no outlet channel) larger than this are
                                      # cut into compact pieces of about this size; None = no split
    "canal_class_field": None,        # e.g. "kelas"; None = every canal is an outlet channel
    "outlet_canal_classes": ["primer", "primary"],
    "use_dem_streams": True,
}


D8 = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}


def _terminal(fdir: np.ndarray, outlet: np.ndarray) -> np.ndarray:
    """Flat index of the first outlet cell reached downstream from every cell (pointer jumping)."""
    h, w = fdir.shape
    idx = np.arange(h * w).reshape(h, w)
    nxt = idx.copy()
    for code, (dr, dc) in D8.items():
        m = fdir == code
        r, c = np.nonzero(m)
        r2, c2 = r + dr, c + dc
        ok = (r2 >= 0) & (r2 < h) & (c2 >= 0) & (c2 < w)
        nxt[r[ok], c[ok]] = idx[r2[ok], c2[ok]]
    nxt = nxt.ravel()
    nxt[outlet.ravel()] = np.flatnonzero(outlet.ravel())      # outlet cells point to themselves
    for _ in range(40):
        n2 = nxt[nxt]
        if np.array_equal(n2, nxt):
            break
        nxt = n2
    return nxt


def _segments(outlet: np.ndarray) -> np.ndarray:
    """Split the outlet channel network into reaches between junctions."""
    k = np.ones((3, 3))
    nb = ndimage.convolve(outlet.astype(int), k, mode="constant") - 1
    junction = outlet & (nb >= 3)
    seg, n = ndimage.label(outlet & ~junction, structure=np.ones((3, 3)))
    if junction.any() and n:
        idx = ndimage.distance_transform_edt(seg == 0, return_distances=False, return_indices=True)
        seg = np.where(junction, seg[idx[0], idx[1]], seg)
    return seg


def _split_large(lab, coastal, max_px):
    """Cut coastal components larger than max_px into compact pieces (k-means on cell coordinates)."""
    from scipy.cluster.vq import kmeans2
    from scipy.spatial import cKDTree
    nxt = int(lab.max()) + 1
    ids, sizes = np.unique(lab[coastal & (lab > 0)], return_counts=True)
    for i, sz in zip(ids, sizes):
        if sz <= max_px:
            continue
        r, c = np.nonzero(lab == i)
        k = int(np.ceil(sz / max_px))
        pts = np.column_stack([r, c]).astype(float)
        sub = pts[np.random.default_rng(0).choice(len(pts), min(len(pts), 20000), replace=False)]
        cent, _ = kmeans2(sub, k, minit="++", seed=0)
        part = cKDTree(cent).query(pts)[1]
        lab[r, c] = np.where(part == 0, i, nxt + part)
        nxt += k
    return lab


def _components(fdir, outlet, khg, split_px=None, with_reach=False):
    """Catchment of each channel reach, split into the two banks, as connected pieces (small ones still in).
    with_reach: also return, per piece, the reach it drains to and the reach links (see _reach_links)."""
    seg = _segments(outlet)
    term = _terminal(fdir, outlet)
    catch = seg.ravel()[term].reshape(seg.shape)
    # flow that leaves the DEM (coast, nodata edge, pits) without meeting an outlet channel is a
    # coastal catchment of its own; split into connected stretches below instead of being dropped
    drains_out = ~outlet.ravel()[term].reshape(seg.shape)
    catch[drains_out] = seg.max() + 1
    catch[~khg] = 0
    land = khg & ~outlet & (catch > 0)
    # barriers between different catchments so components never cross a divide
    barrier = np.zeros_like(land)
    barrier[:, :-1] |= catch[:, :-1] != catch[:, 1:]
    barrier[:-1, :] |= catch[:-1, :] != catch[1:, :]
    lab, n = ndimage.label(land & ~barrier, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    if n and split_px:
        lab = _split_large(lab, catch == seg.max() + 1, split_px)
        n = int(lab.max())
    if not with_reach:
        return lab, int(n)
    reach = np.zeros(int(n) + 1, np.int64)
    reach[lab[lab > 0]] = catch[lab > 0]                 # all cells of a piece drain to one reach (or the coast)
    return lab, int(n), (reach, _reach_links(outlet, seg, int(seg.max()) + 1))


def _reach_links(outlet, seg, coast):
    """Reaches that meet at a junction, as {reach: set of linked reaches, itself included}. The coast is its own link."""
    k = np.ones((3, 3))
    nb = ndimage.convolve(outlet.astype(int), k, mode="constant") - 1
    jl, nj = ndimage.label(outlet & (nb >= 3), structure=np.ones((3, 3)))
    near = ndimage.grey_dilation(jl, size=(3, 3))
    m = (seg > 0) & (near > 0)
    links = {int(r): {int(r)} for r in np.unique(seg[seg > 0])}
    links[coast] = {coast}
    for j in np.unique(near[m]):
        rs = set(np.unique(seg[m & (near == j)]).astype(int).tolist())
        for r in rs:
            links[r] |= rs
    return links


def _fill_labels(lab, khg):
    """Give every cell of the KHG the label of the nearest piece (channel cells and divide lines get filled)."""
    idx = ndimage.distance_transform_edt(lab == 0, return_distances=False, return_indices=True)
    return np.where(khg, lab[idx[0], idx[1]], 0)


def _contacts(filled, chan):
    """Border length (in cell edges) between every pair of adjacent pieces, split in two kinds:
    across a divide (both cells are land) and across an outlet channel (one cell is a channel cell)."""
    n = int(filled.max())
    div, chn = {}, {}
    for sa, sb in (((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
                   ((slice(None, -1), slice(None)), (slice(1, None), slice(None)))):
        a, b = filled[sa], filled[sb]
        m = (a != b) & (a > 0) & (b > 0)
        ch = (chan[sa] | chan[sb])[m]
        key = np.minimum(a[m], b[m]).astype(np.int64) * (n + 1) + np.maximum(a[m], b[m]).astype(np.int64)
        for flag, d in ((False, div), (True, chn)):
            u, c = np.unique(key[ch == flag], return_counts=True)
            for k, v in zip(u.tolist(), c.tolist()):
                pair = (k // (n + 1), k % (n + 1))
                d[pair] = d.get(pair, 0) + v
    return n, div, chn


def _merge_small(sizes, contacts, min_px, flow=None):
    """Merge pieces under min_px, smallest first, into the neighbour with the longest shared border across a divide
    (else across an outlet channel). Returns an array mapping old labels to new consecutive labels (0 stays 0).
    A piece with no neighbour is left alone. flow = (reach per piece, reach links): then a divide neighbour that drains
    to the same channel line (its reaches meet the piece's reaches at a junction) is preferred over a longer border."""
    import heapq
    n, div, chn = contacts
    size = {i: int(sizes[i]) for i in range(1, n + 1) if sizes[i] > 0}
    nbr = {i: {} for i in size}
    for (a, b), c in div.items():
        for x, y in ((a, b), (b, a)):
            nbr[x].setdefault(y, [0, 0])[0] += c
    for (a, b), c in chn.items():
        for x, y in ((a, b), (b, a)):
            nbr[x].setdefault(y, [0, 0])[1] += c
    parent = list(range(n + 1))
    original = list(size)
    if flow is not None:
        reach_of, links = flow
        reaches = {i: {int(reach_of[i])} for i in size}
        linked = lambda i, j: any(links.get(r, {r}) & reaches[j] for r in reaches[i])                                # every piece that exists, including those merged away below
    heap = [(v, i) for i, v in size.items() if v < min_px]
    heapq.heapify(heap)
    while heap:
        sz, i = heapq.heappop(heap)
        if i not in size or size[i] != sz or sz >= min_px or not nbr[i]:
            continue
        cand = nbr[i]
        if flow is None:
            j = max(cand, key=lambda k: (cand[k][0] > 0, cand[k][0] if cand[k][0] > 0 else cand[k][1]))
        else:
            j = max(cand, key=lambda k: (cand[k][0] > 0 and linked(i, k), cand[k][0] > 0, cand[k][0] if cand[k][0] > 0 else cand[k][1]))
            reaches[j] |= reaches.pop(i)
        size[j] += size.pop(i)
        for k, (d, c) in nbr.pop(i).items():
            if k == j:
                continue
            e = nbr[j].setdefault(k, [0, 0])
            e[0] += d
            e[1] += c
            nbr[k].pop(i, None)
            f = nbr[k].setdefault(j, [0, 0])
            f[0] += d
            f[1] += c
        nbr[j].pop(i, None)
        parent[i] = j
        if size[j] < min_px:
            heapq.heappush(heap, (size[j], j))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    roots = {i: root(i) for i in original}
    new = {r: k + 1 for k, r in enumerate(sorted(set(roots.values())))}
    out = np.zeros(n + 1, "int32")
    for i, r in roots.items():
        out[i] = new[r]
    return out


def _units_from_outlets(fdir, outlet, khg, min_px, split_px=None, merge="adjacent"):
    """Catchment of each channel reach, split into the two banks; small pieces merged.
    Returns (labels, number of units, number of pieces before merging). Cells not in a piece stay 0."""
    lab, n, flow = _components(fdir, outlet, khg, split_px, with_reach=True)
    if n == 0:
        return lab, 0, 0
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    pieces = int((sizes[1:] > 0).sum())
    if merge in ("adjacent", "flow"):
        mapping = _merge_small(sizes, _contacts(_fill_labels(lab, khg), outlet), min_px, flow if merge == "flow" else None)
        lab = mapping[lab]
        return lab, int(len(np.unique(mapping[1:][sizes[1:] > 0]))), pieces
    keep = np.zeros(n + 1, bool)
    keep[1:] = sizes[1:] >= min_px
    lab = np.where(keep[lab] & (lab > 0), lab, 0)
    return lab, int(keep.sum()), pieces


def _units_auto(fdir, outlet, khg, split_px, hi, res, merge="adjacent"):
    """Like _units_from_outlets with the adjacent (or flow) merge, but the minimum area is the smallest value (stepping
    10% from 5 km2) at which at most `hi` units remain. Returns (labels, units, pieces, minimum area in km2)."""
    lab, n, flow = _components(fdir, outlet, khg, split_px, with_reach=True)
    if n == 0:
        return lab, 0, 0, 0.0
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    pieces = int((sizes[1:] > 0).sum())
    contacts = _contacts(_fill_labels(lab, khg), outlet)
    a = 5.0
    while True:
        mapping = _merge_small(sizes, contacts, int(a * 1e6 / res ** 2), flow if merge == "flow" else None)
        cnt = int(len(np.unique(mapping[1:][sizes[1:] > 0])))
        if cnt <= hi or a > 20000:
            return mapping[lab], cnt, pieces, round(a, 1)
        a *= 1.1


def _prepare(cfg: dict, khg_id: str, geom_ll, grid: Grid, cache: Path) -> dict:
    """Everything before the units: KHG mask, channel masks, conditioned DTM, D8 flow directions and stream order."""
    import geopandas as gpd
    p = {**SUBKHG_DEFAULTS, **(cfg["units"].get("subkhg") or {})}
    dtm_path = cfg["rasters"].get("dtm")
    if not dtm_path:
        raise ValueError("units.mode = subkhg needs rasters.dtm (e.g. dtm_smooth.tif from demnas_dtm)")
    z = warp_to_grid(dtm_path, grid)
    khg = rasterize([(gpd.GeoSeries([geom_ll], crs=4326).to_crs(grid.crs).iloc[0], 1)],
                    out_shape=grid.shape, transform=grid.transform, fill=0, dtype="uint8").astype(bool)
    nodata_share = float((khg & ~np.isfinite(z)).sum() / max(khg.sum(), 1))
    if nodata_share > 0.01:
        log.warning("KHG %s: %.0f%% of the KHG has no DTM; units are limited to DTM coverage", khg_id, 100 * nodata_share)
    khg &= np.isfinite(z)

    # mapped channels
    rivers = np.zeros(grid.shape, bool)
    for key in ("rivers", "lakes"):
        if cfg["vectors"].get(key):
            rivers |= vector_to_grid_mask(cfg["vectors"][key], grid)
    canals_all = np.zeros(grid.shape, bool)
    canals_out = np.zeros(grid.shape, bool)
    if cfg["vectors"].get("canals"):
        cg = gpd.read_file(cfg["vectors"]["canals"])
        canals_all = vector_to_grid_mask(cg, grid)
        f = p["canal_class_field"]
        if f and f in cg.columns:
            sel = cg[cg[f].astype(str).str.lower().isin([c.lower() for c in p["outlet_canal_classes"]])]
            canals_out = vector_to_grid_mask(sel, grid) if len(sel) else canals_out
        else:
            canals_out = canals_all

    # hydrological conditioning and D8
    zb = np.where(np.isfinite(z), z, np.nan) - p["burn_depth_m"] * (rivers | canals_all)
    order = np.zeros(grid.shape, "int32")
    fdir = None
    if p["use_dem_streams"]:
        try:
            if not hasattr(np, "in1d"):          # pysheds still calls np.in1d, removed in NumPy 2
                np.in1d = np.isin
            from pysheds.grid import Grid as PGrid
            with tempfile.TemporaryDirectory() as td:
                f = Path(td) / "dem.tif"
                write(f, np.where(np.isfinite(zb), zb, -9999).astype("float32"), grid, nodata=-9999)
                pg = PGrid.from_raster(str(f))
                dem = pg.read_raster(str(f))
                dem = pg.resolve_flats(pg.fill_depressions(pg.fill_pits(dem)))
                fdir_r = pg.flowdir(dem)
                acc = pg.accumulation(fdir_r)
                cells = float(p["stream_threshold_km2"]) * 1e6 / (grid.res ** 2)
                streams = acc >= cells                      # keep as a pysheds Raster
                order = np.asarray(pg.stream_order(fdir_r, streams)).astype("int32")
                fdir = np.asarray(fdir_r).astype("int64")
                write(cache / "static" / "stream_order.tif", order.astype("float32"), grid)
                write(cache / "static" / "flow_accumulation.tif", np.asarray(acc).astype("float32"), grid)
        except Exception as e:
            log.warning("KHG %s: DEM stream network failed (%s); using mapped channels only", khg_id, e)
    return {"p": p, "khg": khg, "fdir": fdir, "order": order, "base": rivers | canals_out,
            "auto": str(p["min_area_km2"]).lower() == "auto",
            "min_px": None if str(p["min_area_km2"]).lower() == "auto" else int(float(p["min_area_km2"]) * 1e6 / grid.res ** 2),
            "split_px": int(float(p["coastal_split_km2"]) * 1e6 / grid.res ** 2) if p["coastal_split_km2"] else None}


def sweep(cfg: dict, khg_id: str, geom_ll, grid: Grid, cache: Path, min_areas) -> dict:
    """How many sub-KHG each minimum area gives, with the adjacent merge and with the earlier nearest-unit merge.
    Outlets are the mapped rivers and outlet canals only. The hydrology is computed once."""
    ctx = _prepare(cfg, khg_id, geom_ll, grid, cache)
    if ctx["fdir"] is None:
        raise RuntimeError("no flow directions: the sweep needs the DEM streams (use_dem_streams)")
    lab, n = _components(ctx["fdir"], ctx["base"], ctx["khg"], ctx["split_px"])
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    contacts = _contacts(_fill_labels(lab, ctx["khg"]), ctx["base"])
    out = {"khg_id": khg_id, "min_area_km2": [float(a) for a in min_areas], "pieces": int((sizes[1:] > 0).sum()),
           "units_adjacent": [], "units_nearest": [], "unit_areas_adjacent": []}
    for a in min_areas:
        mp = int(float(a) * 1e6 / grid.res ** 2)
        out["units_nearest"].append(int((sizes[1:] >= mp).sum()))
        mapping = _merge_small(sizes, contacts, mp)
        out["units_adjacent"].append(int(len(np.unique(mapping[1:][sizes[1:] > 0]))))
        tot = np.bincount(mapping[lab[lab > 0]]) if (lab > 0).any() else np.array([0])
        out["unit_areas_adjacent"].append([round(float(v) * grid.res ** 2 / 1e6, 1) for v in sorted(tot[1:][tot[1:] > 0])])
    return out


def delineate(cfg: dict, khg_id: str, geom_ll, grid: Grid, cache: Path):
    import geopandas as gpd
    ctx = _prepare(cfg, khg_id, geom_ll, grid, cache)
    p, khg, fdir, order, base = ctx["p"], ctx["khg"], ctx["fdir"], ctx["order"], ctx["base"]
    min_px, split_px = ctx["min_px"], ctx["split_px"]
    merge = p.get("merge", "adjacent")
    pieces, chosen_area = None, (None if ctx["auto"] else float(p["min_area_km2"]))
    if ctx["auto"] and merge not in ("adjacent", "flow"):
        raise ValueError('min_area_km2: auto needs merge: adjacent or flow')
    if fdir is None:   # no flow directions: fall back to land enclosed by mapped channels
        land = khg & ~ndimage.binary_dilation(base, iterations=1)
        lab, n = ndimage.label(land, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        used, tried, pieces = None, [("mapped channels, enclosed areas", n)], n
    else:
        if p["outlet_order"] is not None:
            orders = [int(p["outlet_order"])]
        else:
            orders = [None] + list(range(int(order.max()), 0, -1))
        lo, hi = p["target_count"]
        tried, best = [], None
        for o in orders:
            out_mask = base if o is None else (base | (order >= o))
            if not out_mask.any():
                continue
            if ctx["auto"]:
                l2, n2, pc, ca = _units_auto(fdir, out_mask, khg, split_px, hi, grid.res, merge)
            else:
                l2, n2, pc = _units_from_outlets(fdir, out_mask, khg, min_px, split_px, merge)
                ca = chosen_area
            tried.append(("mapped channels" if o is None else f"order >= {o}", n2))
            score = 0 if lo <= n2 <= hi else min(abs(n2 - lo), abs(n2 - hi))
            if best is None or score < best[0]:
                best = (score, l2, n2, o, pc, ca)
            if n2 > hi or score == 0:                       # too many even for the largest outlet set, or already in range
                break
        _, lab, n, used, pieces, chosen_area = best
    log.info("KHG %s: %d sub-KHG from %s pieces (merge %s, minimum %s km2%s; outlets: mapped channels%s); tried %s", khg_id, n, pieces,
             merge, chosen_area, " chosen by auto" if ctx["auto"] else "", f" + streams of order >= {used}" if used else "", tried)
    if n == 0:
        raise RuntimeError("sub-KHG delineation produced no units")

    # assign channel cells (and, with merge nearest, dropped small pieces) to the nearest sub-KHG, then clip to the KHG
    filled_idx = ndimage.distance_transform_edt(lab == 0, return_distances=False, return_indices=True)
    full = lab[filled_idx[0], filled_idx[1]]
    full = np.where(khg, full, 0).astype("int32")
    _, full = np.unique(full, return_inverse=True)
    full = full.reshape(grid.shape).astype("int32")
    polys = [(s, int(v)) for s, v in shapes(full, mask=full > 0, transform=grid.transform)]
    from shapely.geometry import shape
    g = gpd.GeoDataFrame({"lab": [v for _, v in polys]}, geometry=[shape(s) for s, _ in polys], crs=grid.crs)
    g = g.dissolve(by="lab").reset_index()
    c = g.geometry.representative_point()
    g = g.assign(_y=-c.y, _x=c.x).sort_values(["_y", "_x"]).reset_index(drop=True)
    g["khg_id"] = khg_id
    g["unit_id"] = [f"{khg_id}_S{i + 1:02d}" for i in range(len(g))]
    g = g[["khg_id", "unit_id", "geometry"]].to_crs(4326)
    g.to_file(cache / "subkhg.gpkg", driver="GPKG")
    return g, {"n_subkhg": int(len(g)), "pieces_before_merge": pieces, "merge": merge, "min_area_km2": chosen_area,
               "min_area_setting": p["min_area_km2"], "target_count": list(p["target_count"]), "outlet_order": used, "tried": tried}
