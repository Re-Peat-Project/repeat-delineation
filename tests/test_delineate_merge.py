"""Checks of how small sub-KHG pieces are merged (repeat_delineation/delineate.py) on small hand-made rasters.

    python -m tests.test_delineate_merge
"""
import numpy as np

from repeat_delineation.delineate import _contacts, _fill_labels, _merge_small, _units_auto, _units_from_outlets


def _merge(lab, chan=None, min_px=5):
    lab = np.asarray(lab)
    chan = np.zeros(lab.shape, bool) if chan is None else np.asarray(chan, bool)
    sizes = np.bincount(lab.ravel(), minlength=int(lab.max()) + 1)
    return _merge_small(sizes, _contacts(lab, chan), min_px)


def test_small_piece_joins_the_neighbour_with_the_longer_border():
    # piece 2 (a strip of 6 cells) touches piece 1 along 6 cell edges, piece 3 along 2 and piece 4 along 4
    lab = np.zeros((6, 10), int)
    lab[:, 0:4] = 1
    lab[:, 4] = 2
    lab[0:2, 5:] = 3
    lab[2:, 5:] = 4
    m = _merge(lab, min_px=8)
    assert m[2] == m[1] and len({m[1], m[3], m[4]}) == 3, m


def test_contacts_are_counted_per_cell_edge_and_split_into_divide_and_channel():
    lab = np.array([[1, 1, 2, 2],
                    [1, 1, 2, 2],
                    [3, 3, 3, 3]])
    chan = np.zeros(lab.shape, bool)
    chan[2, :] = True                                      # the bottom row is an outlet channel
    n, div, chn = _contacts(lab, chan)
    assert n == 3
    assert div == {(1, 2): 2}, div                         # 2 edges across a divide
    assert chn == {(1, 3): 2, (2, 3): 2}, chn              # 2 + 2 edges into the channel row


def test_a_divide_is_preferred_to_an_outlet_channel_even_if_the_channel_border_is_longer():
    sizes = np.array([0, 100, 10, 100])
    m = _merge_small(sizes, (3, {(1, 2): 2}, {(2, 3): 6}), 20)      # piece 2: divide border 2 with piece 1, channel border 6 with 3
    assert m[2] == m[1] and m[3] != m[1], m
    m = _merge_small(sizes, (3, {}, {(2, 3): 6, (1, 2): 1}), 20)    # only channel borders: the longer one wins
    assert m[2] == m[3] and m[1] != m[2], m


def test_chains_of_small_pieces_end_above_the_minimum():
    lab = np.zeros((4, 12), int)
    for i, c in enumerate(range(0, 12, 2)):
        lab[:, c:c + 2] = i + 1                          # six pieces of 8 cells in a row
    m = _merge(lab, min_px=20)
    merged = m[lab]
    sizes = np.bincount(merged.ravel())[1:]
    assert (sizes >= 20).all() and sizes.sum() == lab.size, sizes
    assert sorted(set(m[1:])) == list(range(1, int(m.max()) + 1)), "labels are not consecutive"


def test_a_piece_without_neighbours_is_left_alone():
    lab = np.zeros((6, 6), int)
    lab[0:3, 0:3] = 1
    lab[5, 5] = 2                                        # far away, touches nothing
    m = _merge(lab, min_px=5)
    assert m[2] != 0 and m[2] != m[1]


def test_fill_labels_gives_every_cell_of_the_khg_a_piece():
    lab = np.zeros((5, 9), int)
    lab[:, :4] = 1
    lab[:, 5:] = 2                                       # column 4 is an outlet channel (label 0)
    khg = np.ones(lab.shape, bool)
    f = _fill_labels(lab, khg)
    assert (f > 0).all() and f[0, 4] in (1, 2)


def test_units_from_outlets_end_to_end():
    # a river along row 15; every cell flows to it (south above it, north below it)
    h, w = 31, 40
    fdir = np.full((h, w), 4, "int64")
    fdir[16:, :] = 64
    outlet = np.zeros((h, w), bool)
    outlet[15, :] = True
    outlet[:15, 20] = True                               # a tributary canal joins from the north
    # cells next to the tributary flow towards it so each side has its own catchment
    fdir[:15, :20] = 1                                   # east
    fdir[:15, 21:] = 16                                  # west
    khg = np.ones((h, w), bool)
    _, n_all, pieces = _units_from_outlets(fdir, outlet, khg, min_px=1)
    assert n_all == pieces and pieces >= 3, (n_all, pieces)
    lab, n1, _ = _units_from_outlets(fdir, outlet, khg, min_px=10 ** 6)       # nothing reaches the minimum: all merge
    assert n1 == 1, n1
    lab, n2, _ = _units_from_outlets(fdir, outlet, khg, min_px=10 ** 6, merge="nearest")
    assert n2 == 0                                       # the earlier behaviour drops everything when no piece is large enough


def test_auto_minimum_area_reaches_the_target_count():
    h, w = 31, 40
    fdir = np.full((h, w), 4, "int64")
    fdir[16:, :] = 64
    fdir[:15, :20] = 1
    fdir[:15, 21:] = 16
    outlet = np.zeros((h, w), bool)
    outlet[15, :] = True
    outlet[:15, 20] = True
    khg = np.ones((h, w), bool)
    _, n_many, pieces, a = _units_auto(fdir, outlet, khg, None, hi=100, res=1000.0)   # 1 km cells: 5 km2 = 5 cells
    assert n_many == pieces and a == 5.0, (n_many, pieces, a)
    lab, n_one, _, a2 = _units_auto(fdir, outlet, khg, None, hi=1, res=1000.0)
    assert n_one == 1 and a2 > 5.0 and (lab > 0).sum() > 0, (n_one, a2)
    lab, n_two, _, a3 = _units_auto(fdir, outlet, khg, None, hi=2, res=1000.0)
    assert 1 <= n_two <= 2 and a3 <= a2, (n_two, a3, a2)


def test_flow_merge_prefers_the_neighbour_on_the_same_channel_line():
    # same layout as the first test: piece 2 shares 6 edges with piece 1, 2 with piece 3, 4 with piece 4. Piece 2 drains to
    # reach 20, which meets reach 30 (piece 3) at a junction; pieces 1 and 4 drain to unlinked reaches 10 and 40
    lab = np.zeros((6, 10), int)
    lab[:, 0:4] = 1
    lab[:, 4] = 2
    lab[0:2, 5:] = 3
    lab[2:, 5:] = 4
    sizes = np.bincount(lab.ravel(), minlength=5)
    reach = np.array([0, 10, 20, 30, 40])
    links = {10: {10}, 20: {20, 30}, 30: {20, 30}, 40: {40}}
    contacts = _contacts(lab, np.zeros(lab.shape, bool))
    m = _merge_small(sizes, contacts, 8, (reach, links))
    assert m[2] == m[3] and len({m[1], m[3], m[4]}) == 3, m
    assert _merge_small(sizes, contacts, 8)[2] == _merge_small(sizes, contacts, 8)[1]   # without flow: the longest border


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
