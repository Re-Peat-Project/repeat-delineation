"""Why are cells left unassigned? Classifies them and saves output/unassigned_30.png (run: python diagnose_unassigned.py)"""
import sys
sys.argv = ["harness.py", "30", '{"outlet_order":3}', "_diag"]
from _common import OUT
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import repeat_delineation.delineate as D

orig_units = D._units_from_outlets
def hook(fdir, outlet, khg, min_px, split_px=None):
    lab, n = orig_units(fdir, outlet, khg, min_px, split_px)
    term = D._terminal(fdir, outlet)
    h, w = fdir.shape
    tr, tc = np.divmod(term, w)
    idx = np.arange(h * w).reshape(h, w)
    reaches_outlet = outlet.ravel()[term].reshape(h, w)
    un = khg & (lab == 0) & ~outlet
    print("fdir values in KHG:", dict(zip(*[a.tolist() for a in np.unique(fdir[khg], return_counts=True)])))
    print(f"unassigned {un.sum()/khg.sum():.1%} of KHG;  of which reach no outlet: {(un & ~reaches_outlet).sum()/un.sum():.1%},"
          f" small pieces dropped: {(un & reaches_outlet).sum()/un.sum():.1%}")
    term_in_khg = khg.ravel()[term].reshape(h, w)
    print(f"  unassigned whose terminal cell lies outside KHG footprint: {(un & ~term_in_khg).sum()/un.sum():.1%}")
    D._dbg = (khg, un, reaches_outlet, outlet)
    return lab, n
D._units_from_outlets = hook
exec(open("harness.py").read())
khg, un, reach, outlet = D._dbg
img = np.zeros(khg.shape + (3,)); img[khg] = .85
img[un & ~reach] = (.9, .2, .2); img[un & reach] = (.2, .2, .9); img[outlet & khg] = (0, .6, 0)
plt.figure(figsize=(8, 13)); plt.imshow(img); plt.title("red: drains to non-outlet, blue: small pieces dropped, green: outlet channels")
plt.axis("off"); plt.savefig(OUT / "unassigned_30.png", dpi=80, bbox_inches="tight")
