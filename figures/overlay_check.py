"""Overlay the .su2 node lines on the compiled TikZ grid (fig_mesh_check.pdf) and
measure how well they coincide.

    python3 make_mesh_figure.py --check WORKDIR
    (cd WORKDIR && pdflatex fig_mesh_check.tex && pdftoppm -png -gray -r 600 -singlefile fig_mesh_check.pdf check)
    python3 overlay_check.py . sa_mi_baseline.su2 WORKDIR     # writes WORKDIR/overlay.png
"""
import sys, importlib.util
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import binary_dilation
figdir, su2, work = sys.argv[1], sys.argv[2], sys.argv[3]
spec = importlib.util.spec_from_file_location("cm", f"{figdir}/check_mesh.py"); cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)
quads, pts, marks = cm.read_su2(su2)
NI = 295; P = pts.reshape(-1, NI, 2) * 1000.0          # (NJ, NI, 2) node coordinates, mm
d = np.load(f"{work}/fig_mesh_check_lines.npz")
dpi = 600; px = dpi / 25.4; border = 1.0 / 72 * dpi
img = Image.open(f"{work}/check.png").convert("L"); A = np.array(img); Hpx, Wpx = A.shape
H = float(d["H"]); OA = d["OA"]
def to_px(xy, origin, s, win):
    X = origin[0] + (xy[:, 0] - win[0]) * s; Y = origin[1] + (xy[:, 1] - win[2]) * s
    return np.c_[(X - OA[0]) * px + border, Hpx - ((Y - OA[1]) * px + border)]
red = Image.new("L", (Wpx, Hpx), 0); dr = ImageDraw.Draw(red)
masks = []
for key, origin, s, win in (("a", d["OA"], float(d["SA"]), d["WIN_A"]), ("b", d["OB"], float(d["SB"]), d["WIN_B"])):
    ki, kj = d["ki_" + key], d["kj_" + key]
    jlo, jhi = kj.min(), kj.max()
    lines = [P[jlo:jhi + 1, i] for i in ki] + [P[j, : ki.max() + 1] for j in kj]
    lines += [P[:, 0], P[:, -1], P[0, :], P[-1, :]]           # wall, inlet, outlets
    for L in lines:
        dr.line([tuple(p) for p in to_px(L, origin, s, win)], fill=255, width=1)
    x0 = (origin[0] - OA[0]) * px + border; x1 = x0 + (win[1] - win[0]) * s * px
    y1 = Hpx - ((origin[1] - OA[1]) * px + border); y0 = y1 - (win[3] - win[2]) * s * px
    m = np.zeros_like(A, bool); m[int(y0) + 8:int(y1) - 8, int(x0) + 8:int(x1) - 8] = True   # inside frame, away from it
    masks.append((key, m))
R = np.array(red) > 0; dark = A < 150
for key, m in masks:
    su = R & m; tz = dark & m
    prec = (su & binary_dilation(tz, iterations=2)).sum() / su.sum()
    rec = (tz & binary_dilation(su, iterations=2)).sum() / tz.sum(); rec5 = (tz & binary_dilation(su, iterations=5)).sum() / tz.sum(); print(f"  within 5 px (half the 1 pt wall line): {100*rec5:.2f}%")
    print(f"panel ({key}): .su2 line pixels on a TikZ line {100*prec:.2f}%  |  TikZ line pixels on a .su2 line {100*rec:.2f}%  "
          f"({su.sum()} / {tz.sum()} px)")
rgb = np.stack([A, A, A], -1).astype(np.uint8); rgb[R] = [220, 30, 30]
Image.fromarray(rgb).save(f"{work}/overlay.png")
