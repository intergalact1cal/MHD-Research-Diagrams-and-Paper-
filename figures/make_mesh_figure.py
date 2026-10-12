"""Write fig_mesh.tex: the computational grid of sa_mi_baseline.su2.

    (a) the whole grid
    (b) the near-wall grid at the stagnation point

Every line drawn is a grid line of the mesh: data/mesh_r_i.txt and
data/mesh_phi_j.txt are written by check_mesh.py, which reads the .su2 file
and verifies that the grid is an exact polar O-grid (nodes on circles r_i and
rays phi_j to 1e-12 mm), so circles are drawn as exact TikZ arcs and rays as
straight lines through the nodes.  To keep the figure legible, lines are
skipped where neighbouring lines would print closer than a minimum spacing;
the kept lines are listed in the .tex comments.

    python3 check_mesh.py sa_mi_baseline.su2   # verify the mesh, write the levels
    python3 make_mesh_figure.py                # write fig_mesh.tex
    pdflatex fig_mesh.tex

`python3 make_mesh_figure.py --check` also writes fig_mesh_check.tex, the grid
lines alone with a fixed bounding box, used to overlay the .su2 nodes.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R_I = np.loadtxt(HERE / "data" / "mesh_r_i.txt")
PHI_J = np.loadtxt(HERE / "data" / "mesh_phi_j.txt")
R, RO = R_I[0], R_I[-1]
PHI_OUT = PHI_J[0]

# ------------------------------------------------------------------ layout (mm)
#  (a) whole grid           (b) near-wall grid
WIN_A = (-156.0, -10.0, -142.0, 142.0)
WIN_B = (-41.5, -37.7, -4.0, 4.0)
H = 140.0                                 # frame height of both panels
SA = H / (WIN_A[3] - WIN_A[2])
SB = H / (WIN_B[3] - WIN_B[2])
WA = (WIN_A[1] - WIN_A[0]) * SA
WB = (WIN_B[1] - WIN_B[0]) * SB
OA = (12.0, 10.0)                         # lower-left corner of frame (a)
OB = (OA[0] + WA + 17.0, OA[1])

# line selection
RAY_STRIDE_A = 4                          # every 4th ray (1 deg) in (a)
DMIN_A = 0.42                             # min printed spacing of circles in (a)
RAY_STRIDE_B = 1                          # every ray (0.25 deg) in (b)
DMIN_B = 0.40                             # min printed spacing of circles in (b)


def tk(phi):
    return 180.0 - phi


def f(v, nd=3):
    s = f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def thin(idx_levels, scale, dmin):
    """Indices of levels kept so that kept neighbours print >= dmin apart.

    Starts from the outer end so the coarse (outer) lines are all kept and
    only the tightly packed lines next to the wall are thinned; the wall
    itself (index 0) is always kept.
    """
    idx = list(idx_levels)
    kept = [idx[-1]]
    for k in reversed(idx[:-1]):
        if (R_I[kept[-1]] - R_I[k]) * scale >= dmin:
            kept.append(k)
    if kept[-1] != 0:
        if (R_I[kept[-1]] - R_I[0]) * scale < dmin:
            kept[-1] = 0
        else:
            kept.append(0)
    return sorted(kept)


def window_polar(win):
    x0, x1, y0, y1 = win
    cs = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]])
    ang = np.degrees(np.arctan2(cs[:, 1], -cs[:, 0]))
    rmax = np.hypot(cs[:, 0], cs[:, 1]).max()
    return ang.min(), ang.max(), rmax


def panel(win, origin, s, ray_stride, dmin, letter, ticks, extra=""):
    x0, x1, y0, y1 = win
    W, Hh = (x1 - x0) * s, (y1 - y0) * s
    pmin, pmax, rmax = window_polar(win)
    pmin, pmax = max(pmin - 1.0, -PHI_OUT), min(pmax + 1.0, PHI_OUT)
    ii = [i for i in range(len(R_I)) if R_I[i] <= rmax + 1e-9]
    if ii[-1] + 1 < len(R_I):
        ii.append(ii[-1] + 1)
    keep_i = thin(ii, s, dmin)
    jj = list(range(0, len(PHI_J), ray_stride))
    if jj[-1] != len(PHI_J) - 1:
        jj.append(len(PHI_J) - 1)
    keep_j = [j for j in jj if pmin - 1e-9 <= PHI_J[j] <= pmax + 1e-9]
    r_hi = R_I[keep_i[-1]]
    ox, oy = origin
    out = [f"% kept circles i = {keep_i}\n" if len(keep_i) < 400 else "",
           f"% kept rays j = {keep_j[0]}..{keep_j[-1]} step {ray_stride}\n",
           f"\\begin{{scope}}[shift={{({f(ox - x0 * s)},{f(oy - y0 * s)})}}, x={s:.6f}mm, y={s:.6f}mm]\n",
           "  \\begin{scope}\n",
           f"    \\clip ({f(x0, 4)},{f(y0, 4)}) rectangle ({f(x1, 4)},{f(y1, 4)});\n",
           f"    \\fill[bodygrey] (0,0) circle[radius={f(R, 4)}];\n",
           "    % circles r = r_i\n"]
    a0, a1 = tk(pmax), tk(pmin)
    for i in keep_i[1:]:
        r = R_I[i]
        out.append(f"    \\draw[grid] ({f(a0)}:{f(r, 9)}) arc[start angle={f(a0)}, end angle={f(a1)}, "
                   f"radius={f(r, 9)}];\n")
    out.append("    % rays phi = phi_j\n")
    for j in keep_j:
        a = tk(PHI_J[j])
        out.append(f"    \\draw[grid] ({f(a, 4)}:{f(R, 4)}) -- ({f(a, 4)}:{f(r_hi, 4)});\n")
    out.append("    % boundaries: wall, inlet, outlets\n")
    out.append(f"    \\draw[wall] ({f(tk(PHI_OUT))}:{f(R, 4)}) arc[start angle={f(tk(PHI_OUT))}, "
               f"end angle={f(tk(-PHI_OUT))}, radius={f(R, 4)}];\n")
    out.append(f"    \\draw[edge] ({f(tk(PHI_OUT))}:{f(RO, 4)}) arc[start angle={f(tk(PHI_OUT))}, "
               f"end angle={f(tk(-PHI_OUT))}, radius={f(RO, 4)}];\n")
    for p in (PHI_OUT, -PHI_OUT):
        out.append(f"    \\draw[edge] ({f(tk(p))}:{f(R, 4)}) -- ({f(tk(p))}:{f(RO, 4)});\n")
    out.append(extra)
    out.append("  \\end{scope}\n\\end{scope}\n")
    # frame, ticks, labels
    out.append(f"\\draw[frame] ({f(ox)},{f(oy)}) rectangle ({f(ox + W)},{f(oy + Hh)});\n")
    xt, yt = ticks
    for v in xt:
        X = ox + (v - x0) * s
        out.append(f"\\draw[tick] ({f(X)},{f(oy)}) -- ++(0,1.3) ({f(X)},{f(oy + Hh)}) -- ++(0,-1.3);\n")
        out.append(f"\\node[ticklbl, anchor=north] at ({f(X)},{f(oy - 0.4)}) {{{tick(v)}}};\n")
    for v in yt:
        Y = oy + (v - y0) * s
        out.append(f"\\draw[tick] ({f(ox)},{f(Y)}) -- ++(1.3,0) ({f(ox + W)},{f(Y)}) -- ++(-1.3,0);\n")
        out.append(f"\\node[ticklbl, anchor=east] at ({f(ox - 0.4)},{f(Y)}) {{{tick(v)}}};\n")
    out.append(f"\\node[axlbl, anchor=north] at ({f(ox + W / 2)},{f(oy - 5.4)}) {{$x$ (mm)}};\n")
    out.append(f"\\node[axlbl, anchor=south, rotate=90] at ({f(ox - 10.0)},{f(oy + Hh / 2)}) {{$y$ (mm)}};\n")
    out.append(f"\\node[panel, anchor=south west] at ({f(ox)},{f(oy + Hh + 1.4)}) {{({letter})}};\n")
    return "".join(out), (W, Hh), keep_i, keep_j


def tick(v):
    s = f"{abs(v):g}"
    return ("$-$" + s) if v < 0 else s


PREAMBLE = r"""\documentclass[border=1pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage[scaled=0.92]{helvet}
\renewcommand{\familydefault}{\sfdefault}
\IfFileExists{newtxsf.sty}{\usepackage{newtxsf}}{}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
\definecolor{ink}{HTML}{000000}
\definecolor{bodygrey}{HTML}{D9D9D9}
\definecolor{accent}{HTML}{C0392B}
\definecolor{dimgrey}{HTML}{555555}
\tikzset{
  grid/.style={draw=ink, line width=0.22pt},
  wall/.style={draw=ink, line width=1.0pt},
  edge/.style={draw=ink, line width=0.6pt},
  frame/.style={draw=ink, line width=0.5pt},
  tick/.style={draw=ink, line width=0.45pt},
  ticklbl/.style={font=\sffamily\small, inner sep=1pt},
  axlbl/.style={font=\sffamily\small, inner sep=1pt},
  panel/.style={font=\sffamily\bfseries\small, inner sep=0pt},
  note/.style={font=\sffamily\footnotesize, inner sep=1pt, align=center},
  zoom/.style={draw=accent, line width=0.8pt},
  dim/.style={draw=dimgrey, line width=0.45pt,
              {Stealth[length=1.5mm, width=1mm]}-{Stealth[length=1.5mm, width=1mm]}},
}
"""


def build(check=False):
    zx0, zx1, zy0, zy1 = WIN_B
    zoom_a = (f"    \\draw[zoom] ({f(zx0, 4)},{f(zy0, 4)}) rectangle ({f(zx1, 4)},{f(zy1, 4)});\n"
              if not check else "")
    pa, (wa, ha), ki_a, kj_a = panel(WIN_A, OA, SA, RAY_STRIDE_A, DMIN_A, "a",
                                     ([-150, -125, -100, -75, -50, -25], [-125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125]),
                                     zoom_a)
    pb, (wb, hb), ki_b, kj_b = panel(WIN_B, OB, SB, RAY_STRIDE_B, DMIN_B, "b",
                                     ([-41, -40, -39, -38], [-4, -3, -2, -1, 0, 1, 2, 3, 4]))
    # annotations
    ann = []
    if not check:
        # zoom pointer from (a) to (b): label next to the box in (a)
        bx = OA[0] + (zx0 - WIN_A[0]) * SA
        by = OA[1] + (zy1 - WIN_A[2]) * SA
        ann.append(f"\\node[note, text=accent, anchor=south east, fill=white, inner sep=0.8pt] at ({f(bx - 0.6)},{f(by + 0.6)}) {{(b)}};\n")
        # spacing zones above (b), along y = 0 (where |x| = r)
        def X(v):
            return OB[0] + (v - WIN_B[0]) * SB
        y1 = OB[1] + (3.35 - WIN_B[2]) * SB
        r_wall_end, r_uni0 = R_I[53], R_I[54]
        ann.append(f"\\draw[dim] ({f(X(-r_uni0))},{f(y1)}) -- ({f(X(WIN_B[0]))},{f(y1)});\n")
        ann.append(f"\\node[note, anchor=south, fill=white, inner sep=0.7pt] at ({f(0.5 * (X(-r_uni0) + X(WIN_B[0])))},{f(y1 + 0.3)}) "
                   f"{{uniform $\\Delta r$ = 0.20\\,mm}};\n")
        ann.append(f"\\draw[dim] ({f(X(-R))},{f(y1)}) -- ({f(X(-r_wall_end))},{f(y1)});\n")
        ann.append(f"\\node[note, anchor=south, fill=white, inner sep=0.7pt] at ({f(0.5 * (X(-R) + X(-r_wall_end)))},{f(y1 + 0.3)}) "
                   f"{{ratio 1.15}};\n")
    body = pa + pb + "".join(ann)
    if check:
        # fixed bounding box: page = frames only, so page mm map directly to data
        x_end = OB[0] + wb
        body = (f"\\useasboundingbox ({f(OA[0])},{f(OA[1])}) rectangle ({f(x_end)},{f(OA[1] + H)});\n"
                + body)
    head = ["% " + "=" * 74,
            "%  Computational grid sa_mi_baseline.su2: (a) whole grid, (b) near-wall grid",
            "%  at the stagnation point.",
            "%",
            "%  Standalone LaTeX document: on Overleaf replace ALL of main.tex with this",
            "%  file and compile with pdfLaTeX; include the compiled PDF in the paper.",
            "%  Generated by make_mesh_figure.py from data/mesh_r_i.txt and",
            "%  data/mesh_phi_j.txt (written and verified by check_mesh.py).",
            f"%  Grid: {len(R_I) - 1} x {len(PHI_J) - 1} cells, exact polar O-grid, R = {R:.1f} mm,",
            f"%  inlet r = {RO:.3f} mm, dphi = 0.25 deg, first cell 0.122 um, ratio 1.15 to",
            f"%  r = {R_I[53]:.3f} mm, uniform 0.2005 mm to r = {R_I[233]:.3f} mm, ratio 1.05 outward.",
            f"%  Lines shown: (a) every {RAY_STRIDE_A}th ray, {len(ki_a)} of {len(R_I)} circles;",
            f"%  (b) every ray, {len(ki_b)} circles.  Circles are skipped only where",
            "%  neighbouring lines would print closer than "
            f"{DMIN_A} mm (a) / {DMIN_B} mm (b).",
            "% " + "=" * 74, ""]
    tex = ("\n".join(head) + PREAMBLE + "\n\\begin{document}\n\\sffamily\\small\n"
           "\\begin{tikzpicture}[x=1mm, y=1mm, line cap=round, line join=round]\n"
           + body + "\\end{tikzpicture}\n\\end{document}\n")
    info = dict(ki_a=ki_a, kj_a=kj_a, ki_b=ki_b, kj_b=kj_b, wa=wa, wb=wb)
    return tex, info


if __name__ == "__main__":
    tex, info = build()
    (HERE / "fig_mesh.tex").write_text(tex)
    print(f"wrote fig_mesh.tex: (a) {info['wa']:.1f} x {H:.0f} mm, {len(info['ki_a'])} circles, "
          f"{len(info['kj_a'])} rays; (b) {info['wb']:.1f} x {H:.0f} mm, {len(info['ki_b'])} circles, "
          f"{len(info['kj_b'])} rays; total width {OB[0] + info['wb'] + 1:.1f} mm")
    if "--check" in sys.argv:
        tex, info = build(check=True)
        out = Path(sys.argv[sys.argv.index("--check") + 1]) if len(sys.argv) > sys.argv.index("--check") + 1 else HERE
        (out / "fig_mesh_check.tex").write_text(tex)
        np.savez(out / "fig_mesh_check_lines.npz", **{k: np.array(v) for k, v in info.items()},
                 OA=np.array(OA), OB=np.array(OB), SA=SA, SB=SB, WIN_A=np.array(WIN_A), WIN_B=np.array(WIN_B), H=H)
        print("wrote fig_mesh_check.tex")
