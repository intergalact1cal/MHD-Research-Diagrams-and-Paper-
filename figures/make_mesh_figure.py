"""Write the grid figures of sa_mi_baseline.su2 (standalone TikZ, 170 mm wide).

    fig_mesh.tex        (a) upper half of the grid (the grid is symmetric about
                            y = 0), (b) near-wall grid at the stagnation point
    fig_mesh_full.tex   (a) whole grid, (b) near-wall grid

Every line drawn is a grid line of the mesh.  data/mesh_r_i.txt and
data/mesh_phi_j.txt are written by check_mesh.py, which reads the .su2 file and
verifies that it is an exact polar O-grid (nodes on circles r_i and rays phi_j
to 1e-12 mm), symmetric about y = 0.  Circles are drawn as exact arcs, rays as
straight lines through the nodes.  Panel (a) shows every 4th line in each
direction (a fixed stride, so the drawn spacing varies exactly as the grid
does); panel (b) shows every grid line.

    python3 check_mesh.py sa_mi_baseline.su2   # verify the mesh, write the levels
    python3 make_mesh_figure.py                # write the .tex files
    pdflatex fig_mesh.tex

`python3 make_mesh_figure.py --check DIR` also writes grid-only versions with a
fixed bounding box (DIR/<name>_check.tex) for overlay_check.py.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R_I = np.loadtxt(HERE / "data" / "mesh_r_i.txt")
PHI_J = np.loadtxt(HERE / "data" / "mesh_phi_j.txt")
R, RO = R_I[0], R_I[-1]
PHI_OUT = PHI_J[0]
NI, NJ = len(R_I), len(PHI_J)

TOTAL_W = 170.0
STRIDE_A = 4                         # every 4th line in each direction in (a)
STRIDE_B = 1                         # every line in (b)

VARIANTS = {
    "fig_mesh": dict(
        title="upper half of the grid (symmetric about y = 0)",
        win_a=(-154.0, -12.0, 0.0, 140.0), height=90.0,
        b_x=(-41.6, -37.6), b_y0=0.0,
        ticks_a=([-150, -125, -100, -75, -50, -25], [0, 25, 50, 75, 100, 125]),
        ticks_bx=[-41, -40, -39, -38]),
    "fig_mesh_full": dict(
        title="whole grid",
        win_a=(-154.0, -12.0, -140.0, 140.0), height=118.0,
        b_x=(-41.6, -37.6), b_y0=None,
        ticks_a=([-150, -125, -100, -75, -50, -25], [-100, -50, 0, 50, 100]),
        ticks_bx=[-41, -40, -39, -38]),
}


def tk(phi):
    return 180.0 - phi


def f(v, nd=3):
    s = f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def tick(v):
    s = f"{abs(v):g}"
    return ("$-$" + s) if v < 0 else s


def strided(n, stride):
    """Indices 0, stride, 2 stride, ... plus the last one."""
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return idx


def layout(v):
    x0, x1, y0, y1 = v["win_a"]
    H = v["height"]
    sa = H / (y1 - y0)
    wa = (x1 - x0) * sa
    oa = (12.0, 10.0)
    ob = (oa[0] + wa + 16.0, oa[1])
    wb = TOTAL_W - ob[0] - 1.0
    bx0, bx1 = v["b_x"]
    sb = wb / (bx1 - bx0)
    span = H / sb
    by0 = v["b_y0"] if v["b_y0"] is not None else -span / 2
    win_b = (bx0, bx1, by0, by0 + span)
    return dict(sa=sa, wa=wa, oa=oa, ob=ob, wb=wb, sb=sb, win_b=win_b, H=H)


def window_angles(win):
    x0, x1, y0, y1 = win
    cs = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]])
    ang = np.degrees(np.arctan2(cs[:, 1], -cs[:, 0]))
    return ang.min(), ang.max(), np.hypot(cs[:, 0], cs[:, 1]).max()


def panel(win, origin, s, stride, letter, ticks, extra=""):
    x0, x1, y0, y1 = win
    W, Hh = (x1 - x0) * s, (y1 - y0) * s
    pmin, pmax, rmax = window_angles(win)
    pmin, pmax = max(pmin - 1.0, -PHI_OUT), min(pmax + 1.0, PHI_OUT)
    ki = [i for i in strided(NI, stride) if R_I[i] <= rmax or i == 0]
    nxt = [i for i in strided(NI, stride) if R_I[i] > rmax]
    ki += nxt[:1]
    kj = [j for j in strided(NJ, stride) if pmin - 1e-9 <= PHI_J[j] <= pmax + 1e-9]
    r_hi = R_I[ki[-1]]
    ox, oy = origin
    a0, a1 = tk(pmax), tk(pmin)
    out = [f"% ({letter}) circles i = {ki[0]}..{ki[-1]} step {stride}; rays j = {kj[0]}..{kj[-1]} step {stride}\n",
           f"\\begin{{scope}}[shift={{({f(ox - x0 * s)},{f(oy - y0 * s)})}}, x={s:.6f}mm, y={s:.6f}mm]\n",
           "  \\begin{scope}\n",
           f"    \\clip ({f(x0, 4)},{f(y0, 4)}) rectangle ({f(x1, 4)},{f(y1, 4)});\n",
           f"    \\fill[bodygrey] (0,0) circle[radius={f(R, 4)}];\n",
           "    % circles r = r_i\n"]
    for i in ki[1:]:
        r = R_I[i]
        out.append(f"    \\draw[grid] ({f(a0)}:{f(r, 9)}) arc[start angle={f(a0)}, end angle={f(a1)}, "
                   f"radius={f(r, 9)}];\n")
    out.append("    % rays phi = phi_j\n")
    for j in kj:
        a = tk(PHI_J[j])
        out.append(f"    \\draw[grid] ({f(a, 4)}:{f(R, 4)}) -- ({f(a, 4)}:{f(min(r_hi, RO), 4)});\n")
    out.append("    % boundaries: wall, inlet, outlets\n")
    out.append(f"    \\draw[wall] ({f(tk(PHI_OUT))}:{f(R, 4)}) arc[start angle={f(tk(PHI_OUT))}, "
               f"end angle={f(tk(-PHI_OUT))}, radius={f(R, 4)}];\n")
    out.append(f"    \\draw[edge] ({f(tk(PHI_OUT))}:{f(RO, 4)}) arc[start angle={f(tk(PHI_OUT))}, "
               f"end angle={f(tk(-PHI_OUT))}, radius={f(RO, 4)}];\n")
    for p in (PHI_OUT, -PHI_OUT):
        out.append(f"    \\draw[edge] ({f(tk(p))}:{f(R, 4)}) -- ({f(tk(p))}:{f(RO, 4)});\n")
    out.append(extra)
    out.append("  \\end{scope}\n\\end{scope}\n")
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
    return "".join(out), ki, kj


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
  note/.style={font=\sffamily\footnotesize, inner sep=0.7pt, align=center, fill=white},
  zoom/.style={draw=accent, line width=0.8pt},
  dim/.style={draw=dimgrey, line width=0.45pt,
              {Stealth[length=1.5mm, width=1mm]}-{Stealth[length=1.5mm, width=1mm]}},
}
"""


def build(name, check=False):
    v = VARIANTS[name]
    L = layout(v)
    zx0, zx1, zy0, zy1 = L["win_b"]
    zoom = "" if check else f"    \\draw[zoom] ({f(zx0, 4)},{f(zy0, 4)}) rectangle ({f(zx1, 4)},{f(zy1, 4)});\n"
    pa, ki_a, kj_a = panel(v["win_a"], L["oa"], L["sa"], STRIDE_A, "a", v["ticks_a"], zoom)
    yb = [t for t in range(int(np.ceil(zy0)), int(np.floor(zy1)) + 1)]
    pb, ki_b, kj_b = panel(L["win_b"], L["ob"], L["sb"], STRIDE_B, "b", (v["ticks_bx"], yb))
    ann = []
    if not check:
        x0a, _, y0a, _ = v["win_a"]
        bx = L["oa"][0] + (zx0 - x0a) * L["sa"]
        by = L["oa"][1] + (zy1 - y0a) * L["sa"]
        ann.append(f"\\node[note, text=accent, anchor=south east] at ({f(bx - 0.6)},{f(by + 0.6)}) {{(b)}};\n")
        def X(val):
            return L["ob"][0] + (val - zx0) * L["sb"]
        yb_page = L["ob"][1] + L["H"] - 4.2
        r_wall_end, r_uni0 = R_I[53], R_I[54]
        ann.append(f"\\draw[dim] ({f(X(-r_uni0))},{f(yb_page)}) -- ({f(X(zx0))},{f(yb_page)});\n")
        ann.append(f"\\node[note, anchor=south] at ({f(0.5 * (X(-r_uni0) + X(zx0)))},{f(yb_page + 0.4)}) "
                   f"{{uniform 0.2\\,mm}};\n")
        ann.append(f"\\draw[dim] ({f(X(-R))},{f(yb_page)}) -- ({f(X(-r_wall_end))},{f(yb_page)});\n")
        ann.append(f"\\node[note, anchor=south] at ({f(0.5 * (X(-R) + X(-r_wall_end)))},{f(yb_page + 0.4)}) "
                   f"{{ratio 1.15}};\n")
    body = pa + pb + "".join(ann)
    if check:
        body = (f"\\useasboundingbox ({f(L['oa'][0])},{f(L['oa'][1])}) rectangle "
                f"({f(L['ob'][0] + L['wb'])},{f(L['oa'][1] + L['H'])});\n" + body)
    head = ["% " + "=" * 74,
            f"%  Computational grid sa_mi_baseline.su2: (a) {v['title']},",
            "%  (b) near-wall grid at the stagnation point.",
            "%",
            "%  Standalone LaTeX document: on Overleaf replace ALL of main.tex with this",
            "%  file and compile with pdfLaTeX; include the compiled PDF in the paper.",
            "%  Generated by make_mesh_figure.py from data/mesh_r_i.txt and",
            "%  data/mesh_phi_j.txt (written and verified by check_mesh.py).",
            f"%  Grid: {NI - 1} x {NJ - 1} cells, exact polar O-grid, symmetric about y = 0;",
            f"%  R = {R:.1f} mm, inlet r = {RO:.3f} mm, dphi = 0.25 deg; radial spacing: first cell",
            f"%  0.122 um, ratio 1.15 to r = {R_I[53]:.3f} mm, uniform 0.2005 mm to r = {R_I[233]:.3f} mm,",
            "%  ratio 1.05 to the inlet.",
            f"%  Lines shown: (a) every {STRIDE_A}th line in each direction; (b) every grid line.",
            "% " + "=" * 74, ""]
    tex = ("\n".join(head) + PREAMBLE + "\n\\begin{document}\n\\sffamily\\small\n"
           "\\begin{tikzpicture}[x=1mm, y=1mm, line cap=round, line join=round]\n"
           + body + "\\end{tikzpicture}\n\\end{document}\n")
    info = dict(ki_a=ki_a, kj_a=kj_a, ki_b=ki_b, kj_b=kj_b, OA=L["oa"], OB=L["ob"], SA=L["sa"],
                SB=L["sb"], WIN_A=v["win_a"], WIN_B=L["win_b"], H=L["H"], WB=L["wb"])
    return tex, info


if __name__ == "__main__":
    check_dir = Path(sys.argv[sys.argv.index("--check") + 1]) if "--check" in sys.argv else None
    for name in VARIANTS:
        tex, info = build(name)
        (HERE / f"{name}.tex").write_text(tex)
        print(f"wrote {name}.tex: (a) {len(info['ki_a'])} circles, {len(info['kj_a'])} rays; "
              f"(b) {len(info['ki_b'])} circles, {len(info['kj_b'])} rays; window (b) "
              f"{tuple(round(x, 2) for x in info['WIN_B'])}; size {TOTAL_W:.0f} x {info['H'] + 17:.0f} mm")
        if check_dir:
            tex, info = build(name, check=True)
            (check_dir / f"{name}_check.tex").write_text(tex)
            np.savez(check_dir / f"{name}_check_lines.npz",
                     **{k: np.array(val) for k, val in info.items()})
