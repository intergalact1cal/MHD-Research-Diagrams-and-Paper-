"""Write the combined set-up figure (standalone TikZ, 170 mm wide).

    fig_setup.tex          (a) domain, boundary conditions and grid
                           (b) Edney type IV interaction, with the jet shock cells
                           (c) grid in the shock layer at the stagnation point
    fig_setup_clean.tex    same, with the simpler (b) (no shock cells)

Layout: (a) is a tall panel on the left; (b) and (c) are stacked on the
right.  The figure is drawn at its print size (170 mm wide) so the 9 pt / 8 pt
labels appear at that size in the paper.  Each .tex is a complete standalone
LaTeX document (pdfLaTeX, e.g. Overleaf); all coordinates are written out.

Grid lines are the real grid lines of sa_mi_baseline.su2 (an exact polar
O-grid: circles at r_i, rays at phi_j).  Where neighbouring lines would be
closer than a minimum spacing on paper, lines are skipped, so the grading is
visible without lines merging.

Units: mm, origin at the cylinder centre, x downstream, y up; phi is measured
from the stagnation point, positive up.
"""
from pathlib import Path

import numpy as np

import fig_setup_typeiv as g
import jet_geometry as jg

HERE = Path(__file__).resolve().parent

R, RO = g.R, g.R_OUT
D = np.degrees
PHI_OUT, PHI_NS = D(g.PHI_OUT), D(g.PHI_NOSLIP)
PHI_SPLIT, BETA, THETA = D(g.PHI_SPLIT), D(g.BETA), D(g.THETA)
T1, T2, SPLIT = g.T1, g.T2, g.SPLIT

R_I = np.loadtxt(HERE / "data" / "mesh_r_i.txt")
PHI_J = np.loadtxt(HERE / "data" / "mesh_phi_j.txt")
R_WALL_END, R_UNI0, R_UNI1 = R_I[53], R_I[54], R_I[233]

# ---------------------------------------------------------------- layout (mm)
TOTAL_W = 170.0
WIN_A = (-158.0, 10.0, -145.0, 145.0)      # panel (a) data window
WA = 62.0                                   # width of (a) on the page
GAP = 7.0
X_RIGHT = WA + GAP                          # left edge of (b), (c)
W_RIGHT = TOTAL_W - X_RIGHT
SA = WA / (WIN_A[1] - WIN_A[0])
HA = (WIN_A[3] - WIN_A[2]) * SA             # total height
WIN_B_X = (-86.0, -30.0)
SB = W_RIGHT / (WIN_B_X[1] - WIN_B_X[0])
H_C = 30.0                                  # height of (c)
GAP_BC = 7.0
H_B = HA - H_C - GAP_BC
WIN_B = (*WIN_B_X, -32.0, -32.0 + H_B / SB)
WIN_C_X = (-80.0, -36.0)
SC = W_RIGHT / (WIN_C_X[1] - WIN_C_X[0])
WIN_C = (*WIN_C_X, -0.5 * H_C / SC, 0.5 * H_C / SC)


# ======================================================================
#  helpers
# ======================================================================
def tk(phi):
    return 180.0 - phi


def f(v, nd=2):
    s = f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def pt(p, nd=2):
    return f"({f(p[0], nd)},{f(p[1], nd)})"


def poly(xy, step=1, indent="      "):
    xy = np.asarray(xy, dtype=float)
    idx = list(range(0, len(xy), step))
    if idx[-1] != len(xy) - 1:
        idx.append(len(xy) - 1)
    pts = [pt(xy[i]) for i in idx]
    lines, cur = [], ""
    for p in pts:
        if cur and len(cur) + len(p) > 84:
            lines.append(cur.rstrip())
            cur = ""
        cur += p + " -- "
    lines.append(cur[:-4])
    return ("\n" + indent).join(lines)


def wall(phi_deg, r=R):
    p = np.radians(phi_deg)
    return np.array([-r * np.cos(p), r * np.sin(p)])


def arc_cmd(r, phi0, phi1, nd=4):
    a0, a1 = tk(phi0), tk(phi1)
    return (f"({f(a0, 3)}:{f(r, nd)}) arc[start angle={f(a0, 3)}, "
            f"end angle={f(a1, 3)}, radius={f(r, nd)}]")


_names = {}


def leader(text, at, target, anchor, frm=None, style="lbl", extra=()):
    frm = frm or anchor
    key = "".join(c for c in text if c.isalnum())[:10]
    _names[key] = _names.get(key, 0) + 1
    name = f"L{key}{_names[key]}"
    s = f"  \\node[{style}, anchor={anchor}] ({name}) at {pt(at)} {{{text}}};\n"
    for t in (target, *extra):
        s += f"  \\draw[leader] ({name}.{frm}) -- {pt(t)};\n"
    return s


def grid_lines(stride, phi_lo, phi_hi, r_hi, style, ray_stride=None):
    """TikZ commands for every `stride`-th circle and every `ray_stride`-th ray
    (default: the same), exactly at the node levels of the mesh."""
    ray_stride = ray_stride or stride
    ii = list(range(0, len(R_I), stride))
    if ii[-1] != len(R_I) - 1:
        ii.append(len(R_I) - 1)
    ri = R_I[ii]
    ri = ri[ri <= r_hi + 1e-9]
    jj = list(range(0, len(PHI_J), ray_stride))
    if jj[-1] != len(PHI_J) - 1:
        jj.append(len(PHI_J) - 1)
    pj = PHI_J[jj]
    pj = pj[(pj >= phi_lo - 1e-9) & (pj <= phi_hi + 1e-9)]
    out = [f"    \\draw[{style}] {arc_cmd(r, phi_hi, phi_lo, nd=9)};\n" for r in ri[1:]]
    rr = min(r_hi, R_I[-1])
    out += [f"    \\draw[{style}] ({f(tk(p), 4)}:{R}) -- ({f(tk(p), 4)}:{f(rr, 4)});\n"
            for p in pj]
    return "".join(out), len(ri), len(pj)


PREAMBLE = r"""\documentclass[border=1pt]{standalone}
% ---- fonts: Helvetica text; sans-serif math (newtxsf) when available
\usepackage[T1]{fontenc}
\usepackage[scaled=0.92]{helvet}
\renewcommand{\familydefault}{\sfdefault}
\IfFileExists{newtxsf.sty}{\usepackage{newtxsf}}{}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}

% ---- colours
\definecolor{ink}{HTML}{111111}
\definecolor{shockgrey}{HTML}{3C3C3C}
\definecolor{dimgrey}{HTML}{7A7A7A}
\definecolor{leadergrey}{HTML}{9A9A9A}
\definecolor{gridgrey}{HTML}{BDBDBD}
\definecolor{bodygrey}{HTML}{E2E2E2}
\definecolor{freeblue}{HTML}{1A5FA8}
\definecolor{postgreen}{HTML}{0E8A62}
\definecolor{slipochre}{HTML}{B06A10}
\definecolor{jetfill}{HTML}{FBEBD3}
\definecolor{accent}{HTML}{C0392B}
\definecolor{shockref}{HTML}{7B5AA6}

% ---- styles
\tikzset{
  lbl/.style={font=\sffamily\small, text=ink, inner sep=1.1pt, align=center},
  sml/.style={font=\sffamily\footnotesize, text=ink, inner sep=1.1pt, align=center},
  halo/.style={fill=white, fill opacity=0.85, text opacity=1},
  region/.style={font=\sffamily\small, text=dimgrey, inner sep=1pt},
  leader/.style={draw=leadergrey, line width=0.4pt},
  freeinlet/.style={draw=freeblue, line width=1.5pt},
  postinlet/.style={draw=postgreen, line width=1.5pt},
  outlet/.style={draw=dimgrey, line width=1.1pt, dash pattern=on 3pt off 2pt, line cap=butt},
  noslip/.style={draw=ink, line width=1.7pt},
  slip/.style={draw=slipochre, line width=1.7pt, dash pattern=on 2pt off 1.1pt, line cap=butt},
  shockA/.style={draw=shockref, line width=0.75pt, draw opacity=0.9},
  shockB/.style={draw=ink, line width=1.1pt},
  shear/.style={draw=ink, line width=0.7pt, dash pattern=on 2.4pt off 1.4pt, line cap=butt},
  cellshock/.style={draw=ink, line width=0.45pt},
  fan/.style={draw=ink, line width=0.45pt, dash pattern=on 1.4pt off 0.9pt, line cap=butt},
  hatch/.style={draw=dimgrey, line width=0.35pt, line cap=butt},
  flow/.style={line width=0.9pt, -{Stealth[length=2.3mm, width=1.7mm]}},
  jetflow/.style={draw=ink, line width=0.75pt, -{Stealth[length=2mm, width=1.5mm]}},
  dim/.style={draw=dimgrey, line width=0.45pt,
              {Stealth[length=1.5mm, width=1mm]}-{Stealth[length=1.5mm, width=1mm]}},
  triple/.style={circle, fill=ink, inner sep=0pt, minimum size=1.6mm},
  dot/.style={circle, fill=leadergrey, inner sep=0pt, minimum size=1mm},
  panel/.style={font=\sffamily\bfseries\small, inner sep=0pt},
  frame/.style={draw=ink, line width=0.45pt},
  gridA/.style={draw=gridgrey, line width=0.25pt},
  gridC/.style={draw=ink, line width=0.3pt},
  gridedge/.style={draw=ink, line width=0.6pt},
  window/.style={draw=ink, line width=0.5pt, dash pattern=on 1.6pt off 1.1pt, line cap=butt},
  windowC/.style={draw=accent, line width=0.6pt},
}
"""


# ======================================================================
#  shared flow geometry
# ======================================================================
xt, yt = g.transmitted(40)
TRANS = np.c_[xt, yt]
JSHOCK = jg.jet_shock()
AFTER_U, AFTER_L = jg.after_shock(+1, 30), jg.after_shock(-1, 30)
SHOCKS, FANS = jg.wave_pattern()
JET_FILL = np.r_[jg.UPPER_FULL, jg.LOWER[::-1], TRANS[::-1][1:]]
TRACE = np.loadtxt(HERE / "data" / "bow_shock_no_field.csv", delimiter=",")


def jet_mid(frac):
    return 0.5 * (jg.at(frac, +1) + jg.at(frac, -1))


def wedge_point(x):
    yu = np.interp(x, jg.UPPER_FULL[:, 0], jg.UPPER_FULL[:, 1])
    yv = np.interp(x, TRANS[:, 0], TRANS[:, 1])
    return np.array([x, 0.5 * (yu + yv)])


def shift_scope(origin, win, s):
    ox, oy = origin
    return (f"\\begin{{scope}}[shift={{({f(ox - win[0] * s, 3)},{f(oy - win[2] * s, 3)})}}, "
            f"x={s:.5f}mm, y={s:.5f}mm]\n")


# ======================================================================
#  (a) domain, boundary conditions, grid
# ======================================================================
def panel_a():
    s = SA
    up_a = g.upper_bow(np.linspace(g.y_t1, 165, 140)).T
    lo_a = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -175, 140)])
    inc_up = np.array([WIN_A[0] + 1, g.incident(WIN_A[0] + 1)])
    grid, nr, nphi = grid_lines(4, -PHI_OUT, PHI_OUT, R_I[-1], "gridA", ray_stride=8)
    # 150 mm dimension along the lower outlet, offset outward
    phi_d = -PHI_OUT
    u = wall(phi_d, 1.0)
    nrm = np.array([u[1], -u[0]])
    if nrm[0] < 0:
        nrm = -nrm
    off = 9.0 * nrm
    d0, d1 = off, wall(phi_d, RO) + off
    rot = np.degrees(np.arctan2(d1[1] - d0[1], d1[0] - d0[0]))
    rot = rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot
    lbl150 = 0.55 * (d0 + d1) + 4.5 * nrm
    bx0, bx1, by0, by1 = WIN_B
    cx0, cx1, cy0, cy1 = WIN_C

    return rf"""
% ======================================================================
%  (a) computational domain, boundary conditions and grid
%      grid: every 4th circle and every 8th ray (2 deg), at the exact node levels
%      shocks: no-field solution, drawn for reference (violet), not geometry
% ======================================================================
{shift_scope((0, 0), WIN_A, s)}  \begin{{scope}}
    \clip ({f(WIN_A[0])},{f(WIN_A[2])}) rectangle ({f(WIN_A[1])},{f(WIN_A[3])});
    \fill[bodygrey] (0,0) circle[radius={R}];
    \draw[dimgrey, line width=0.45pt, dash pattern=on 0.4pt off 1.3pt]
      {arc_cmd(R, PHI_OUT, -PHI_OUT + 360)};
  \end{{scope}}

  % grid ({nr} circles, {nphi} rays shown)
{grid}
  % shock system (no-field solution)
  \begin{{scope}}
    \clip {arc_cmd(RO, PHI_OUT, -PHI_OUT)} -- ({f(tk(-PHI_OUT), 3)}:{R})
      arc[start angle={f(tk(-PHI_OUT), 3)}, end angle={f(tk(PHI_OUT), 3)}, radius={R}] -- cycle;
    \draw[shockA]
      {poly(up_a, 2)};
    \draw[shockA]
      {poly(lo_a, 2)};
    \draw[shockA]
      {poly(TRANS, 2)};
    \draw[shockA] {pt(SPLIT)} -- {pt(T1)};
  \end{{scope}}
  \draw[shockA, dash pattern=on 2.4pt off 1.6pt, line cap=butt] {pt(inc_up)} -- {pt(SPLIT)};

  % boundaries
  \draw[freeinlet] {arc_cmd(RO, PHI_SPLIT, PHI_OUT)};
  \draw[postinlet] {arc_cmd(RO, -PHI_OUT, PHI_SPLIT)};
  \draw[outlet] ({f(tk(PHI_OUT))}:{R}) -- ({f(tk(PHI_OUT))}:{f(RO, 3)});
  \draw[outlet] ({f(tk(-PHI_OUT))}:{R}) -- ({f(tk(-PHI_OUT))}:{f(RO, 3)});
  \draw[slip]   {arc_cmd(R, PHI_OUT, PHI_NS)};
  \draw[slip]   {arc_cmd(R, -PHI_NS, -PHI_OUT)};
  \draw[noslip] {arc_cmd(R, PHI_NS, -PHI_NS)};
  \draw[ink, line width=0.5pt] (-4,0) -- (4,0) (0,-4) -- (0,4);

  % windows of (b) and (c)
  \draw[window] ({f(bx0)},{f(by0)}) rectangle ({f(bx1)},{f(by1)});
  \node[sml, anchor=south west, inner sep=0.8pt, halo] at ({f(bx0)},{f(by1)}) {{(b)}};
  \draw[windowC] ({f(cx0)},{f(cy0, 3)}) rectangle ({f(cx1)},{f(cy1, 3)});
  \node[sml, text=accent, anchor=south east, inner sep=0.8pt, halo] at ({f(cx1)},{f(by1 + 0.3, 3)}) {{(c)}};

  % inlet split and shock angle
  \draw[ink, line width=0.45pt] {pt(SPLIT + [40, 0])} arc[start angle=0, end angle={f(BETA)}, radius=40];
  \draw[dimgrey, line width=0.4pt, dash pattern=on 0.8pt off 1pt] {pt(SPLIT)} -- ++(52,0);
  \node[lbl, halo, anchor=north west] at {pt(SPLIT + [9, -3])} {{$\beta$ = 18.1$^\circ$}};
  \node[sml, halo, text=dimgrey, anchor=north west] at {pt(SPLIT + [9, -17])} {{$\phi$ = $-$12.25$^\circ$}};
  \filldraw[fill=white, draw=ink, line width=0.75pt] {pt(SPLIT)} circle[radius=4.6];

  % inflow
  \foreach \yy in {{116, 104}} {{\draw[flow, freeblue] (-156,\yy) -- ++(28,0);}}
  \node[lbl, text=freeblue, anchor=north west] at (-157,143) {{$M_\infty$ = 8.03}};
  \foreach \yy in {{-103, -115}} {{\draw[flow, postgreen] (-156,\yy) -- ++({THETA:.2f}:28.7);}}
  \node[lbl, text=postgreen, anchor=south west] at (-157,-143) {{$M_2$ = 5.25}};

  % boundary labels
  \node[lbl, text=dimgrey, anchor=west] at (-30,84) {{outlet}};
  \node[lbl, text=dimgrey, anchor=west] at (-30,-84) {{outlet}};
  \node[lbl, anchor=west] (sw) at (-24,66) {{slip wall}};
  \draw[leader] (sw.west) -- {pt(wall(57))};
  \node[sml, anchor=center] (ns) at (-11.5,12) {{no-slip wall\\$T_w$ = 294\,K}};
  \draw[leader] (ns.west) -- {pt(wall(12))};
  \node[lbl, halo, text=shockref, anchor=west, align=left] (bs) at (-108,66) {{bow shock\\[-0.4ex]{{\footnotesize (reference)}}}};
  \draw[leader] (bs.east) -- {pt(g.upper_bow(64.0))};

  % dimensions
  \draw[dim] (0,0) -- (270:{R});
  \node[sml, anchor=north] at (0,{f(-R - 2.5)}) {{$R$ = {R}\,mm}};
  \draw[dim] (0,0) -- ({f(tk(-55))}:{f(RO, 3)});
  \node[sml, text=dimgrey, halo, rotate=55] at {pt(wall(-55, 0.70 * RO) + [-3.6, 2.5])} {{150\,mm}};
\end{{scope}}
\node[panel, anchor=south west] at (0,{f(HA + 1.2)}) {{(a)}};
"""


# ======================================================================
#  (b) type IV interaction
# ======================================================================
def panel_b(detailed):
    zx0, zx1, zy0, zy1 = WIN_B
    s = SB
    origin = (X_RIGHT, H_C + GAP_BC)
    up_b = g.upper_bow(np.linspace(g.y_t1, zy1 + 6, 90)).T
    lo_b = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -80, 160)])
    lo_b = lo_b[lo_b[:, 1] > zy0 - 6]
    inc0 = np.array([zx0 - 4, g.incident(zx0 - 4)])

    if detailed:
        bodyfill = [f"    \\fill[white] (0,0) circle[radius={R}];\n",
                    "    \\begin{scope}\n", f"      \\clip (0,0) circle[radius={R}];\n"]
        dd = np.array([np.cos(np.radians(-50)), np.sin(np.radians(-50))])
        for phi in np.arange(-50, 50.01, 2.4):
            p = wall(phi)
            bodyfill.append(f"      \\draw[hatch] {pt(p - 0.2 * dd)} -- {pt(p + 2.6 * dd)};\n")
        bodyfill.append("    \\end{scope}\n")
        bodyfill = "".join(bodyfill)
        jetfill = ""
    else:
        bodyfill = f"    \\fill[bodygrey] (0,0) circle[radius={R}];\n"
        jetfill = "    \\fill[jetfill]\n      " + poly(JET_FILL, 2) + " -- cycle;\n"

    out = [rf"""
% ======================================================================
%  (b) Edney type IV interaction
% ======================================================================
{shift_scope(origin, WIN_B, s)}  \begin{{scope}}
    \clip ({f(zx0)},{f(zy0, 3)}) rectangle ({f(zx1)},{f(zy1, 3)});
{bodyfill}    \draw[noslip, line width=2pt] {arc_cmd(R, PHI_NS, -PHI_NS)};
{jetfill}    \draw[shockB]
      {poly(up_b)};
    \draw[shockB]
      {poly(lo_b, 2)};
    \draw[shockB] {pt(inc0)} -- {pt(T1)};
"""]
    if detailed:
        out.append("    % jet shock cells: compressions (solid), expansion fans (dashed)\n")
        out += [f"    \\draw[cellshock] {pt(a)} -- {pt(b)};\n" for a, b in SHOCKS]
        out += [f"    \\draw[fan] {pt(a)} -- {pt(b)};\n" for a, b in FANS]
    a0 = wedge_point(-61.0)
    out.append(rf"""    \draw[shockB]
      {poly(TRANS)};
    \draw[shear]
      {poly(jg.UPPER_FULL, 4)};
    \draw[shear]
      {poly(jg.LOWER, 4)};
    \draw[shear]
      {poly(AFTER_U, 3)};
    \draw[shear]
      {poly(AFTER_L, 3)};
    \draw[shockB, line width=1.1pt] {pt(JSHOCK[0])} -- {pt(JSHOCK[1])};
  \end{{scope}}
  \draw[frame] ({f(zx0)},{f(zy0, 3)}) rectangle ({f(zx1)},{f(zy1, 3)});
  \node[triple] at {pt(T1)} {{}};
  \node[triple] at {pt(T2)} {{}};
  \draw[jetflow] {pt(a0)} -- {pt(a0 + 4.4 * jg.D0)};

  % scale bar
  \draw[ink, line width=0.8pt, line cap=butt] (-84.5,{f(zy0 + 2.4)}) -- (-74.5,{f(zy0 + 2.4)});
  \draw[ink, line width=0.6pt] (-84.5,{f(zy0 + 1.9)}) -- (-84.5,{f(zy0 + 2.9)}) (-74.5,{f(zy0 + 1.9)}) -- (-74.5,{f(zy0 + 2.9)});
  \node[sml, anchor=south] at (-79.5,{f(zy0 + 3.0)}) {{10\,mm}};

  % labels
""")
    sl_a = 0.5 * (T1 + jg.UPPER[0])
    lab = ""
    lab += f"  \\node[lbl, text=freeblue, anchor=north west] at ({f(zx0 + 0.8)},{f(zy1 - 0.6)}) {{$M_\\infty$}};\n"
    lab += f"  \\node[lbl, text=postgreen, anchor=north west] at ({f(zx0 + 0.8)},{f(g.incident(zx0) - 2.6)}) {{$M_2$}};\n"
    lab += leader("bow shock", (-74.2, 0.5), g.upper_bow(0.5), "east")
    lab += leader(r"upper\\[-0.3ex]triple point", (-69.5, -2.2), T1, "west")
    lab += leader("incident shock", (-73.5, -17.5), (-77.5, g.incident(-77.5)), "north")
    lab += leader(r"transmitted\\[-0.3ex]shock", (-68.5, -21.8), TRANS[18], "north")
    lab += leader("lower triple point", (-56, -28.2), T2, "north")
    lab += leader("shear layers", (-52.5, 1.2), sl_a, "south", extra=[jg.at(0.12, -1)])
    lab += leader(r"jet bow\\[-0.3ex]shock", (-43.5, -5.8), 0.75 * jg.UPPER_END + 0.25 * jg.LOWER_END, "south")
    if detailed:
        lab += f"  \\node[dot] at {pt(jet_mid(0.52))} {{}};\n"
        lab += leader(r"jet, $M > 1$", (-36, -25.5), jet_mid(0.52), "north")
        lab += f"  \\node[region] at (-62,{f(zy1 - 2.6)}) {{$M < 1$}};\n"
        lab += f"  \\node[region] at (-41.5,{f(zy0 + 2.2)}) {{$M < 1$}};\n"
    else:
        lab += leader(r"supersonic jet", (-36, -25.5), jet_mid(0.52), "north")
    out.append(lab)
    out.append(f"\\end{{scope}}\n\\node[panel, anchor=south west] at ({f(X_RIGHT)},{f(HA + 1.2)}) {{(b)}};\n")
    return "".join(out)


# ======================================================================
#  (c) grid in the shock layer at the stagnation point
# ======================================================================
def panel_c():
    zx0, zx1, zy0, zy1 = WIN_C
    s = SC
    phi_lo, phi_hi = -10.0, 10.0
    grid, nr, nphi = grid_lines(2, phi_lo, phi_hi, 85.0, "gridC")
    yb = zy1 + 1.2
    return rf"""
% ======================================================================
%  (c) grid in the shock layer at the stagnation point
%      every 2nd grid line in each direction (exact node levels)
% ======================================================================
{shift_scope((X_RIGHT, 0), WIN_C, s)}  \begin{{scope}}
    \clip ({f(zx0)},{f(zy0, 4)}) rectangle ({f(zx1)},{f(zy1, 4)});
    \fill[bodygrey] (0,0) circle[radius={R}];
    % grid ({nr} circles, {nphi} rays shown)
{grid}    \draw[gridedge, line width=0.9pt] {arc_cmd(R, phi_hi, phi_lo)};
    % no-field bow shock
    \draw[shockref, line width=0.9pt, dash pattern=on 2.4pt off 1.3pt, line cap=butt]
      {poly(TRACE)};
  \end{{scope}}
  \draw[windowC, line width=0.7pt] ({f(zx0)},{f(zy0, 4)}) rectangle ({f(zx1)},{f(zy1, 4)});

  % uniform-spacing band
  \draw[dim] ({f(-R_UNI1)},{f(yb, 3)}) -- ({f(-R_UNI0)},{f(yb, 3)});
  \node[sml, anchor=south, fill=white, inner sep=0.6pt] at ({f(-0.5 * (R_UNI0 + R_UNI1))},{f(yb + 0.1, 3)})
    {{uniform $\Delta r$ = 0.2\,mm}};
  \node[sml, text=shockref, anchor=south west, fill=white, inner sep=0.8pt] at ({f(zx0 + 0.6)},{f(zy0 + 0.5, 3)})
    {{bow shock (reference)}};
  \node[sml, anchor=south east, fill=white, inner sep=0.8pt] at ({f(-R - 0.6)},{f(zy0 + 0.5, 3)})
    {{$\Delta r_w$ = 0.12\,$\mu$m}};
\end{{scope}}
\node[panel, anchor=south west] at ({f(X_RIGHT)},{f(H_C + 1.2)}) {{(c)}};
"""


def build(detailed):
    notes = [
        "Combined set-up figure: (a) domain, boundary conditions and grid;",
        "(b) Edney type IV interaction; (c) grid in the shock layer.",
        "",
        "Standalone LaTeX document.  On Overleaf: New Project > Blank Project, replace",
        "ALL of main.tex with this file (or upload it and set it as the main",
        "document) and compile with pdfLaTeX.  In the paper, include the compiled",
        "PDF with \\includegraphics[width=\\textwidth]{...}.",
        "",
        "Generated by make_figure.py.  Coordinates in mm, origin at the cylinder",
        "centre, x downstream, y up.",
        f"Case: M_inf = 8.03, beta = 18.1 deg -> theta = {THETA:.2f} deg, M_2 = {g.M_2:.3f} (gamma = 1.4);",
        f"R = {R} mm, inlet radius {RO:.2f} mm, outlets at phi = +/-{PHI_OUT:.0f} deg, no-slip wall",
        f"|phi| <= {PHI_NS:.0f} deg, slip wall beyond, inlet split phi = {PHI_SPLIT} deg.",
        f"Grid sa_mi_baseline.su2: {len(R_I) - 1} x {len(PHI_J) - 1} cells, polar O-grid, dphi = 0.25 deg,",
        f"first cell 0.12 um, ratio 1.15 to r = {R_WALL_END:.2f} mm, uniform 0.20 mm for",
        f"{R_UNI0:.2f} <= r <= {R_UNI1:.2f} mm, ratio 1.05 to the inlet.  Grid lines shown: (a) every",
        "4th circle and 8th ray, (c) every 2nd line; exact node levels.",
        f"Jet: transmitted shock {D(jg.SIGMA):.1f} deg -> jet direction {D(jg.DELTA0):.1f} deg, "
        f"M_jet = {jg.M_JET:.2f}, width {jg.WIDTH:.1f} mm;",
        "jet arch, jet bow shock stand-off and shock cells are schematic"
        + (" (solid: compressions, dashed: expansion fans)." if detailed else "."),
    ]
    head = "% " + "=" * 74 + "\n" + "".join(f"%  {n}\n" for n in notes) + "% " + "=" * 74 + "\n"
    body = panel_a() + panel_b(detailed) + panel_c()
    return (head + PREAMBLE + "\n\\begin{document}\n\\sffamily\\small\n"
            "\\begin{tikzpicture}[x=1mm, y=1mm, line cap=round, line join=round]\n"
            + body + "\\end{tikzpicture}\n\\end{document}\n")


if __name__ == "__main__":
    print(f"(a) {WA:.1f} x {HA:.1f} mm, scale {SA:.4f};  (b) {W_RIGHT:.1f} x {H_B:.1f} mm, "
          f"scale {SB:.3f}, window {tuple(round(v, 2) for v in WIN_B)};  (c) {W_RIGHT:.1f} x {H_C:.1f} mm, "
          f"scale {SC:.3f}, window {tuple(round(v, 2) for v in WIN_C)}")
    for name, det in (("fig_setup", True), ("fig_setup_clean", False)):
        tex = build(det)
        (HERE / f"{name}.tex").write_text(tex)
        print(f"wrote {name}.tex ({len(tex.splitlines())} lines)")
