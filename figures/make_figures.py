"""Write the standalone TikZ figures for the paper.

    fig_domain.tex             computational domain and boundary conditions
    fig_typeiv.tex             Edney type IV interaction, close-up (clean)
    fig_typeiv_detailed.tex    same, with the jet shock-cell structure
    fig_mesh.tex               computational grid (sa_mi_baseline.su2)

fig_domain and fig_typeiv* are drawn at their final print size (about 82 mm
wide) so two of them sit side by side in a 170 mm text width without
shrinking the text.  Every file is a complete standalone LaTeX document
(pdfLaTeX, e.g. Overleaf) with all coordinates written out, so no data files
are needed to compile it.  Units: mm, origin at the cylinder centre, x
downstream, y up; phi is measured from the stagnation point, positive up.

    python3 make_figures.py && for f in fig_*.tex; do pdflatex $f; done
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

FIG_W = 82.0                     # print width of fig_domain / fig_typeiv* (mm)


# ======================================================================
#  helpers
# ======================================================================
def tk(phi):
    """Polar angle from the stagnation point -> TikZ angle from +x."""
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
    """TikZ arc of radius r from phi0 to phi1 (angles from the stagnation point)."""
    a0, a1 = tk(phi0), tk(phi1)
    return f"({f(a0, 3)}:{f(r, nd)}) arc[start angle={f(a0, 3)}, end angle={f(a1, 3)}, radius={f(r, nd)}]"


_names = {}


def leader(text, at, target, anchor, frm=None, style="lbl", name=None, extra=()):
    """Label node plus a thin leader line from `frm` (default: the anchor) to target."""
    frm = frm or anchor
    if name is None:
        _names[text] = _names.get(text, 0) + 1
        name = "L" + "".join(c for c in text if c.isalnum())[:10] + str(_names[text])
    s = f"  \\node[{style}, anchor={anchor}] ({name}) at {pt(at)} {{{text}}};\n"
    s += f"  \\draw[leader] ({name}.{frm}) -- {pt(target)};\n"
    for t in extra:
        s += f"  \\draw[leader] ({name}.{frm}) -- {pt(t)};\n"
    return s


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
\definecolor{shockgrey}{HTML}{4A4A4A}
\definecolor{dimgrey}{HTML}{808080}
\definecolor{leadergrey}{HTML}{9A9A9A}
\definecolor{bodygrey}{HTML}{E4E4E4}
\definecolor{freeblue}{HTML}{1A5FA8}
\definecolor{postgreen}{HTML}{0E8A62}
\definecolor{slipochre}{HTML}{B06A10}
\definecolor{jetfill}{HTML}{FBEBD3}
\definecolor{zoomred}{HTML}{C0392B}

% ---- styles
\tikzset{
  lbl/.style={font=\sffamily\small, text=ink, inner sep=1.2pt, align=center},
  sml/.style={font=\sffamily\footnotesize, text=ink, inner sep=1.2pt, align=center},
  region/.style={font=\sffamily\small, text=dimgrey, inner sep=1pt},
  leader/.style={draw=leadergrey, line width=0.4pt},
  freeinlet/.style={draw=freeblue, line width=1.5pt},
  postinlet/.style={draw=postgreen, line width=1.5pt},
  outlet/.style={draw=dimgrey, line width=1.1pt, dash pattern=on 3.2pt off 2pt, line cap=butt},
  noslip/.style={draw=ink, line width=1.7pt},
  slip/.style={draw=slipochre, line width=1.7pt, dash pattern=on 2pt off 1.1pt, line cap=butt},
  shockA/.style={draw=shockgrey, line width=0.7pt},
  shockB/.style={draw=ink, line width=1.1pt},
  jetshock/.style={draw=ink, line width=1.1pt},
  shear/.style={draw=ink, line width=0.7pt, dash pattern=on 2.4pt off 1.4pt, line cap=butt},
  cellshock/.style={draw=ink, line width=0.45pt},
  fan/.style={draw=ink, line width=0.45pt, dash pattern=on 1.4pt off 0.9pt, line cap=butt},
  hatch/.style={draw=dimgrey, line width=0.35pt, line cap=butt},
  flow/.style={line width=0.9pt, -{Stealth[length=2.4mm, width=1.8mm]}},
  jetflow/.style={draw=ink, line width=0.75pt, -{Stealth[length=2mm, width=1.5mm]}},
  dim/.style={draw=dimgrey, line width=0.45pt,
              {Stealth[length=1.6mm, width=1.1mm]}-{Stealth[length=1.6mm, width=1.1mm]}},
  triple/.style={circle, fill=ink, inner sep=0pt, minimum size=1.6mm},
  dot/.style={circle, fill=leadergrey, inner sep=0pt, minimum size=1mm},
  panel/.style={font=\sffamily\bfseries\small, inner sep=0pt},
  % grid
  gridline/.style={draw=ink, line width=0.12pt},
  gridedge/.style={draw=ink, line width=0.6pt},
  axis/.style={draw=ink, line width=0.5pt},
  tick/.style={draw=ink, line width=0.4pt},
  ticklbl/.style={font=\sffamily\footnotesize, text=ink, inner sep=1.5pt},
  zoombox/.style={draw=zoomred, line width=0.6pt},
}
"""


def header(title, notes):
    return (f"% {'=' * 74}\n%  {title}\n%\n"
            "%  Standalone LaTeX document.  On Overleaf: New Project > Blank Project,\n"
            "%  replace ALL of main.tex with this file (or upload it and set it as the\n"
            "%  main document) and compile with pdfLaTeX.  In the paper, include the\n"
            "%  compiled PDF with \\includegraphics.\n%\n"
            "%  Generated by make_figures.py; coordinates in mm, origin at the cylinder\n"
            "%  centre, x downstream, y up.\n"
            + "".join(f"%  {n}\n" for n in notes) + f"% {'=' * 74}\n")


def document(title, notes, body, panel_label):
    lab = (f"% panel letter: change to \\newcommand{{\\panellabel}}{{}} to remove it\n"
           f"\\newcommand{{\\panellabel}}{{{panel_label}}}\n") if panel_label is not None else ""
    return (header(title, notes) + PREAMBLE + lab + "\n\\begin{document}\n\\sffamily\\small\n"
            "\\begin{tikzpicture}[x=1mm, y=1mm, line cap=round, line join=round]\n"
            + body + "\\end{tikzpicture}\n\\end{document}\n")


CASE_NOTES = [
    f"Case: M_inf = 8.03, beta = 18.1 deg -> theta = {THETA:.2f} deg, M_2 = {g.M_2:.3f} (gamma = 1.4).",
    f"R = {R} mm, inlet radius {RO:.2f} mm, outlets at phi = +/-{PHI_OUT:.0f} deg, no-slip wall",
    f"|phi| <= {PHI_NS:.0f} deg, slip wall beyond, inlet split phi = {PHI_SPLIT} deg.",
]


# ======================================================================
#  geometry shared by the figures
# ======================================================================
xt, yt = g.transmitted(40)
TRANS = np.c_[xt, yt]
JSHOCK = jg.jet_shock()
AFTER_U, AFTER_L = jg.after_shock(+1, 30), jg.after_shock(-1, 30)
SHOCKS, FANS = jg.wave_pattern()
JET_FILL = np.r_[jg.UPPER_FULL, jg.LOWER[::-1], TRANS[::-1][1:]]
ZOOM = (-85.0, -30.0, -42.0, 15.0)          # close-up window (fig_typeiv*)


def jet_mid(frac):
    return 0.5 * (jg.at(frac, +1) + jg.at(frac, -1))


def wedge_point(x):
    yu = np.interp(x, jg.UPPER_FULL[:, 0], jg.UPPER_FULL[:, 1])
    yv = np.interp(x, TRANS[:, 0], TRANS[:, 1])
    return np.array([x, 0.5 * (yu + yv)])


# ======================================================================
#  fig_domain
# ======================================================================
def fig_domain():
    x0, x1, y0, y1 = -226.0, 58.0, -153.0, 153.0
    s = FIG_W / (x1 - x0)
    up_a = g.upper_bow(np.linspace(g.y_t1, 165, 140)).T
    lo_a = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -175, 140)])
    inc_up = np.array([-222.0, g.incident(-222.0)])
    pd = -37.0
    q1 = wall(pd, RO)
    d = q1 / np.linalg.norm(q1)
    n = np.array([-d[1], d[0]])
    n = n if n[1] > 0 else -n
    lbl150 = 0.66 * q1 + 7.0 * n
    rot = np.degrees(np.arctan2(d[1], d[0]))
    rot = rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot
    zx0, zx1, zy0, zy1 = ZOOM

    body = rf"""
\begin{{scope}}[shift={{({f(-x0 * s)},{f(-y0 * s)})}}, x={s:.5f}mm, y={s:.5f}mm]

  % cylinder (dotted where it lies outside the domain)
  \fill[bodygrey] (0,0) circle[radius={R}];
  \draw[dimgrey, line width=0.45pt, dash pattern=on 0.4pt off 1.3pt]
    {arc_cmd(R, PHI_OUT, -PHI_OUT + 360)};

  % shock system (no-field solution), clipped to the domain
  \begin{{scope}}
    \clip {arc_cmd(RO, PHI_OUT, -PHI_OUT)} -- {arc_cmd(R, -PHI_OUT, PHI_OUT)[0:0]}({f(tk(-PHI_OUT), 3)}:{R})
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

  % inlet split and imposed shock angle
  \draw[dimgrey, line width=0.4pt, dash pattern=on 0.8pt off 1pt] {pt(SPLIT)} -- ++(50,0);
  \draw[ink, line width=0.45pt] {pt(SPLIT + [38, 0])} arc[start angle=0, end angle={f(BETA)}, radius=38];
  \node[lbl, anchor=north west] at {pt(SPLIT + [7, -3])} {{$\beta$ = 18.1$^\circ$}};
  \filldraw[fill=white, draw=ink, line width=0.75pt] {pt(SPLIT)} circle[radius=4.4];
  \node[sml, text=dimgrey, anchor=south east, align=right] at {pt(SPLIT + [-6, 8])}
    {{inlet split\\[-0.2ex]$\phi$ = $-$12.25$^\circ$}};

  % inflow
  \foreach \yy in {{106, 80}} {{\draw[flow, freeblue] (-222,\yy) -- ++(32,0);}}
  \foreach \yy in {{-82, -108}} {{\draw[flow, postgreen] (-222,\yy) -- ++({THETA:.2f}:32.8);}}
  \node[lbl, text=freeblue, anchor=south west, align=left] at (-223,120)
    {{freestream\\$M_\infty$ = 8.03}};
  \node[lbl, text=postgreen, anchor=north west, align=left] at (-223,-118)
    {{post-shock\\$M_2$ = 5.25}};

  % boundary labels
  \node[sml] (ns) at (2,16) {{no-slip wall\\$T_w$ = 294\,K}};
  \draw[leader] (ns.west) -- {pt(wall(15))};
  \node[lbl, anchor=west] (sw) at (16,66) {{slip wall}};
  \draw[leader] (sw.west) -- {pt(wall(57))};
  \node[lbl, text=dimgrey, anchor=south west] at {pt(wall(PHI_OUT, RO) + [3, -6])} {{outlet}};
  \node[lbl, text=dimgrey, anchor=north west] at {pt(wall(-PHI_OUT, RO) + [3, 6])} {{outlet}};
  \node[lbl, anchor=west] (bs) at (-100,98) {{bow shock}};
  \draw[leader] ([xshift=-1mm]bs.south) -- {pt(g.upper_bow(74.0))};

  % dimensions
  \draw[dim] (0,0) -- (-45:{R});
  \node[lbl, anchor=north] at {pt(wall(-135) + [9, -8])} {{$R$ = {R}\,mm}};
  \draw[dim] (0,0) -- ({f(tk(pd))}:{f(RO, 3)});
  \node[sml, text=dimgrey, rotate={rot:.2f}] at {pt(lbl150)} {{150\,mm}};

  % close-up window (other figure)
  \draw[ink, line width=0.45pt, dash pattern=on 1.6pt off 1.2pt, line cap=butt]
    ({f(zx0)},{f(zy0)}) rectangle ({f(zx1)},{f(zy1)});
\end{{scope}}

% panel letter
\node[panel, anchor=north west] at (0,{f((y1 - y0) * s + 2)}) {{\panellabel}};
"""
    return document("Computational domain and boundary conditions",
                    CASE_NOTES + ["Dashed box: window of fig_typeiv.tex / fig_typeiv_detailed.tex."],
                    body, "(a)"), s


# ======================================================================
#  fig_typeiv (clean and detailed)
# ======================================================================
def fig_typeiv(detailed):
    zx0, zx1, zy0, zy1 = ZOOM
    s = FIG_W / (zx1 - zx0)
    up_b = g.upper_bow(np.linspace(g.y_t1, zy1 + 6, 90)).T
    lo_b = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -80, 160)])
    lo_b = lo_b[lo_b[:, 1] > zy0 - 6]
    inc0 = np.array([zx0 - 4, g.incident(zx0 - 4)])

    if detailed:
        body_fill = [f"    \\fill[white] (0,0) circle[radius={R}];\n",
                     "    \\begin{scope}\n", f"      \\clip (0,0) circle[radius={R}];\n"]
        dd = np.array([np.cos(np.radians(-50)), np.sin(np.radians(-50))])
        for phi in np.arange(-60, 60.01, 2.6):
            p = wall(phi)
            body_fill.append(f"      \\draw[hatch] {pt(p - 0.2 * dd)} -- {pt(p + 3.0 * dd)};\n")
        body_fill.append("    \\end{scope}\n")
        body_fill = "".join(body_fill)
        jetfill = ""
    else:
        body_fill = f"    \\fill[bodygrey] (0,0) circle[radius={R}];\n"
        jetfill = ("    % supersonic jet\n    \\fill[jetfill]\n      "
                   + poly(JET_FILL, 2) + " -- cycle;\n")

    out = [rf"""
\begin{{scope}}[shift={{({f(-zx0 * s)},{f(-zy0 * s)})}}, x={s:.5f}mm, y={s:.5f}mm]
  \begin{{scope}}
    \clip ({f(zx0)},{f(zy0)}) rectangle ({f(zx1)},{f(zy1)});
    % body and wall
{body_fill}    \draw[noslip, line width=2pt] {arc_cmd(R, PHI_NS, -PHI_NS)};
    \draw[slip, line width=2pt] {arc_cmd(R, PHI_OUT, PHI_NS)};
    \draw[slip, line width=2pt] {arc_cmd(R, -PHI_NS, -PHI_OUT)};
{jetfill}
    % bow shock (upper and lower branches) and incident shock
    \draw[shockB]
      {poly(up_b)};
    \draw[shockB]
      {poly(lo_b, 2)};
    \draw[shockB] {pt(inc0)} -- {pt(T1)};
  \end{{scope}}
"""]
    if detailed:
        out.append("\n  % jet shock cells: compressions (solid), expansion fans (dashed)\n")
        out += [f"  \\draw[cellshock] {pt(a)} -- {pt(b)};\n" for a, b in SHOCKS]
        out += [f"  \\draw[fan] {pt(a)} -- {pt(b)};\n" for a, b in FANS]
    a0 = wedge_point(-61.0)
    out.append(rf"""
  % transmitted shock
  \draw[shockB]
    {poly(TRANS)};

  % shear layers, continued past the jet bow shock along the wall
  \draw[shear]
    {poly(jg.UPPER_FULL, 4)};
  \draw[shear]
    {poly(jg.LOWER, 4)};
  \draw[shear]
    {poly(AFTER_U, 3)};
  \draw[shear]
    {poly(AFTER_L, 3)};

  % jet bow shock
  \draw[jetshock] {pt(JSHOCK[0])} -- {pt(JSHOCK[1])};

  % triple points
  \node[triple] at {pt(T1)} {{}};
  \node[triple] at {pt(T2)} {{}};

  % jet direction
  \draw[jetflow] {pt(a0)} -- {pt(a0 + 4.4 * jg.D0)};

  % oncoming flow
  \draw[flow, freeblue] (-84,9.5) -- ++(7,0);
  \draw[flow, postgreen] (-84,-37.5) -- ++({THETA:.2f}:7.1);

  % scale bar
  \draw[ink, line width=0.8pt, line cap=butt] (-50,11.5) -- (-40,11.5);
  \draw[ink, line width=0.6pt] (-50,11) -- (-50,12) (-40,11) -- (-40,12);
  \node[sml, anchor=north] at (-45,10.8) {{10\,mm}};

  % labels
""")
    sl_a = 0.5 * (T1 + jg.UPPER[0])
    lab = ""
    lab += f"  \\node[lbl, text=freeblue, anchor=south west] at (-84.6,10.6) {{$M_\\infty$ = 8.03}};\n"
    lab += f"  \\node[lbl, text=postgreen, anchor=south west] at (-84.6,-35.6) {{$M_2$ = 5.25}};\n"
    lab += leader("bow shock", (-67.5, 2.5), g.upper_bow(2.5), "west",
                  extra=[] if not detailed else [])
    lab += leader(r"upper\\[-0.3ex]triple point", (-66, -4.2), T1, "west")
    lab += leader("incident\\\\[-0.3ex]shock", (-79.5, -16.5), (-81, g.incident(-81)), "north")
    lab += leader("transmitted\\\\[-0.3ex]shock", (-64.5, -22), TRANS[20], "north")
    lab += leader(r"lower\\[-0.3ex]triple point", (-58.5, -29.5), T2, "north")
    lab += leader("shear layers", (-50.5, -3.2), sl_a, "south",
                  extra=[jg.at(0.12, -1)])
    lab += leader(r"jet bow\\[-0.3ex]shock", (-39.5, -26), 0.4 * jg.UPPER_END + 0.6 * jg.LOWER_END, "north")
    if detailed:
        lab += r"""
  % Mach-number regions
  \node[region] at (-60,8) {$M < 1$};
  \node[region] at (-34,-37.5) {$M < 1$};
"""
        lab += f"  \\node[dot] at {pt(jet_mid(0.52))} {{}};\n"
        lab += leader(r"jet, $M > 1$", (-44.5, 3.2), jet_mid(0.52), "south")
    else:
        lab += leader(r"supersonic\\[-0.3ex]jet", (-44.5, 3.2), jet_mid(0.52), "south")
    out.append(lab)
    out.append(f"""\\end{{scope}}

% panel letter
\\node[panel, anchor=north west] at (0,{f((zy1 - zy0) * s + 2)}) {{\\panellabel}};
""")
    notes = CASE_NOTES + [
        f"Jet: transmitted shock {D(jg.SIGMA):.1f} deg -> jet direction {D(jg.DELTA0):.1f} deg,",
        f"M_jet = {jg.M_JET:.2f}, width {jg.WIDTH:.1f} mm (oblique-shock relations).",
        "Jet arch, jet bow shock stand-off and shock cells are schematic.",
    ]
    if detailed:
        notes.append("Shock cells: compressions solid, expansion fans dashed (state this in the caption).")
    title = "Edney type IV interaction, close-up" + (" with jet shock-cell structure" if detailed else "")
    return document(title, notes, "".join(out), "(b)"), s


# ======================================================================
#  fig_mesh
# ======================================================================
def load_grid():
    r = np.loadtxt(HERE / "data" / "mesh_r_i.txt")
    phi = np.loadtxt(HERE / "data" / "mesh_phi_j.txt")
    return r, phi


def mesh_panel(name, win, origin, scale, r, phi, stride_r, stride_phi, ticks, letter,
               zoom=None, extra=""):
    """One mesh panel: grid lines clipped to the window, frame, ticks, labels."""
    x0, x1, y0, y1 = win
    ox, oy = origin
    W, Hh = (x1 - x0) * scale, (y1 - y0) * scale
    # angular range of the window (from the stagnation point), padded
    corners = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]])
    ang = np.degrees(np.arctan2(corners[:, 1], -corners[:, 0]))
    rmin_w = max(R, np.hypot(np.clip(0, x0, x1), np.clip(0, y0, y1)))
    rmax_w = np.hypot(corners[:, 0], corners[:, 1]).max()
    pmin, pmax = max(ang.min() - 2, -PHI_OUT), min(ang.max() + 2, PHI_OUT)
    ri = np.unique(np.r_[r[::stride_r], r[-1]])
    pj = np.unique(np.r_[phi[::stride_phi], phi[0], phi[-1]])
    ri_in = ri[(ri >= rmin_w - 1e-9) & (ri <= rmax_w + 1e-9)]
    pj_in = pj[(pj >= pmin - 1e-9) & (pj <= pmax + 1e-9)]
    rlo, rhi = R, min(rmax_w, r[-1])
    nd = 5 if scale > 10 else 3

    s = [f"\n% ---------------- ({letter}) {name}\n",
         f"\\begin{{scope}}[shift={{({f(ox - x0 * scale, 3)},{f(oy - y0 * scale, 3)})}}, "
         f"x={scale:.5f}mm, y={scale:.5f}mm]\n",
         "  \\begin{scope}\n",
         f"    \\clip ({f(x0, nd)},{f(y0, nd)}) rectangle ({f(x1, nd)},{f(y1, nd)});\n",
         f"    \\fill[bodygrey] (0,0) circle[radius={R}];\n",
         "    % circumferential grid lines (constant r)\n"]
    for rr in ri_in:
        s.append(f"    \\draw[gridline] {arc_cmd(rr, pmin, pmax, nd=6)};\n")
    s.append("    % radial grid lines (constant phi)\n")
    for p in pj_in:
        a = tk(p)
        s.append(f"    \\draw[gridline] ({f(a, 3)}:{f(rlo, 4)}) -- ({f(a, 3)}:{f(rhi, 4)});\n")
    s.append("    % domain boundaries\n")
    s.append(f"    \\draw[gridedge, line width=0.8pt] {arc_cmd(R, PHI_OUT, -PHI_OUT)};\n")
    s.append(f"    \\draw[gridedge] {arc_cmd(r[-1], PHI_OUT, -PHI_OUT, nd=4)};\n")
    for p in (PHI_OUT, -PHI_OUT):
        s.append(f"    \\draw[gridedge] ({f(tk(p), 3)}:{R}) -- ({f(tk(p), 3)}:{f(r[-1], 4)});\n")
    s.append(extra)
    s.append("  \\end{scope}\n")
    if zoom is not None:
        zx0, zx1, zy0, zy1, zl = zoom
        s.append(f"  \\draw[zoombox] ({f(zx0, nd)},{f(zy0, nd)}) rectangle ({f(zx1, nd)},{f(zy1, nd)});\n")
        s.append(f"  \\node[sml, text=zoomred, fill=white, anchor=south west, inner sep=0.8pt, "
                 f"font=\\sffamily\\footnotesize\\bfseries] at ({f(zx0, nd)},{f(zy1, nd)}) {{{zl}}};\n")
    s.append("\\end{scope}\n")

    # frame, ticks and tick labels (in page coordinates)
    s.append(f"\\draw[axis] ({f(ox, 3)},{f(oy, 3)}) rectangle ({f(ox + W, 3)},{f(oy + Hh, 3)});\n")
    xt_, yt_ = ticks
    for v in xt_:
        X = ox + (v - x0) * scale
        s.append(f"\\draw[tick] ({f(X, 3)},{f(oy, 3)}) -- ++(0,1.2) ({f(X, 3)},{f(oy + Hh, 3)}) -- ++(0,-1.2);\n")
        s.append(f"\\node[ticklbl, anchor=north] at ({f(X, 3)},{f(oy, 3)}) {{{fmt_tick(v)}}};\n")
    for v in yt_:
        Y = oy + (v - y0) * scale
        s.append(f"\\draw[tick] ({f(ox, 3)},{f(Y, 3)}) -- ++(1.2,0) ({f(ox + W, 3)},{f(Y, 3)}) -- ++(-1.2,0);\n")
        s.append(f"\\node[ticklbl, anchor=east] at ({f(ox, 3)},{f(Y, 3)}) {{{fmt_tick(v)}}};\n")
    s.append(f"\\node[ticklbl, anchor=north] at ({f(ox + W / 2, 3)},{f(oy - 4.2, 3)}) {{$x$ (mm)}};\n")
    s.append(f"\\node[ticklbl, anchor=south, rotate=90] at ({f(ox - 8.2, 3)},{f(oy + Hh / 2, 3)}) {{$y$ (mm)}};\n")
    s.append(f"\\node[panel, anchor=south west] at ({f(ox, 3)},{f(oy + Hh + 1.2, 3)}) {{({letter})}};\n")
    return "".join(s), (ox, oy, W, Hh)


def fmt_tick(v):
    """Tick label in the text font, with a proper minus sign."""
    s = f"{abs(v):g}"
    return ("$-$" + s) if v < 0 else s


def fig_mesh():
    r, phi = load_grid()
    r_wall_end, r_uni0, r_uni1 = r[53], r[54], r[233]
    HT = 100.0                                  # height of panel (a) (mm)

    # (a) full grid, every 4th line
    win_a = (-158.0, -10.0, -142.0, 142.0)
    sa = HT / (win_a[3] - win_a[2])
    win_b = (-82.0, -34.0, -20.0, 4.0)
    pa, (oxa, oya, Wa, Ha) = mesh_panel(
        "full grid (every 4th line)", win_a, (11.0, 9.0), sa, r, phi, 4, 4,
        ([-150, -100, -50], [-100, -50, 0, 50, 100]), "a", zoom=(*win_b, "(b)"))

    # (b) shock layer around the stagnation point, every 2nd line
    xb0 = oxa + Wa + 19.0
    WB = 168.0 - xb0
    sb = WB / (win_b[1] - win_b[0])
    HB = (win_b[3] - win_b[2]) * sb
    HC = (oya + HT) - (WB / (win_b[1] - win_b[0])) * (win_b[3] - win_b[2]) - 15.0 - oya
    sc = HC / 1.6
    span_c = (WB - 44.0) / sc                    # leave room for the notes beside (c)
    win_c = (-37.85 - span_c, -37.85, -0.8, 0.8)
    trace = np.loadtxt(HERE / "data" / "bow_shock_no_field.csv", delimiter=",")
    extra_b = ("    % no-field bow shock (for reference)\n"
               "    \\draw[zoomred, line width=0.7pt, dash pattern=on 2.2pt off 1.2pt, line cap=butt]\n      "
               + poly(trace) + ";\n")
    pb, (oxb, oyb, Wb, Hb) = mesh_panel(
        "shock layer at the stagnation point (every 2nd line)", win_b,
        (xb0, oya + Ha - HB), sb, r, phi, 2, 2,
        ([-80, -70, -60, -50, -40], [-20, -10, 0]), "b",
        zoom=(*win_c, "(c)"), extra=extra_b)

    # (c) near-wall cells at the stagnation point, every line
    pc, (oxc, oyc, Wc, Hc) = mesh_panel(
        "near-wall grid at the stagnation point (every line)", win_c,
        (xb0, oya), sc, r, phi, 1, 1,
        ([-39.5, -39, -38.5, -38], [-0.5, 0, 0.5]), "c")

    # spacing zones: brackets above (b), measured along y = 0 (r = -x)
    def X(v):
        return oxb + (v - win_b[0]) * sb
    yb = oyb + Hb + 2.2
    ann = rf"""
% ---------------- annotations
\draw[dim] ({f(X(-r_uni1), 3)},{f(yb, 3)}) -- ({f(X(-r_uni0), 3)},{f(yb, 3)});
\node[sml, anchor=south] at ({f(X(-0.5 * (r_uni0 + r_uni1)), 3)},{f(yb, 3)})
  {{uniform $\Delta r$ = 0.20\,mm ($r$ = {r_uni0:.1f}--{r_uni1:.1f}\,mm)}};
\node[panel, anchor=south west] at ({f(oxb, 3)},{f(yb + 0.4, 3)}) {{\phantom{{(b)}}}};
\node[sml, anchor=north west, align=left] at ({f(oxc + Wc + 4, 3)},{f(oyc + Hc, 3)})
  {{{len(r) - 1} $\times$ {len(phi) - 1} cells\\[0.6ex]
    $\Delta\phi$ = 0.25$^\circ$\\[0.6ex]
    first cell 0.12\,$\mu$m,\\ ratio 1.15 to $r$ = {r_wall_end:.2f}\,mm\\[0.6ex]
    uniform 0.20\,mm to\\ $r$ = {r_uni1:.1f}\,mm\\[0.6ex]
    ratio 1.05 to the inlet\\[0.6ex]
    {{\color{{zoomred}}-\,-\,-}} bow shock (no field)}};
"""
    body = pa + pb + pc + ann
    notes = [f"Grid: {len(r) - 1} x {len(phi) - 1} cells (radial x circumferential), polar O-grid,",
             f"Delta phi = 0.25 deg; radial: first cell 0.12 um, ratio 1.15 to r = {r_wall_end:.2f} mm,",
             f"uniform 0.20 mm for {r_uni0:.2f} <= r <= {r_uni1:.2f} mm, ratio 1.05 to the inlet.",
             "Lines drawn: (a) every 4th, (b) every 2nd, (c) every grid line.",
             "Red dashed line in (b): no-field bow shock (data/bow_shock_no_field.csv)."]
    return document("Computational grid (sa_mi_baseline.su2)", notes, body, None)


# ======================================================================
if __name__ == "__main__":
    for name, (tex, *_) in {
        "fig_domain": fig_domain(),
        "fig_typeiv": fig_typeiv(False),
        "fig_typeiv_detailed": fig_typeiv(True),
    }.items():
        (HERE / f"{name}.tex").write_text(tex)
        print(f"wrote {name}.tex ({len(tex.splitlines())} lines)")
    tex = fig_mesh()
    (HERE / "fig_mesh.tex").write_text(tex)
    print(f"wrote fig_mesh.tex ({len(tex.splitlines())} lines)")
