"""Write the two TikZ versions of the set-up / type IV figure.

    fig_setup_typeiv_tikz.tex           (a) domain + (b) type IV close-up
    fig_setup_typeiv_detailed_tikz.tex  same, with the full jet structure
                                        (shock cells, expansion fans, Mach regions)

Both are self-contained standalone LaTeX files (pdfLaTeX, e.g. Overleaf); this
script only converts the geometry in fig_setup_typeiv.py / jet_geometry.py into
literal TikZ coordinates, so the .tex files need no data.  Units: mm, origin
at the cylinder centre, x downstream, y up.

    python3 make_tikz.py
    pdflatex fig_setup_typeiv_tikz.tex
    pdflatex fig_setup_typeiv_detailed_tikz.tex
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

# panel (a) layout: data x in [-218, 52], y in [-152, 152]
SA = 0.30
A_ORIGIN = (218 * SA, 152 * SA)
A_HEIGHT = 304 * SA
B_LEFT, B_WIDTH = 94.0, 77.0          # panel (b) frame on the page (mm)


# ----------------------------------------------------------------- helpers
def tk(phi):
    """Polar angle from the stagnation point -> TikZ angle from +x."""
    return 180.0 - phi


def f(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def pt(p):
    return f"({f(p[0])},{f(p[1])})"


def poly(xy, step=1, indent="      "):
    xy = np.asarray(xy, dtype=float)
    if xy.ndim == 2 and xy.shape[0] == 2 and xy.shape[1] != 2:
        xy = xy.T
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


# ----------------------------------------------------------------- geometry
up_a = g.upper_bow(np.linspace(g.y_t1, 160, 140)).T
lo_a = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -170, 140)])
xt, yt = g.transmitted(40)
TRANS = np.c_[xt, yt]
INC_UP = np.array([-214.0, g.incident(-214.0)])

# 150 mm dimension
PD = -37.0
Q1 = wall(PD, RO)
_d = Q1 / np.linalg.norm(Q1)
_n = np.array([-_d[1], _d[0]])
_n = _n if _n[1] > 0 else -_n
LBL150 = 0.70 * Q1 + 5.5 * _n
ROT150 = np.degrees(np.arctan2(_d[1], _d[0]))
ROT150 = ROT150 - 180 if ROT150 > 90 else ROT150 + 180 if ROT150 < -90 else ROT150

JET_FILL = np.r_[jg.UPPER_FULL, jg.LOWER[::-1], TRANS[::-1][1:]]
JSHOCK = jg.jet_shock()                     # straight, nearly normal to the jet
AFTER_U, AFTER_L = jg.after_shock(+1, 30), jg.after_shock(-1, 30)
SHOCKS, FANS = jg.wave_pattern()


def jet_mid(frac):
    return 0.5 * (jg.at(frac, +1) + jg.at(frac, -1))


# ----------------------------------------------------------------- preamble
PREAMBLE = r"""\documentclass[border=1.5pt]{standalone}
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
\definecolor{dimgrey}{HTML}{8A8A8A}
\definecolor{leadergrey}{HTML}{A0A0A0}
\definecolor{bodygrey}{HTML}{E4E4E4}
\definecolor{freeblue}{HTML}{1A5FA8}
\definecolor{postgreen}{HTML}{0E8A62}
\definecolor{slipochre}{HTML}{B06A10}
\definecolor{jetfill}{HTML}{FBEBD3}
\definecolor{jetedge}{HTML}{8C5A14}

% ---- styles
\tikzset{
  lbl/.style={font=\sffamily\footnotesize, text=ink, inner sep=1.2pt, align=center},
  small/.style={font=\sffamily\scriptsize},
  region/.style={font=\sffamily\footnotesize, text=dimgrey, inner sep=1pt},
  leader/.style={draw=leadergrey, line width=0.35pt},
  freeinlet/.style={draw=freeblue, line width=1.5pt},
  postinlet/.style={draw=postgreen, line width=1.5pt},
  outlet/.style={draw=dimgrey, line width=1.1pt, dash pattern=on 3.2pt off 2pt, line cap=butt},
  noslip/.style={draw=ink, line width=1.7pt},
  slip/.style={draw=slipochre, line width=1.7pt, dash pattern=on 2pt off 1.1pt, line cap=butt},
  shockA/.style={draw=shockgrey, line width=0.7pt},
  shockB/.style={draw=ink, line width=1.15pt},
  jetshock/.style={draw=ink, line width=1.0pt},
  shear/.style={draw=ink, line width=0.7pt, dash pattern=on 2.4pt off 1.4pt, line cap=butt},
  cellshock/.style={draw=ink, line width=0.4pt},
  fan/.style={draw=ink, line width=0.4pt, dash pattern=on 1.3pt off 0.8pt, line cap=butt},
  hatch/.style={draw=dimgrey, line width=0.35pt, line cap=butt},
  flow/.style={line width=0.9pt, -{Stealth[length=2.6mm, width=1.9mm]}},
  jetflow/.style={draw=ink, line width=0.7pt, -{Stealth[length=1.9mm, width=1.4mm]}},
  dim/.style={draw=dimgrey, line width=0.45pt,
              {Stealth[length=1.6mm, width=1.1mm]}-{Stealth[length=1.6mm, width=1.1mm]}},
  triple/.style={circle, fill=ink, inner sep=0pt, minimum size=1.5mm},
  dot/.style={circle, fill=leadergrey, inner sep=0pt, minimum size=0.9mm},
}
"""


def header(title, extra):
    return rf"""% ==========================================================================
%  {title}
%  (a) computational domain and boundary conditions,
%  (b) close-up of the Edney type IV interaction (dashed box in (a)).
%
%  Standalone LaTeX file.  On Overleaf: New Project > Blank Project, replace
%  the whole of main.tex with this file (or upload it and set it as the main
%  document), compiler pdfLaTeX.  To use the figure in a paper, compile this
%  file and include the PDF with \includegraphics.
%
%  Generated by make_tikz.py.  Coordinates in mm, origin at the cylinder
%  centre, x downstream, y up; phi measured from the stagnation point.
%  Shock positions: no-field shock trace (data/bow_shock_no_field.csv) and its
%  Billig-type far-field fit.
%  Case: M_inf = 8.03, beta = 18.1 deg -> theta = {THETA:.2f} deg, M_2 = {g.M_2:.3f};
%        transmitted shock {D(jg.SIGMA):.1f} deg -> jet direction {D(jg.DELTA0):.1f} deg,
%        M_jet = {jg.M_JET:.2f}, jet width {jg.WIDTH:.1f} mm  (oblique-shock relations, gamma = 1.4).
%        R = {R} mm, inlet radius {RO:.0f} mm, outlets at phi = +/-{PHI_OUT:.0f} deg,
%        no-slip wall |phi| <= {PHI_NS:.0f} deg, slip wall beyond, inlet split phi = {PHI_SPLIT} deg.
{extra}% ==========================================================================
"""


# ----------------------------------------------------------------- panel (a)
def panel_a(zoom):
    zx0, zx1, zy0, zy1 = zoom
    return rf"""
% ==========================================================================
%  (a)  computational domain              scale 1 : {1 / SA:.2f}
% ==========================================================================
\begin{{scope}}[shift={{({f(A_ORIGIN[0])},{f(A_ORIGIN[1])})}}, x={SA}mm, y={SA}mm]

  % cylinder (dotted outline where it lies outside the domain)
  \fill[bodygrey] (0,0) circle[radius={R}];
  \draw[dimgrey, line width=0.45pt, dash pattern=on 0.4pt off 1.3pt]
    ({tk(PHI_OUT)}:{R}) arc[start angle={tk(PHI_OUT)}, end angle={tk(-PHI_OUT) - 360}, radius={R}];

  % shock system, clipped to the domain
  \begin{{scope}}
    \clip ({tk(PHI_OUT)}:{RO}) arc[start angle={tk(PHI_OUT)}, end angle={tk(-PHI_OUT)}, radius={RO}]
      -- ({tk(-PHI_OUT)}:{R}) arc[start angle={tk(-PHI_OUT)}, end angle={tk(PHI_OUT)}, radius={R}] -- cycle;
    \draw[shockA]
      {poly(up_a, 2)};
    \draw[shockA]
      {poly(lo_a, 2)};
    \draw[shockA]
      {poly(TRANS, 2)};
    \draw[shockA] {pt(SPLIT)} -- {pt(T1)};
  \end{{scope}}

  % imposed oblique shock upstream of the inlet
  \draw[shockA, dash pattern=on 2.4pt off 1.6pt, line cap=butt] {pt(INC_UP)} -- {pt(SPLIT)};

  % inlets and outlets
  \draw[freeinlet] ({tk(PHI_SPLIT)}:{RO}) arc[start angle={tk(PHI_SPLIT)}, end angle={tk(PHI_OUT)}, radius={RO}];
  \draw[postinlet] ({tk(-PHI_OUT)}:{RO}) arc[start angle={tk(-PHI_OUT)}, end angle={tk(PHI_SPLIT)}, radius={RO}];
  \draw[outlet] ({tk(PHI_OUT)}:{R}) -- ({tk(PHI_OUT)}:{RO});
  \draw[outlet] ({tk(-PHI_OUT)}:{R}) -- ({tk(-PHI_OUT)}:{RO});

  % walls: no-slip for |phi| <= {PHI_NS:.0f} deg, slip for {PHI_NS:.0f} < |phi| <= {PHI_OUT:.0f} deg
  \draw[slip]   ({tk(PHI_OUT)}:{R}) arc[start angle={tk(PHI_OUT)}, end angle={tk(PHI_NS)}, radius={R}];
  \draw[slip]   ({tk(-PHI_NS)}:{R}) arc[start angle={tk(-PHI_NS)}, end angle={tk(-PHI_OUT)}, radius={R}];
  \draw[noslip] ({tk(PHI_NS)}:{R}) arc[start angle={tk(PHI_NS)}, end angle={tk(-PHI_NS)}, radius={R}];
  \draw[ink, line width=0.5pt] (-3.5,0) -- (3.5,0) (0,-3.5) -- (0,3.5);

  % inlet split and shock angle
  \draw[dimgrey, line width=0.4pt, dash pattern=on 0.8pt off 1pt] {pt(SPLIT)} -- ++(48,0);
  \draw[ink, line width=0.45pt] {pt(SPLIT + [36, 0])} arc[start angle=0, end angle={BETA}, radius=36];
  \node[lbl, anchor=north west] at {pt(SPLIT + [6, -2.5])} {{$\beta$ = 18.1$^\circ$}};
  \filldraw[fill=white, draw=ink, line width=0.75pt] {pt(SPLIT)} circle[radius=4.2];
  \node[lbl, small, text=dimgrey, anchor=south east, align=right] at {pt(SPLIT + [-6, 7])}
    {{inlet split\\[-0.3ex]$\phi$ = $-$12.25$^\circ$}};

  % inflow
  \foreach \yy in {{116, 94, 72}} {{\draw[flow, freeblue] (-214,\yy) -- ++(30,0);}}
  \foreach \yy in {{-72, -94, -116}} {{\draw[flow, postgreen] (-214,\yy) -- ++({THETA:.2f}:30.6);}}
  \node[lbl, text=freeblue, anchor=south west, align=left] at (-215,127)
    {{freestream inlet\\$M_\infty$ = 8.03}};
  \node[lbl, text=postgreen, anchor=north west, align=left] at (-215,-122)
    {{post-shock inlet\\$M_2$ = 5.25}};

  % boundary labels
  \node[lbl] (ns) at (5,15) {{no-slip wall\\$T_w$ = 294\,K}};
  \draw[leader] (ns.west) -- {pt(wall(15))};
  \node[lbl, anchor=west] (sw) at (24,64) {{slip wall}};
  \draw[leader] (sw.west) -- {pt(wall(56))};
  \node[lbl, text=dimgrey, anchor=west] (o1) at (-2,128) {{outlet}};
  \draw[leader] (o1.west) -- ({tk(PHI_OUT)}:{0.8 * RO:.0f});
  \node[lbl, text=dimgrey, anchor=west] (o2) at (-2,-128) {{outlet}};
  \draw[leader] (o2.west) -- ({tk(-PHI_OUT)}:{0.8 * RO:.0f});
  \node[lbl, anchor=west] (bs) at (-98,100) {{bow shock}};
  \draw[leader] ([xshift=-6mm]bs.south east) -- {pt(g.upper_bow(72.0))};

  % dimensions
  \draw[dim] (0,0) -- (-45:{R});
  \node[lbl, anchor=north] at {pt(wall(-135) + [6, -4.5])} {{$R$ = {R}\,mm}};
  \draw[dim] (0,0) -- ({tk(PD)}:{RO});
  \node[lbl, small, text=dimgrey, rotate={ROT150:.2f}] at {pt(LBL150)} {{{RO:.0f}\,mm}};

  % window shown in (b)
  \draw[ink, line width=0.4pt, dash pattern=on 1.6pt off 1.2pt, line cap=butt]
    ({f(zx0)},{f(zy0)}) rectangle ({f(zx1)},{f(zy1)});
  \node[lbl, small, anchor=south west, inner sep=1.2pt] at ({f(zx0)},{f(zy1)}) {{(b)}};
\end{{scope}}
"""


# ----------------------------------------------------------------- panel (b)
def panel_b(zoom, labels, detailed):
    zx0, zx1, zy0, zy1 = zoom
    sb = min(B_WIDTH / (zx1 - zx0), A_HEIGHT / (zy1 - zy0))
    ox = B_LEFT - zx0 * sb
    oy = -zy0 * sb + (A_HEIGHT - (zy1 - zy0) * sb) / 2
    up_b = g.upper_bow(np.linspace(g.y_t1, zy1 + 6, 90)).T
    lo_b = np.array([g.lower_bow(v) for v in np.linspace(g.v_t2, -80, 160)])
    lo_b = lo_b[lo_b[:, 1] > zy0 - 6]
    inc0 = np.array([zx0 - 4, g.incident(zx0 - 4)])

    out = [rf"""
% ==========================================================================
%  (b)  type IV interaction               scale {sb:.3f} : 1
% ==========================================================================
\begin{{scope}}[shift={{({f(ox)},{f(oy)})}}, x={sb:.4f}mm, y={sb:.4f}mm]
  \begin{{scope}}
    \clip ({f(zx0)},{f(zy0)}) rectangle ({f(zx1)},{f(zy1)});

    % body and wall
{body_b(detailed)}    \draw[noslip, line width=2pt] ({tk(PHI_NS)}:{R}) arc[start angle={tk(PHI_NS)}, end angle={tk(-PHI_NS)}, radius={R}];
    \draw[slip, line width=2pt] ({tk(PHI_OUT)}:{R}) arc[start angle={tk(PHI_OUT)}, end angle={tk(PHI_NS)}, radius={R}];
    \draw[slip, line width=2pt] ({tk(-PHI_NS)}:{R}) arc[start angle={tk(-PHI_NS)}, end angle={tk(-PHI_OUT)}, radius={R}];

{jet_fill(detailed)}
    % bow shock above and below the interaction, incident shock
    \draw[shockB]
      {poly(up_b)};
    \draw[shockB]
      {poly(lo_b, 2)};
    \draw[shockB] {pt(inc0)} -- {pt(T1)};
  \end{{scope}}
"""]
    if detailed:
        out.append("\n  % shock cells inside the jet: compressions (solid), reflected expansion fans (dashed)\n")
        for a, b in SHOCKS:
            out.append(f"  \\draw[cellshock] {pt(a)} -- {pt(b)};\n")
        for a, b in FANS:
            out.append(f"  \\draw[fan] {pt(a)} -- {pt(b)};\n")
    out.append(rf"""
  % transmitted shock
  \draw[shockB]
    {poly(TRANS)};

  % shear layers bounding the jet, continued past the jet bow shock to the wall
  \draw[shear]
    {poly(jg.UPPER_FULL, 4)};
  \draw[shear]
    {poly(jg.LOWER, 4)};
  \draw[shear]
    {poly(AFTER_U, 3)};
  \draw[shear]
    {poly(AFTER_L, 3)};

  % jet bow shock
  \draw[jetshock]
    {poly(JSHOCK, 2)};

  % triple points
  \node[triple] (T1) at {pt(T1)} {{}};
  \node[triple] (T2) at {pt(T2)} {{}};
""")
    out.append(labels(sb))
    out.append("\\end{scope}\n")
    return "".join(out), sb


def body_b(detailed):
    if not detailed:
        return f"    \\fill[bodygrey] (0,0) circle[radius={R}];\n"
    # white body with hatching along the surface, as in the classical sketch
    out = [f"    \\fill[white] (0,0) circle[radius={R}];\n",
           "    \\begin{scope}\n",
           f"      \\clip (0,0) circle[radius={R}];\n"]
    d = np.array([np.cos(np.radians(-50)), np.sin(np.radians(-50))])
    for phi in np.arange(-60, 60.01, 2.4):
        p = wall(phi)
        out.append(f"      \\draw[hatch] {pt(p - 0.2 * d)} -- {pt(p + 3.2 * d)};\n")
    out.append("    \\end{scope}\n")
    return "".join(out)


def jet_fill(detailed):
    if detailed:
        return ""
    return ("    % supersonic jet (light fill)\n    \\fill[jetfill]\n      "
            + poly(JET_FILL, 2) + " -- cycle;\n")


def leader(text, at, target, anchor, frm=None, style="lbl", name=None):
    frm = frm or {"east": "east", "west": "west", "north": "north", "south": "south",
                  "north east": "north east", "north west": "north west",
                  "south east": "south east", "south west": "south west",
                  "center": "center"}[anchor]
    nm = name or f"n{abs(hash((text, at[0], at[1]))) % 10**8}"
    return (f"  \\node[{style}, anchor={anchor}] ({nm}) at {pt(at)} {{{text}}};\n"
            f"  \\draw[leader] ({nm}.{frm}) -- {pt(target)};\n")


# ---- labels (both versions use the same close-up window)
ZOOM = (-82.0, -30.0, -40.0, 12.0)


def wedge_point(x):
    """Point midway between the upper shear layer and the transmitted shock at x."""
    yu = np.interp(x, jg.UPPER_FULL[:, 0], jg.UPPER_FULL[:, 1])
    yt = np.interp(x, TRANS[:, 0], TRANS[:, 1])
    return np.array([x, 0.5 * (yu + yt)])


def common_labels(detailed):
    a0 = wedge_point(-61.5)
    s = rf"""
  % jet direction (upstream of the shock cells)
  \draw[jetflow] {pt(a0)} -- {pt(a0 + 4.2 * jg.D0)};

  % inflow
  \foreach \yy in {{8.5, 4.8}} {{\draw[flow, freeblue] (-81.2,\yy) -- ++(6,0);}}
  \node[lbl, text=freeblue, anchor=south west] at (-81.8,9.9) {{$M_\infty$ = 8.03}};
  \foreach \yy in {{-35, -38.6}} {{\draw[flow, postgreen] (-81.2,\yy) -- ++({THETA:.2f}:6.1);}}
  \node[lbl, text=postgreen, anchor=south west] at (-81.8,-33.4) {{$M_2$ = 5.25}};

  % scale bar
  \draw[ink, line width=0.8pt, line cap=butt] (-51,10) -- (-41,10);
  \draw[ink, line width=0.6pt] (-51,9.5) -- (-51,10.5) (-41,9.5) -- (-41,10.5);
  \node[lbl, small, anchor=north] at (-46,9.3) {{10\,mm}};

  % labels
"""
    s += leader("bow shock", (-68.5, 2.5), g.upper_bow(2.5), "west")
    s += leader(r"upper\\[-0.4ex]triple point", (-66.5, -3.5), T1, "west", frm="west")
    s += leader("incident shock", (-75.5, -14.3), (-78, g.incident(-78)), "north")
    s += leader("transmitted shock", (-64, -21.5), TRANS[20], "north")
    s += leader("lower triple point", (-59, -28.5), T2, "north")
    sl_a = 0.5 * (T1 + jg.UPPER[0])
    s += leader("shear layers", (-51.5, -3), sl_a, "south", name="lsl")
    s += f"  \\draw[leader] (lsl.south) -- {pt(jg.at(0.62, +1))};\n"
    s += f"  \\draw[leader] (lsl.south) -- {pt(jg.at(0.1, -1))};\n"
    s += leader(r"jet bow\\[-0.4ex]shock", (-40.5, -25.5), 0.4 * jg.UPPER_END + 0.6 * jg.LOWER_END, "north")
    if detailed:
        s += r"""
  % Mach-number regions
  \node[region] at (-78,1.5) {$M > 1$};
  \node[region] at (-77.5,-23) {$M > 1$};
  \node[region] at (-61,8) {$M < 1$};
  \node[region] at (-33.5,-31.5) {$M < 1$};

  % key for the wave pattern inside the jet
  \draw[cellshock, line width=0.5pt] (-72,-34.6) -- ++(4,0);
  \node[lbl, small, anchor=west] at (-67.6,-34.6) {compression};
  \draw[fan, line width=0.5pt] (-72,-38) -- ++(4,0);
  \node[lbl, small, anchor=west] at (-67.6,-38) {expansion fan};
"""
        s += leader("bow shock", (-40.5, -38), g.lower_bow(-24.5), "west")
        s += f"  \\node[dot] at {pt(jet_mid(0.52))} {{}};\n"
        s += leader(r"jet, $M > 1$", (-44.5, 2.5), jet_mid(0.52), "south")
    else:
        s += leader(r"supersonic\\[-0.4ex]jet", (-44.5, 2.5), jet_mid(0.52), "south")
    return s


def labels_clean(sb):
    return common_labels(False)


def labels_detailed(sb):
    return common_labels(True)


ZOOM_CLEAN = ZOOM_DET = ZOOM


# ----------------------------------------------------------------- assemble
def build(path, title, zoom, labels, detailed, extra=""):
    pb, sb = panel_b(zoom, labels, detailed)
    tex = (header(title, extra) + PREAMBLE + "\n\\begin{document}\n\\sffamily\\footnotesize\n"
           "\\begin{tikzpicture}[x=1mm, y=1mm, line cap=round, line join=round]\n"
           + panel_a(zoom) + pb +
           f"""
% panel letters
\\node[anchor=north west, inner sep=0pt, font=\\sffamily\\bfseries\\small] at (0,{A_HEIGHT + 5:.1f}) {{(a)}};
\\node[anchor=north west, inner sep=0pt, font=\\sffamily\\bfseries\\small] at ({B_LEFT - 2:.1f},{A_HEIGHT + 5:.1f}) {{(b)}};

\\end{{tikzpicture}}
\\end{{document}}
""")
    tex = tex.replace("\n\n\n", "\n\n")
    path.write_text(tex)
    print(f"wrote {path.name}: {len(tex.splitlines())} lines, panel (b) scale {sb:.3f}")


if __name__ == "__main__":
    build(HERE / "fig_setup_typeiv_tikz.tex",
          "Computational set-up and Edney type IV shock structure",
          ZOOM_CLEAN, labels_clean, detailed=False)
    build(HERE / "fig_setup_typeiv_detailed_tikz.tex",
          "Computational set-up and Edney type IV shock structure (detailed)",
          ZOOM_DET, labels_detailed, detailed=True,
          extra="%  Jet shock cells and expansion fans are schematic (drawn steeper than\n"
                "%  the jet Mach angle for legibility).\n")
