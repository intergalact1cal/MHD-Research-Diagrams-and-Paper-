"""Computational set-up and Edney type IV shock structure.

Panel (a): full computational domain with boundary conditions.
Panel (b): close-up of the type IV interaction (region boxed in (a)).

Geometry is in millimetres with the origin at the cylinder centre, x pointing
downstream and y up.  phi is the polar angle measured from the stagnation
point, positive upward.  The shock positions are the no-field solution
(data/bow_shock_no_field.csv); the far-field shock is a Billig-type
hyperbola fitted to that trace with its asymptote fixed at the local Mach
angle.  The incident-shock deflection and M2 follow from the oblique-shock
relations (gamma = 1.4).

Run:  python3 fig_setup_typeiv.py   ->  fig_setup_typeiv.{pdf,svg,png}
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Rectangle
from matplotlib.path import Path as MPath
from scipy.optimize import brentq, least_squares

HERE = Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# Case parameters
# --------------------------------------------------------------------------
GAMMA = 1.4
M_INF = 8.03
BETA = np.radians(18.1)       # imposed oblique-shock angle
R = 38.0                      # cylinder radius
R_OUT = 150.0                 # inlet-arc radius
PHI_OUT = np.radians(65.0)    # outlet half-angle (from stagnation line)
PHI_NOSLIP = np.radians(45.0) # no-slip / slip wall junction
PHI_SPLIT = np.radians(-12.25)
T_WALL = 294

M1n = M_INF * np.sin(BETA)
THETA = np.arctan(2 / np.tan(BETA) * (M1n**2 - 1)
                  / (M_INF**2 * (GAMMA + np.cos(2 * BETA)) + 2))
M2n = np.sqrt((1 + 0.5 * (GAMMA - 1) * M1n**2) / (GAMMA * M1n**2 - 0.5 * (GAMMA - 1)))
M_2 = M2n / np.sin(BETA - THETA)

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
INK = "#111111"
GREY = "#8a8a8a"
LEADER = "#9a9a9a"
BODY = "#e4e4e4"
C_FREE = "#1a5fa8"   # freestream inlet  (same palette as fig_mhd_structure)
C_POST = "#0e8a62"   # post-shock inlet
C_SLIP = "#b06a10"   # slip wall
FS = 8.0
DEG = r"$\!$°"   # FreeSans° carries a wide left bearing

mpl.rcParams.update({
    "font.family": "FreeSans",
    "font.size": FS,
    "mathtext.fontset": "custom",
    "mathtext.rm": "FreeSans",
    "mathtext.it": "FreeSans:italic",
    "mathtext.bf": "FreeSans:bold",
    "axes.unicode_minus": True,
    "pdf.fonttype": 42,
    "svg.fonttype": "path",
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "butt",
})


def pol(r, phi):
    """Polar (r, phi from stagnation point) -> Cartesian."""
    return np.array([-r * np.cos(phi), r * np.sin(phi)])


def arc(r, p0, p1, n=200):
    phi = np.linspace(p0, p1, n)
    return -r * np.cos(phi), r * np.sin(phi)


# --------------------------------------------------------------------------
# Shock system
# --------------------------------------------------------------------------
trace = np.loadtxt(HERE / "data" / "bow_shock_no_field.csv", delimiter=",")


def hyperbola(p, s, k):
    x0, a, s0 = p
    return x0 + np.sqrt(a**2 + (k * (s - s0))**2) - a


rot = np.array([[np.cos(THETA), np.sin(THETA)], [-np.sin(THETA), np.cos(THETA)]])
k_up = 1 / np.tan(np.arcsin(1 / M_INF))
k_lo = 1 / np.tan(np.arcsin(1 / M_2))
up_pts = trace[trace[:, 1] > -3]
lo_pts = trace[trace[:, 1] < -19.5] @ rot.T
p_up = least_squares(lambda p: hyperbola(p, up_pts[:, 1], k_up) - up_pts[:, 0],
                     [-72, 500, 0]).x
p_lo = least_squares(lambda p: hyperbola(p, lo_pts[:, 1], k_lo) - lo_pts[:, 0],
                     [-55, 500, 0]).x


def upper_bow(y):
    return np.array([hyperbola(p_up, y, k_up), y])


def lower_bow(v):
    """Lower bow shock parametrised by the cross-stream coordinate of the M2 flow."""
    return rot.T @ np.array([hyperbola(p_lo, v, k_lo), v])


SPLIT = pol(R_OUT, PHI_SPLIT)


def incident(x):
    return SPLIT[1] + np.tan(BETA) * (x - SPLIT[0])


# upper triple point: incident shock meets the upper bow shock
y_t1 = brentq(lambda y: upper_bow(y)[0] - (SPLIT[0] + (y - SPLIT[1]) / np.tan(BETA)), -20, 10)
T1 = upper_bow(y_t1)

# transmitted shock: straight part of the trace, extended to the lower bow shock
mid = trace[(trace[:, 1] < -10) & (trace[:, 1] > -16.2)]
slope, icpt = np.polyfit(mid[:, 0], mid[:, 1], 1)
v_t2 = brentq(lambda v: lower_bow(v)[1] - (slope * lower_bow(v)[0] + icpt), -40, 0)
T2 = lower_bow(v_t2)


def transmitted(n=60):
    """Quadratic through T1 and T2 that follows the straight part of the trace."""
    xs = np.r_[T1[0], mid[:, 0], T2[0]]
    ys = np.r_[T1[1], mid[:, 1], T2[1]]
    w = np.r_[50, np.ones(len(mid)), 50]
    c = np.polyfit(xs, ys, 2, w=w)
    x = np.linspace(T1[0], T2[0], n)
    return x, np.polyval(c, x)


# supersonic jet: shear layers from the triple points to the jet bow shock
PHI_JET = np.radians(-19.5)          # impingement location
HALF_JET = np.radians(6.5)
R_JS = 44.0                          # jet bow-shock stand-off radius
JS_SPAN = 1.45 * HALF_JET


def r_jet_shock(phi):
    return R_JS + 0.35 - 2.2 * ((phi - PHI_JET) / JS_SPAN)**2


J1 = pol(r_jet_shock(PHI_JET + HALF_JET), PHI_JET + HALF_JET)
J2 = pol(r_jet_shock(PHI_JET - HALF_JET), PHI_JET - HALF_JET)


def bezier(p0, p1, p2, n=80):
    t = np.linspace(0, 1, n)[:, None]
    return ((1 - t)**2 * p0 + 2 * (1 - t) * t * p1 + t**2 * p2).T


def shear_layer(p0, p1, bend):
    c = 0.5 * (p0 + p1) + np.array([0.0, bend])
    return bezier(p0, c, p1)


SL1 = shear_layer(T1, J1, 0.9)
SL2 = shear_layer(T2, J2, 0.4)


def jet_shock():
    phi = np.linspace(PHI_JET - JS_SPAN, PHI_JET + JS_SPAN, 80)
    r = r_jet_shock(phi)
    return -r * np.cos(phi), r * np.sin(phi)


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def arrow(ax, p0, p1, color=INK, lw=1.2, ms=7, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=color, shrinkA=0, shrinkB=0, **kw))


def dim_arrow(ax, p0, p1, color=GREY, lw=0.7):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<|-|>", mutation_scale=6,
                                 lw=lw, color=color, shrinkA=0, shrinkB=0))


def label(ax, text, xy_text, xy_target, ha="left", va="center", color=INK,
          gap=None, **kw):
    """Text with a thin grey leader line ending on the feature."""
    ax.annotate(text, xy=xy_target, xytext=xy_text, ha=ha, va=va, color=color,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=LEADER,
                                shrinkA=2.5, shrinkB=0.6),
                **kw)


def flow_arrows(ax, start, ys, angle, length, color):
    d = length * np.array([np.cos(angle), np.sin(angle)])
    for y in ys:
        p = np.array([start, y])
        arrow(ax, p, p + d, color=color, lw=1.3, ms=8)


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(7.0, 4.15))
ax_a = fig.add_axes([0.005, 0.02, 0.455, 0.93])
ax_b = fig.add_axes([0.505, 0.02, 0.49, 0.93])
for ax in (ax_a, ax_b):
    ax.set_aspect("equal")
    ax.axis("off")

ZOOM = (-100.0, -26.0, -50.0, 38.0)   # x0, x1, y0, y1 of panel (b)

# ---------------------------- panel (a) -----------------------------------
ax = ax_a
ax.set_xlim(-218, 52)
ax.set_ylim(-152, 152)

# domain outline (used as clip path for the shocks)
xi, yi = arc(R_OUT, -PHI_OUT, PHI_OUT)
xw, yw = arc(R, PHI_OUT, -PHI_OUT)
domain = MPath(np.c_[np.r_[xi, xw], np.r_[yi, yw]])
clip = PathPatch(domain, transform=ax.transData, fc="none", ec="none")
ax.add_patch(clip)

# cylinder
ax.add_patch(Circle((0, 0), R, fc=BODY, ec=GREY, lw=0.7, ls=(0, (1, 1.6)), zorder=1))
ax.plot(*arc(R, -PHI_NOSLIP, PHI_NOSLIP), color=INK, lw=2.4, zorder=4)
for s in (1, -1):
    ax.plot(*arc(R, s * PHI_NOSLIP, s * PHI_OUT), color=C_SLIP, lw=2.4,
            ls=(0, (2.2, 1.2)), zorder=4)
ax.plot(0, 0, marker="+", ms=6, mew=0.9, color=INK, zorder=5)

# inlets and outlets
ax.plot(*arc(R_OUT, PHI_SPLIT, PHI_OUT), color=C_FREE, lw=2.2, zorder=4)
ax.plot(*arc(R_OUT, -PHI_OUT, PHI_SPLIT), color=C_POST, lw=2.2, zorder=4)
for s in (1, -1):
    p0, p1 = pol(R, s * PHI_OUT), pol(R_OUT, s * PHI_OUT)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=GREY, lw=1.6,
            ls=(0, (4, 2.4)), zorder=3)

# shock system (no-field solution), clipped to the domain
sh = dict(color="#444444", lw=1.0, zorder=2, clip_path=clip)
y = np.linspace(y_t1, 160, 300)
ax.plot(*upper_bow(y), **sh)
v = np.linspace(v_t2, -170, 300)
ax.plot(*np.array([lower_bow(vv) for vv in v]).T, **sh)
ax.plot(*transmitted(), **sh)
ax.plot([SPLIT[0], T1[0]], [SPLIT[1], T1[1]], **sh)

# imposed oblique shock upstream of the inlet (outside the domain)
x_up = -214
ax.plot([x_up, SPLIT[0]], [incident(x_up), SPLIT[1]], color="#444444", lw=1.0,
        ls=(0, (3, 2)), zorder=3)
ax.plot(*SPLIT, marker="o", ms=4.6, mfc="white", mec=INK, mew=1.0, zorder=6)

# beta
ax.plot([SPLIT[0], SPLIT[0] + 46], [SPLIT[1]] * 2, color=GREY, lw=0.6,
        ls=(0, (1.5, 1.5)), zorder=3)
a = np.linspace(0, BETA, 30)
ax.plot(SPLIT[0] + 34 * np.cos(a), SPLIT[1] + 34 * np.sin(a), color=INK, lw=0.6)
ax.text(SPLIT[0] + 7, SPLIT[1] - 3, r"$\beta$ = 18.1" + DEG, ha="left", va="top")
ax.text(SPLIT[0] - 7, SPLIT[1] + 9, "inlet split\n" r"$\phi$ = −12.25" + DEG,
        ha="right", va="bottom", color=GREY, fontsize=7, linespacing=1.15)

# inflow arrows
flow_arrows(ax, -214, [116, 94, 72], 0.0, 30, C_FREE)
flow_arrows(ax, -214, [-72, -94, -116], THETA, 30, C_POST)
ax.text(-214, 128, "freestream inlet\n" r"$M_\infty$ = 8.03", color=C_FREE,
        ha="left", va="bottom", linespacing=1.25)
ax.text(-214, -122, "post-shock inlet\n" r"$M_2$ = 5.25", color=C_POST,
        ha="left", va="top", linespacing=1.25)

# walls
ax.text(6, 18, "no-slip wall\n" r"$T_w$ = 294 K", ha="center", va="center",
        linespacing=1.25, fontsize=7.5)
ax.annotate("", xy=pol(R, np.radians(18)), xytext=(-16, 21),
            arrowprops=dict(arrowstyle="-", lw=0.6, color=LEADER, shrinkA=0, shrinkB=0.6))
label(ax, "slip wall", (26, 64), pol(R, np.radians(56)), ha="left")
for s, yy in ((1, 128), (-1, -128)):
    label(ax, "outlet", (-4, yy), pol(0.8 * R_OUT, s * PHI_OUT), ha="left",
          color=GREY)
label(ax, "bow shock", (-96, 100), upper_bow(70.0), ha="left")

# dimensions
dim_arrow(ax, (0, 0), pol(R, np.radians(-150)))
ax.text(40, -27, r"$R$ = 38 mm", ha="left", va="center", fontsize=7.5)
pd = np.radians(-37)
q1 = pol(R_OUT, pd)
dim_arrow(ax, (0, 0), q1)
d = q1 / np.hypot(*q1)
nrm = np.array([-d[1], d[0]])
if nrm[1] < 0:
    nrm = -nrm
ang = np.degrees(np.arctan2(d[1], d[0]))
ang = ang - 180 if ang > 90 else ang + 180 if ang < -90 else ang
ax.text(*(0.72 * q1 + 5.0 * nrm), "150 mm", rotation=ang, ha="center",
        va="center", fontsize=7, color=GREY, rotation_mode="anchor")

# zoom box
zx0, zx1, zy0, zy1 = ZOOM
ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fc="none", ec=INK,
                       lw=0.6, ls=(0, (2, 1.5)), zorder=7))
ax.text(zx0 + 2, zy1 - 2.5, "(b)", ha="left", va="top", fontsize=7.5)

# ---------------------------- panel (b) -----------------------------------
ax = ax_b
ax.set_xlim(zx0, zx1)
ax.set_ylim(zy0, zy1)
frame = Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, transform=ax.transData,
                  fc="none", ec="none")
ax.add_patch(frame)
cl = dict(clip_path=frame)

ax.add_patch(Circle((0, 0), R, fc=BODY, ec="none", zorder=1, **cl))
ax.plot(*arc(R, -PHI_NOSLIP, PHI_NOSLIP), color=INK, lw=2.6, zorder=4, **cl)

lw_s = 1.6
y = np.linspace(y_t1, zy1 + 5, 200)
ax.plot(*upper_bow(y), color=INK, lw=lw_s, zorder=3, **cl)
v = np.linspace(v_t2, -80, 300)
ax.plot(*np.array([lower_bow(vv) for vv in v]).T, color=INK, lw=lw_s, zorder=3, **cl)
ax.plot(*transmitted(), color=INK, lw=lw_s, zorder=3)
xx = np.array([zx0 - 5, T1[0]])
ax.plot(xx, incident(xx), color=INK, lw=lw_s, zorder=3, **cl)
dash = (0, (3.2, 1.8))
ax.plot(*SL1, color=INK, lw=1.0, ls=dash, zorder=3)
ax.plot(*SL2, color=INK, lw=1.0, ls=dash, zorder=3)
ax.plot(*jet_shock(), color=INK, lw=lw_s, zorder=3)
for T in (T1, T2):
    ax.plot(*T, marker="o", ms=4.2, color=INK, zorder=6)

# jet direction (left part of the jet, clear of the leader lines)
jc0 = 0.5 * (SL1[:, 54] + SL2[:, 2]) + [-0.4, 0]
arrow(ax, jc0, jc0 + 4.2 * np.array([np.cos(-0.07), np.sin(-0.07)]), lw=1.1, ms=7)

# inflow
flow_arrows(ax, -98, [22, 15], 0.0, 12, C_FREE)
ax.text(-98, 27, r"$M_\infty$ = 8.03", color=C_FREE, ha="left", va="bottom")
flow_arrows(ax, -98, [-38, -45], THETA, 12, C_POST)
ax.text(-98, -32.5, r"$M_2$ = 5.25", color=C_POST, ha="left", va="bottom")

# labels
js = jet_shock()
label(ax, "bow shock", (-82, 9), upper_bow(9.0), ha="right")
label(ax, "upper triple point", (-89, -3.5), T1, ha="center", va="bottom")
label(ax, "incident shock", (-86, -22), (-90, incident(-90)), ha="center", va="top")
xt, yt = transmitted()
label(ax, "transmitted shock", (-68, -27), (xt[30], yt[30]), ha="center", va="top")
label(ax, "lower triple point", (-56, -40), T2, ha="center", va="top")
label(ax, "shear layers", (-55, 4), SL1[:, 30], ha="center", va="bottom")
ax.annotate("", xy=SL2[:, 60], xytext=(-55, 4),
            arrowprops=dict(arrowstyle="-", lw=0.6, color=LEADER, shrinkA=2.5,
                            shrinkB=0.6))
label(ax, "jet bow\nshock", (-40.5, 13), (js[0][72], js[1][72]), ha="center", va="bottom",
      linespacing=1.15)
label(ax, "supersonic\njet", (-33.5, -33), jc0 + [4.2, -2.6],
      ha="center", va="top", linespacing=1.15)

# scale bar
sb_x, sb_y = -64, 34
ax.plot([sb_x, sb_x + 10], [sb_y, sb_y], color=INK, lw=1.2, solid_capstyle="butt")
for xb in (sb_x, sb_x + 10):
    ax.plot([xb, xb], [sb_y - 0.9, sb_y + 0.9], color=INK, lw=1.0)
ax.text(sb_x + 5, sb_y - 2.2, "10 mm", ha="center", va="top", fontsize=7.5)

# panel letters
fig.text(0.008, 0.975, "(a)", fontweight="bold", fontsize=10, va="top")
fig.text(0.505, 0.975, "(b)", fontweight="bold", fontsize=10, va="top")

if __name__ == "__main__":
    print(f"theta = {np.degrees(THETA):.2f} deg, M2 = {M_2:.3f}")
    print(f"T1 = ({T1[0]:.1f}, {T1[1]:.1f}) mm, phi = {np.degrees(np.arctan2(T1[1], -T1[0])):.1f} deg")
    print(f"T2 = ({T2[0]:.1f}, {T2[1]:.1f}) mm, phi = {np.degrees(np.arctan2(T2[1], -T2[0])):.1f} deg")
    for ext in ("pdf", "svg"):
        fig.savefig(HERE / f"fig_setup_typeiv.{ext}")
    fig.savefig(HERE / "fig_setup_typeiv.png", dpi=600)
