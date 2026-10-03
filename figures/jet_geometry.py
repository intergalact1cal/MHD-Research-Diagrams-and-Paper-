"""Supersonic-jet geometry for the Edney type IV figure.

The jet is the post-shock (M2) flow that has also crossed the transmitted
shock T1-T2.  Its initial direction and Mach number follow from the
oblique-shock relations with the transmitted-shock inclination taken from the
shock trace; its width is the separation of T1 and T2 normal to that
direction.  Downstream the jet is drawn curving toward the wall-normal
direction and terminating in a jet bow shock a short stand-off from the wall,
as in the classical type IV structure (Edney 1968; Wieting & Holden 1989).
That downstream part (curvature, impingement point, stand-off, internal wave
pattern) is schematic.

Coordinates: mm, origin at the cylinder centre, x downstream, y up.
"""
import numpy as np
from scipy.optimize import brentq

import fig_setup_typeiv as g

GAMMA = g.GAMMA


def theta_beta_m(M, beta):
    return np.arctan(2 / np.tan(beta) * (M**2 * np.sin(beta)**2 - 1)
                     / (M**2 * (GAMMA + np.cos(2 * beta)) + 2))


def mach_after(M, beta):
    th = theta_beta_m(M, beta)
    Mn = M * np.sin(beta)
    Mn2 = np.sqrt((1 + 0.5 * (GAMMA - 1) * Mn**2) / (GAMMA * Mn**2 - 0.5 * (GAMMA - 1)))
    return Mn2 / np.sin(beta - th)


T1, T2, R = g.T1, g.T2, g.R

# transmitted shock (chord T1 -> T2) and the flow behind it
SIGMA = np.arctan2(T2[1] - T1[1], T2[0] - T1[0])
BETA_T = g.THETA - SIGMA                       # wave angle seen by the M2 flow
DELTA0 = g.THETA - theta_beta_m(g.M_2, BETA_T)  # jet direction at the triple points
M_JET = mach_after(g.M_2, BETA_T)

D0 = np.array([np.cos(DELTA0), np.sin(DELTA0)])
N0 = np.array([-np.sin(DELTA0), np.cos(DELTA0)])
WIDTH = float((T1 - T2) @ N0)
H = 0.5 * WIDTH

# schematic downstream part (shape after the classical type IV sketches):
# the jet follows a circular arc, turning by TURN, and ends in the jet bow
# shock STANDOFF from the wall.
TURN = np.radians(36.0)
STANDOFF = 0.9 * WIDTH
DELTA1 = DELTA0 + TURN
D1 = np.array([np.cos(DELTA1), np.sin(DELTA1)])
N1 = np.array([-np.sin(DELTA1), np.cos(DELTA1)])

C0 = T2 + H * N0                                   # jet centreline start


def _arc_point(rho, delta):
    centre = C0 + rho * N0
    return centre + rho * np.array([np.sin(delta), -np.cos(delta)])


def _end_gap(rho):
    return np.linalg.norm(_arc_point(rho, DELTA1)) - (R + STANDOFF)


RHO = brentq(_end_gap, 2.0, 80.0)                  # radius of curvature of the jet
JC = _arc_point(RHO, DELTA1)                       # jet bow-shock centre
PHI_IMP = np.arctan2(JC[1], -JC[0])                # impingement direction
WALL_IMP = R * np.array([-np.cos(PHI_IMP), np.sin(PHI_IMP)])


def centreline(n=240):
    d = np.linspace(DELTA0, DELTA1, n)
    c = np.array([_arc_point(RHO, x) for x in d])
    nrm = np.c_[-np.sin(d), np.cos(d)]
    s = RHO * (d - DELTA0)
    return c, nrm, s


C, NRM, S = centreline()

# boundaries run from the triple points to the jet bow shock
UPPER = C + H * NRM            # starts on the line through T1
LOWER = C - H * NRM            # starts at T2
UPPER_FULL = np.r_[T1[None, :], UPPER]
UPPER_END, LOWER_END = UPPER[-1], LOWER[-1]


def at(s_frac, side):
    """Point on a jet boundary (side=+1 upper, -1 lower) at fraction of length."""
    s = np.clip(s_frac, 0, 1) * S[-1]
    i = np.searchsorted(S, s).clip(1, len(S) - 1)
    a = (s - S[i - 1]) / (S[i] - S[i - 1])
    c = C[i - 1] + a * (C[i] - C[i - 1])
    n = NRM[i - 1] + a * (NRM[i] - NRM[i - 1])
    return c + side * H * n / np.linalg.norm(n)


def jet_shock():
    """Jet bow shock: a short, nearly normal shock closing the end of the jet."""
    return np.array([UPPER_END, LOWER_END])


jet_shock_between = jet_shock


def _polar(p):
    return np.hypot(*p), np.arctan2(p[1], -p[0])


def after_shock(side, n=30, sweep_deg=17.0, approach=0.45):
    """Shear layer beyond the jet bow shock.

    It turns sharply at the end of the jet bow shock and runs along the wall,
    away from the impingement point (upper boundary upward, lower boundary
    downward), closing in on the wall from the jet-shock stand-off to
    `approach` times that distance.
    """
    p0 = UPPER_END if side > 0 else LOWER_END
    r0, phi0 = _polar(p0)
    phi = phi0 + side * np.linspace(0, np.radians(sweep_deg), n)
    r = r0 + (R + approach * (r0 - R) - r0) * np.linspace(0, 1, n)
    return np.c_[-r * np.cos(phi), r * np.sin(phi)]


def wave_pattern(n_cells=2, first=0.2, last=0.86,
                 fan=(0.22, 0.42, 0.62, 0.82), comp=(0.12, 0.32, 0.52, 0.72)):
    """Shock-cell structure of the jet, as in the classical type IV sketch.

    The wave starting at the lower triple point T2 crosses the jet to a point
    F0 on the upper (T1-side) boundary.  From each focus F_k an expansion fan
    (dashed) spreads to the lower boundary; compressions (solid) leave the
    lower boundary and converge on the next focus F_k+1, and so on to the jet
    bow shock.  Positions are fractions of the jet length; the cells are
    schematic (drawn steeper than the jet Mach angle, about 25 deg).
    Returns (shocks, fans) as lists of 2-point segments.
    """
    shocks = [(T2.copy(), at(first, +1))]
    fans = []
    cell = (last - first) / n_cells
    for k in range(n_cells + 1):
        fk = first + k * cell
        if k < n_cells:
            for a in fan:
                fans.append((at(fk, +1), at(fk + a * cell, -1)))
            for a in comp:
                shocks.append((at(fk + a * cell, -1), at(fk + cell, +1)))
        else:
            # last focus: a short fan toward the jet bow shock
            for a in fan[:2]:
                if fk + a * cell <= 0.995:
                    fans.append((at(fk, +1), at(fk + a * cell, -1)))
    return shocks, fans


if __name__ == "__main__":
    print(f"transmitted shock {np.degrees(SIGMA):.2f} deg, wave angle {np.degrees(BETA_T):.2f} deg")
    print(f"jet direction {np.degrees(DELTA0):.2f} deg, M_jet = {M_JET:.2f}, width = {WIDTH:.2f} mm")
    print("jet shock centre", JC, "r =", np.linalg.norm(JC), "length", S[-1],
          "rho", RHO, "phi_imp", np.degrees(PHI_IMP))
