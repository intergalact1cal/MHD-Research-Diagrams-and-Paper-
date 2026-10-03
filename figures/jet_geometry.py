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

# schematic downstream part
TURN = np.radians(24.0)          # total turning of the jet toward the wall normal
PHI_IMP = np.radians(-27.5)      # impingement point on the wall (from stagnation)
STANDOFF = 0.85 * WIDTH          # jet bow shock stand-off from the wall
DELTA1 = DELTA0 + TURN
D1 = np.array([np.cos(DELTA1), np.sin(DELTA1)])
N1 = np.array([-np.sin(DELTA1), np.cos(DELTA1)])

WALL_IMP = np.array([-R * np.cos(PHI_IMP), R * np.sin(PHI_IMP)])


def _wall_hit(p, d):
    """First intersection of the ray p + s d with the cylinder."""
    b = p @ d
    c = p @ p - R**2
    return p + (-b - np.sqrt(b * b - c)) * d


C0 = T2 + H * N0                                   # jet centreline start
JC = _wall_hit(WALL_IMP - 30 * D1, D1) - STANDOFF * D1   # jet bow-shock centre
# quadratic centreline: tangent D0 at C0, tangent D1 at JC (uniform turning)
_A = np.c_[D0, -D1]
_u = np.linalg.solve(_A, JC - C0)
CTRL = (C0, C0 + _u[0] * D0, JC)


def centreline(n=200):
    t = np.linspace(0, 1, n)[:, None]
    P0, P1, P2 = CTRL
    c = (1 - t)**2 * P0 + 2 * (1 - t) * t * P1 + t**2 * P2
    dc = 2 * (1 - t) * (P1 - P0) + 2 * t * (P2 - P1)
    tang = dc / np.linalg.norm(dc, axis=1, keepdims=True)
    nrm = np.c_[-tang[:, 1], tang[:, 0]]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))]
    return c, nrm, s


C, NRM, S = centreline()
JS_EXT, JS_BULGE = 1.18, 0.45


def _shock_point(side):
    """Where the jet boundary on `side` meets the jet bow shock."""
    return JC + side * H * N1 + JS_BULGE * H * (1 / JS_EXT)**2 * D1


# boundaries run from the triple points to the jet bow shock
UPPER = np.r_[C + H * NRM, _shock_point(+1)[None, :]]   # starts on the line through T1
LOWER = np.r_[C - H * NRM, _shock_point(-1)[None, :]]   # starts at T2
UPPER_FULL = np.r_[T1[None, :], UPPER]


def at(s_frac, side):
    """Point on a jet boundary (side=+1 upper, -1 lower) at fraction of length."""
    s = s_frac * S[-1]
    i = np.searchsorted(S, s).clip(1, len(S) - 1)
    a = (s - S[i - 1]) / (S[i] - S[i - 1])
    c = C[i - 1] + a * (C[i] - C[i - 1])
    n = NRM[i - 1] + a * (NRM[i] - NRM[i - 1])
    return c + side * H * n / np.linalg.norm(n)


def jet_shock(n=30, ext=JS_EXT, bulge=JS_BULGE):
    """Jet bow shock: short, slightly curved, convex toward the oncoming jet."""
    s = np.linspace(-ext * H, ext * H, n)[:, None]
    return JC + s * N1 + bulge * H * (s / (ext * H))**2 * D1


def after_shock(side, n=40, spread_deg=15.0, gap=0.9):
    """Shear layer past the jet bow shock, turning onto the wall as a wall jet.

    Leaves the jet bow shock along the jet direction and ends `gap` mm off the
    wall, `spread_deg` from the impingement point, running parallel to it.
    """
    p0 = _shock_point(side)
    phi_end = PHI_IMP + side * np.radians(spread_deg)
    p3 = (R + gap) * np.array([-np.cos(phi_end), np.sin(phi_end)])
    tw = side * np.array([np.sin(phi_end), np.cos(phi_end)])   # along the wall, away
    dist = np.linalg.norm(p3 - p0)
    d_start = D1 + 0.8 * tw
    d_start /= np.linalg.norm(d_start)
    p1 = p0 + 0.35 * dist * d_start
    p2 = p3 - 0.45 * dist * tw
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t)**3 * p0 + 3 * (1 - t)**2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3


def jet_shock_between(n=20):
    """Part of the jet bow shock between the two jet boundaries (upper -> lower)."""
    s = np.linspace(H, -H, n)[:, None]
    return JC + s * N1 + JS_BULGE * H * (s / (JS_EXT * H))**2 * D1


def wave_pattern(shock_deg=52.0, fan_deg=(43.0, 52.0, 63.0), end=0.97):
    """Shock (solid) / expansion-fan (dashed) cells reflecting off the jet boundaries.

    A compression starting at the lower triple point crosses the jet, reflects
    from the constant-pressure upper boundary as an expansion fan, which
    reflects back as compressions, and so on until the jet bow shock.
    Waves are inclined at shock_deg / fan_deg to the local jet axis (the Mach
    angle of the jet is about 25 deg; the cells are drawn shorter (steeper)
    for legibility).  Returns (shocks, fans) as lists of 2-point segments.
    """
    shocks, fans = [], []
    L = S[-1]
    s, side = 0.0, -1
    while True:
        s1 = s + WIDTH / np.tan(np.radians(shock_deg))
        if s1 > end * L:
            break
        shocks.append((at(s / L, side), at(s1 / L, -side)))
        land = [s1 + WIDTH / np.tan(np.radians(a)) for a in fan_deg]
        for sl in land:
            if sl <= end * L:
                fans.append((at(s1 / L, -side), at(sl / L, side)))
        s = land[0]
        if s > end * L:
            break
    return shocks, fans


if __name__ == "__main__":
    print(f"transmitted shock {np.degrees(SIGMA):.2f} deg, wave angle {np.degrees(BETA_T):.2f} deg")
    print(f"jet direction {np.degrees(DELTA0):.2f} deg, M_jet = {M_JET:.2f}, width = {WIDTH:.2f} mm")
    print("jet shock centre", JC, "r =", np.linalg.norm(JC), "length", S[-1])
