"""Read sa_mi_baseline.su2, verify its structure and write the grid-line data.

    python3 check_mesh.py path/to/sa_mi_baseline.su2

Checks (each one stops with an error if it fails):
  * 2-D mesh of quadrilaterals only; node and element counts;
  * the nodes form a structured NI x NJ array (i radial, j circumferential)
    and every element is exactly the cell (i, j)-(i+1, j+1) of that array;
  * the grid is polar about the cylinder centre: every node lies on its
    circle r_i and its ray phi_j (max deviations reported);
  * every boundary marker lies where the figure draws it
    (wall i = 0, inlet i = NI-1, outlets j = 0 and j = NJ-1) and spans the
    stated angles.
Then writes data/mesh_r_i.txt and data/mesh_phi_j.txt (all node levels, mm and
deg, 9 decimals), which make_mesh_figure.py draws.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read_su2(path):
    lines = Path(path).read_text().split("\n")
    assert lines[0].split("=")[1].strip() == "2", "not a 2-D mesh"
    ne = int(lines[1].split("=")[1])
    el = np.array([l.split() for l in lines[2:2 + ne]], dtype=np.int64)
    assert (el[:, 0] == 9).all(), "non-quadrilateral elements present"
    k = 2 + ne
    npn = int(lines[k].split("=")[1].split()[0])
    pts = np.array([l.split()[:2] for l in lines[k + 1:k + 1 + npn]], dtype=float)
    k += 1 + npn
    nm = int(lines[k].split("=")[1])
    k += 1
    marks = {}
    for _ in range(nm):
        tag = lines[k].split("=")[1].strip()
        n = int(lines[k + 1].split("=")[1])
        k += 2
        marks[tag] = np.array([l.split()[1:3] for l in lines[k:k + n]], dtype=np.int64)
        k += n
    return el[:, 1:5], pts, marks


def main(path):
    quads, pts, marks = read_su2(path)
    P = pts * 1000.0                                           # m -> mm
    print(f"elements {len(quads)} (all quadrilateral), nodes {len(P)}, markers {list(marks)}")

    # ---- structured layout: first element is (0, 1, NI+1, NI)
    NI = int(quads[0, 3])
    NJ = len(P) // NI
    assert NI * NJ == len(P), "node count is not NI x NJ"
    assert len(quads) == (NI - 1) * (NJ - 1), "element count is not (NI-1)(NJ-1)"
    j, i = np.divmod(np.arange(len(quads)), NI - 1)
    expect = np.c_[j * NI + i, j * NI + i + 1, (j + 1) * NI + i + 1, (j + 1) * NI + i]
    assert (np.sort(quads, 1) == np.sort(expect, 1)).all(), "elements are not the structured cells"
    print(f"structured: NI = {NI} (radial), NJ = {NJ} (circumferential); "
          f"every element is cell (i, j)-(i+1, j+1)")

    # ---- polar about the cylinder centre
    wall_ids = np.unique(np.r_[marks["cylinder"].ravel(), marks["cyl_lee"].ravel()])
    q = P[wall_ids]
    A = np.c_[2 * q, np.ones(len(q))]
    c = np.linalg.lstsq(A, (q**2).sum(1), rcond=None)[0]
    centre, R = c[:2], np.sqrt(c[2] + c[0]**2 + c[1]**2)
    G = (P - centre).reshape(NJ, NI, 2)
    r = np.hypot(G[..., 0], G[..., 1])
    phi = np.degrees(np.arctan2(G[..., 1], -G[..., 0]))
    r_i, phi_j = r.mean(0), phi.mean(1)
    dr = np.abs(r - r_i[None, :]).max()
    dp = np.abs(phi - phi_j[:, None]).max()
    print(f"cylinder centre ({centre[0]:.2e}, {centre[1]:.2e}) mm, R = {R:.9f} mm")
    print(f"polar grid: max |r - r_i| = {dr:.1e} mm, max |phi - phi_j| = {dp:.1e} deg")
    assert np.abs(centre).max() < 1e-9 and dr < 1e-9 and dp < 1e-9, "grid is not exactly polar"
    assert np.all(np.diff(r_i) > 0) and np.all(np.diff(phi_j) < 0), "levels not monotonic"

    # ---- symmetry about y = 0 (same circles; rays at +phi and -phi)
    sym = np.abs(phi_j + phi_j[::-1]).max()
    print(f"symmetric about y = 0: max |phi_j + phi_(NJ-1-j)| = {sym:.1e} deg")
    assert sym < 1e-9, "grid is not symmetric about y = 0"

    # ---- markers
    def nodes(tag):
        ids = np.unique(marks[tag])
        jj, ii = np.divmod(ids, NI)
        return ii, jj, phi_j[jj]

    checks = {
        "cylinder": (0, None, (-45.0, 45.0)),
        "cyl_lee": (0, None, (-65.0, 65.0)),
        "inflow_upper": (NI - 1, None, None),
        "inflow_lower": (NI - 1, None, None),
        "cut_upper": (None, 0, None),
        "cut_lower": (None, NJ - 1, None),
    }
    for tag, (ii_req, jj_req, _) in checks.items():
        ii, jj, ph = nodes(tag)
        if ii_req is not None:
            assert (ii == ii_req).all(), f"{tag} not on i = {ii_req}"
        if jj_req is not None:
            assert (jj == jj_req).all(), f"{tag} not on j = {jj_req}"
        rr = r_i[ii]
        print(f"  {tag:13s} {len(marks[tag]):3d} edges  i {ii.min()}..{ii.max()}  j {jj.min()}..{jj.max()}  "
              f"r {rr.min():.3f}..{rr.max():.3f} mm  phi {ph.min():+.2f}..{ph.max():+.2f} deg")
    _, _, pc = nodes("cylinder")
    assert np.isclose(pc.min(), -45) and np.isclose(pc.max(), 45), "no-slip wall not |phi| <= 45"
    _, _, pu = nodes("inflow_upper")
    _, _, pl = nodes("inflow_lower")
    print(f"  inlet split at phi = {pu.min():.4f} deg (= {pl.max():.4f} deg)")

    # ---- spacing summary
    d = np.diff(r_i)
    g = d[1:] / d[:-1]
    print(f"radial: {NI - 1} cells, first {d[0] * 1e3:.4f} um, wall ratio {g[20]:.4f}, "
          f"uniform {d[100]:.5f} mm for r = {r_i[54]:.4f}..{r_i[233]:.4f}, outer ratio {g[260]:.4f}, "
          f"last {d[-1]:.4f} mm")
    print(f"circumferential: {NJ - 1} cells, dphi = {np.ptp(np.diff(phi_j)):.1e} spread, "
          f"{-np.diff(phi_j)[0]:.6f} deg")

    out = HERE / "data"
    np.savetxt(out / "mesh_r_i.txt", r_i, fmt="%.9f",
               header=f"radial node levels r_i (mm) of sa_mi_baseline.su2, i = 0..{NI - 1} (wall to inlet);\n"
                      "written by check_mesh.py")
    np.savetxt(out / "mesh_phi_j.txt", phi_j, fmt="%.9f",
               header=f"circumferential node angles phi_j (deg from the stagnation point) of "
                      f"sa_mi_baseline.su2, j = 0..{NJ - 1};\nwritten by check_mesh.py")
    print("wrote data/mesh_r_i.txt, data/mesh_phi_j.txt")
    return G


if __name__ == "__main__":
    main(sys.argv[1])
