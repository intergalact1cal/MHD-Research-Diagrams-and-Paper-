"""Check that every grid line drawn in fig_setup*.tex sits on a node level of the .su2 mesh.

    python3 verify_setup_grid.py . sa_mi_baseline.su2

Parses the circle radii and ray angles written in panels (a) and (c) and
compares them with radii and angles computed directly from the .su2 nodes.
"""
import re, sys, importlib.util, numpy as np
figdir, su2 = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("cm", f"{figdir}/check_mesh.py"); cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)
q, pts, marks = cm.read_su2(su2)
G = (pts * 1000).reshape(-1, 295, 2)
r_su2 = np.hypot(G[0, :, 0], G[0, :, 1])                        # radii of node row j = 0
phi_su2 = np.degrees(np.arctan2(G[:, 0, 1], -G[:, 0, 0]))         # angles of wall nodes
for name in ("fig_setup", "fig_setup_clean"):
    tex = open(f"{figdir}/{name}.tex").read()
    for panel, style, cs, rs in (("a", "gridA", 4, 8), ("c", "gridC", 2, 2)):
        arcs = [float(x) for x in re.findall(rf"\\draw\[{style}\] \([-\d.]+:[-\d.]+\) arc\[[^]]*radius=([-\d.]+)\]", tex)]
        rays = [180 - float(a) for a in re.findall(rf"\\draw\[{style}\] \(([-\d.]+):[\d.]+\) --", tex)]
        ii = np.unique(np.r_[np.arange(0, len(r_su2), cs), len(r_su2) - 1]); exp_r = r_su2[ii][1:]; exp_r = exp_r[exp_r <= max(arcs) + 1e-6]
        er = np.abs(np.array(arcs) - exp_r).max() if len(arcs) == len(exp_r) else None
        jj = np.arange(0, len(phi_su2), rs); jj = np.unique(np.r_[jj, len(phi_su2) - 1])
        exp_p = phi_su2[jj]; exp_p = exp_p[(exp_p >= min(rays) - 1e-6) & (exp_p <= max(rays) + 1e-6)]
        ep = np.abs(np.sort(rays) - np.sort(exp_p)).max() if len(rays) == len(exp_p) else None
        print(f"{name} ({panel}): {len(arcs)} circles, max |r_tex - r_su2| = {er:.1e} mm;  "
              f"{len(rays)} rays, max |phi_tex - phi_su2| = {ep:.1e} deg")
