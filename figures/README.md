# Figures

| file | content |
|---|---|
| `fig_setup.tex` / `.pdf` / `.png` | **(a)** domain, boundary conditions and grid; **(b)** Edney type IV interaction with the jet shock cells; **(c)** grid in the shock layer at the stagnation point |
| `fig_setup_clean.*` | same, with a simpler (b): no shock cells, shaded jet |

The figure is drawn at its print size, **170 × 108 mm**, with 9 pt and 8 pt
labels. Include it at full text width:

```latex
\begin{figure*}
  \centering
  \includegraphics[width=\textwidth]{fig_setup.pdf}
  \caption{(a) Computational domain, boundary conditions and grid (every 8th
  radial line shown); dashed and red boxes mark (b) and (c).
  (b) Edney type IV interaction; inside the jet, solid lines are
  compressions and dashed lines expansion fans.
  (c) Grid in the shock layer at the stagnation point (every 2nd radial line),
  with the no-field bow shock.}
\end{figure*}
```

Grid lines are the real lines of `sa_mi_baseline.su2`. Where two lines would
be closer than about 0.85 mm (a) or 0.55 mm (c) on paper, lines are skipped,
so the wall clustering shows as grading rather than a solid band. The full
spacing is in the notes below.

## `fig_mesh` — computational grid (standalone)

`fig_mesh.tex` / `.pdf` / `.png` show the grid of `sa_mi_baseline.su2` at
168 × 147 mm, full text width: **(a)** the whole grid, **(b)** the near-wall
grid at the stagnation point. The red box in (a) marks (b).

Suggested caption: *Computational grid (294 × 520 cells). (a) Whole domain,
every 4th radial line shown. (b) Near-wall region at the stagnation point,
every radial line shown. In both panels, circumferential lines closer than
0.4 mm on the page are omitted.*

**How it is made and checked** (no hand drawing; all from the `.su2` file):

1. `check_mesh.py sa_mi_baseline.su2` reads the mesh and stops if any check
   fails:
   * all 152 880 elements are quadrilaterals;
   * every element is cell (i, j)–(i+1, j+1) of a structured 295 × 521 node
     array;
   * every node lies on its circle r_i and ray φ_j about the cylinder
     centre, to 1e−12 mm / 4e−14°;
   * the wall, inlet and outlet markers lie on i = 0, i = 294 and
     j = 0 / 520.

   It then writes the exact node levels (`data/mesh_r_i.txt`,
   `data/mesh_phi_j.txt`, 9 decimals).
2. `make_mesh_figure.py` draws those levels: circles as exact arcs, rays as
   straight lines through the nodes.
3. `overlay_check.py` redraws the same grid lines as polylines through the raw
   `.su2` node coordinates and overlays them on the compiled figure at 600 dpi.
   Result: 100 % of the `.su2` line pixels fall on a drawn line in both panels.

Mesh facts confirmed by the check: R = 38.1 mm; inlet r = 150.387 mm;
outlets at φ = ±65°; no-slip wall (`cylinder`) |φ| ≤ 45°, slip wall
(`cyl_lee`) beyond; inlet split at **φ = −12.50°**. Radially, the first cell
is 0.122 µm with ratio 1.15 to r = 39.44 mm, then a uniform 0.2005 mm out to
r = 75.53 mm, then ratio 1.05 to the inlet. Circumferentially, Δφ = 0.25°.

### Compiling on Overleaf

Don't paste a file below another document's preamble. Use either:

* **New Project → Blank Project**, replace all of `main.tex` with the file,
  and compile with **pdfLaTeX**, or
* **Upload** the `.tex` and set it as the main document.

The files need only `tikz`, `helvet` and optionally `newtxsf`. If `newtxsf`
is missing it is skipped, and Greek letters fall back to the default math
font.

### What is drawn to scale

Units are mm, with the origin at the cylinder centre.

* Geometry is from the mesh `sa_mi_baseline.su2`: R = 38.1 mm, inlet arc
  150.39 mm, outlets at φ = ±65°, no-slip wall |φ| ≤ 45° (`cylinder`), slip
  wall 45° < |φ| ≤ 65° (`cyl_lee`).
* Shock positions come from the no-field shock trace
  (`data/bow_shock_no_field.csv`), with a Billig-type far-field fit.
* The incident shock leaves the inlet split at β = 18.1°. Oblique-shock
  relations (γ = 1.4) give θ = 12.49° and M₂ = 5.252.
* Supersonic jet (`jet_geometry.py`): the M₂ flow crosses the transmitted
  shock (−26.1°), so the jet leaves at −15.1° with M ≈ 2.4 and is 4.5 mm wide.
* Grid: the mesh is an exact polar O-grid, so it is drawn with
  arcs and straight rays at the node positions in `data/mesh_r_i.txt` and
  `data/mesh_phi_j.txt`. It has 294 × 520 cells and Δφ = 0.25°. Radially:
  first cell 0.12 µm, ratio 1.15 to r = 39.44 mm, uniform 0.20 mm for
  39.6 ≤ r ≤ 75.5 mm (the shock layer), then ratio 1.05 to the inlet.

Schematic only: the jet's arch, the jet bow shock stand-off and the
shock-cell pattern. The cells are drawn steeper than the jet Mach angle
(≈ 25°) for legibility.

### Regenerating

    cd figures
    pip install numpy scipy matplotlib
    python3 check_mesh.py sa_mi_baseline.su2   # verifies the mesh, writes data/mesh_*.txt
    python3 make_figure.py
    python3 make_mesh_figure.py
    pdflatex fig_mesh.tex
    pdflatex fig_setup.tex
    pdflatex fig_setup_clean.tex

Case parameters and shock fits are in `fig_setup_typeiv.py`, and the jet is
in `jet_geometry.py`. Edit those, not the numbers in the `.tex` files.
