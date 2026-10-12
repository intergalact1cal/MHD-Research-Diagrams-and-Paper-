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
    python3 make_figure.py
    pdflatex fig_setup.tex
    pdflatex fig_setup_clean.tex

Case parameters and shock fits are in `fig_setup_typeiv.py`, and the jet is
in `jet_geometry.py`. Edit those, not the numbers in the `.tex` files.
