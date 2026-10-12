# Figures

Each `.tex` is a complete standalone LaTeX document drawn at its final print
size. Each `.pdf` is the compiled vector figure, and each `.png` is a 600 dpi
preview.

| file | content | print size |
|---|---|---|
| `fig_domain` | computational domain and boundary conditions, panel (a) | 84 × 88 mm |
| `fig_typeiv` | Edney type IV interaction close-up, panel (b) | 83 × 88 mm |
| `fig_typeiv_detailed` | same, with the jet's shock cells and Mach-number regions, panel (b) | 83 × 88 mm |
| `fig_mesh` | computational grid: (a) full grid, (b) shock layer, (c) near-wall cells | 170 × 114 mm |

`fig_domain` and `fig_typeiv*` are each about 82 mm wide, so they sit side by
side in a 170 mm text width at full size, with 9 pt labels:

```latex
\begin{figure*}
  \centering
  \includegraphics{fig_domain.pdf}\hfill
  \includegraphics{fig_typeiv_detailed.pdf}
  \caption{(a) Computational domain and boundary conditions; the dashed box is
  the region shown in (b). (b) Edney type IV interaction: solid lines inside
  the jet are compressions, dashed lines expansion fans.}
\end{figure*}
```

The panel letters (a)/(b) are drawn in the figures. To remove one, change the
line `\newcommand{\panellabel}{(a)}` near the top of that `.tex` to
`\newcommand{\panellabel}{}`.

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
* Grid (`fig_mesh`): the mesh is an exact polar O-grid, so it is drawn with
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
    python3 make_figures.py
    for f in fig_*.tex; do pdflatex $f; done

Case parameters and shock fits are in `fig_setup_typeiv.py`, and the jet is
in `jet_geometry.py`. Edit those, not the numbers in the `.tex` files.
