# Figures

## Set-up and Edney type IV structure (TikZ)

| file | content |
|---|---|
| `fig_setup_typeiv_tikz.tex` / `.pdf` / `.png` | (a) domain and boundary conditions, (b) type IV close-up |
| `fig_setup_typeiv_detailed_tikz.tex` / `.pdf` / `.png` | same, with the jet's shock cells and expansion fans, Mach-number regions and a line key |

### Compiling on Overleaf

Each `.tex` is a complete standalone document. Don't paste it below another
document's preamble. Use either:

* **New Project → Blank Project**, delete everything in `main.tex`, paste the
  whole file, and compile with **pdfLaTeX** (Menu → Compiler), or
* **Upload** the `.tex` and set it as the main document (Menu → Main document).

It needs only standard packages: `tikz`, `helvet`, and optionally `newtxsf`.
If `newtxsf` is missing it is skipped, and Greek letters fall back to the
default math font. To put the figure in a paper, compile it once and use the
PDF: `\includegraphics[width=\textwidth]{fig_setup_typeiv_tikz.pdf}`.

### What is drawn to scale

Units are mm, with the origin at the cylinder centre.

* Cylinder R = 38.1 mm, inlet arc 150 mm, outlets at φ = ±65°, no-slip wall
  |φ| ≤ 45° with slip wall beyond, inlet split at φ = −12.25°.
* Bow, transmitted and lower shocks come from the no-field shock trace
  (`data/bow_shock_no_field.csv`, digitised from `fig_mhd_structure.pdf`),
  with a Billig-type far-field fit.
* The incident shock leaves the inlet split at β = 18.1°. Oblique-shock
  relations (γ = 1.4) give θ = 12.49° and M₂ = 5.252.
* **Supersonic jet** (`jet_geometry.py`): the M₂ flow crosses the transmitted
  shock (inclined at −26.1°, i.e. a 38.6° wave angle). That gives a 27.6°
  deflection, so the jet leaves the triple points at −15.1° with M ≈ 2.4. Its
  width, 4.5 mm, is the triple-point separation normal to that direction.
  This is why the jet is long and thin, as in the classical type IV sketches.

Schematic only, drawn after the classical type IV sketch (Edney 1968): the
jet's arch (a circular arc turning 36°), the jet bow shock (a short, nearly
normal shock closing the jet about 4 mm off the wall), and the shear layers,
which turn sharply at its ends and run along the wall either side of the
impingement point. The detailed version also shows the jet's shock cells:
compressions converging on one boundary and an expansion fan spreading from
that point, repeated, starting at the lower triple point. The cells are drawn
steeper than the jet Mach angle (≈ 25°) for legibility.

### Regenerating

    cd figures
    pip install numpy scipy matplotlib
    python3 make_tikz.py          # rewrites both .tex files
    pdflatex fig_setup_typeiv_tikz.tex
    pdflatex fig_setup_typeiv_detailed_tikz.tex

`fig_setup_typeiv.py` holds the case parameters and shock fits, and
`jet_geometry.py` holds the jet. Edit those, not the numbers in the `.tex`.
