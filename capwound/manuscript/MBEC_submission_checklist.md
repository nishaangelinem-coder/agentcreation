# Medical & Biological Engineering & Computing (MBEC): submission checklist

Applied in `main.tex`:

- [x] Original Article with separate Introduction, Methods, Results, Discussion and Conclusion sections
- [x] Unstructured abstract of no more than 200 words (currently 196) covering rationale, objective, methods, findings and conclusions, with few abbreviations
- [x] Up to 5 keywords, all MeSH terms: Plasma gases, Electric impedance, Dielectric spectroscopy, Biofilms, Wound infection
- [x] Title page shows the main-text and abstract word counts and the number of figures and tables (produced by `wordcount.py` during `build.sh`)
- [x] Graphical abstract (mandatory): `figures/graphical_abstract.pdf` / `.png`. Also upload it as a separate file.
- [x] Author biography of no more than 200 characters (currently 174). Add a photograph if requested.
- [x] Numbered references in square brackets (`sn-mathphys-num`)
- [x] Declarations: funding, competing interests, ethics, data/code availability, AI-tool use, author contribution
- [x] Figures are vector PDFs; 600-dpi PNG copies are in `../figures/`

To do before upload:

- [ ] Copy the official `sn-jnl.cls` and `sn-mathphys-num.bst` from the Springer Nature LaTeX template into this folder and rebuild. Then check that the page count is 18 or fewer (the preview is 15 pages).
- [ ] Re-check the current MBEC guidelines (https://link.springer.com/journal/11517/submission-guidelines), especially the graphical-abstract size and file format and the reference style. These could not be read directly from this environment.
- [ ] Add co-authors, funding and the code-archive DOI (Zenodo).
- [ ] Verify the DOIs of the references cited from memory (Cole 1941, Asami 2002, Schönleber 2014, Mafart 2002, Graves 2012, Lu 2016, Brug 1984, Igreja 2004, Wold 2001, Markx 1999, Virtanen 2020).
- [ ] Cover letter: state clearly that this is a computational study and propose experimental validation.
