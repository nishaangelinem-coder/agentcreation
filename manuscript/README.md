# Closed-loop CAP review: Q1 submission package

**Manuscript:** *Closed-loop cold atmospheric plasma delivery with embedded sensing: a critical review of sensing, control and integration, with an in-silico analysis of where dose error originates.*
**Author:** Dr. M. Nisha Angeline, Velalar College of Engineering and Technology.
**Format:** Elsevier review article (`elsarticle`, author–year). Set up for *Biosensors and Bioelectronics*; change `\journal{}` in `main.tex` for another venue.

| File | What it is |
|---|---|
| `main.pdf` | Typeset manuscript (double-spaced review format with line numbers) |
| `main.tex`, `references.bib` | LaTeX source and BibTeX (50 references) |
| `main.docx` | Word version for journals that require an editable file (`python build_docx.py`) |
| `figures/` | Fig. 1–3 as vector PDF and 600 dpi TIFF/PNG (`python figures/make_figures.py`) |
| `supplementary/` | S1 simulation code and its JSON output |
| `submission/` | Cover letter and highlights as separate files |
| `provenance/source_manuscript_export.md` | Export of the earlier draft, including the integrity audit (Appendix A). Internal only; do not submit. |

## Rebuild

```
python supplementary/S1_closed_loop_cap_simulation.py --out supplementary/S1_simulation_results.json
python figures/make_figures.py
latexmk -pdf main.tex
python build_docx.py
```

## Changes from the earlier draft

- **S1 re-implemented.** The original script was not available, so it was rewritten from the documented specification. The new script is deterministic, but it uses a different random stream, so the Section 7 numbers changed slightly. All text now quotes the new output: open-loop error 90% → 89%, closed-loop error 16% → 17%, degraded sensor 51% → 50%, open-loop safety exceedance 29.6% → 31.5%, degraded-sensor exceedance 11.5% → 13.9%, open-loop duration 52.4 s → 52.6 s. The conclusions are unchanged. If you still have the original script, either use it instead or keep these numbers. Do not mix the two sets.
- **New Fig. 1** (four-layer architecture). The simulation figure is now Fig. 2 and the roadmap is Fig. 3.
- **Appendix A (integrity audit) removed** from the manuscript. A summary of the screening counts is now in Section 2, and the full audit is kept in `provenance/`.

## Before submission (still open)

- [ ] **Complete the references.** 42 `.bib` entries are marked `% TODO-complete`: their author lists end in "and others", or their volume/pages are missing. Crossref was not reachable from the build environment, so nothing was guessed. In Zotero, use *Add item by identifier* with each DOI, export to BibTeX, and replace the entries.
- [ ] **Boehm & Bourke (2018).** Confirm DOI `10.1515/hsz-2018-0222`.
- [ ] Confirm the venue's current quartile (JCR/SJR), word and reference limits, and whether it needs a graphical abstract (Fig. 1 can serve as one).
- [ ] Add your ORCID, any co-authors who meet the CRediT/ICMJE criteria, and funding/acknowledgements.
- [ ] Run iThenticate.
- [ ] Deposit `supplementary/` on Zenodo and add the DOI to the Data availability statement.
- [ ] Optional: expand to 80–120 references through a documented Scopus/WoS search. Call it a systematic review only if you add a PRISMA flow and dual screening.
