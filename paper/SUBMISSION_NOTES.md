# Submission checklist (ICISPES 2027)

- Build: `pdflatex main && bibtex main && pdflatex main && pdflatex main` (numbers are auto-filled by `python -m fdia_dt.build_paper` from `results/digest.json`).
- Edit `abstract.tex`, `sec_method.tex`, `sec_results.tex`, `sens_text.tex`; do not edit the `*_filled.tex` files by hand.
- Confirm the page limit and template on the ICISPES 2027 call for papers (draft is 6 pages, IEEEtran conference).
- Add co-authors, ORCID and funding acknowledgement if required.
- Two references need page numbers that could not be verified online from this session: [Esmalifalak et al., IEEE Syst. J. 2017] and [Boyaci et al., IEEE Syst. J. 2022, DOI 10.1109/JSYST.2021.3109082]. Please complete them from IEEE Xplore.
- Run the institutional plagiarism/similarity check before submission.
