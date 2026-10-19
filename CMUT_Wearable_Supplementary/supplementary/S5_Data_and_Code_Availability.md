# S5. Data and code availability, and ethics statements (ready to paste)

## Data availability statement
All data supporting this study are simulation outputs. They include
acoustic field maps, profiles, thermal and displacement fields, controller
traces, and the Monte Carlo digital-twin outputs. They are available in
the Supplementary Information and in [REPOSITORY, e.g. Zenodo] at
https://doi.org/[DOI to be assigned]. No experimental or human/animal data
were generated.

## Code availability statement
The simulation code (`cmutsim` Python package and `run_all.py`)
reproduces every figure and table with a fixed random seed. It is
available under the MIT licence at [GitHub URL] and archived at
https://doi.org/[DOI]. Dependencies are listed in `requirements.txt`.

## Statement on emulated hardware results
"Results labelled *emulated* (Supplementary Note S4) are Monte Carlo
predictions from a hardware digital twin. They are not experimental
measurements. No device was fabricated in this study."

## Ethics
Not applicable. The study involves no human participants, animals or
biological material. The proposed Phase III spheroid work will use
commercially sourced cell lines and will be submitted for institutional
biosafety approval before it starts.

## Author contributions (CRediT) — edit as appropriate
* **M. N. A.:** conceptualisation, methodology, software, formal
  analysis, writing (original draft).
* **S. K. M.:** biomedical device requirements, validation planning,
  writing (review and editing).

## Competing interests
The authors declare no competing interests.

## Before deposit: checklist
- [ ] Create a GitHub repository; push `CMUT_Wearable_Supplementary/`.
- [ ] Connect it to Zenodo and make a release to obtain a DOI. Replace the
      placeholders above.
- [ ] Run `python run_all.py` once more from a clean environment and
      confirm that `results/summary.json` is unchanged.
- [ ] Apply the fixes listed in S2 (items A1–A8) to the manuscript.
- [ ] Upload `S1_Supplementary_Methods_and_Results.docx` as the
      Supplementary Information file.
