# High Bandwidth Flash (HBF) for LLM Inference — manuscript

- `HBF_LLM_Inference_Paper.pdf` / `.docx` / `.doc` — the manuscript (30 pages, 16 figures, 15 tables)
- `hbf_model.py` — analytical model (Eqs. 2–7); `make_figures.py` — Figs. 1–6 + `results.json`
- `hw/` — near-flash multiplier study (Section 7): RTL, testbenches, synthesis/power flow, system model; see `hw/README.md`
- `build_paper.js` + `hw_section.js` — regenerate the .docx (`node build_paper.js`), then
  `soffice --headless --convert-to pdf|doc HBF_LLM_Inference_Paper.docx`
