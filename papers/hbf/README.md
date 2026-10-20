# High Bandwidth Flash (HBF) for LLM Inference — manuscript

- `HBF_LLM_Inference_Paper.pdf` / `.docx` / `.doc` — the manuscript
- `hbf_model.py` — analytical model (Eqs. 2–7); `make_figures.py` — figures + `results.json`
- `build_paper.js` — regenerates the .docx (`node build_paper.js`), then
  `soffice --headless --convert-to pdf|doc HBF_LLM_Inference_Paper.docx`
