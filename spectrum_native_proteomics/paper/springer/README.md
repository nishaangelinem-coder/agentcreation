# Springer Nature submission package

| File | Purpose |
|---|---|
| `sn-article.tex` | **Submission manuscript.** Uses the official Springer Nature LaTeX template class (`\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}`) |
| `references.bib` | BibTeX database (33 references) |
| `figures/fig1–fig7.png` | Figures, 300 dpi |
| `sn-article-preview.pdf` | Local preview typeset with `sn-preview.cls` |
| `sn-preview.cls` | Stand-in class; **preview only, do not submit** |
| `build_springer.py` | Regenerates both `.tex` files from the Markdown sources in `../` |

## Compile for submission

1. Download the official template from Springer Nature (Overleaf: "Springer Nature LaTeX Template"; or the template zip linked from the journal's submission guidelines).
2. Copy `sn-jnl.cls` and `sn-mathphys-num.bst` next to `sn-article.tex`.
3. Build:
   ```bash
   pdflatex sn-article && bibtex sn-article && pdflatex sn-article && pdflatex sn-article
   ```

For a journal that uses another reference style (e.g. Nature Portfolio: `sn-nature`), change the class option.

Upload `sn-article.tex`, `references.bib` (or the generated `.bbl`) and the `figures/` folder. Most Springer Nature journals accept LaTeX source and also accept a PDF for initial submission.

## Preview build (no sn-jnl.cls needed)

```bash
pdflatex sn-article-preview && bibtex sn-article-preview && pdflatex sn-article-preview && pdflatex sn-article-preview
```
