"""Build main.docx from main.tex with pandoc.

Pandoc's LaTeX reader drops tabularx tables, siunitx and mhchem macros and
table cross-references, so this script rewrites them into plain LaTeX first.
Usage:  python build_docx.py   (needs pandoc on PATH)
"""

import os
import re
import subprocess

UNITS = {r"\celsius": "°C", r"\second": "s", r"\watt": "W", r"\minute": "min"}
CHEM = {"H2O2": "H₂O₂", "NO2-": "NO₂⁻", "N2": "N₂", "Ca^2+": "Ca²⁺", "O3": "O₃", "NO_x": "NOₓ"}


def unit(u):
    parts = [p for p in re.split(r"(\\per)", u) if p]
    out = ""
    for p in parts:
        out += "/" if p == r"\per" else UNITS.get(p.strip(), p)
    return out


def main():
    src = open("main.tex", encoding="utf-8").read()
    body = src[src.index(r"\section{Introduction}"):]

    # Number tables in order of appearance and resolve their references.
    labels = re.findall(r"\\label\{(tab:[^}]+)\}", body)
    for i, lab in enumerate(labels, 1):
        body = body.replace(r"\ref{" + lab + "}", str(i))

    body = re.sub(r"\\SI\{([^}]*)\}\{([^}]*)\}", lambda m: f"{m[1]} {unit(m[2])}", body)
    body = re.sub(r"\\si\{([^}]*)\}", lambda m: unit(m[1]), body)
    body = re.sub(r"\\ce\{([^}]*)\}", lambda m: CHEM.get(m[1], m[1]), body)
    def colspec(m):
        spec = re.sub(r"p\{[^}]*\}", "l", m[1]).replace("@{}", "")
        return r"\begin{tabular}{" + "l" * len(spec) + "}"

    body = re.sub(r"\\begin\{tabularx\}\{\\textwidth\}\{(.*)\}\s*$", colspec, body, flags=re.M)
    body = body.replace(".pdf}", ".png}")
    body = body.replace(r"\end{tabularx}", r"\end{tabular}")
    body = body.replace(r"\small", "")

    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", src, re.S)[1].strip()
    highlights = re.findall(r"\\item (.*)", re.search(r"\\begin\{highlights\}(.*?)\\end\{highlights\}", src, re.S)[1])
    keywords = re.search(r"\\begin\{keyword\}(.*?)\\end\{keyword\}", src, re.S)[1].replace(r"\sep", ";").strip()
    title = re.search(r"\\title\{(.*?)\}\n", src)[1]

    front = (
        "\\documentclass{article}\n\\begin{document}\n"
        f"\\section*{{{title}}}\n\n"
        "M. Nisha Angeline\\textsuperscript{*}\n\n"
        "Department of Electronics and Communication Engineering, Velalar College of Engineering and "
        "Technology, Thindal, Erode 638012, Tamil Nadu, India\n\n"
        "\\textsuperscript{*}Corresponding author. E-mail: nishavlsidesign@gmail.com. ORCID: [to be added]\n\n"
        "\\section*{Highlights}\n\\begin{itemize}\n" + "".join(f"\\item {h}\n" for h in highlights) + "\\end{itemize}\n\n"
        f"\\section*{{Abstract}}\n{abstract}\n\n"
        f"\\textbf{{Keywords:}} {keywords}\n\n"
    )
    tmp = front + body.replace(r"\bibliographystyle{elsarticle-harv}", "").replace(
        r"\bibliography{references}", "\\section*{References}\n")
    open("_docx_src.tex", "w", encoding="utf-8").write(tmp)
    subprocess.run(["pandoc", "_docx_src.tex", "--citeproc", "--bibliography=references.bib",
                    "--resource-path=.:figures", "-o", "main.docx"], check=True)
    os.remove("_docx_src.tex")


if __name__ == "__main__":
    main()
