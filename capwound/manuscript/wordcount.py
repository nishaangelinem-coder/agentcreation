"""Approximate word counts (abstract, main text) for the MBEC title page."""
import re
from pathlib import Path

tex = Path("main.tex").read_text()


def words(t):
    t = re.sub(r"(?<!\\)%.*", "", t)
    t = re.sub(r"\\begin\{(figure\*?|table|sidewaystable|equation|center)\}.*?\\end\{\1\}", " ", t, flags=re.S)
    t = re.sub(r"\$[^$]*\$", " x ", t)
    t = re.sub(r"\\(cite|ref|label|eqref|includegraphics|inputtable|graphicspath)\{[^}]*\}", " ", t)
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", t)
    t = re.sub(r"[{}~\\]", " ", t)
    return len(re.findall(r"[A-Za-z0-9][\w\-–.,%/]*", t))


a = tex[tex.index(r"\abstract{") + 10: tex.index(r"\keywords")]
m = tex[tex.index(r"\section{Introduction}"): tex.index(r"\backmatter")]
Path("tables").mkdir(exist_ok=True)
Path("tables/wordcount_abstract.txt").write_text(str(words(a)))
Path("tables/wordcount_main.txt").write_text(str(words(m)))
print("abstract", words(a), "main", words(m))
