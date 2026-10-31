"""Convert the Markdown manuscript into Springer Nature LaTeX (sn-jnl template)."""
import re, subprocess, os, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
P = '..'
p1 = open(f'{P}/part1_intro_methods.md').read()
p2 = open(f'{P}/part2_results.md').read()
intro, methods = p1.split('# Materials and methods', 1)
intro = intro.replace('# Introduction', '', 1)
results, rest = p2.split('# Discussion', 1)
results = results.replace('# Results', '', 1)
discussion, rest = rest.split('# Conclusions', 1)
conclusion = rest.split('# Data and code availability', 1)[0]

# tables from tables.md
tabs, cur = {}, None
for block in open(f'{P}/tables.md').read().split('\n\n'):
    m = re.match(r'\*\*Table (\d+)\.', block)
    if m: cur = m.group(1); tabs[cur] = block
    elif cur and block.strip(): tabs[cur] += '\n\n' + block
for k, v in tabs.items(): results = results.replace('{{TABLE%s}}' % k, v)

KEYS = ['aebersold2016', 'cox2008', 'elias2007', 'lazar2016', 'mann2021', 'perez2022', 'wang2018', 'vaswani2017', 'yilmaz2024',
        'eloff2025', 'yang2024', 'lee2024', 'lapin2025', 'ebrahimi2023', 'tran2017', 'qiao2021', 'krull2026', 'qiu2025',
        'bittremieux2022', 'gessulat2019', 'zeng2022', 'gorishniy2021', 'devlin2019', 'he2022', 'dancik1999', 'ma2003',
        'shrikumar2017', 'maaten2008', 'paszke2019', 'pedregosa2011', 'dai2024', 'dai2021', 'cox2011']


def cite(m):
    nums = []
    for part in m.group(1).split(','):
        if re.search('[–-]', part):
            a, b = re.split('[–-]', part); nums += list(range(int(a), int(b) + 1))
        else: nums.append(int(part))
    if any(n < 1 or n > len(KEYS) for n in nums): return m.group(0)
    return '\\cite{' + ','.join(KEYS[n - 1] for n in nums) + '}'


def inline(s):
    """Markdown inline -> LaTeX for table cells and captions."""
    s = s.replace('\\', '\\textbackslash{}')
    for a, b in [('&', '\\&'), ('%', '\\%'), ('#', '\\#'), ('_', '\\_'), ('{', '\\{'), ('}', '\\}'), ('$', '\\$')]:
        s = s.replace(a, b)
    s = s.replace('\\textbackslash\\{\\}', '\\textbackslash{}')
    s = re.sub(r'`([^`]+)`', r'\\texttt{\1}', s)
    s = re.sub(r'\*\*\*([^*]+)\*\*\*', r'\\textbf{\\emph{\1}}', s)
    s = re.sub(r'\*\*([^*]+)\*\*', r'\\textbf{\1}', s)
    s = re.sub(r'\*([^*]+)\*', r'\\emph{\1}', s)
    s = re.sub(r'~([^~\s]+)~', r'\\textsubscript{\1}', s)
    s = re.sub(r'\^([^\^\s]+)\^', r'\\textsuperscript{\1}', s)
    s = re.sub(r'\[(\d{1,2}(?:[,–-]\d{1,2})*)\]', cite, s)
    return s


def table_tex(lines, caption=None, label=None):
    rows = [[c.strip() for c in l.strip().strip('|').split('|')] for l in lines if not re.match(r'^\|[-|: ]+\|$', l.strip())]
    nc = len(rows[0])
    small = '\\footnotesize' if nc >= 5 else '\\small'
    lens = [max(min(max(len(r[j]) for r in rows), 45), 6) for j in range(nc)]
    w = [nc * x / sum(lens) for x in lens]
    spec = '@{}' + ''.join('>{\\raggedright\\arraybackslash\\hsize=%.3f\\hsize}X' % x for x in w) + '@{}'
    body = [' & '.join(inline(c) for c in rows[0]) + ' \\\\', '\\midrule']
    body += [' & '.join(inline(c) for c in r) + ' \\\\' for r in rows[1:]]
    width = '\\textwidth' if caption is not None else '0.7\\textwidth'
    tab = '\\setlength{\\tabcolsep}{%s}\n\\begin{tabularx}{%s}{%s}\n\\toprule\n%s\n\\botrule\n\\end{tabularx}' % (
        '3pt' if nc >= 6 else '4pt', width, spec, '\n'.join(body))
    if caption is None:
        return '\n\\begin{center}%s\n%s\n\\end{center}\n' % (small, tab)
    return '\n\\begin{table}[htbp]\n%s\n\\caption{%s}\\label{%s}\n%s\n\\end{table}\n' % (small, caption, label, tab)


def convert_tables(md):
    out, lines, i = [], md.split('\n'), 0
    pending_cap = None
    while i < len(lines):
        l = lines[i]
        m = re.match(r'^\*\*Table (\d+)\. (.*)\*\*\s*$', l)
        if m:
            pending_cap = (m.group(1), m.group(2)); i += 1; continue
        if l.strip().startswith('|'):
            blk = []
            while i < len(lines) and lines[i].strip().startswith('|'): blk.append(lines[i]); i += 1
            if pending_cap:
                n, cap = pending_cap
                cap = inline(cap.rstrip('.')) + '.'
                out.append('```{=latex}' + table_tex(blk, cap, 'tab:%s' % n) + '```')
                pending_cap = None
            else:
                out.append('```{=latex}' + table_tex(blk) + '```')
            continue
        if pending_cap and l.strip():
            # caption paragraph not followed directly by a table
            pass
        out.append(l); i += 1
    return '\n'.join(out)


def convert_figs(md):
    def f(m):
        cap, path = m.group(1), m.group(2)
        n = re.search(r'Fig\. (\d+)\.', cap).group(1)
        cap = re.sub(r'^\*\*Fig\. \d+\. ([^*]+)\*\*\s*', lambda q: '**' + q.group(1).strip() + '** ', cap)
        fn = os.path.basename(path)
        return ('```{=latex}\n\\begin{figure}[htbp]\n\\centering\n\\includegraphics[width=\\textwidth]{figures/%s}\n'
                '\\caption{%s}\\label{fig:%s}\n\\end{figure}\n```' % (fn, inline(cap), n))
    return re.sub(r'!\[(.*?)\]\(([^)]+)\)\{width=100%\}', f, md, flags=re.S)


def refs(md):
    md = re.sub(r'\bFig\. (\d+)', r'Fig.~\\ref{fig:\1}', md)
    md = re.sub(r'\bTable (\d+)', r'Table~\\ref{tab:\1}', md)
    md = re.sub(r'\bTables (\d+)–(\d+)', r'Tables~\\ref{tab:\1}--\\ref{tab:\2}', md)
    md = re.sub(r'\[(\d{1,2}(?:[,–-]\d{1,2})*)\]', cite, md)
    return md


def to_latex(md, shift=0):
    md = convert_figs(md)
    md = convert_tables(md)
    md = refs(md)
    if shift: md = re.sub(r'^(#+) ', lambda m: '#' * (len(m.group(1)) - shift) + ' ', md, flags=re.M)
    r = subprocess.run(['pandoc', '-f', 'markdown+raw_tex+raw_attribute', '-t', 'latex', '--wrap=none'],
                       input=md, capture_output=True, text=True, check=True)
    tex = r.stdout
    tex = re.sub(r'\\hypertarget\{[^}]*\}\{%\n', '', tex)
    tex = re.sub(r'\}\\label\{[^}]*\}\}', '}', tex)
    tex = tex.replace('\\tightlist\n', '')
    return tex


body = {
    'intro': to_latex(intro),
    'results': to_latex(results),
    'discussion': to_latex(discussion),
    'conclusion': to_latex(conclusion),
    'methods': to_latex(methods),
}

TITLE = 'From identification-first to spectrum-native proteomics: a ground-truth benchmark with real-data validation'
ABSTRACT = (
    "Mass-spectrometry proteomics is shifting from identification-first workflows, in which tandem (MS2) spectra are matched "
    "to a sequence database before any statistics are done, to spectrum-native models that learn directly from MS1 and MS2 signal. "
    "Because the ground truth of real spectra is itself defined by identification, how much information identification discards has been hard to quantify. "
    "We built an open, CPU-trainable virtual instrument with exact ground truth and used it to benchmark database search, a spectrum-graph algorithm and "
    "SpecNovo, a transformer for database-free sequencing, together with identification-first and identification-free phenotype classifiers, "
    "including a ProFormer-style transformer with self-supervised pretraining. A differentiable fragment-ion-evidence term raised SpecNovo peptide recall "
    "from 0.5\\% to 41.6\\%. Unlike database search, SpecNovo recovered variant and novel-proteome peptides (38.9\\% and 38.7\\%) and gave the "
    "best-calibrated confidence. Removing identification was the dominant effect for phenotype prediction (AUROC 0.975 vs 0.648 when half of the "
    "biomarkers were unidentifiable), whereas pretraining mainly improved transfer to unseen acquisition sites (0.728 vs 0.651). Validation on "
    "PRIDE-derived LC-MS/MS runs and MassIVE-KB spectra confirmed the bottleneck: only 0.19\\% of 25,613 real MS1 features were identified at 1\\% FDR. "
    "Calibrating the simulator to measured fragmentation statistics doubled SpecNovo's real-spectrum peptide recall from 7.6\\% to 15.1\\%, matching a "
    "spectrum-graph baseline in precision--coverage AUC. Spectrum-native learning recovers signal that identification-first pipelines lose, but at "
    "current data scales it depends on physics-informed inductive biases and on calibration to real instruments.")
KEYWORDS = 'De novo peptide sequencing, Transformer, Spectrum-to-phenotype, Self-supervised learning, Mass spectrometry, Benchmark'

UNI = {'≥': '\\ensuremath{\\geq}', '≤': '\\ensuremath{\\leq}', '±': '\\ensuremath{\\pm}', '×': '\\ensuremath{\\times}',
       '−': '\\ensuremath{-}', '→': '\\ensuremath{\\rightarrow}', '≈': '\\ensuremath{\\approx}', '∼': '\\ensuremath{\\sim}',
       '∈': '\\ensuremath{\\in}', 'λ': '\\ensuremath{\\lambda}', 'σ': '\\ensuremath{\\sigma}', 'π': '\\ensuremath{\\pi}',
       'Σ': '\\ensuremath{\\Sigma}', '·': '\\ensuremath{\\cdot}', '₂': '\\textsubscript{2}', '₃': '\\textsubscript{3}',
       '⁻': '\\textsuperscript{\\ensuremath{-}}', '¹': '\\textsuperscript{1}', '²': '\\textsuperscript{2}', '³': '\\textsuperscript{3}',
       '⁴': '\\textsuperscript{4}', '⁺': '\\textsuperscript{+}', '–': '--', '—': '---', '‘': '`', '’': "'", '“': '``', '”': "''",
       'ä': '\\"a', 'ö': '\\"o', 'ü': '\\"u', 'é': "\\'e", 'č': '\\v{c}', 'í': "\\'\\i{}", 'α': '\\ensuremath{\\alpha}',
       'μ': '\\ensuremath{\\mu}', 'Δ': '\\ensuremath{\\Delta}', '√': '\\ensuremath{\\surd}', '∞': '\\ensuremath{\\infty}',
       '⁰': '\\textsuperscript{0}', '₀': '\\textsubscript{0}', '₁': '\\textsubscript{1}', '\u2009': '\\,', '\u00a0': '~', '…': '\\ldots{}', 'Å': '\\AA{}'}


def preamble(cls):
    opts = '[pdflatex,sn-mathphys-num]' if cls == 'sn-jnl' else ''
    s = '\\documentclass%s{%s}\n' % (opts, cls)
    s += ('%% Springer Nature LaTeX template (sn-jnl). Reference style: sn-mathphys-num (numbered).\n'
          '\\usepackage{graphicx}\n\\usepackage{multirow}\n\\usepackage{amsmath,amssymb,amsfonts}\n\\usepackage{booktabs}\n\\usepackage{array,tabularx}\n\\usepackage{xurl}\n'
          '\\usepackage{textcomp}\n\\usepackage{xcolor}\n\\usepackage{newunicodechar}\n'
          '\\providecommand{\\pandocbounded}[1]{#1}\n')
    for k, v in UNI.items(): s += '\\newunicodechar{%s}{%s}\n' % (k, v)
    return s


def document(cls):
    t = preamble(cls)
    t += r"""
\raggedbottom

\begin{document}

\title[Spectrum-native proteomics benchmark]{%s}

\author*[1]{\fnm{M.} \sur{Nisha Angeline}}\email{nishavlsidesign@gmail.com}

\affil*[1]{\orgdiv{Department of Electronics and Communication Engineering}, \orgname{Velalar College of Engineering and Technology}, \orgaddress{\street{Thindal}, \city{Erode}, \postcode{638012}, \state{Tamil Nadu}, \country{India}}}

\abstract{%s}

\keywords{%s}

\maketitle

\section{Introduction}\label{sec:intro}
%s
\section{Results}\label{sec:results}
%s
\section{Discussion}\label{sec:discussion}
%s
\section{Conclusion}\label{sec:conclusion}
%s
\section{Methods}\label{sec:methods}
%s
\backmatter

\bmhead{Acknowledgements}
The author thanks the PRIDE/bigbio and Casanovo developers for making the test spectra used in the real-data validation publicly available.

\section*{Declarations}

\begin{itemize}
\item \textbf{Funding:} Not applicable.
\item \textbf{Conflict of interest/Competing interests:} The author declares no competing interests.
\item \textbf{Ethics approval and consent to participate:} Not applicable. The study uses simulated data and publicly available proteomics data; no human participants or animals were involved.
\item \textbf{Consent for publication:} Not applicable.
\item \textbf{Data availability:} Simulated corpora are regenerated deterministically by the released code. Real data: PRIDE-derived BSA LC-MS/MS runs with SDRF metadata (\url{https://github.com/bigbio/quantms-test-datasets}, \texttt{testdata/lfq\_ci/BSA}) and MassIVE-KB annotated spectra (\url{https://github.com/Noble-Lab/casanovo}, \texttt{sample\_data/}). All result files (JSON/pickle) are included in the code repository.
\item \textbf{Materials availability:} Not applicable.
\item \textbf{Code availability:} All source code (virtual instrument, models, baselines, real-data validation, evaluation, figure and manuscript scripts) and trained checkpoints are available in the project repository (\texttt{spectrum\_native\_proteomics/}); every experiment runs on a CPU with Python $\geq$ 3.11 and PyTorch $\geq$ 2.0.
\item \textbf{Author contribution:} M.N.A.: conceptualisation, methodology, software, investigation, formal analysis, visualisation, writing -- original draft, review and editing.
\end{itemize}

\bibliography{references}

\end{document}
""" % (TITLE, ABSTRACT, KEYWORDS, body['intro'], body['results'], body['discussion'], body['conclusion'], body['methods'])
    return t


open('sn-article.tex', 'w').write(document('sn-jnl'))
prev = document('sn-preview').replace('\\bibliography{references}', '\\bibliographystyle{unsrtnat}\n\\bibliography{references}')
open('sn-article-preview.tex', 'w').write(prev)
left = sorted({c for c in open('sn-article.tex').read() if ord(c) > 127} - set(UNI))
print('unmapped unicode:', left)
