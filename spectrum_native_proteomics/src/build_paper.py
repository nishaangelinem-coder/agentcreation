"""Assemble manuscript Markdown, then DOCX (pandoc) / PDF (LibreOffice) / LaTeX."""
import re, subprocess, os
os.chdir('paper')
tabs = {}
cur = None
for block in open('tables.md').read().split('\n\n'):
    m = re.match(r'\*\*Table (\d+)\.', block)
    if m: cur = m.group(1); tabs[cur] = block
    elif cur and block.strip(): tabs[cur] += '\n\n' + block
body = open('part2_results.md').read()
for k, v in tabs.items():
    body = body.replace('{{TABLE%s}}' % k, v)
md = open('part0_front.md').read() + '\n\n' + open('part1_intro_methods.md').read() + '\n\n' + body + '\n\n' + open('part3_refs.md').read()
open('manuscript.md', 'w').write(md)
subprocess.run(['pandoc', 'manuscript.md', '-o', 'manuscript.docx', '--resource-path=.:..', '--reference-doc=reference.docx'] if os.path.exists('reference.docx')
               else ['pandoc', 'manuscript.md', '-o', 'manuscript.docx', '--resource-path=.:..'], check=True)
subprocess.run(['pandoc', 'manuscript.md', '-s', '-o', 'manuscript.tex', '--resource-path=.:..'], check=True)
print('built')
