% sn-jnl-compat.cls -- LOCAL PREVIEW SHIM ONLY. NOT the official Springer Nature class.
% It emulates the sn-jnl front-matter macros on top of `article` so that main.tex
% compiles when the official template is absent. For submission, download the
% Springer Nature LaTeX template (sn-jnl.cls + sn-mathphys-num.bst) and place
% them in this folder; main.tex then uses them automatically.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{sn-jnl-compat}[2026/10/26 preview shim emulating sn-jnl]
\LoadClass[11pt]{article}
\RequirePackage[a4paper,margin=2.5cm]{geometry}
\RequirePackage[numbers,sort&compress]{natbib}
\RequirePackage{graphicx,amsmath,amssymb,booktabs,xcolor}
\RequirePackage[hidelinks]{hyperref}
\RequirePackage{lineno}
\linenumbers
\def\@affils{}
\newcommand{\affil}{\@ifstar\snc@affil\snc@affil}
\newcommand{\snc@affil}[2][]{\g@addto@macro\@affils{\par\noindent\textsuperscript{#1}#2}}
\newcommand{\orgdiv}[1]{#1}\newcommand{\orgname}[1]{#1}\newcommand{\orgaddress}[1]{#1}
\newcommand{\street}[1]{#1}\newcommand{\city}[1]{#1}\newcommand{\postcode}[1]{#1}
\newcommand{\state}[1]{#1}\newcommand{\country}[1]{#1}
\def\@emails{}
\newcommand{\email}[1]{\g@addto@macro\@emails{\par\noindent Corresponding author: \texttt{#1}}}
\renewcommand{\author}{\@ifstar\snc@author\snc@author}
\newcommand{\snc@author}[2][]{\g@addto@macro\@authlist{#2\textsuperscript{#1} }}
\def\@authlist{}
\newcommand{\fnm}[1]{#1}\newcommand{\sur}[1]{#1}
\newcommand{\equalcont}[1]{}
\renewcommand{\title}[2][]{\gdef\@title{#2}}
\def\@abstr{}\def\@kw{}
\renewcommand{\abstract}[1]{\gdef\@abstr{#1}}
\newcommand{\keywords}[1]{\gdef\@kw{#1}}
\renewcommand{\maketitle}{%
  \begin{center}{\LARGE\bfseries\@title\par}\vskip 1em{\large\@authlist\par}\end{center}
  {\small\@affils\par\@emails\par}\vskip 1em
  \noindent\textbf{Abstract}\par\noindent\@abstr\par\vskip .5em
  \noindent\textbf{Keywords:} \@kw\par\vskip 1em}
\newcommand{\bmhead}[1]{\paragraph*{#1}}
\newcommand{\backmatter}{}
\endinput
