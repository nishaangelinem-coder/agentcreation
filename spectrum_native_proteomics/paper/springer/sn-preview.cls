%% sn-preview.cls -- LOCAL PREVIEW ONLY. A minimal stand-in that implements the
%% front-matter macros of the official Springer Nature LaTeX template (sn-jnl.cls)
%% so that sn-article.tex can be typeset where sn-jnl.cls is unavailable.
%% For submission compile sn-article.tex with the official sn-jnl.cls.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{sn-preview}[2026/10/31 preview stand-in for sn-jnl]
\DeclareOption*{}
\ProcessOptions\relax
\LoadClass[10pt,a4paper]{article}
\RequirePackage[a4paper,top=2.6cm,bottom=2.6cm,left=2.7cm,right=2.7cm]{geometry}
\RequirePackage{graphicx,amsmath,amssymb,booktabs,xcolor,url}
\RequirePackage[numbers,sort&compress]{natbib}
\RequirePackage[hidelinks]{hyperref}
\RequirePackage[font=small,labelfont=bf,labelsep=space]{caption}
\RequirePackage{titlesec}
\RequirePackage{mathptmx}
\RequirePackage[scaled=0.9]{helvet}
\captionsetup[figure]{name=Fig.,labelsep=period}
\captionsetup[table]{labelsep=newline,justification=raggedright,singlelinecheck=false}
\titleformat{\section}{\normalfont\large\bfseries}{\thesection}{1em}{}
\titleformat{\subsection}{\normalfont\normalsize\bfseries}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalfont\normalsize\itshape}{\thesubsubsection}{1em}{}
\setlength{\parindent}{1em}\setlength{\parskip}{0pt}
\newcommand\botrule{\bottomrule}
\newcommand\fnm[1]{#1}\newcommand\sur[1]{#1}\newcommand\spfx[1]{#1}\newcommand\pfx[1]{#1}\newcommand\sfx[1]{#1}
\newcommand\orgdiv[1]{#1}\newcommand\orgname[1]{#1}\newcommand\orgaddress[1]{#1}
\newcommand\street[1]{#1}\newcommand\city[1]{#1}\newcommand\postcode[1]{#1}\newcommand\state[1]{#1}\newcommand\country[1]{#1}
\newcommand\equalcont[1]{}
\newcommand\jyear[1]{}\newcommand\theyear{}
\def\@snauthors{}\def\@snaffils{}\def\@snemails{}\def\@snabstract{}\def\@snkeywords{}
\renewcommand\title[2][]{\gdef\@title{#2}}
\newcommand\email[1]{\g@addto@macro\@snemails{#1}}
\renewcommand\author{\@ifstar{\@snauthorS}{\@snauthorN}}
\newcommand\@snauthorS[2][]{\g@addto@macro\@snauthors{#2\textsuperscript{#1*}}}
\newcommand\@snauthorN[2][]{\g@addto@macro\@snauthors{#2\textsuperscript{#1}}}
\newcommand\affil{\@ifstar{\@snaffil}{\@snaffil}}
\newcommand\@snaffil[2][]{\g@addto@macro\@snaffils{\par\noindent\textsuperscript{#1}#2}}
\renewcommand\abstract[1]{\gdef\@snabstract{#1}}
\newcommand\keywords[1]{\gdef\@snkeywords{#1}}
\newcommand\bmhead[1]{\par\medskip\noindent\textbf{#1}\par\smallskip}
\newcommand\backmatter{}
\renewcommand\maketitle{%
  \thispagestyle{plain}%
  \noindent{\color[HTML]{1F2A44}\rule{\textwidth}{1.2pt}}\par\vspace{6pt}
  \noindent{\sffamily\fontsize{17}{21}\selectfont\bfseries \@title\par}\vspace{10pt}
  \noindent{\large\@snauthors\par}\vspace{6pt}
  {\small\raggedright\@snaffils\par\medskip\noindent *Corresponding author(s). E-mail(s): \@snemails\par}\vspace{10pt}
  \noindent{\color[HTML]{1F2A44}\rule{\textwidth}{0.6pt}}\par\vspace{4pt}
  {\small\noindent\textbf{Abstract}\par\smallskip\noindent\@snabstract\par\medskip
  \noindent\textbf{Keywords:} \@snkeywords\par}
  \vspace{4pt}\noindent{\color[HTML]{1F2A44}\rule{\textwidth}{0.6pt}}\par\vspace{12pt}}
\endinput
