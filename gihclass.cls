%
% This file is copyright (C) 2003 Sam Noble. It may be modified so long
% as my name is not removed and the modifier adds his name to the file.
% Redistribution permitted.
%
% 27 Jan 2004	Sam Noble	Removed tocbibind dependency.
% 04 Dec 2001 	Sam Noble	Class file
% 03 Sep 1995	David Perkinson	Title Page
% Acknowledgements Page, David Perkinson & Sam Noble
% May 2005 Patrick Carlisle Table of contents chapter definition
% 2004-2005 Ben Salzberg (BTS) a variety of tweaks here and in the template
% 2015 Chester Ismay (CII) modified to create R Markdown template
% 2019 Daniel Hammarström (DH) modified to fit GIH standards as part
% of the thesisdowngih package. 
% ----- The 2019 update borrows heavily from 
% PhDthesis class file by Andreas Solders 2011 based on
% PhDthesis version 2 style file by Jakob Suckale, 2007
% and CUEDthesis version 1 by Harish Bhanderi 2002


% Oddities:
        %
%   We *ASSUME* that any time \cleardoublepage is called
%  we actually want a blank back side with NO page number/heading
%
%   Minor bug -- seems to be a more general LaTeX thing:
        %  If you use \frontmatter \mainmatter without any chapters inbetween
%  be prepared to have the page numbering messed up. Not a big deal,
%  but I'm not sure how to fix it.
%
%

\NeedsTeXFormat{LaTeX2e}
%% Specify new class --- gihclass %% 
\ProvidesClass{gihclass}[2019/01/05 GIH thesis class]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions\relax
\LoadClass{book}
\RequirePackage{fancyhdr}
% Added by CII to produce nice looking code in BeginVerbatimEnvironment
\RequirePackage{fancyvrb}


%:--------- Borrowed from PhDthesisSU.cls ----------------------------- 


%:-------------------------- packages for fancy things -----------------------
\usepackage{mathptmx}
\usepackage[scaled=.90]{helvet}
\usepackage{courier}
% Default font for dissertations is Times.
%\usepackage{fourier} % If mathematics don't display well using Times, then use Fourier.
\usepackage{blindtext}
\usepackage[scaled=.90]{helvet}% Helvetica, served as a model for arial
\usepackage{titlesec}
\usepackage[letterspace=150]{microtype}
\usepackage{geometry}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{amsbsy}
%\usepackage{graphics} 
\usepackage[pdftex]{graphicx} % for improved inclusion of graphics
\usepackage{epstopdf} % converts eps figures to pdf
%\usepackage{wrapfig} % to include figure with text wrapping around it
\usepackage[margin=8pt,font=small,labelfont=bf]{caption} % for improved layout of figure captions with extra margin, smaller font than text
\usepackage{fancyhdr} % for better header layout
\usepackage{eucal}
\usepackage[english]{babel}
\usepackage[usenames, dvipsnames]{color}
\usepackage[perpage]{footmisc}
\usepackage{enumerate} % enumerated list for List of Papers
\usepackage[sort,compress]{natbib}
\usepackage{ifthen}
\usepackage{multicol} % for pages with multiple text columns, e.g. References
\setlength{\columnsep}{20pt} % space between columns; default 10pt quite narrow
\usepackage[nottoc]{tocbibind} % correct page numbers for bib in TOC, nottoc suppresses an entry for TOC itself
%\usepackage{nextpage}
\usepackage[T1]{fontenc}
\usepackage{calligra}
% For text positioning 
\usepackage{textpos}
\usepackage{changepage}




\renewcommand\boldsymbol[1]{\pmb{#1}}

 \DeclareFixedFont{\timesEightteen}{\encodingdefault}{\sfdefault}{\seriesdefault}{\shapedefault}{17\p@}

\graphicspath{{artwork/}}




%
%
%:------ Removed content compared to reedthesis.cls to clean up file (DH)  
%
%

%:----------------From PhDThesisSU:  Glossary/Abbrev./Symbols -----------------------
% Not sure what this does ...

\usepackage[intoc]{nomencl} % load nomencl extension; include in TOC
%\nomrefpage % to include page numbers after abbrevations
\renewcommand{\nomlabel}[1]{\textbf{#1}} % make abbreviations bold
\makenomenclature % used to be \makeglossary
\newcommand{\g}{\footnote{For all abbreviations see the glossary on page \pageref{nom}.}} % type "\g" to refer to glossary




\setlength{\emergencystretch}{3em}  % prevent overfull lines
\providecommand{\tightlist}{%
\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
% End of CII add

% This gives us rules below the headers
\AtBeginDocument{%
 \fancyhf{}
 \fancyfoot[LE,RO]{\thepage}
%  \fancyhead[RE]{\slshape \leftmark}
%  \fancyhead[LO]{\slshape \rightmark}
%	The above makes your headers in all caps. If you would like different headers, choose one of the following options (be sure to remove the % symbol from both the right and left headers):
% \fancyfoot[RE]{\slshape \nouppercase} % This makes the headers on the RIGHT side pages be italic and use lowercase With Capitals When Specified.
% \fancyfoot[LO]{\slshape \nouppercase} % This does the same thing to the LEFT side pages

\renewcommand{\headrulewidth}{0pt}
%	or
%	\fancyhead[RE]{\scshape \leftmark} % The RIGHT headers will be in small caps.
%	\fancyhead[LO]{\scshape \rightmark} %  And so will the LEFT headers
 \pagestyle{fancy}






% Psych majors: You do not need the following six lines, as it conflicts with apacite, so comment them out.
\let\oldthebibliography=\thebibliography
\let\endoldthebibliography=\endthebibliography
\renewenvironment{thebibliography}[1]{
\oldthebibliography{#1}
\addcontentsline{toc}{chapter}{\bibname}
}{\endoldthebibliography}
%%%%%% end of things for psych majors to comment out

\let\oldtheindex=\theindex
\let\endoldtheindex=\endtheindex
\renewenvironment{theindex}{
\oldtheindex
\addcontentsline{toc}{chapter}{\indexname}
}{\endoldtheindex}
}

% Stolen from book.cls and modified
\let\RToldchapter\chapter
\renewcommand{\chapter}{\if@openright\RTcleardoublepage
\else\clearpage\fi
\thispagestyle{empty}%
\global\@topnum\z@
\@afterindentfalse
\secdef\@chapter\@schapter}

% Stolen from book.cls PBC 5/12/05
% % Using this to actually show "Chapter 1" in TOC instead of "1"
% \def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
% \if@mainmatter
% \refstepcounter{chapter}%
% \typeout{\@chapapp\space\thechapter.}%
% \addcontentsline{toc}{chapter}%
% {\@chapapp\space\thechapter:\space#1}%
% \else
% \addcontentsline{toc}{chapter}{#1}%
% \fi
% \else
% \addcontentsline{toc}{chapter}{#1}%
% \fi
% \chaptermark{#1}%
% \addtocontents{lof}{\protect\addvspace{10\p@}}%
% \addtocontents{lot}{\protect\addvspace{10\p@}}%
% \if@twocolumn
% \@topnewpage[\@makechapterhead{#2}]%
% \else
% \@makechapterhead{#2}%
% \@afterheading
% \fi}

\newcommand{\RTcleardoublepage}{
\clearpage\if@twoside \ifodd\c@page\else
\thispagestyle{empty}\hbox{}\newpage
\if@twocolumn\hbox{}\newpage\fi\fi\fi}

\let\RToldcleardoublepage\cleardoublepage
\renewcommand{\cleardoublepage}{\RTcleardoublepage}

% adjust margins for binding  (changed 2007-04-24 tgp)
\setlength{\oddsidemargin}{.5in}
\setlength{\evensidemargin}{0in}
\setlength{\textwidth}{6.0in}
\setlength{\textheight}{9.0in}
\setlength\topmargin{0in}
\addtolength\topmargin{-\headheight}
\addtolength\topmargin{-\headsep}

%:-------------------- From SU template -----------------

%S5 settings (default) 

\geometry{paperheight=240mm,paperwidth=170mm}
\geometry{left=22.5mm,right=22.5mm}
\geometry{top=22.5mm,bottom=22.5mm}

\setlength{\headheight}{22.49pt} % fixes warning?


%Make al headlines non-fat
\def\section{\@startsection {section}{1}{\z@}{-3.5ex plus -1ex minus 
    -.2ex}{2.3ex plus .2ex}{\Large\sffamily}}
\def\subsection{\@startsection{subsection}{2}{\z@}{-3.25ex plus -1ex minus 
   -.2ex}{1.5ex plus .2ex}{\large\sffamily}}
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}{-3.25ex plus
-1ex minus -.2ex}{1.5ex plus .2ex}{\normalsize\sffamily}}


% Change font for text and titles (DH 2019-02-05)
% This is in conflict with SU template code below that removes "Chapter n"

% \titleformat{\chapter}[display]
%   {\normalfont\sffamily\huge}
%   {\chaptertitlename\ \thechapter}{16pt}{\Huge}
% \titleformat{\section}
%   {\normalfont\sffamily\Large\bfseries}
%   {\thesection}{1em}{}
% \titleformat{\subsection}{\normalsize\bfseries\itshape}{}{0pt}{}

% Removing the "Chapter n" text an replacing it by a single number.
% From SU template
\def\@makechapterhead#1{%
  \vspace{0\p@}%
  {\parindent \z@ \raggedright \normalfont\sffamily
    \ifnum \c@secnumdepth >\m@ne
      % Skipping the "Chapter n" text 
      \if@mainmatter
        \baselineskip=20pt\rightskip 2cm plus 3cm\huge\thechapter\hspace{0.5em} 
        %If a period is needed after chapters add . to after thechapter ex: "\thechapter.\ .
      \fi
    \fi
    \huge #1\par\nobreak
    \vskip 84\p@
  }}
\def\@makeschapterhead#1{%
  \vspace*{0\p@}%
  {\parindent \z@ \raggedright
    \normalfont\sffamily
    \interlinepenalty\@M
    \huge  #1\par\nobreak
    \vskip 84\p@
  }}





% Line spacing 
\renewcommand{\baselinestretch}{1.25}



%\setlength{\oddsidemargin}{.6in}
%\setlength{\evensidemargin}{0in}
%\setlength{\textwidth}{5.9in}
%\setlength\topmargin{0in}
%\addtolength\headheight{2.5pt}
%\addtolength\topmargin{-\headheight}
%\addtolength\topmargin{-\headsep}
%\addtolength\textheight{1in}
%\addtolength\textheight{\headheight}
%\addtolength\textheight{\headsep}


% Defining the artwork (crest)
\def\crest#1{\gdef\@crest{#1}}



\def\isbn#1{\gdef \@isbn{#1}}
\def\@isbn{\@latex@warning@no@line{No \noexpand\isbn given}}

\def\place#1{\gdef \@place{#1}}
\def\@place{\@latex@warning@no@line{No \noexpand\place given}}

\def\printeby#1{\gdef \@printeby{#1}}
\def\@printeby{\@latex@warning@no@line{No \noexpand\printeby given}}

\def\coverinfo#1{\gdef \@coverinfo{#1}}
\def\@coverinfo{\@latex@warning@no@line{No \noexpand\coverinfo given}}

\def\year#1{\gdef \@year{#1}}
\def\@year{\@latex@warning@no@line{No \noexpand\year given}}


\def\sernr#1{\gdef \@sernr{#1}}
\def\@sernr{\@latex@warning@no@line{No \noexpand\sernr given}}



\renewcommand\l@chapter[2]{%
\ifnum \c@tocdepth >\m@ne
\addpenalty{-\@highpenalty}%
\vskip 1.0em \@plus\p@
\setlength\@tempdima{1.5em}%
\begingroup
\parindent \z@ \rightskip \@pnumwidth
\parfillskip -\@pnumwidth
\leavevmode \bfseries
\advance\leftskip\@tempdima
\hskip -\leftskip
#1\nobreak\normalfont
\leaders\hbox{$\m@th \mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill
\nobreak\hb@xt@\@pnumwidth{\bfseries \hss #2}\par
\penalty\@highpenalty
\endgroup
\fi}

%% Added "defence information" use the Rmd file to write it up.
\newenvironment{defence}{%
\thispagestyle{empty}
\if@twocolumn
\@restonecoltrue\onecolumn
\else
\@restonecolfalse
\fi
\fontsize{8}{10}\selectfont
}{\clearpage \if@restonecol\twocolumn\fi}%



%%% Added by CII

% Removed acknowledgements in GIH version, goes as a separate chapter
% \newenvironment{acknowledgements}{%
% \if@twocolumn
% \@restonecoltrue\onecolumn
% \else
% \@restonecolfalse
% \fi
% \chapter*{Acknowledgements}
% \fontsize{12}{14}\selectfont
% }{\clearpage \if@restonecol\twocolumn\fi}%

\newenvironment{preface}{%
\if@twocolumn
\@restonecoltrue\onecolumn
\else
\@restonecolfalse
\fi
\chapter*{Preface}
\fontsize{10}{12}\selectfont
}{\clearpage \if@restonecol\twocolumn\fi}%
% End of CII add

\newenvironment{abstract}{%
\if@twocolumn
\@restonecoltrue\onecolumn
\else
\@restonecolfalse
\fi
\chapter*{Abstract}
\fontsize{10}{12}\selectfont
}{\clearpage \if@restonecol\twocolumn\fi}%


\newenvironment{listofpapers}{%
\if@twocolumn
\@restonecoltrue\onecolumn
\else
\@restonecolfalse
\fi
\chapter*{List of scientific papers}
\fontsize{10}{12}\selectfont
}{\clearpage \if@restonecol\twocolumn\fi}%




\newenvironment{dedication}{%
\thispagestyle{empty}
\if@twocolumn
\@restonecoltrue\onecolumn
\else
\@restonecolfalse
\flushright

\phantom{.}

\vspace{\stretch{1}}
\itshape
% \chapter*{Dedication} % removed dedication header 
\fontsize{12}{14}\selectfont
}{\clearpage \if@restonecol\twocolumn\fi}%


\ifx\@pdfoutput\@undefined
\newcommand{\RTpercent}{\@percentchar\space}
\AtBeginDvi{\special{!\RTpercent Reed College LaTeX Thesis Class 2001/12/04 SN}}
\AtBeginDvi{\special{rawpostscript \RTpercent Reed College LaTeX Thesis Class 2001/12/04 SN}}
\else
\AtEndDocument{\pdfinfo{/Creator	(Reed College LaTeX Thesis Class 2001/12/04 SN)}}
\fi

% I hacked the title page to all be the same font size
% as requested by the library, BTS 2005



\newcommand{\fromtop}[1]{%
 \dimexpr-1in-\topskip-\topmargin-\headheight-\headsep+#1\relax
}
\newcommand{\fromleft}[1]{%
 \dimexpr-1in-\oddsidemargin+#1\relax
}


\crest{\includegraphics[width=80mm]{eng_logoGIH_tryck}}

\renewcommand{\maketitle}{%
        {\pagenumbering{alph}
        \setcounter{page}{1}
        \thispagestyle{empty}
        \fontsize{12}{14}\selectfont
        \begin{titlepage}
       \frontmatter
       \pagestyle{plain}
% \thispagestyle{empty}

\begin{textblock*}{\paperwidth}(\fromleft{0cm},\fromtop{4cm})
\centering
\fontfamily{ptm}\lsstyle\normalsize{Avhandlingsserie f\"{o}r \\ Gymnastik- och idrottsh\"{o}gskolan}
\end{textblock*}


\begin{textblock*}{\paperwidth}(\fromleft{0cm},\fromtop{6cm})
\fontfamily{ptm}\centering
Nr~\@sernr
\end{textblock*}



\begin{textblock*}{\paperwidth}(\fromleft{0cm},\fromtop{8cm})

\begin{adjustwidth}{1.5cm}{1.5cm}
\centering
\fontfamily{ptm}\large\MakeUppercase{\@title}
\end{adjustwidth}
\end{textblock*}



\cleardoublepage
        
    \renewcommand{\footnotesize}{\small}
    \renewcommand{\footnoterule}{\relax}
   % \thispagestyle{empty}
{\large \vspace*{0mm} \hspace*{-5mm} {{\@crest} \par} } %Move the logo
    { \fontsize{20pt}{16pt}\fontfamily{ptm}\selectfont \vspace*{30mm}{\noindent}{\@title} \par \vspace*{8mm}}
    \bigskip
    \bigskip
    {\fontsize{16pt}{16pt}\fontfamily{ptm}\selectfont{\LARGE {\noindent}\@author} \par}
  \null\vfill


\pagebreak

%:------ Copyright information ----------- %%

\phantom{.}

\vspace{\stretch{1}}


{\scriptsize
\noindent
\copyright~ \@author % \@place~ \@year % Name of author, location year
\\*
\noindent
Gymnastik- och idrottsh\"{o}gskolan~\@year
\\*
\noindent
ISBN \@isbn % Provided by the library

\vspace{1mm}
\noindent
\@coverinfo % 
\\*
\noindent
Printed by: \@printeby %name of printing company
\\*
Distributor: Gymnastik- och idrottsh\"{o}gskolan  %name of department

\vspace{2mm}
}


        
        \end{titlepage}

   }
}

% From JSS
% Added by CII
\renewcommand{\title}[1]{\def\@title{#1}}
%\newcommand{\Shorttitle}[1]{\def\@Shorttitle{#1}}
\newcommand{\Plainauthor}[1]{\def\@Plainauthor{#1}}
\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
\newcommand{\Listofpapers}[1]{\def\@Listofpapers{#1}}
%\newcommand{\Keywords}[1]{\def\@Keywords{#1}}
%\newcommand{\Plainkeywords}[1]{\def\@Plainkeywords{#1}}
\newcommand{\Acknowledgements}[1]{\def\@Acknowledgements{#1}}
\newcommand{\Dedication}[1]{\def\@Dedication{#1}}
\newcommand{\Preface}[1]{\def\@Preface{#1}}
% End of CII add




