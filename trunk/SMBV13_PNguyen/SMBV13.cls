%
% LLNCS DOCUMENT CLASS -- version 2.10
% for LaTeX2e
%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{SMBV13}[2008/28/05 Seminar medizinische Bildverarbeitung an der RWTH-Aachen]

\typeout{**********************************************}
\typeout{* Seminar Medizinische Bildverarbeitung 2013 *}
\typeout{**********************************************}

%\let\@langsel\false
%
%\DeclareOption{deutsch}{
%	\let\@langsel\true
%	\PassOptionsToPackage{english,ngerman}{babel}}
%\DeclareOption{englisch}{
%	\let\@langsel\true
%	\PassOptionsToPackage{english}{babel}}
%
%\DeclareOption*{\typeout{\CurrentOption\ not defined for class SemMbv09}}
%
%\ExecuteOptions{deutsch}
%
%\ifx\@langsel\false
%	\ExecuteOptions{englisch}
%\fi
%
%\ProcessOptions

\LoadClass[twoside,a4paper,10pt]{article}

\RequirePackage[latin1]{inputenc}
\RequirePackage[english, ngerman]{babel}
\RequirePackage{graphicx}
\RequirePackage[T1]{fontenc}

\renewcommand\section{\@startsection{section}{1}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {12\p@ \@plus 4\p@ \@minus 4\p@}%
                       {\normalfont\large\bfseries\boldmath
                        \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {8\p@ \@plus 4\p@ \@minus 4\p@}%
                       {\normalfont\normalsize\bfseries\boldmath
                        \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {-0.5em \@plus -0.22em \@minus -0.1em}%
                       {\normalfont\normalsize\bfseries\boldmath}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                       {-12\p@ \@plus -4\p@ \@minus -4\p@}%
                       {-0.5em \@plus -0.22em \@minus -0.1em}%
                       {\normalfont\normalsize\itshape}}
\renewcommand\subparagraph[1]{\typeout{SemMbv12 warning: You should not use
                  \string\subparagraph\space with this class}\vskip0.5cm
You should not use \verb|\subparagraph| with this class.\vskip0.5cm}

\DeclareMathSymbol{\Gamma}{\mathalpha}{letters}{"00}
\DeclareMathSymbol{\Delta}{\mathalpha}{letters}{"01}
\DeclareMathSymbol{\Theta}{\mathalpha}{letters}{"02}
\DeclareMathSymbol{\Lambda}{\mathalpha}{letters}{"03}
\DeclareMathSymbol{\Xi}{\mathalpha}{letters}{"04}
\DeclareMathSymbol{\Pi}{\mathalpha}{letters}{"05}
\DeclareMathSymbol{\Sigma}{\mathalpha}{letters}{"06}
\DeclareMathSymbol{\Upsilon}{\mathalpha}{letters}{"07}
\DeclareMathSymbol{\Phi}{\mathalpha}{letters}{"08}
\DeclareMathSymbol{\Psi}{\mathalpha}{letters}{"09}
\DeclareMathSymbol{\Omega}{\mathalpha}{letters}{"0A}

\let\footnotesize\small

\AtBeginDocument{\RequirePackage{multicol}}		% needed for the list of participants, index
\AtBeginDocument{\RequirePackage{ae}}
%\AtBeginDocument{\RequirePackage[T1]{fontenc}}% T1 font for correct use of umlauts
\AtBeginDocument{\RequirePackage{eurosym}}		% Euro Symbol
\AtBeginDocument{\RequirePackage{amsmath}}
\AtBeginDocument{\RequirePackage{amsfonts}}
\AtBeginDocument{\RequirePackage{amssymb}}
\AtBeginDocument{\RequirePackage{amsxtra}}
\AtBeginDocument{\RequirePackage{epsfig}}
%\AtBeginDocument{\RequirePackage{graphics}}
%\AtBeginDocument{\RequirePackage{graphicx}}
\AtBeginDocument{\RequirePackage{subfigure}}
\AtBeginDocument{\RequirePackage{algorithm}}
\AtBeginDocument{\RequirePackage{algorithmic}}

\AtBeginDocument{\raggedbottom}

\pagestyle{empty}
%
\frenchspacing
\widowpenalty=10000
\clubpenalty=10000

\setlength\textheight{242mm}
\setlength\textwidth{160mm}
\setlength\oddsidemargin{0mm}
\setlength\evensidemargin{0mm}
\setlength\unitlength{1mm}
\setlength\topmargin{-15mm}

\AtBeginDocument{
\def\abbrsymbname{Abk\"urzungs- und Symbolverzeichnis}%
\def\abstractname{Zusammenfassung}%
\def\ackname{Danksagung.}%
\def\andname{und}%
\def\lastandname{ und}%
\def\appendixname{Anhang}%
\def\bibname{Literaturverzeichnis}%
\def\chaptername{Kapitel}%
\def\claimname{Behauptung}%
\def\conjecturename{Hypothese}%
\def\contentsname{Inhalt}%
\def\corollaryname{Korollar}%
%\def\definitionname{Definition}%
\def\examplename{Beispiel}%
\def\exercisename{\"Ubung}%
\def\figurename{Abb.}%
\def\keywordname{Keywords:}%
\def\indexname{Index}%
%\def\lemmaname{Lemma}%
\def\contriblistname{Mitarbeiter}%
\def\listfigurename{Abbildungsverzeichnis}%
\def\listtablename{Tabellenverzeichnis}%
\def\mailname{{\it Correspondence to\/}:}%
\def\noteaddname{Nachtrag}%
\def\notename{Anmerkung}%
\def\partname{Teil}%
\def\prefacename{Vorwort}%
%\def\problemname{Problem}%
\def\proofname{Beweis}%
\def\propertyname{Eigenschaft}%
%\def\propositionname{Proposition}%
\def\questionname{Frage}%
\def\refname{Literaturverzeichnis}%
\def\remarkname{Anmerkung}%
\def\seename{siehe}%
\def\solutionname{L\"osung}%
\def\subclassname{{\it Subject Classifications\/}:}%
\def\tablename{Tab.}%
}

\def\keywords#1{\vskip5mm\noindent{\bf\keywordname\ }#1\par}

\def\fnum@figure{\textbf{\figurename\thinspace\thefigure}}
\def\fnum@table{\textbf{\tablename\thinspace\thetable}}
\renewenvironment{figure}
               {\setcounter{subfigure}{0}
                \@float{figure}}
               {\end@float}
\renewenvironment{table}
               {\setcounter{subtable}{0}
                \@float{table}}
               {\end@float}


\renewenvironment{abstract}
	{\noindent{\Large\bf\abstractname} 
	\vskip10pt%
	\noindent\it}%
	{\par}

\renewcommand\maketitle{%
	\thispagestyle{empty}
  \begin{flushleft}%
 		\noindent
    {\Large \bf \@title \par}%
    \vskip 5mm%
    \noindent
    {\Large von\ \@author \par}%
  \end{flushleft}\par
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
  \tableofcontents
  \newpage
}