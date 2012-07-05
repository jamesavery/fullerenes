%%
%% This is file `standalone.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% standalone.dtx  (with options: `cls')
%% 
%% Copyright (c) 2010 by Martin Scharrer <martin@scharrer-online.de>
%% 
%% This work may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License, either version 1.3c
%% of this license or (at your option) any later version.
%% The latest version of this license is in
%% 
%%    http://www.latex-project.org/lppl.txt
%% 
%% and version 1.3c or later is part of all distributions of LaTeX
%% version 2008/05/04 or later.
%% 
%% This work has the LPPL maintenance status `maintained'.
%% 
%% The Current Maintainer of this work is Martin Scharrer.
%% 
%% This work consists of the files standalone.dtx, standalone.ins
%% and the derived file standalone.sty.
%% 
% $Id: standalone.dtx 2451 2011-04-07 11:40:38Z martin $
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{standalone}   [2011/04/07 v0.4a Class to compile TeX sub-files standalone]
\newif\ifstandalone
\standalonetrue
\newif\ifstandalonebeamer
\standalonebeamerfalse
\let\onlyifstandalone\@firstofone
\RequirePackage{kvoptions}
\SetupKeyvalOptions{prefix=sa@}
\DeclareBoolOption[true]{preview}
\DeclareStringOption{border}
\let\sa@border\relax
\DeclareVoidOption{beamer}{%
  \def\sa@class{beamer}%
  \sa@previewfalse
  \standalonebeamertrue
}
\DeclareStringOption[article]{class}
\def\sa@classoptionslist{}
\DeclareDefaultOption{%
  \xdef\sa@classoptionslist{\sa@classoptionslist,\CurrentOption}%
}
\ProcessKeyvalOptions*\relax
\let\@classoptionslist\sa@classoptionslist
\begingroup
\xdef\@tempa{[\sa@classoptionslist]{\sa@class}}
\expandafter
\endgroup
\expandafter\LoadClass\@tempa
\let\standalone\empty
\let\endstandalone\relax
\input{standalone.cfg}
\ifsa@preview
  \@ifundefined{endstandalone}{%
    \renewenvironment{standalone}
      {\preview   }
      {\endpreview}
  }{}%
  \RequirePackage{preview}
  \ifx\sa@border\relax\else
    \begingroup
    \def\rem@bp#1bp\relax#2\@nnil{#1}%
    \def\default@bp#1{\expandafter\rem@bp\the\dimexpr#1 bp\relax bp\relax\@nnil}%
    \def\sa@readborder#1 #2 #3 #4 {%
       \ifx\\#2#3#4\\%
         \@defaultunits\global\PreviewBorder=\dimexpr#1 bp\relax\@nnil%
       \else
       \ifx\\#4\\%
        \xdef\PreviewBbAdjust{-\default@bp{#1} -\default@bp{#2} \default@bp{#1} \default@bp{#2}}%
       \else
        \xdef\PreviewBbAdjust{-\default@bp{#1} -\default@bp{#2} \default@bp{#3} \default@bp{#4}}%
       \fi\fi
    }%
    \@firstofone{\expandafter\sa@readborder\sa@border} {} {} {} {} \relax
    \endgroup
  \fi
\fi
\ifstandalonebeamer
\newenvironment{standaloneframe}{%
  \@ifnextchar<%
    {\@standaloneframe}%
    {\@@standaloneframe{}}%
}{\end{frame}}%
\def\@standaloneframe<#1>{%
    \@@standaloneframe{<#1>}%
}
\def\@@standaloneframe#1{%
  \@ifnextchar[%]
    {\@@@standaloneframe{#1}}%
    {\@@@standaloneframe{#1}[]}%
}%
\def\@@@standaloneframe#1[{%
  \@ifnextchar<%
    {\@@@@standaloneframe{#1}[}%
    {\@@@@@@standaloneframe{#1}[}%
}%
\def\@@@@standaloneframe#1[#2]{%
  \@ifnextchar[%]
    {\@@@@@standaloneframe{#1}{#2}}%
    {\begin{frame}#1[#2][environment=standaloneframe]}%
}%
\def\@@@@@standaloneframe#1#2[#3]{%
    \begin{frame}#1[#2][environment=standaloneframe,#3]%
}%
\def\@@@@@@standaloneframe#1[#2]{%
    \begin{frame}#1[environment=standaloneframe,#2]%
}%
\fi
\let\sa@cls@orig@document\document
\let\sa@cls@orig@enddocument\enddocument
\def\document{%
  \sa@cls@orig@document
  \let\documentclass\sa@documentclass  % TODO: really required?
  \sa@cls@afterbegindocument
}
\def\enddocument{%
  \sa@cls@beforeenddocument
  \sa@cls@orig@enddocument
}
\def\sa@cls@afterbegindocument{\expandafter\standalone\ignorespaces}
\def\sa@cls@beforeenddocument{\unskip\endstandalone}
\endinput
%%
%% End of file `standalone.cls'.
