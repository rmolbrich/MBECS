\name{NEWS}
\title{News for Package \pkg{MBECS}}

\section{Changes in version 0.99.8, Bioconductor Pre-Release}{
    \itemize{
        \item Included NEWS file.
        \item Added Bioconductor installation instructions in vignette.
        \item Added package man-page.
        \item Set lazyData to \code{false}.
        \item Added URL and BugReports fields to Description.
        \item Remove bapred and permute packages from suggest because they are
        no longer required.
        \item Also removed pals package, but included reference because I use
        the tableau color-scheme.
        \item All packages that are required for execution are now 'Depends'
        instead of 'Suggests'.
        \item Revised vignette content, formatting and spelling for more
        convenient user experience.
        \item Fixed test for percentile normalization to actually test
        PN-function.
        \item Prelim report now checks if clr/tss transformed values are present
        prior to calculating them.
        \item Report-names/directories can now be changed from default by user.
        \item Reporting functions included in testing.
        \item Formatted code to adhere to 4-space indentation and 80 characters
        width requirements. For the most part.
        \item Included code for generation of dummy-data.
    }
}


\section{Changes in version 0.99.10, Bioconductor Pre-Release}{

    \itemize{
        \item Added BiocViews: ReportWriting, Visualization, Normalization and
        QualityControl.
        \item Inlcuded 'importFrom' for all ggplot2 functionality.
        \item Required packages now live in Imports section.
        \item Switched to 'aes_string' to remove visible binding NOTEs in plot
        functions.
        \item Included local variables to shut up the remaining NOTEs.
        \item Moved 'match.arg' choices to function heads.
        \item Moved code for generation of mock-up dummy-data from data-raw
        into the data.R file.
        \item Added testing for dummy data.
        \item Fixed highlighting in vignette.
        \item Fixed gray-area issue in box-plots.

    }
}


\section{Changes in version 0.99.11, Bioconductor Pre-Release}{

    \itemize{
        \item Cleaned repository.
        \item GitHub Readme includes Installation instructions and brief
        workflow tutorial.

    }
}







