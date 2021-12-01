% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbecs_corrections.R
\name{mbecSVD}
\alias{mbecSVD}
\title{Singular Value Decomposition (SVD)}
\usage{
mbecSVD(input.obj, model.vars)
}
\arguments{
\item{input.obj}{phyloseq object or numeric matrix (correct orientation is
handeled internally)}

\item{model.vars}{two covariates of interest to select by first variable
selects panels and second one determines coloring}
}
\value{
A matrix of batch-effect corrected counts
}
\description{
Basically perform matrix factorization and compute singular eigenvectors
(SEV). Assume that the first SEV captures the batch-effect and remove this
effect from the data. The interesting thing is that this works pretty well.
But since the SEVs are latent factors that are (most likely) confounded with
other effects it is not obvious to me that this is the optimal approach to
solve this issue.
}
\details{
ToDo: IF I find the time to works on "my-own-approach" then this is the
point to start from!!!

The function returns a matrix of corrected counts Correct orientation of
counts will be handled internally.
}
\keyword{BECA}
\keyword{Duvallet}
\keyword{Normalisation}
\keyword{Percentile}