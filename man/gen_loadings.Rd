% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataGeneration.R
\name{gen_loadings}
\alias{gen_loadings}
\title{Construct loading matrix}
\usage{
gen_loadings(FIXED = NULL, CONSTRMAT, SEED = 123, LB = 0, UB = 1)
}
\arguments{
\item{FIXED}{Fixed value to assign to all free loadings. If 'NULL' it draws
them randomly from Unif(0,1)}

\item{CONSTRMAT}{Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
that the corresponding element in the loading matrix is free to be estimated.
A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
to \eqn{0}.}

\item{SEED}{Random seed.}

\item{LB}{Lower bound for uniform random generator. Default set to 0.}

\item{UB}{Upper bound for uniform random generator. Default set to 1.}
}
\description{
gen_loadings() construct a (possibly) random loading matrix following the
constraints passed.
}
