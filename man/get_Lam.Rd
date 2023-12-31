% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_Lam}
\alias{get_Lam}
\title{Get loading matrix from theta}
\usage{
get_Lam(A, C, THETA)
}
\arguments{
\item{A}{Binary matrix of dimension \eqn{p*q} where \eqn{p} is the number
of items and \eqn{q} the number of latent variables. Entries equal to
\eqn{1} refer to free loadings, while entries equal to \eqn{0} indicate
loadings constrained to be null.}

\item{C}{Sum of the number of categories for each item.}

\item{THETA}{Numerical vector of parameters.}
}
\description{
get_Lam() constructs the loading matrix from the parameter vector
}
