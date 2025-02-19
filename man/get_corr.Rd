% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_corr}
\alias{get_corr}
\title{Extract vector of latent correlations from unconstrained parameter vector}
\usage{
get_corr(THETA, NTHR, NLOAD, NCORR, NVAR, Q)
}
\arguments{
\item{THETA}{Unconstrained parameter vector}

\item{NTHR}{Number of thresholds parameters.}

\item{NLOAD}{Number of free loadings parameters}

\item{NCORR}{Number of free latent correlations parameters.}

\item{NVAR}{Number of free latent variance parameters.}

\item{Q}{Number of latent variables}
}
\value{
Vector of latent correlations
}
\description{
Extract vector of latent correlations from unconstrained parameter vector
}
