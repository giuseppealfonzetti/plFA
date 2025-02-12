% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_R}
\alias{get_R}
\title{Extract latent correlation matrix from unconstrained parameter vector}
\usage{
get_R(THETA, NTHR, NLOAD, NCORR, NVAR, Q)
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
Latent correlation matrix
}
\description{
Extract latent correlation matrix from unconstrained parameter vector
}
