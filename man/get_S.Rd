% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_S}
\alias{get_S}
\title{Get latent correlation matrix from theta}
\usage{
get_S(THETA, Q)
}
\arguments{
\item{THETA}{Numerical vector of parameters.}

\item{Q}{Number of latent variables.}
}
\description{
get_S() extracts the latent correlation matrix from theta assuming
theta elements to be reparametrised following the
Lewandowski-Kurowicka-Joe (2009) transform.
}
