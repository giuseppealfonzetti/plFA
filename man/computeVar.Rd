% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allMethods.R
\name{computeVar}
\alias{computeVar}
\title{Compute variance of the estimates}
\usage{
computeVar(OBJ, DATA, NUMDERIV = F, OPTION = "transformed")
}
\arguments{
\item{OBJ}{Object of class plFaFit}

\item{DATA}{Original data}

\item{NUMDERIV}{TRUE if the hessian must be computed using \link[numDeriv]{jacobian}}

\item{OPTION}{\code{transformed} if correlations are of interest. \code{raw} for inference on the reparametrised level.}
}
\description{
Compute variance of the estimates
}
