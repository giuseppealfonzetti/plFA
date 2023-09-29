% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allMethods.R
\name{getThetaPath}
\alias{getThetaPath}
\title{Extract parameters trajectory along the optimisation}
\usage{
getThetaPath(OBJ, LAB)
}
\arguments{
\item{OBJ}{Object of class StoFit or PlFaFit}

\item{LAB}{Can take values \code{'pathTheta'} for raw trajectories or \code{'pathAvTheta'} for averaged ones. Set by default at \code{pathAvTheta}.}
}
\description{
Extract parameters trajectory along the optimisation
}