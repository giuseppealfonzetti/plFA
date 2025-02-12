% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfa.R
\name{cfa}
\alias{cfa}
\title{Fit Confirmatory Factor Analysis Models}
\usage{
cfa(
  model,
  data,
  std.lv = FALSE,
  estimator = "PML",
  estimator.args = list(method = c("ucminf", "SA"), init_method = c("SA", "custom",
    "standard"), cpp_control_init = NULL, ncores = 1, valdata = NULL, computevar_numderiv
    = FALSE),
  start = NULL,
  control = list(),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
is described using the lavaan model syntax. See \code{\link[lavaan:model.syntax]{lavaan::model.syntax}}
for more information. Alternatively, a parameter table (eg. the output of
the \code{\link[lavaan:model.syntax]{lavaan::lavaanify()}} function) is also accepted.}

\item{data}{A data frame containing the observed variables used in the model.
Variables must be declared as ordered factors.}

\item{std.lv}{If \code{TRUE}, the metric of each latent variable is determined by
fixing their (residual) variances to 1.0. If \code{FALSE}, the metric of each
latent variable is determined by fixing the factor loading of the first
indicator to 1.0.}

\item{estimator}{The estimator is \code{PML}. If any other estimator is provided,
then \code{\link[lavaan:cfa]{lavaan::cfa()}} is used.}

\item{estimator.args}{A list of arguments for \code{\link[=fit_plFA]{fit_plFA()}}--see the help file
for more details. Possible options are:
\describe{
\item{\code{method}}{One of \code{"ucminf"} (default) for the quasi-Newton BFGS optimiser, or \code{"SA"} for stochastic approximation.}
\item{\code{init_method}}{One of \code{"SA"} (default) for stochastic approximation, \code{"custom"} for user-defined starting values, or \code{"standard"} for standard starting values.}
\item{\code{cpp_control_init}}{A list of control parameters for the initialisation algorithm.}
\item{\code{ncores}}{The number of cores to use for parallel computation.}
\item{\code{valdata}}{Validation data.}
\item{\code{computevar_numderiv}}{If \code{TRUE}, the asymptotic variance-covariance matrix is computed using numerical derivatives.}
}}

\item{start}{A vector of starting values to use (in the order of free
loadings, thresholds, and then factor correlations). If not provided, the
starting point is computed according to \code{\link[=fit_plFA]{fit_plFA()}}'s \code{INIT_METHOD}.}

\item{control}{A list of control parameters for the estimation algorithm. See
\code{\link[=fit_plFA]{fit_plFA()}} for more information.}

\item{verbose}{If \code{TRUE}, print additional information during the estimation
process.}

\item{...}{Additional arguments to be passed to \code{\link[=lavaan]{lavaan()}}.}
}
\value{
A \code{plFAlavaan} object, which is a subclass of the \code{lavaan} class.
Therefore, all methods available for \code{lavaan} objects are expected to be
compatible with \code{plFAlavaan} objects.
}
\description{
Fit a Confirmatory Factor Analysis (CFA) model using the \code{lavaan.pl()} engine
using the \code{\link{lavaan}} framework. By default, a quasi-Newton BFGS optimiser
from the \code{\link{ucminf}} package is used. For large data sets, the stochastic
approximation algorithm can be used.
}
\details{
Not all \code{\link{lavaan}} options can be used at present. Some options of
interest are:
\describe{
\item{\code{information}}{The information matrix to use. Only \code{"observed"} is currently supported.}
\item{\code{se}}{The standard error method to use. Only \code{"robust.huber.white"} is currently supported.}
\item{\code{test}}{No GOF tests are available as of now, so this is set to \code{"none"}.}
}
}
\examples{

# A simple binary factor model using the LSAT data
fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
summary(fit)

}
\references{
Katsikatsou, M., Moustaki, I., Yang-Wallentin, F., & Jöreskog, K. G. (2012).
Pairwise likelihood estimation for factor analysis models with ordinal data.
\emph{Computational Statistics & Data Analysis}, \emph{56}(12), 4243--4258.
\url{https://doi.org/10.1016/j.csda.2012.04.010}

Alfonzetti, G., Bellio, R., Chen, Y., & Moustaki, I. (2025). Pairwise
stochastic approximation for confirmatory factor analysis of categorical
data. \emph{British Journal of Mathematical and Statistical Psychology}, \emph{78}(1),
22--43. \url{https://doi.org/10.1111/bmsp.12347}
}
