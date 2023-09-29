
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plFA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/plFA)](https://CRAN.R-project.org/package=plFA)
<!-- badges: end -->

The plFA package allows the estimation of confirmatory factor models for
ordinal data using stochastic and numeric pairwise likelihood
optimisation.

## Installation

You can install the development version of plFAfrom
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/plFA")
```

## Example

First, generate a synthetic dataset

``` r
library(plFA)
set.seed(123)

# p = number of items, q = number of latent variables, n = number of observations
p <- 30; q <- 5; n <- 5000L

# Simple loading matrix constraints
A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')

# Draw some random loadings according to constraints
Load <- gen_loadings(CONSTRMAT = A)

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)

# Generate random latent correlation matrix
S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q)

# simulate the data
D <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)

# categories per item
cat <- apply(D, 2, max) + 1

# store true parameter vector (interpretable parameterisation)
theta <- getPar(get_theta(rep(thr, p), Load, S, cat, A), C = sum(cat), P = p, Q = q, CONSTRMAT = A)
```

Estimate the model with the stochastic optimiser

``` r

# Set options for cpp stochastic optimiser
cpp_ctrl <- list(
  MAXT = 1000,
  PAIRS_PER_ITERATION = 8
)

# Fit the model with stochastic optimiser
stFit <- fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'hyper',
  CONTROL = cpp_ctrl,
  ITERATIONS_SUBSET = seq(0, cpp_ctrl$MAXT, 100)
)
#> 1. Initialising at default values
#> 2. Optimising with hyper...
#> 
#> 3. Rearranging output...
#> 4. Done! (0.9 secs)

# extract estimated parameter vector
stPar <- getPar(stFit)

# extract list of parameter estimates
stParList <- getPar(stFit, OPTION = 'list')
matrixcalc::is.positive.definite(stParList$latent_correlations)
#> [1] TRUE

# mean square error
mean((stPar-theta)^2)
#> [1] 0.0006180528
```

Numerical estimation as comparison

``` r
numFit <- fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'=1),
  METHOD = 'ucminf')
#> 1. Initialising at default values
#> 2. Optimising with ucminf...
#> 3. Done! (20.63 secs)

# extract estimated parameter vector
numPar <- getPar(numFit)

# extract list of parameter estimates
numParList <- getPar(numFit, OPTION = 'list')

# mean square error
mean((numPar-theta)^2)
#> [1] 0.0005450474
```
