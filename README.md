
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
p <- 30; q <- 8; n <- 1000L

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
  MAXT = 2*n,
  PAIRS_PER_ITERATION = 8
)

# Fit the model with stochastic optimiser
stFit <- fit_plFA(
  DATA = D,
  VALDATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'hyper',
  CONTROL = cpp_ctrl,
  ITERATIONS_SUBSET = seq(0, cpp_ctrl$MAXT, 100)
)
#> 1. Initialising at default values
#> 2. Computing frequencies...
#> 3. Optimising with hyper...
#> 4. Done! (1.7 secs)
stFit
#> - Stochastic estimate
#> 
#>      Sampling scheme: hyper 
#>      Pairs per iteration: 8  out of  435 
#>      Iterations: 2000 
#>      Total time: 1.7 s (Data reduction: 0.02 s)
#>      Cores used: 1 
#>  
#> - Use getPar() to extract parameter estimates.
# extract estimated parameter vector
stPar <- getPar(stFit)

# extract list of parameter estimates
stParList <- getPar(stFit, OPTION = 'list')
matrixcalc::is.positive.definite(stParList$latent_correlations)
#> [1] TRUE

# mean square error
mean((stPar-theta)^2)
#> [1] 0.004670866
```

Numerical estimation as comparison

``` r
numFit <- fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'=1),
  METHOD = 'ucminf')
#> 1. Initialising at default values
#> 2. Computing frequencies...
#> 3. Optimising with ucminf...
#> 4. Done! (31.7 secs)

# extract estimated parameter vector
numPar <- getPar(numFit)

# extract list of parameter estimates
numParList <- getPar(numFit, OPTION = 'list')

# mean square error
mean((numPar-theta)^2)
#> [1] 0.003324634
```

Compute standard errors

``` r
stovar <- computeVar(OBJ = stFit, DATA = D, NUMDERIV = T)
#> 1. Computing H numerically...
#> 2. Estimating J...
#> 3. Inverting H...
#> 3. Computing the variances...
#> Done!
stose <- sqrt(stovar$optimisation_noise+stovar$asymptotic_variance)
stose[1:10]
#>  [1] 0.06241961 0.03966074 0.06040693 0.06736675 0.03962925 0.06303031
#>  [7] 0.06558621 0.03968211 0.06241817 0.06179535
```

``` r

numvar <- computeVar(OBJ = numFit, DATA = D, NUMDERIV = T)
#> 1. Computing H numerically...
#> 2. Estimating J...
#> 3. Inverting H...
#> 3. Computing the variances...
#> Done!
numse <- sqrt(numvar$asymptotic_variance)
numse[1:10]
#>  [1] 0.06236862 0.03962693 0.06035836 0.06734567 0.03960554 0.06298093
#>  [7] 0.06553732 0.03964214 0.06236992 0.06175864
```

``` r
sessionInfo()
#> R version 4.4.2 (2024-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=it_IT.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=it_IT.UTF-8        LC_COLLATE=it_IT.UTF-8    
#>  [5] LC_MONETARY=it_IT.UTF-8    LC_MESSAGES=it_IT.UTF-8   
#>  [7] LC_PAPER=it_IT.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=it_IT.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Rome
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] plFA_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5         cli_3.6.3           knitr_1.48         
#>  [4] rlang_1.1.4         xfun_0.47           generics_0.1.3     
#>  [7] RcppParallel_5.1.8  glue_1.8.0          colorspace_2.1-1   
#> [10] htmltools_0.5.8.1   matrixcalc_1.0-6    RcppClock_1.1      
#> [13] scales_1.3.0        fansi_1.0.6         rmarkdown_2.28     
#> [16] grid_4.4.2          evaluate_1.0.1      munsell_0.5.1      
#> [19] tibble_3.2.1        fastmap_1.2.0       numDeriv_2016.8-1.1
#> [22] mvtnorm_1.3-1       yaml_2.3.10         lifecycle_1.0.4    
#> [25] compiler_4.4.2      dplyr_1.1.4         ucminf_1.2.2       
#> [28] pkgconfig_2.0.3     Rcpp_1.0.13-1       RcppEigen_0.3.4.0.2
#> [31] rstudioapi_0.16.0   digest_0.6.37       R6_2.5.1           
#> [34] tidyselect_1.2.1    utf8_1.2.4          pillar_1.9.0       
#> [37] magrittr_2.0.3      tools_4.4.2         gtable_0.3.6       
#> [40] ggplot2_3.5.1
```
