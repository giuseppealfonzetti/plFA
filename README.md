
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plFA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/plFA)](https://CRAN.R-project.org/package=plFA)
[![R-CMD-check](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check.yaml)
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
#> Loading required package: lavaan
#> This is lavaan 0.6-19
#> lavaan is FREE software! Please report any bugs.
#> 
#> Attaching package: 'plFA'
#> The following object is masked from 'package:lavaan':
#> 
#>     cfa
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
S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q, CORRFLAG = 1)

# simulate the data
D <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)

# categories per item
cat <- apply(D, 2, max) + 1

# store true parameter vector (interpretable parameterisation)
theta <- get_theta(rep(thr, p), Load, S, cat, A, CORRFLAG = 1)
```

<!-- Estimate the model with the stochastic optimiser -->
<!-- ```{r sto} -->
<!-- # Set options for cpp stochastic optimiser -->
<!-- cpp_ctrl <- list( -->
<!--   MAXT = 2*n, -->
<!--   PAIRS_PER_ITERATION = 8 -->
<!-- ) -->
<!-- # Fit the model with stochastic optimiser -->
<!-- stFit <- fit_plFA( -->
<!--   DATA = D, -->
<!--   VALDATA = D, -->
<!--   CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1), -->
<!--   METHOD = 'hyper', -->
<!--   CONTROL = cpp_ctrl, -->
<!--   ITERATIONS_SUBSET = seq(0, cpp_ctrl$MAXT, 100) -->
<!-- ) -->
<!-- stFit -->
<!-- # extract estimated parameter vector -->
<!-- stPar <- getPar(stFit) -->
<!-- # extract list of parameter estimates -->
<!-- stParList <- getPar(stFit, OPTION = 'list') -->
<!-- matrixcalc::is.positive.definite(stParList$latent_correlations) -->
<!-- # mean square error -->
<!-- mean((stPar-theta)^2) -->
<!-- ``` -->

Numerical estimation as comparison

``` r
numFit <- fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'=1),
  METHOD = 'ucminf')
#> 1. Initialising at default values
#> 2. Computing frequencies...
#> 3. Optimising with ucminf...
#> 4. Done! (22.59 secs)

# extract estimated parameter vector
numPar <- getPar(numFit, 'raw')

# extract list of parameter estimates
numParList <- getPar(numFit, OPTION = 'list')

# mean square error
mean((numPar-theta)^2)
#> [1] 0.09474404
```

Compute standard errors

<!-- ```{r stovar} -->
<!-- stovar <- computeVar(OBJ = stFit, DATA = D, NUMDERIV = T) -->
<!-- stose <- sqrt(stovar$optimisation_noise+stovar$asymptotic_variance) -->
<!-- stose[1:10] -->
<!-- ``` -->

``` r

numvar <- computeVar(OBJ = numFit, DATA = D)
#> 2. Estimating J...
#> 3. Computing the variances...
#> Done!
numse <- sqrt(numvar$asymptotic_variance)
numse[1:10]
#>  [1] 0.06290281 0.03945790 0.06169743 0.06636763 0.03894489 0.06398055
#>  [7] 0.06212953 0.03703314 0.06324520 0.06125982
```

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 15.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Asia/Brunei
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] plFA_0.0.0.9001 lavaan_0.6-19  
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5         cli_3.6.3           knitr_1.49         
#>  [4] rlang_1.1.4         xfun_0.49           generics_0.1.3     
#>  [7] RcppParallel_5.1.9  glue_1.8.0          colorspace_2.1-1   
#> [10] pbivnorm_0.6.0      htmltools_0.5.8.1   stats4_4.4.1       
#> [13] RcppClock_1.1       scales_1.3.0        rmarkdown_2.29     
#> [16] quadprog_1.5-8      grid_4.4.1          evaluate_1.0.1     
#> [19] munsell_0.5.1       tibble_3.2.1        fastmap_1.2.0      
#> [22] numDeriv_2016.8-1.1 mvtnorm_1.3-3       yaml_2.3.10        
#> [25] lifecycle_1.0.4     compiler_4.4.1      dplyr_1.1.4        
#> [28] ucminf_1.2.2        pkgconfig_2.0.3     Rcpp_1.0.14        
#> [31] RcppEigen_0.3.4.0.2 rstudioapi_0.17.1   digest_0.6.37      
#> [34] R6_2.5.1            tidyselect_1.2.1    pillar_1.10.1      
#> [37] mnormt_2.1.1        magrittr_2.0.3      tools_4.4.1        
#> [40] gtable_0.3.6        ggplot2_3.5.1
```
