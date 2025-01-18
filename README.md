
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plFA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/plFA)](https://CRAN.R-project.org/package=plFA)
[![R-CMD-check](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check/badge.svg)](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check.yaml)
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
#> 4. Done! (23.88 secs)

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

## `{lavaan}` wrapper

The package provides a user-friendly `{lavaan}`-style interface to
estimate plFA models via the `cfa()` function.

``` r
dat <- as.data.frame(lapply(as.data.frame(D), ordered))
names(dat) <- paste0("y", 1:p)

# Model syntax
mod <- "
  eta1 =~ y1 + y2 + y3
  eta2 =~ y4 + y5 + y6
  eta3 =~ y7 + y8 + y9
"

# Fit the model
fit <- cfa(model = mod, data = dat, std.lv = TRUE)
#> 1. Initialising at default values
#> 2. Computing frequencies...
#> 3. Optimising with ucminf...
#> 4. Done! (0.07 secs)
#> 2. Estimating J...
#> 3. Computing the variances...
#> Done!
summary(fit)
#> plFA x
#> lavaan 0.0.0.9001 ended normally after 49 iterations
#> 
#>   Estimator                                        PML
#>   Optimization method                           UCMINF
#>   Number of model parameters                        39
#> 
#>   Number of observations                          1000
#> 
#> 
#> Parameter Estimates:
#> 
#>   Parameterization                               Delta
#>   Standard errors                             Sandwich
#>   Information bread                           Expected
#>   Information bread saturated (h1) model  Unstructured
#>   Information meat saturated (h1) model     Structured
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   eta1 =~                                             
#>     y1                0.291    0.119    2.440    0.015
#>     y2                0.782    0.050   15.525    0.000
#>     y3                0.405    0.120    3.373    0.001
#>   eta2 =~                                             
#>     y4                0.873    0.103    8.487    0.000
#>     y5                0.921    0.047   19.386    0.000
#>     y6                0.104    0.093    1.125    0.260
#>   eta3 =~                                             
#>     y7                0.569    0.091    6.263    0.000
#>     y8                0.874    0.046   18.843    0.000
#>     y9                0.623    0.086    7.214    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   eta1 ~~                                             
#>     eta2              0.931    0.044   21.017    0.000
#>     eta3              0.357    0.060    5.996    0.000
#>   eta2 ~~                                             
#>     eta3              0.073    0.052    1.391    0.164
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     y1|t1            -1.538    0.095  -16.272    0.000
#>     y1|t2             0.015    0.039    0.375    0.707
#>     y1|t3             1.484    0.090   16.449    0.000
#>     y2|t1            -1.661    0.076  -21.728    0.000
#>     y2|t2             0.029    0.044    0.657    0.511
#>     y2|t3             1.557    0.072   21.734    0.000
#>     y3|t1            -1.617    0.106  -15.309    0.000
#>     y3|t2             0.058    0.049    1.184    0.237
#>     y3|t3             1.538    0.114   13.521    0.000
#>     y4|t1            -1.528    0.080  -19.208    0.000
#>     y4|t2             0.024    0.050    0.479    0.632
#>     y4|t3             1.583    0.098   16.188    0.000
#>     y5|t1            -1.460    0.084  -17.437    0.000
#>     y5|t2             0.049    0.048    1.025    0.305
#>     y5|t3             1.492    0.086   17.436    0.000
#>     y6|t1            -1.506    0.083  -18.157    0.000
#>     y6|t2             0.028    0.048    0.575    0.565
#>     y6|t3             1.530    0.083   18.394    0.000
#>     y7|t1            -1.411    0.059  -23.826    0.000
#>     y7|t2             0.034    0.040    0.858    0.391
#>     y7|t3             1.484    0.051   29.192    0.000
#>     y8|t1            -1.424    0.021  -68.134    0.000
#>     y8|t2             0.028    0.022    1.315    0.189
#>     y8|t3             1.463    0.062   23.589    0.000
#>     y9|t1            -1.514    0.036  -41.696    0.000
#>     y9|t2             0.047    0.032    1.476    0.140
#>     y9|t3             1.454    0.035   41.572    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .y1                1.000                           
#>    .y2                1.000                           
#>    .y3                1.000                           
#>    .y4                1.000                           
#>    .y5                1.000                           
#>    .y6                1.000                           
#>    .y7                1.000                           
#>    .y8                1.000                           
#>    .y9                1.000                           
#>     eta1              1.000                           
#>     eta2              1.000                           
#>     eta3              1.000
```

## Session Info

``` r
sessioninfo::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.4.1 (2024-06-14)
#>  os       macOS 15.2
#>  system   aarch64, darwin20
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Asia/Brunei
#>  date     2025-01-18
#>  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package      * version    date (UTC) lib source
#>  cli            3.6.3      2024-06-21 [1] CRAN (R 4.4.0)
#>  colorspace     2.1-1      2024-07-26 [1] CRAN (R 4.4.0)
#>  digest         0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
#>  dplyr          1.1.4      2023-11-17 [1] CRAN (R 4.4.0)
#>  evaluate       1.0.1      2024-10-10 [1] CRAN (R 4.4.1)
#>  fastmap        1.2.0      2024-05-15 [1] CRAN (R 4.4.0)
#>  generics       0.1.3      2022-07-05 [1] CRAN (R 4.4.0)
#>  ggplot2        3.5.1      2024-04-23 [1] CRAN (R 4.4.0)
#>  glue           1.8.0      2024-09-30 [1] CRAN (R 4.4.1)
#>  gtable         0.3.6      2024-10-25 [1] CRAN (R 4.4.1)
#>  htmltools      0.5.8.1    2024-04-04 [1] CRAN (R 4.4.0)
#>  knitr          1.49       2024-11-08 [1] CRAN (R 4.4.1)
#>  lavaan       * 0.6-19     2024-09-26 [1] CRAN (R 4.4.1)
#>  lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.4.0)
#>  magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.4.0)
#>  mnormt         2.1.1      2022-09-26 [1] CRAN (R 4.4.0)
#>  munsell        0.5.1      2024-04-01 [1] CRAN (R 4.4.0)
#>  mvtnorm        1.3-3      2025-01-10 [1] CRAN (R 4.4.1)
#>  numDeriv       2016.8-1.1 2019-06-06 [1] CRAN (R 4.4.0)
#>  pbivnorm       0.6.0      2015-01-23 [1] CRAN (R 4.4.0)
#>  pillar         1.10.1     2025-01-07 [1] CRAN (R 4.4.1)
#>  pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.4.0)
#>  plFA         * 0.0.0.9001 2025-01-18 [1] local
#>  quadprog       1.5-8      2019-11-20 [1] CRAN (R 4.4.0)
#>  R6             2.5.1      2021-08-19 [1] CRAN (R 4.4.0)
#>  Rcpp           1.0.14     2025-01-12 [1] CRAN (R 4.4.1)
#>  RcppClock      1.1        2021-11-06 [1] CRAN (R 4.4.0)
#>  RcppEigen      0.3.4.0.2  2024-08-24 [1] CRAN (R 4.4.1)
#>  RcppParallel   5.1.9      2024-08-19 [1] CRAN (R 4.4.1)
#>  rlang          1.1.4      2024-06-04 [1] CRAN (R 4.4.0)
#>  rmarkdown      2.29       2024-11-04 [1] CRAN (R 4.4.1)
#>  rstudioapi     0.17.1     2024-10-22 [1] CRAN (R 4.4.1)
#>  scales         1.3.0      2023-11-28 [1] CRAN (R 4.4.0)
#>  sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.4.0)
#>  tibble         3.2.1      2023-03-20 [1] CRAN (R 4.4.0)
#>  tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.4.0)
#>  ucminf         1.2.2      2024-06-24 [1] CRAN (R 4.4.0)
#>  vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.4.0)
#>  xfun           0.49       2024-10-31 [1] CRAN (R 4.4.1)
#>  yaml           2.3.10     2024-07-26 [1] CRAN (R 4.4.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
