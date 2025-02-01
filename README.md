
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plFA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/plFA)](https://CRAN.R-project.org/package=plFA)
[![R-CMD-check](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giuseppealfonzetti/plFA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `{plFA}` package allows the estimation of confirmatory factor models
for ordinal data using stochastic and numeric pairwise likelihood
optimisation.

## Installation

You can install the development version of plFA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/plFA")
```

When the package is installed, you can load it with:

``` r
library(plFA)
#> ℹ Loading required package: lavaan
#> This is lavaan 0.6-19
#> lavaan is FREE software! Please report any bugs.
#> 
#> ── Conflicts ───────────────────────────────────────────────────── plFA 0.1.0 ──
#> ✖ plFA::cfa() masks lavaan::cfa()
```

The startup message indicates that the `{plFA}` package provides its own
version of the `cfa()` function, which can be used to estimate plFA
models using the `lavaan` syntax. See below for an example.

## Example

Consider the data containing 25 personality items measured on 6-point
scale from the International Personality Item Pool (IPIP) dataset. The
data set is available in the `{plFA}` package and can be loaded with the
following code:

``` r
data("bfi")
head(bfi)
#> # A tibble: 6 × 28
#>   A1    A2    A3    A4    A5    C1    C2    C3    C4    C5    E1    E2    E3   
#>   <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord> <ord>
#> 1 6     6     5     6     5     6     6     6     1     3     2     1     6    
#> 2 4     3     1     5     1     3     2     4     2     4     3     6     4    
#> 3 4     4     5     6     5     4     3     5     3     2     1     3     2    
#> 4 4     5     2     2     1     5     5     5     2     2     3     4     3    
#> 5 1     5     6     5     6     4     3     2     4     5     2     1     2    
#> 6 2     6     5     6     5     3     5     6     3     6     2     2     4    
#> # ℹ 15 more variables: E4 <ord>, E5 <ord>, N1 <ord>, N2 <ord>, N3 <ord>,
#> #   N4 <ord>, N5 <ord>, O1 <ord>, O2 <ord>, O3 <ord>, O4 <ord>, O5 <ord>,
#> #   gender <int>, education <int>, age <int>
```

To fit a latent factor analysis model with 5 factors, where the first 5
items load on the first factor and so on, we can use the `fit_plFA()`
function. In preparation, we must build the following objects to be
supplied as a list to the `CONSTR_LIST` argument:

- `CONSTRMAT`: A matrix of constraints on the factor loadings.
- `CORRFLAG`: A flag indicating whether the factor correlations are
  estimated.
- `CONSTRVAR`: A vector of constraints on the factor variances, if any.
- `STDLV`: A flag indicating whether the factor loadings are
  standardised.
- `LLC`: A list indicating linear constraints on factor loading
  estimations, if any.

The above ingredients are prepared using the code below.

``` r
p <- 25         # number of items
q <- 5          # number of factors
n <- nrow(bfi)  # number of observations

# Convert data frame to matrix, with each item on a 3-point scale 0,1 and 2
# (for simplicity)
D <- as.matrix(data.frame(lapply(bfi[, 1:p], as.numeric)))
D <- apply(head(D / 2), 2, floor)
head(D)
#>      A1 A2 A3 A4 A5 C1 C2 C3 C4 C5 E1 E2 E3 E4 E5 N1 N2 N3 N4 N5 O1 O2 O3 O4 O5
#> [1,]  3  3  2  3  2  3  3  3  0  1  1  0  3  2  3  1  2  1  1  1  2  1  2  3  0
#> [2,]  2  1  0  2  0  1  1  2  1  2  1  3  2  1  0  3  1  1  3  2  1  1  2  2  1
#> [3,]  2  2  2  3  2  2  1  2  1  1  0  1  1  2  2  1  1  2  1  1  2  1  2  3  1
#> [4,]  2  2  1  1  0  2  2  2  1  1  1  2  1  3  2  1  2  1  1  1  2  1  2  2  2
#> [5,]  0  2  3  2  3  2  1  1  2  2  1  0  1  2  1  1  1  1  1  1  3  0  2  2  1
#> [6,]  1  3  2  3  2  1  2  3  1  3  1  1  2  3  3  2  2  2  3  3  3  0  2  3  0

# Constraint matrix: NA indicate free loadings
constrmat <- build_constrMat(P = p, Q = q, STRUCT = "simple")
head(constrmat, 12)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]   NA    0    0    0    0
#>  [2,]   NA    0    0    0    0
#>  [3,]   NA    0    0    0    0
#>  [4,]   NA    0    0    0    0
#>  [5,]   NA    0    0    0    0
#>  [6,]    0   NA    0    0    0
#>  [7,]    0   NA    0    0    0
#>  [8,]    0   NA    0    0    0
#>  [9,]    0   NA    0    0    0
#> [10,]    0   NA    0    0    0
#> [11,]    0    0   NA    0    0
#> [12,]    0    0   NA    0    0

# Build constraint list
constr_list <- list(
  CONSTRMAT = constrmat,
  CORRFLAG = TRUE,
  CONSTRVAR = rep(1, q),  # these options standardise the  
  STDLV = TRUE,           # factor variances
  LLC = NULL
)
```

We can now fit the factor model:

``` r
fit <- fit_plFA(
  DATA = D, 
  CONSTR_LIST = constr_list,
  METHOD = "ucminf",
  VERBOSE = TRUE
)
#> - Computing frequencies...
#> - Initialising with SA ... Done.
#> - Optimising with ucminf...
#> Done! (0.22 secs)
print(fit)
#> - Dimensions:
#>    - Sample size: 6 
#>    - Items: 25 ( 300  pairs)
#>    - Latent traits: 5 
#> 
#>  - Free parameters:
#>    - Thresholds: 66 
#>    - Loadings: 25 
#>    - Latent correlations: 10 
#>    - Latent variances: 0 
#>    - Total: 101 
#> 
#> - Numerical estimate obtained via ucminf.
#> 
#>     Total time: 0.22 s (Data reduction: 0 s)
#>     Cores used: 1 
#>  
#> - Use getPar() to extract parameter estimates.
```

The default method is using the optimiser from the `{ucminf}` package,
which works well for small data sizes. Alternatively, the user may
choose `METHOD = "SA"` to turn on the stochastic approximation
algorithm. Use the function `getPar()` to extract the parameter
estimates:

``` r
getPar(fit)
#> $thresholds
#>  [1] -0.89753150 -0.26343002  1.04187367 -1.15313352 -0.31877232  0.84434997
#>  [7] -0.89892502 -0.26737989  1.04361055 -1.14395827 -0.23543114  0.58338248
#> [13] -0.48004429 -0.49743283  1.01772946 -1.24247107 -0.14067645  1.06252840
#> [19] -1.31732791  0.10768040  1.01856898 -1.16408048 -0.33623210  0.83345818
#> [25] -0.98028776  0.98583212 -1.27252490  0.09565382  1.00597095 -0.41625753
#> [31] -0.74005975  0.18778446  0.88587696 -1.31908569  0.10353940  1.01329123
#> [37] -1.14380838 -0.30368808  0.84786351 -0.88583131 -0.17008387  0.77214151
#> [43] -1.35736909  0.36693824  0.92432458 -1.19496465  0.53912704 -1.23981375
#> [49]  0.77326184 -1.31028262  0.59477780  0.49861553 -1.34420641  0.35872234
#> [55]  0.92809062 -1.15981733 -0.32439685  0.83422903 -0.21621479 -1.00000000
#> [61] -0.25431263 -1.00000000 -0.52545042  0.65278002 -0.73769272  0.94544565
#> 
#> $loadings
#>            [,1]        [,2]        [,3]      [,4]       [,5]
#>  [1,] 0.2690915  0.00000000  0.00000000 0.0000000  0.0000000
#>  [2,] 0.6816494  0.00000000  0.00000000 0.0000000  0.0000000
#>  [3,] 0.3292058  0.00000000  0.00000000 0.0000000  0.0000000
#>  [4,] 0.6446377  0.00000000  0.00000000 0.0000000  0.0000000
#>  [5,] 0.2342125  0.00000000  0.00000000 0.0000000  0.0000000
#>  [6,] 0.0000000  0.37733899  0.00000000 0.0000000  0.0000000
#>  [7,] 0.0000000  0.55843264  0.00000000 0.0000000  0.0000000
#>  [8,] 0.0000000  0.69156240  0.00000000 0.0000000  0.0000000
#>  [9,] 0.0000000 -0.04276024  0.00000000 0.0000000  0.0000000
#> [10,] 0.0000000  0.32387956  0.00000000 0.0000000  0.0000000
#> [11,] 0.0000000  0.00000000  0.40019739 0.0000000  0.0000000
#> [12,] 0.0000000  0.00000000 -0.03668167 0.0000000  0.0000000
#> [13,] 0.0000000  0.00000000  0.45663712 0.0000000  0.0000000
#> [14,] 0.0000000  0.00000000  0.52605244 0.0000000  0.0000000
#> [15,] 0.0000000  0.00000000  0.53747648 0.0000000  0.0000000
#> [16,] 0.0000000  0.00000000  0.00000000 0.2169872  0.0000000
#> [17,] 0.0000000  0.00000000  0.00000000 0.5895011  0.0000000
#> [18,] 0.0000000  0.00000000  0.00000000 0.4632519  0.0000000
#> [19,] 0.0000000  0.00000000  0.00000000 0.3628013  0.0000000
#> [20,] 0.0000000  0.00000000  0.00000000 0.3774107  0.0000000
#> [21,] 0.0000000  0.00000000  0.00000000 0.0000000  0.5396677
#> [22,] 0.0000000  0.00000000  0.00000000 0.0000000  0.1428574
#> [23,] 0.0000000  0.00000000  0.00000000 0.0000000  0.5725192
#> [24,] 0.0000000  0.00000000  0.00000000 0.0000000  0.5862221
#> [25,] 0.0000000  0.00000000  0.00000000 0.0000000 -0.1031503
#> 
#> $latent_correlations
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9253955 0.8979487 0.8705361 0.9214402
#> [2,] 0.9253955 1.0000000 0.9551117 0.9354997 0.9489682
#> [3,] 0.8979487 0.9551117 1.0000000 0.9581621 0.9618791
#> [4,] 0.8705361 0.9354997 0.9581621 1.0000000 0.9560167
#> [5,] 0.9214402 0.9489682 0.9618791 0.9560167 1.0000000
```

Computing the standard errors is done separately using the
`computeVar()` function:

``` r
var <- computeVar(OBJ = fit, DATA = D)
str(var)
#> List of 8
#>  $ trJacob            : num [1:101, 1:101] 1 0 0 0 0 ...
#>  $ H                  : logi NA
#>  $ invH               : num [1:101, 1:101] 1 0 0 0 0 0 0 0 0 0 ...
#>  $ J                  : num [1:101, 1:101] 348 -220 0 0 -125 ...
#>  $ vcov               : num [1:101, 1:101] 348 -220 0 0 -125 ...
#>  $ asymptotic_variance: num [1:101] 348.1 545 292.3 62.5 446.1 ...
#>  $ optimisation_noise : logi NA
#>  $ RTime              : num 0.0771
head(sqrt(diag(var$vcov)))
#> [1] 18.657817 23.344184 17.096949  7.904421 21.120533 21.164977
```

## `{lavaan}` wrapper

This package also provides a user-friendly `{lavaan}`-style interface to
estimate factor analysis models using the `cfa()` function. Let’s rerun
the above example, but this time using the stochastic approximation
algorithm and the original 6-point scale data.

``` r
# {lavaan} model syntax
mod <- "
  Agreeableness     =~ A1 + A2 + A3 + A4 + A5
  Conscientiousness =~ C1 + C2 + C3 + C4 + C5
  Extraversion      =~ E1 + E2 + E3 + E4 + E5
  Neuroticism       =~ N1 + N2 + N3 + N4 + N5
  Openness          =~ O1 + O2 + O3 + O4 + O5
"

# Fit the model (note, this cfa() function is from the plFA package)
fit <- cfa(
  model = mod, 
  data = bfi, 
  std.lv = TRUE,
  estimator.args = list(
    method = "SA",
    ncores = 4
  )
)
#> Warning in fit_plFA(DATA = D, CONSTR_LIST = constr_list, METHOD = method, :
#> Possible divergent trajectories detected. Try decreasing STEP0.
summary(fit)
#> plFA 0.1.0 
#>   ⨉
#> lavaan 0.6-19 ended normally after 4000 iterations
#> 
#>   Estimator                                        PML
#>   Optimization method                               SA
#>   Number of model parameters                       160
#> 
#>   Number of observations                          2236
#> 
#> 
#> Parameter Estimates:
#> 
#>   Parameterization                               Delta
#>   Standard errors                             Sandwich
#>   Information bread                           Observed
#>   Observed information based on                Hessian
#> 
#> Latent Variables:
#>                        Estimate  Std.Err   z-value  P(>|z|)
#>   Agreeableness =~                                         
#>     A1                   -0.290     1.282   -0.227    0.821
#>     A2                    0.548     1.273    0.430    0.667
#>     A3                    0.639     1.392    0.459    0.646
#>     A4                    0.480     1.681    0.285    0.775
#>     A5                    0.676     2.490    0.271    0.786
#>   Conscientiousness =~                                     
#>     C1                    0.539     3.128    0.172    0.863
#>     C2                    0.541     1.965    0.275    0.783
#>     C3                    0.512     1.639    0.312    0.755
#>     C4                   -0.720     1.301   -0.553    0.580
#>     C5                   -0.690     1.286   -0.536    0.592
#>   Extraversion =~                                          
#>     E1                   -0.520     2.404   -0.216    0.829
#>     E2                   -0.730     1.717   -0.425    0.671
#>     E3                    0.654     1.482    0.441    0.659
#>     E4                    0.713     1.268    0.563    0.574
#>     E5                    0.600     1.319    0.455    0.649
#>   Neuroticism =~                                           
#>     N1                    0.419     2.198    0.191    0.849
#>     N2                    0.388     1.636    0.237    0.813
#>     N3                    0.408     1.464    0.279    0.780
#>     N4                    0.543     1.284    0.423    0.672
#>     N5                    0.358     1.252    0.286    0.775
#>   Openness =~                                              
#>     O1                    0.600     2.757    0.218    0.828
#>     O2                   -0.408     1.761   -0.232    0.817
#>     O3                    0.721     1.467    0.491    0.623
#>     O4                    0.082     1.262    0.065    0.948
#>     O5                   -0.418     1.337   -0.313    0.755
#> 
#> Covariances:
#>                        Estimate  Std.Err   z-value  P(>|z|)
#>   Agreeableness ~~                                         
#>     Conscientisnss        0.525     1.902    0.276    0.783
#>     Extraversion          0.753     2.191    0.344    0.731
#>     Neuroticism          -1.000     2.142   -0.467    0.641
#>     Openness              0.382  4280.654    0.000    1.000
#>   Conscientiousness ~~                                     
#>     Extraversion          0.451 78955.461    0.000    1.000
#>     Neuroticism          -0.545 32258.112   -0.000    1.000
#>     Openness              0.395     1.832    0.216    0.829
#>   Extraversion ~~                                          
#>     Neuroticism          -0.755     1.921   -0.393    0.694
#>     Openness              0.575     2.262    0.254    0.799
#>   Neuroticism ~~                                           
#>     Openness             -0.387 33110.041   -0.000    1.000
#> 
#> Thresholds:
#>                    Estimate  Std.Err   z-value  P(>|z|)
#>     A1|t1            -0.404     2.748   -0.147    0.883
#>     A1|t2             0.363     1.840    0.197    0.844
#>     A1|t3             0.780     1.480    0.527    0.598
#>     A1|t4             1.270     1.258    1.009    0.313
#>     A1|t5             1.889     1.366    1.383    0.167
#>     A2|t1            -2.153     2.487   -0.866    0.387
#>     A2|t2            -1.557     1.632   -0.954    0.340
#>     A2|t3            -1.222     1.396   -0.875    0.381
#>     A2|t4            -0.512     1.245   -0.412    0.681
#>     A2|t5             0.455     1.402    0.325    0.745
#>     A3|t1            -1.848     2.490   -0.742    0.458
#>     A3|t2            -1.331     1.627   -0.818    0.413
#>     A3|t3            -0.995     1.394   -0.714    0.475
#>     A3|t4            -0.349     1.248   -0.280    0.780
#>     A3|t5             0.593     1.448    0.410    0.682
#>     A4|t1            -1.716     1.313   -1.307    0.191
#>     A4|t2            -1.206     1.253   -0.962    0.336
#>     A4|t3            -0.927     1.345   -0.689    0.491
#>     A4|t4            -0.423     1.697   -0.249    0.803
#>     A4|t5             0.202     2.826    0.071    0.943
#>     A5|t1            -2.014     1.434   -1.405    0.160
#>     A5|t2            -1.363     1.259   -1.083    0.279
#>     A5|t3            -0.945     1.243   -0.760    0.447
#>     A5|t4            -0.276     1.337   -0.207    0.836
#>     A5|t5             0.652     1.706    0.382    0.702
#>     C1|t1            -2.011     1.372   -1.466    0.143
#>     C1|t2            -1.432     1.243   -1.152    0.249
#>     C1|t3            -0.960     1.261   -0.762    0.446
#>     C1|t4            -0.269     1.391   -0.193    0.847
#>     C1|t5             0.736     1.775    0.415    0.678
#>     C2|t1            -1.878     1.431   -1.312    0.189
#>     C2|t2            -1.202     1.249   -0.963    0.336
#>     C2|t3            -0.784     1.246   -0.629    0.529
#>     C2|t4            -0.146     1.375   -0.106    0.915
#>     C2|t5             0.827     1.716    0.482    0.630
#>     C3|t1            -1.887     2.063   -0.915    0.360
#>     C3|t2            -1.191     1.500   -0.794    0.427
#>     C3|t3            -0.778     1.305   -0.596    0.551
#>     C3|t4            -0.022     1.254   -0.018    0.986
#>     C3|t5             0.921     1.579    0.583    0.560
#>     C4|t1            -0.550     2.049   -0.269    0.788
#>     C4|t2             0.208     1.518    0.137    0.891
#>     C4|t3             0.664     1.362    0.487    0.626
#>     C4|t4             1.289     1.258    1.024    0.306
#>     C4|t5             2.033     1.341    1.516    0.129
#>     C5|t1            -0.886     2.350   -0.377    0.706
#>     C5|t2            -0.254     1.623   -0.156    0.876
#>     C5|t3             0.047     1.399    0.034    0.973
#>     C5|t4             0.630     1.248    0.504    0.614
#>     C5|t5             1.302     1.384    0.941    0.347
#>     E1|t1            -0.728     1.344   -0.542    0.588
#>     E1|t2            -0.063     1.225   -0.051    0.959
#>     E1|t3             0.311     1.244    0.250    0.802
#>     E1|t4             0.789     1.415    0.558    0.577
#>     E1|t5             1.373     1.857    0.740    0.460
#>     E2|t1            -0.867     1.582   -0.548    0.584
#>     E2|t2            -0.153     1.280   -0.119    0.905
#>     E2|t3             0.169     1.225    0.138    0.891
#>     E2|t4             0.754     1.295    0.582    0.560
#>     E2|t5             1.319     1.637    0.806    0.420
#>     E3|t1            -1.625     1.437   -1.130    0.258
#>     E3|t2            -1.010     1.232   -0.820    0.412
#>     E3|t3            -0.518     1.224   -0.423    0.672
#>     E3|t4             0.267     1.334    0.200    0.841
#>     E3|t5             1.147     1.733    0.662    0.508
#>     E4|t1            -1.627     1.466   -1.110    0.267
#>     E4|t2            -1.061     1.248   -0.850    0.396
#>     E4|t3            -0.715     1.234   -0.580    0.562
#>     E4|t4            -0.263     1.369   -0.192    0.848
#>     E4|t5             0.640     1.727    0.371    0.711
#>     E5|t1            -1.813     1.355   -1.338    0.181
#>     E5|t2            -1.198     1.236   -0.970    0.332
#>     E5|t3            -0.798     1.252   -0.637    0.524
#>     E5|t4            -0.168     1.396   -0.120    0.904
#>     E5|t5             0.790     1.748    0.452    0.651
#>     N1|t1            -0.712     4.243   -0.168    0.867
#>     N1|t2            -0.059     2.182   -0.027    0.978
#>     N1|t3             0.337     1.632    0.206    0.837
#>     N1|t4             0.894     1.281    0.698    0.485
#>     N1|t5             1.490     1.278    1.165    0.244
#>     N2|t1            -1.165     1.313   -0.887    0.375
#>     N2|t2            -0.473     1.243   -0.381    0.704
#>     N2|t3            -0.094     1.297   -0.072    0.942
#>     N2|t4             0.572     1.505    0.380    0.704
#>     N2|t5             1.261     1.977    0.638    0.524
#>     N3|t1            -0.924     2.708   -0.341    0.733
#>     N3|t2            -0.210     1.863   -0.113    0.910
#>     N3|t3             0.106     1.478    0.072    0.943
#>     N3|t4             0.686     1.241    0.552    0.581
#>     N3|t5             1.352     1.408    0.960    0.337
#>     N4|t1            -0.960     3.256   -0.295    0.768
#>     N4|t2            -0.235     2.000   -0.117    0.907
#>     N4|t3             0.155     1.671    0.093    0.926
#>     N4|t4             0.768     1.325    0.580    0.562
#>     N4|t5             1.340     1.238    1.082    0.279
#>     N5|t1            -0.714     1.311   -0.544    0.586
#>     N5|t2            -0.035     1.260   -0.027    0.978
#>     N5|t3             0.298     1.391    0.214    0.830
#>     N5|t4             0.825     1.748    0.472    0.637
#>     N5|t5             1.357     2.596    0.523    0.601
#>     O1|t1            -2.447     1.017   -2.405    0.016
#>     O1|t2            -1.695     0.903   -1.878    0.060
#>     O1|t3            -1.190     0.803   -1.482    0.138
#>     O1|t4            -0.429     0.976   -0.440    0.660
#>     O1|t5             0.455     0.783    0.580    0.562
#>     O2|t1            -0.565     0.932   -0.606    0.544
#>     O2|t2             0.131     0.938    0.139    0.889
#>     O2|t3             0.494     0.940    0.526    0.599
#>     O2|t4             1.011     0.803    1.258    0.208
#>     O2|t5             1.557     0.762    2.045    0.041
#>     O3|t1            -1.990     0.879   -2.264    0.024
#>     O3|t2            -1.451     0.642   -2.258    0.024
#>     O3|t3            -0.955     0.706   -1.353    0.176
#>     O3|t4            -0.122     0.675   -0.180    0.857
#>     O3|t5             0.838     0.751    1.117    0.264
#>     O4|t1            -2.188     0.903   -2.423    0.015
#>     O4|t2            -1.567     0.893   -1.755    0.079
#>     O4|t3            -1.246     0.890   -1.400    0.162
#>     O4|t4            -0.607     1.069   -0.568    0.570
#>     O4|t5             0.255     0.941    0.271    0.787
#>     O5|t1            -0.584     1.037   -0.563    0.573
#>     O5|t2             0.261     1.166    0.224    0.823
#>     O5|t3             0.770     0.996    0.774    0.439
#>     O5|t4             1.341     1.401    0.957    0.338
#>     O5|t5             1.932     1.162    1.663    0.096
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .A1                1.000                           
#>    .A2                1.000                           
#>    .A3                1.000                           
#>    .A4                1.000                           
#>    .A5                1.000                           
#>    .C1                1.000                           
#>    .C2                1.000                           
#>    .C3                1.000                           
#>    .C4                1.000                           
#>    .C5                1.000                           
#>    .E1                1.000                           
#>    .E2                1.000                           
#>    .E3                1.000                           
#>    .E4                1.000                           
#>    .E5                1.000                           
#>    .N1                1.000                           
#>    .N2                1.000                           
#>    .N3                1.000                           
#>    .N4                1.000                           
#>    .N5                1.000                           
#>    .O1                1.000                           
#>    .O2                1.000                           
#>    .O3                1.000                           
#>    .O4                1.000                           
#>    .O5                1.000                           
#>     Agreeableness     1.000                           
#>     Conscientisnss    1.000                           
#>     Extraversion      1.000                           
#>     Neuroticism       1.000                           
#>     Openness          1.000
```

All[^1] `lavaan` methods should carry over, such as:

``` r
coef(fit)
#>                Agreeableness=~A1                Agreeableness=~A2 
#>                           -0.290                            0.548 
#>                Agreeableness=~A3                Agreeableness=~A4 
#>                            0.639                            0.480 
#>                Agreeableness=~A5            Conscientiousness=~C1 
#>                            0.676                            0.539 
#>            Conscientiousness=~C2            Conscientiousness=~C3 
#>                            0.541                            0.512 
#>            Conscientiousness=~C4            Conscientiousness=~C5 
#>                           -0.720                           -0.690 
#>                 Extraversion=~E1                 Extraversion=~E2 
#>                           -0.520                           -0.730 
#>                 Extraversion=~E3                 Extraversion=~E4 
#>                            0.654                            0.713 
#>                 Extraversion=~E5                  Neuroticism=~N1 
#>                            0.600                            0.419 
#>                  Neuroticism=~N2                  Neuroticism=~N3 
#>                            0.388                            0.408 
#>                  Neuroticism=~N4                  Neuroticism=~N5 
#>                            0.543                            0.358 
#>                     Openness=~O1                     Openness=~O2 
#>                            0.600                           -0.408 
#>                     Openness=~O3                     Openness=~O4 
#>                            0.721                            0.082 
#>                     Openness=~O5                            A1|t1 
#>                           -0.418                           -0.404 
#>                            A1|t2                            A1|t3 
#>                            0.363                            0.780 
#>                            A1|t4                            A1|t5 
#>                            1.270                            1.889 
#>                            A2|t1                            A2|t2 
#>                           -2.153                           -1.557 
#>                            A2|t3                            A2|t4 
#>                           -1.222                           -0.512 
#>                            A2|t5                            A3|t1 
#>                            0.455                           -1.848 
#>                            A3|t2                            A3|t3 
#>                           -1.331                           -0.995 
#>                            A3|t4                            A3|t5 
#>                           -0.349                            0.593 
#>                            A4|t1                            A4|t2 
#>                           -1.716                           -1.206 
#>                            A4|t3                            A4|t4 
#>                           -0.927                           -0.423 
#>                            A4|t5                            A5|t1 
#>                            0.202                           -2.014 
#>                            A5|t2                            A5|t3 
#>                           -1.363                           -0.945 
#>                            A5|t4                            A5|t5 
#>                           -0.276                            0.652 
#>                            C1|t1                            C1|t2 
#>                           -2.011                           -1.432 
#>                            C1|t3                            C1|t4 
#>                           -0.960                           -0.269 
#>                            C1|t5                            C2|t1 
#>                            0.736                           -1.878 
#>                            C2|t2                            C2|t3 
#>                           -1.202                           -0.784 
#>                            C2|t4                            C2|t5 
#>                           -0.146                            0.827 
#>                            C3|t1                            C3|t2 
#>                           -1.887                           -1.191 
#>                            C3|t3                            C3|t4 
#>                           -0.778                           -0.022 
#>                            C3|t5                            C4|t1 
#>                            0.921                           -0.550 
#>                            C4|t2                            C4|t3 
#>                            0.208                            0.664 
#>                            C4|t4                            C4|t5 
#>                            1.289                            2.033 
#>                            C5|t1                            C5|t2 
#>                           -0.886                           -0.254 
#>                            C5|t3                            C5|t4 
#>                            0.047                            0.630 
#>                            C5|t5                            E1|t1 
#>                            1.302                           -0.728 
#>                            E1|t2                            E1|t3 
#>                           -0.063                            0.311 
#>                            E1|t4                            E1|t5 
#>                            0.789                            1.373 
#>                            E2|t1                            E2|t2 
#>                           -0.867                           -0.153 
#>                            E2|t3                            E2|t4 
#>                            0.169                            0.754 
#>                            E2|t5                            E3|t1 
#>                            1.319                           -1.625 
#>                            E3|t2                            E3|t3 
#>                           -1.010                           -0.518 
#>                            E3|t4                            E3|t5 
#>                            0.267                            1.147 
#>                            E4|t1                            E4|t2 
#>                           -1.627                           -1.061 
#>                            E4|t3                            E4|t4 
#>                           -0.715                           -0.263 
#>                            E4|t5                            E5|t1 
#>                            0.640                           -1.813 
#>                            E5|t2                            E5|t3 
#>                           -1.198                           -0.798 
#>                            E5|t4                            E5|t5 
#>                           -0.168                            0.790 
#>                            N1|t1                            N1|t2 
#>                           -0.712                           -0.059 
#>                            N1|t3                            N1|t4 
#>                            0.337                            0.894 
#>                            N1|t5                            N2|t1 
#>                            1.490                           -1.165 
#>                            N2|t2                            N2|t3 
#>                           -0.473                           -0.094 
#>                            N2|t4                            N2|t5 
#>                            0.572                            1.261 
#>                            N3|t1                            N3|t2 
#>                           -0.924                           -0.210 
#>                            N3|t3                            N3|t4 
#>                            0.106                            0.686 
#>                            N3|t5                            N4|t1 
#>                            1.352                           -0.960 
#>                            N4|t2                            N4|t3 
#>                           -0.235                            0.155 
#>                            N4|t4                            N4|t5 
#>                            0.768                            1.340 
#>                            N5|t1                            N5|t2 
#>                           -0.714                           -0.035 
#>                            N5|t3                            N5|t4 
#>                            0.298                            0.825 
#>                            N5|t5                            O1|t1 
#>                            1.357                           -2.447 
#>                            O1|t2                            O1|t3 
#>                           -1.695                           -1.190 
#>                            O1|t4                            O1|t5 
#>                           -0.429                            0.455 
#>                            O2|t1                            O2|t2 
#>                           -0.565                            0.131 
#>                            O2|t3                            O2|t4 
#>                            0.494                            1.011 
#>                            O2|t5                            O3|t1 
#>                            1.557                           -1.990 
#>                            O3|t2                            O3|t3 
#>                           -1.451                           -0.955 
#>                            O3|t4                            O3|t5 
#>                           -0.122                            0.838 
#>                            O4|t1                            O4|t2 
#>                           -2.188                           -1.567 
#>                            O4|t3                            O4|t4 
#>                           -1.246                           -0.607 
#>                            O4|t5                            O5|t1 
#>                            0.255                           -0.584 
#>                            O5|t2                            O5|t3 
#>                            0.261                            0.770 
#>                            O5|t4                            O5|t5 
#>                            1.341                            1.932 
#> Agreeableness~~Conscientiousness      Agreeableness~~Extraversion 
#>                            0.525                            0.753 
#>       Agreeableness~~Neuroticism          Agreeableness~~Openness 
#>                           -1.000                            0.382 
#>  Conscientiousness~~Extraversion   Conscientiousness~~Neuroticism 
#>                            0.451                           -0.545 
#>      Conscientiousness~~Openness        Extraversion~~Neuroticism 
#>                            0.395                           -0.755 
#>           Extraversion~~Openness            Neuroticism~~Openness 
#>                            0.575                           -0.387
head(partable(fit))
#>   id               lhs op rhs user block group free ustart exo label plabel
#> 1  1     Agreeableness =~  A1    1     1     1    1     NA   0         .p1.
#> 2  2     Agreeableness =~  A2    1     1     1    2     NA   0         .p2.
#> 3  3     Agreeableness =~  A3    1     1     1    3     NA   0         .p3.
#> 4  4     Agreeableness =~  A4    1     1     1    4     NA   0         .p4.
#> 5  5     Agreeableness =~  A5    1     1     1    5     NA   0         .p5.
#> 6  6 Conscientiousness =~  C1    1     1     1    6     NA   0         .p6.
#>    start    est    se
#> 1  0.483 -0.290 1.282
#> 2 -0.855  0.548 1.273
#> 3 -0.702  0.639 1.392
#> 4 -0.514  0.480 1.681
#> 5 -0.595  0.676 2.490
#> 6  0.638  0.539 3.128
```

## Session information

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
#>  date     2025-02-01
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
#>  numDeriv       2016.8-1.1 2019-06-06 [1] CRAN (R 4.4.0)
#>  pbivnorm       0.6.0      2015-01-23 [1] CRAN (R 4.4.0)
#>  pillar         1.10.1     2025-01-07 [1] CRAN (R 4.4.1)
#>  pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.4.0)
#>  plFA         * 0.1.0      2025-02-01 [1] local
#>  quadprog       1.5-8      2019-11-20 [1] CRAN (R 4.4.0)
#>  R6             2.5.1      2021-08-19 [1] CRAN (R 4.4.0)
#>  Rcpp           1.0.14     2025-01-12 [1] CRAN (R 4.4.1)
#>  RcppClock      1.1        2021-11-06 [1] CRAN (R 4.4.0)
#>  RcppEigen      0.3.4.0.2  2024-08-24 [1] CRAN (R 4.4.1)
#>  RcppParallel   5.1.10     2025-01-24 [1] CRAN (R 4.4.1)
#>  rlang          1.1.5      2025-01-17 [1] CRAN (R 4.4.1)
#>  rmarkdown      2.29       2024-11-04 [1] CRAN (R 4.4.1)
#>  rstudioapi     0.17.1     2024-10-22 [1] CRAN (R 4.4.1)
#>  scales         1.3.0      2023-11-28 [1] CRAN (R 4.4.0)
#>  sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.4.0)
#>  stringi        1.8.4      2024-05-06 [1] CRAN (R 4.4.0)
#>  stringr        1.5.1      2023-11-14 [1] CRAN (R 4.4.0)
#>  tibble         3.2.1      2023-03-20 [1] CRAN (R 4.4.0)
#>  tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.4.0)
#>  ucminf         1.2.2      2025-01-23 [1] Github (hdakpo/ucminf@a3a411f)
#>  utf8           1.2.4      2023-10-22 [1] CRAN (R 4.4.0)
#>  vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.4.0)
#>  xfun           0.49       2024-10-31 [1] CRAN (R 4.4.1)
#>  yaml           2.3.10     2024-07-26 [1] CRAN (R 4.4.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```

[^1]: Most :)
