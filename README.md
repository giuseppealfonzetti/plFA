
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

# Convert data frame to matrix, with each item taking values 0,1,2,...,5
D <- as.matrix(data.frame(lapply(bfi[, 1:p], as.numeric))) - 1
head(D)
#>      A1 A2 A3 A4 A5 C1 C2 C3 C4 C5 E1 E2 E3 E4 E5 N1 N2 N3 N4 N5 O1 O2 O3 O4 O5
#> [1,]  5  5  4  5  4  5  5  5  0  2  1  0  5  4  5  2  4  1  1  2  3  2  4  5  0
#> [2,]  3  2  0  4  0  2  1  3  1  3  2  5  3  1  0  5  2  1  5  3  2  1  3  4  2
#> [3,]  3  3  4  5  4  3  2  4  2  1  0  2  1  4  3  2  2  3  1  2  4  2  4  5  2
#> [4,]  3  4  1  1  0  4  4  4  1  1  2  3  2  5  4  1  3  1  1  2  4  1  4  4  4
#> [5,]  0  4  5  4  5  3  2  1  3  4  1  0  1  4  1  1  1  1  1  1  5  0  4  4  1
#> [6,]  1  5  4  5  4  2  4  5  2  5  1  1  3  5  5  3  3  3  5  5  5  0  4  5  0

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
#> Done! (4.69 secs)
print(fit)
#> - Dimensions:
#>    - Sample size: 2236 
#>    - Items: 25 ( 300  pairs)
#>    - Latent traits: 5 
#> 
#>  - Free parameters:
#>    - Thresholds: 125 
#>    - Loadings: 25 
#>    - Latent correlations: 10 
#>    - Latent variances: 0 
#>    - Total: 160 
#> 
#> - Numerical estimate obtained via ucminf.
#> 
#>     Total time: 4.69 s (Data reduction: 0.05 s)
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
#>   [1] -0.40342975  0.36296158  0.77979396  1.27012725  1.88847654 -2.15645578
#>   [7] -1.55843441 -1.22269675 -0.51268052  0.45551687 -1.84890939 -1.33210658
#>  [13] -0.99519211 -0.34954618  0.59330668 -1.71624967 -1.20615048 -0.92738336
#>  [19] -0.42301003  0.20174764 -2.01475684 -1.36362928 -0.94614150 -0.27696867
#>  [25]  0.65201808 -2.01132631 -1.43219461 -0.96036696 -0.26857279  0.73619150
#>  [31] -1.87965681 -1.20296001 -0.78457218 -0.14617313  0.82755714 -1.88675414
#>  [37] -1.19146229 -0.77840530 -0.02248862  0.92135042 -0.55048487  0.20834125
#>  [43]  0.66399977  1.28873806  2.03317057 -0.88613864 -0.25386319  0.04740555
#>  [49]  0.63010147  1.30278349 -0.72803433 -0.06273648  0.31148908  0.78947419
#>  [55]  1.37340971 -0.86693287 -0.15288018  0.16864075  0.75418319  1.31961507
#>  [61] -1.62489332 -1.01028743 -0.51842151  0.26713987  1.14770838 -1.62791281
#>  [67] -1.06094033 -0.71552961 -0.26268359  0.64044010 -1.81316372 -1.19864360
#>  [73] -0.79811409 -0.16798431  0.78993071 -0.71063704 -0.06018644  0.33492797
#>  [79]  0.89235623  1.48922166 -1.16284396 -0.47258028 -0.09397253  0.57016151
#>  [85]  1.25928746 -0.92354249 -0.21037911  0.10580162  0.68502246  1.35122727
#>  [91] -0.96082085 -0.23512431  0.15468532  0.76811843  1.34146250 -0.71406194
#>  [97] -0.03538097  0.29728057  0.82453916  1.35765031 -2.44741937 -1.69521428
#> [103] -1.19002335 -0.42977506  0.45472814 -0.56518561  0.13085747  0.49458507
#> [109]  1.01070361  1.55717203 -1.98983391 -1.45112391 -0.95609874 -0.12195569
#> [115]  0.83877632 -2.18710334 -1.56666348 -1.24590767 -0.60730957  0.25489866
#> [121] -0.58430883  0.26112088  0.77087279  1.34085101  1.93118790
#> 
#> $loadings
#>             [,1]       [,2]       [,3]      [,4]       [,5]
#>  [1,] -0.3353209  0.0000000  0.0000000 0.0000000  0.0000000
#>  [2,]  0.6627841  0.0000000  0.0000000 0.0000000  0.0000000
#>  [3,]  0.7705674  0.0000000  0.0000000 0.0000000  0.0000000
#>  [4,]  0.5450794  0.0000000  0.0000000 0.0000000  0.0000000
#>  [5,]  0.7737092  0.0000000  0.0000000 0.0000000  0.0000000
#>  [6,]  0.0000000  0.5644095  0.0000000 0.0000000  0.0000000
#>  [7,]  0.0000000  0.5707375  0.0000000 0.0000000  0.0000000
#>  [8,]  0.0000000  0.5430038  0.0000000 0.0000000  0.0000000
#>  [9,]  0.0000000 -0.7602739  0.0000000 0.0000000  0.0000000
#> [10,]  0.0000000 -0.7055654  0.0000000 0.0000000  0.0000000
#> [11,]  0.0000000  0.0000000 -0.5401615 0.0000000  0.0000000
#> [12,]  0.0000000  0.0000000 -0.7270841 0.0000000  0.0000000
#> [13,]  0.0000000  0.0000000  0.6902125 0.0000000  0.0000000
#> [14,]  0.0000000  0.0000000  0.7419254 0.0000000  0.0000000
#> [15,]  0.0000000  0.0000000  0.6202669 0.0000000  0.0000000
#> [16,]  0.0000000  0.0000000  0.0000000 0.8365348  0.0000000
#> [17,]  0.0000000  0.0000000  0.0000000 0.8007918  0.0000000
#> [18,]  0.0000000  0.0000000  0.0000000 0.7554746  0.0000000
#> [19,]  0.0000000  0.0000000  0.0000000 0.6664113  0.0000000
#> [20,]  0.0000000  0.0000000  0.0000000 0.5696504  0.0000000
#> [21,]  0.0000000  0.0000000  0.0000000 0.0000000  0.6474240
#> [22,]  0.0000000  0.0000000  0.0000000 0.0000000 -0.4443222
#> [23,]  0.0000000  0.0000000  0.0000000 0.0000000  0.8077955
#> [24,]  0.0000000  0.0000000  0.0000000 0.0000000  0.1355670
#> [25,]  0.0000000  0.0000000  0.0000000 0.0000000 -0.4862491
#> 
#> $latent_correlations
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  1.0000000  0.3697391  0.7022438 -0.2317069  0.2693949
#> [2,]  0.3697391  1.0000000  0.3901962 -0.3012799  0.3221120
#> [3,]  0.7022438  0.3901962  1.0000000 -0.2674171  0.4659849
#> [4,] -0.2317069 -0.3012799 -0.2674171  1.0000000 -0.1463728
#> [5,]  0.2693949  0.3221120  0.4659849 -0.1463728  1.0000000
```

Computing the standard errors is done separately using the
`computeVar()` function:

``` r
var <- computeVar(OBJ = fit, DATA = D, VERBOSE = TRUE)
#> - Estimating J...
#> - Computing the variances...
#> Done! (6.7 secs)
str(var)
#> List of 8
#>  $ trJacob            : num [1:160, 1:160] 1 0 0 0 0 ...
#>  $ H                  : logi NA
#>  $ invH               : num [1:160, 1:160] 0.1043 0.0793 0.0736 0.0667 0.0582 ...
#>  $ J                  : num [1:160, 1:160] 484 -265 0 0 0 ...
#>  $ vcov               : num [1:160, 1:160] 7.02 6.15 5.33 3.74 4.08 ...
#>  $ asymptotic_variance: num [1:160] 7.02 8.48 10.16 16.8 43.41 ...
#>  $ optimisation_noise : logi NA
#>  $ RTime              : num 6.7

# Standard errors
head(sqrt(diag(var$vcov) / n))
#> [1] 0.05603846 0.06158027 0.06739821 0.08669167 0.13933213 0.13607636
```

## `{lavaan}` wrapper

This package also provides a user-friendly `{lavaan}`-style interface to
estimate factor analysis models using the `cfa()` function. Let’s rerun
the above example, but this time using the stochastic approximation
algorithm.

``` r
# {lavaan} model syntax, not in the same order as the data
mod <- "
  opn =~ O1 + O2 + O3 + O4 + O5  # Openness
  con =~ C1 + C2 + C3 + C4 + C5  # Conscientiousness
  ext =~ E1 + E2 + E3 + E4 + E5  # Extraversion
  agr =~ A1 + A2 + A3 + A4 + A5  # Agreeableness
  neu =~ N1 + N2 + N3 + N4 + N5  # Neuroticism
"

# Fit the model (note, this cfa() function is from the plFA package)
fit <- cfa(
  model = mod, 
  data = bfi, 
  std.lv = TRUE,
  estimator.args = list(method = "SA", ncores = 4)
)
#> Warning in fit_plFA(DATA = D, CONSTR_LIST = constr_list, METHOD = method, :
#> Possible divergent trajectories detected. Try decreasing STEP0.

fit@timing$total  # total time in seconds
#> [1] 9.26392

summary(fit)
#> plFA 0.1.0 
#>   ⨉
#> lavaan 0.6-19 ended normally after 3500 iterations
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
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn =~                                              
#>     O1                0.359    0.090    3.999    0.000
#>     O2               -0.267    0.046   -5.809    0.000
#>     O3                0.422    0.034   12.271    0.000
#>     O4               -0.075    0.027   -2.776    0.005
#>     O5               -0.255    0.027   -9.508    0.000
#>   con =~                                              
#>     C1                0.515    0.028   18.642    0.000
#>     C2                0.528    0.026   20.138    0.000
#>     C3                0.497    0.027   18.146    0.000
#>     C4               -0.719    0.032  -22.607    0.000
#>     C5               -0.682    0.042  -16.297    0.000
#>   ext =~                                              
#>     E1               -0.517    0.057   -9.035    0.000
#>     E2               -0.705    0.039  -18.020    0.000
#>     E3                0.637    0.031   20.508    0.000
#>     E4                0.719    0.026   27.742    0.000
#>     E5                0.575    0.029   19.540    0.000
#>   agr =~                                              
#>     A1               -0.304    0.069   -4.378    0.000
#>     A2                0.611    0.042   14.413    0.000
#>     A3                0.714    0.035   20.175    0.000
#>     A4                0.529    0.028   18.887    0.000
#>     A5                0.741    0.026   28.276    0.000
#>   neu =~                                              
#>     N1                0.722    0.028   26.227    0.000
#>     N2                0.669    0.026   25.232    0.000
#>     N3                0.672    0.029   22.886    0.000
#>     N4                0.669    0.037   18.097    0.000
#>     N5                0.546    0.055    9.936    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn ~~                                              
#>     con               0.599    0.082    7.256    0.000
#>     ext               0.784    0.133    5.893    0.000
#>     agr               0.569    0.127    4.476    0.000
#>     neu              -0.882    0.082  -10.799    0.000
#>   con ~~                                              
#>     ext               0.445    0.071    6.290    0.000
#>     agr               0.453    0.121    3.751    0.000
#>     neu              -0.571    0.286   -1.992    0.046
#>   ext ~~                                              
#>     agr               0.833    0.167    4.984    0.000
#>     neu              -0.527    0.659   -0.800    0.424
#>   agr ~~                                              
#>     neu              -0.449    0.891   -0.503    0.615
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     O1|t1            -2.449    0.058  -42.201    0.000
#>     O1|t2            -1.695    0.039  -43.647    0.000
#>     O1|t3            -1.190    0.031  -38.100    0.000
#>     O1|t4            -0.429    0.027  -16.168    0.000
#>     O1|t5             0.455    0.029   15.777    0.000
#>     O2|t1            -0.565    0.052  -10.783    0.000
#>     O2|t2             0.130    0.034    3.778    0.000
#>     O2|t3             0.494    0.029   16.782    0.000
#>     O2|t4             1.011    0.026   38.542    0.000
#>     O2|t5             1.558    0.030   52.813    0.000
#>     O3|t1            -1.992    0.053  -37.807    0.000
#>     O3|t2            -1.451    0.034  -42.216    0.000
#>     O3|t3            -0.955    0.029  -32.431    0.000
#>     O3|t4            -0.121    0.026   -4.595    0.000
#>     O3|t5             0.838    0.031   27.420    0.000
#>     O4|t1            -2.189    0.028  -78.922    0.000
#>     O4|t2            -1.566    0.027  -59.101    0.000
#>     O4|t3            -1.245    0.028  -43.738    0.000
#>     O4|t4            -0.607    0.036  -16.896    0.000
#>     O4|t5             0.255    0.060    4.265    0.000
#>     O5|t1            -0.584    0.030  -19.262    0.000
#>     O5|t2             0.260    0.027    9.751    0.000
#>     O5|t3             0.770    0.026   29.231    0.000
#>     O5|t4             1.341    0.028   47.339    0.000
#>     O5|t5             1.934    0.036   53.506    0.000
#>     C1|t1            -2.012    0.029  -69.405    0.000
#>     C1|t2            -1.432    0.026  -54.520    0.000
#>     C1|t3            -0.960    0.027  -36.058    0.000
#>     C1|t4            -0.268    0.029   -9.132    0.000
#>     C1|t5             0.736    0.037   19.632    0.000
#>     C2|t1            -1.880    0.030  -62.168    0.000
#>     C2|t2            -1.203    0.026  -45.558    0.000
#>     C2|t3            -0.784    0.026  -29.770    0.000
#>     C2|t4            -0.146    0.029   -5.019    0.000
#>     C2|t5             0.827    0.036   22.805    0.000
#>     C3|t1            -1.887    0.044  -43.306    0.000
#>     C3|t2            -1.191    0.032  -37.571    0.000
#>     C3|t3            -0.778    0.028  -28.237    0.000
#>     C3|t4            -0.022    0.026   -0.845    0.398
#>     C3|t5             0.921    0.033   27.623    0.000
#>     C4|t1            -0.550    0.043  -12.690    0.000
#>     C4|t2             0.208    0.032    6.463    0.000
#>     C4|t3             0.663    0.029   22.995    0.000
#>     C4|t4             1.289    0.027   48.408    0.000
#>     C4|t5             2.034    0.028   71.768    0.000
#>     C5|t1            -0.886    0.050  -17.844    0.000
#>     C5|t2            -0.254    0.034   -7.405    0.000
#>     C5|t3             0.047    0.030    1.595    0.111
#>     C5|t4             0.630    0.026   23.894    0.000
#>     C5|t5             1.302    0.029   44.621    0.000
#>     E1|t1            -0.728    0.027  -26.826    0.000
#>     E1|t2            -0.063    0.027   -2.333    0.020
#>     E1|t3             0.311    0.029   10.577    0.000
#>     E1|t4             0.789    0.036   22.195    0.000
#>     E1|t5             1.373    0.053   26.085    0.000
#>     E2|t1            -0.867    0.066  -13.114    0.000
#>     E2|t2            -0.153    0.042   -3.678    0.000
#>     E2|t3             0.168    0.035    4.843    0.000
#>     E2|t4             0.754    0.028   27.330    0.000
#>     E2|t5             1.320    0.027   48.419    0.000
#>     E3|t1            -1.625    0.051  -31.954    0.000
#>     E3|t2            -1.010    0.036  -27.737    0.000
#>     E3|t3            -0.518    0.031  -16.475    0.000
#>     E3|t4             0.267    0.027    9.934    0.000
#>     E3|t5             1.147    0.028   41.000    0.000
#>     E4|t1            -1.628    0.046  -35.038    0.000
#>     E4|t2            -1.061    0.035  -30.655    0.000
#>     E4|t3            -0.715    0.031  -23.091    0.000
#>     E4|t4            -0.262    0.027   -9.659    0.000
#>     E4|t5             0.640    0.026   24.175    0.000
#>     E5|t1            -1.813    0.058  -31.145    0.000
#>     E5|t2            -1.198    0.037  -32.168    0.000
#>     E5|t3            -0.798    0.031  -25.666    0.000
#>     E5|t4            -0.168    0.027   -6.272    0.000
#>     E5|t5             0.790    0.028   27.848    0.000
#>     A1|t1            -0.403    0.029  -14.103    0.000
#>     A1|t2             0.363    0.026   13.948    0.000
#>     A1|t3             0.780    0.026   29.516    0.000
#>     A1|t4             1.270    0.030   42.199    0.000
#>     A1|t5             1.889    0.039   47.840    0.000
#>     A2|t1            -2.157    0.034  -64.008    0.000
#>     A2|t2            -1.558    0.027  -57.251    0.000
#>     A2|t3            -1.222    0.026  -46.985    0.000
#>     A2|t4            -0.512    0.027  -18.648    0.000
#>     A2|t5             0.456    0.035   13.091    0.000
#>     A3|t1            -1.849    0.031  -60.284    0.000
#>     A3|t2            -1.332    0.026  -50.800    0.000
#>     A3|t3            -0.995    0.026  -38.165    0.000
#>     A3|t4            -0.349    0.028  -12.304    0.000
#>     A3|t5             0.593    0.037   16.082    0.000
#>     A4|t1            -1.716    0.031  -54.930    0.000
#>     A4|t2            -1.206    0.027  -45.336    0.000
#>     A4|t3            -0.927    0.026  -35.346    0.000
#>     A4|t4            -0.423    0.029  -14.551    0.000
#>     A4|t5             0.202    0.037    5.503    0.000
#>     A5|t1            -2.014    0.029  -69.801    0.000
#>     A5|t2            -1.363    0.026  -51.849    0.000
#>     A5|t3            -0.945    0.027  -35.505    0.000
#>     A5|t4            -0.277    0.030   -9.313    0.000
#>     A5|t5             0.652    0.037   17.563    0.000
#>     N1|t1            -0.710    0.025  -28.770    0.000
#>     N1|t2            -0.059    0.023   -2.569    0.010
#>     N1|t3             0.336    0.026   13.099    0.000
#>     N1|t4             0.892    0.024   37.313    0.000
#>     N1|t5             1.488    0.023   63.701    0.000
#>     N2|t1            -1.163    0.020  -59.336    0.000
#>     N2|t2            -0.472    0.020  -24.002    0.000
#>     N2|t3            -0.094    0.019   -4.818    0.000
#>     N2|t4             0.571    0.016   34.667    0.000
#>     N2|t5             1.259    0.015   81.596    0.000
#>     N3|t1            -0.923    0.018  -51.121    0.000
#>     N3|t2            -0.210    0.013  -15.612    0.000
#>     N3|t3             0.106    0.015    7.078    0.000
#>     N3|t4             0.685    0.014   49.377    0.000
#>     N3|t5             1.350    0.016   84.750    0.000
#>     N4|t1            -0.960    0.021  -44.915    0.000
#>     N4|t2            -0.234    0.017  -13.490    0.000
#>     N4|t3             0.155    0.015   10.479    0.000
#>     N4|t4             0.767    0.020   39.032    0.000
#>     N4|t5             1.339    0.014   93.697    0.000
#>     N5|t1            -0.714    0.015  -48.961    0.000
#>     N5|t2            -0.035    0.015   -2.269    0.023
#>     N5|t3             0.298    0.015   19.716    0.000
#>     N5|t4             0.824    0.015   55.371    0.000
#>     N5|t5             1.357    0.017   78.962    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .O1                0.871                           
#>    .O2                0.928                           
#>    .O3                0.822                           
#>    .O4                0.994                           
#>    .O5                0.935                           
#>    .C1                0.735                           
#>    .C2                0.721                           
#>    .C3                0.753                           
#>    .C4                0.483                           
#>    .C5                0.534                           
#>    .E1                0.733                           
#>    .E2                0.502                           
#>    .E3                0.595                           
#>    .E4                0.483                           
#>    .E5                0.669                           
#>    .A1                0.908                           
#>    .A2                0.627                           
#>    .A3                0.491                           
#>    .A4                0.720                           
#>    .A5                0.451                           
#>    .N1                0.479                           
#>    .N2                0.553                           
#>    .N3                0.549                           
#>    .N4                0.552                           
#>    .N5                0.701                           
#>     opn               1.000                           
#>     con               1.000                           
#>     ext               1.000                           
#>     agr               1.000                           
#>     neu               1.000
```

All[^1] `lavaan` methods should carry over, such as:

``` r
coef(fit)
#>  opn=~O1  opn=~O2  opn=~O3  opn=~O4  opn=~O5  con=~C1  con=~C2  con=~C3 
#>    0.359   -0.267    0.422   -0.075   -0.255    0.515    0.528    0.497 
#>  con=~C4  con=~C5  ext=~E1  ext=~E2  ext=~E3  ext=~E4  ext=~E5  agr=~A1 
#>   -0.719   -0.682   -0.517   -0.705    0.637    0.719    0.575   -0.304 
#>  agr=~A2  agr=~A3  agr=~A4  agr=~A5  neu=~N1  neu=~N2  neu=~N3  neu=~N4 
#>    0.611    0.714    0.529    0.741    0.722    0.669    0.672    0.669 
#>  neu=~N5    O1|t1    O1|t2    O1|t3    O1|t4    O1|t5    O2|t1    O2|t2 
#>    0.546   -2.449   -1.695   -1.190   -0.429    0.455   -0.565    0.130 
#>    O2|t3    O2|t4    O2|t5    O3|t1    O3|t2    O3|t3    O3|t4    O3|t5 
#>    0.494    1.011    1.558   -1.992   -1.451   -0.955   -0.121    0.838 
#>    O4|t1    O4|t2    O4|t3    O4|t4    O4|t5    O5|t1    O5|t2    O5|t3 
#>   -2.189   -1.566   -1.245   -0.607    0.255   -0.584    0.260    0.770 
#>    O5|t4    O5|t5    C1|t1    C1|t2    C1|t3    C1|t4    C1|t5    C2|t1 
#>    1.341    1.934   -2.012   -1.432   -0.960   -0.268    0.736   -1.880 
#>    C2|t2    C2|t3    C2|t4    C2|t5    C3|t1    C3|t2    C3|t3    C3|t4 
#>   -1.203   -0.784   -0.146    0.827   -1.887   -1.191   -0.778   -0.022 
#>    C3|t5    C4|t1    C4|t2    C4|t3    C4|t4    C4|t5    C5|t1    C5|t2 
#>    0.921   -0.550    0.208    0.663    1.289    2.034   -0.886   -0.254 
#>    C5|t3    C5|t4    C5|t5    E1|t1    E1|t2    E1|t3    E1|t4    E1|t5 
#>    0.047    0.630    1.302   -0.728   -0.063    0.311    0.789    1.373 
#>    E2|t1    E2|t2    E2|t3    E2|t4    E2|t5    E3|t1    E3|t2    E3|t3 
#>   -0.867   -0.153    0.168    0.754    1.320   -1.625   -1.010   -0.518 
#>    E3|t4    E3|t5    E4|t1    E4|t2    E4|t3    E4|t4    E4|t5    E5|t1 
#>    0.267    1.147   -1.628   -1.061   -0.715   -0.262    0.640   -1.813 
#>    E5|t2    E5|t3    E5|t4    E5|t5    A1|t1    A1|t2    A1|t3    A1|t4 
#>   -1.198   -0.798   -0.168    0.790   -0.403    0.363    0.780    1.270 
#>    A1|t5    A2|t1    A2|t2    A2|t3    A2|t4    A2|t5    A3|t1    A3|t2 
#>    1.889   -2.157   -1.558   -1.222   -0.512    0.456   -1.849   -1.332 
#>    A3|t3    A3|t4    A3|t5    A4|t1    A4|t2    A4|t3    A4|t4    A4|t5 
#>   -0.995   -0.349    0.593   -1.716   -1.206   -0.927   -0.423    0.202 
#>    A5|t1    A5|t2    A5|t3    A5|t4    A5|t5    N1|t1    N1|t2    N1|t3 
#>   -2.014   -1.363   -0.945   -0.277    0.652   -0.710   -0.059    0.336 
#>    N1|t4    N1|t5    N2|t1    N2|t2    N2|t3    N2|t4    N2|t5    N3|t1 
#>    0.892    1.488   -1.163   -0.472   -0.094    0.571    1.259   -0.923 
#>    N3|t2    N3|t3    N3|t4    N3|t5    N4|t1    N4|t2    N4|t3    N4|t4 
#>   -0.210    0.106    0.685    1.350   -0.960   -0.234    0.155    0.767 
#>    N4|t5    N5|t1    N5|t2    N5|t3    N5|t4    N5|t5 opn~~con opn~~ext 
#>    1.339   -0.714   -0.035    0.298    0.824    1.357    0.599    0.784 
#> opn~~agr opn~~neu con~~ext con~~agr con~~neu ext~~agr ext~~neu agr~~neu 
#>    0.569   -0.882    0.445    0.453   -0.571    0.833   -0.527   -0.449
head(partable(fit))
#>   id lhs op rhs user block group free ustart exo label plabel  start    est
#> 1  1 opn =~  O1    1     1     1    1     NA   0         .p1.  0.614  0.359
#> 2  2 opn =~  O2    1     1     1    2     NA   0         .p2. -0.504 -0.267
#> 3  3 opn =~  O3    1     1     1    3     NA   0         .p3.  0.709  0.422
#> 4  4 opn =~  O4    1     1     1    4     NA   0         .p4.  0.346 -0.075
#> 5  5 opn =~  O5    1     1     1    5     NA   0         .p5. -0.590 -0.255
#> 6  6 con =~  C1    1     1     1    6     NA   0         .p6.  0.638  0.515
#>      se
#> 1 0.090
#> 2 0.046
#> 3 0.034
#> 4 0.027
#> 5 0.027
#> 6 0.028
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

[^1]: Most :)–we’re working on it!
