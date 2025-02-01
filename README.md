
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
#> Done! (4.64 secs)
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
#>     Total time: 4.64 s (Data reduction: 0.05 s)
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
#> Done! (6.71 secs)
str(var)
#> List of 8
#>  $ trJacob            : num [1:160, 1:160] 1 0 0 0 0 ...
#>  $ H                  : logi NA
#>  $ invH               : num [1:160, 1:160] 0.1043 0.0793 0.0736 0.0667 0.0582 ...
#>  $ J                  : num [1:160, 1:160] 484 -265 0 0 0 ...
#>  $ vcov               : num [1:160, 1:160] 7.02 6.15 5.33 3.74 4.08 ...
#>  $ asymptotic_variance: num [1:160] 7.02 8.48 10.16 16.8 43.41 ...
#>  $ optimisation_noise : logi NA
#>  $ RTime              : num 6.71
head(sqrt(diag(var$vcov)))
#> [1] 2.649855 2.911907 3.187016 4.099334 6.588509 6.434556
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
#> [1] 9.406044

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
#>     O1                0.359    4.250    0.085    0.933
#>     O2               -0.267    2.178   -0.123    0.902
#>     O3                0.422    1.626    0.260    0.795
#>     O4               -0.075    1.273   -0.059    0.953
#>     O5               -0.255    1.268   -0.201    0.841
#>   con =~                                              
#>     C1                0.515    1.306    0.394    0.693
#>     C2                0.528    1.240    0.426    0.670
#>     C3                0.497    1.296    0.384    0.701
#>     C4               -0.719    1.505   -0.478    0.633
#>     C5               -0.682    1.980   -0.345    0.730
#>   ext =~                                              
#>     E1               -0.517    2.705   -0.191    0.848
#>     E2               -0.705    1.851   -0.381    0.703
#>     E3                0.637    1.468    0.434    0.665
#>     E4                0.719    1.226    0.587    0.557
#>     E5                0.575    1.392    0.413    0.679
#>   agr =~                                              
#>     A1               -0.304    3.280   -0.093    0.926
#>     A2                0.611    2.004    0.305    0.761
#>     A3                0.714    1.672    0.427    0.670
#>     A4                0.529    1.325    0.399    0.690
#>     A5                0.741    1.239    0.598    0.550
#>   neu =~                                              
#>     N1                0.722    1.302    0.555    0.579
#>     N2                0.669    1.253    0.534    0.594
#>     N3                0.672    1.388    0.484    0.628
#>     N4                0.669    1.749    0.383    0.702
#>     N5                0.546    2.600    0.210    0.834
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn ~~                                              
#>     con               0.599    3.901    0.153    0.878
#>     ext               0.784    6.294    0.125    0.901
#>     agr               0.569    6.013    0.095    0.925
#>     neu              -0.882    3.862   -0.228    0.819
#>   con ~~                                              
#>     ext               0.445    3.344    0.133    0.894
#>     agr               0.453    5.715    0.079    0.937
#>     neu              -0.571   13.545   -0.042    0.966
#>   ext ~~                                              
#>     agr               0.833    7.902    0.105    0.916
#>     neu              -0.527   31.155   -0.017    0.987
#>   agr ~~                                              
#>     neu              -0.449   42.127   -0.011    0.992
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     O1|t1            -2.449    2.744   -0.892    0.372
#>     O1|t2            -1.695    1.836   -0.923    0.356
#>     O1|t3            -1.190    1.476   -0.806    0.420
#>     O1|t4            -0.429    1.255   -0.342    0.732
#>     O1|t5             0.455    1.363    0.334    0.739
#>     O2|t1            -0.565    2.479   -0.228    0.820
#>     O2|t2             0.130    1.626    0.080    0.936
#>     O2|t3             0.494    1.391    0.355    0.723
#>     O2|t4             1.011    1.240    0.815    0.415
#>     O2|t5             1.558    1.395    1.117    0.264
#>     O3|t1            -1.992    2.491   -0.800    0.424
#>     O3|t2            -1.451    1.626   -0.893    0.372
#>     O3|t3            -0.955    1.393   -0.686    0.493
#>     O3|t4            -0.121    1.245   -0.097    0.923
#>     O3|t5             0.838    1.445    0.580    0.562
#>     O4|t1            -2.189    1.312   -1.669    0.095
#>     O4|t2            -1.566    1.253   -1.250    0.211
#>     O4|t3            -1.245    1.346   -0.925    0.355
#>     O4|t4            -0.607    1.698   -0.357    0.721
#>     O4|t5             0.255    2.826    0.090    0.928
#>     O5|t1            -0.584    1.435   -0.407    0.684
#>     O5|t2             0.260    1.261    0.206    0.837
#>     O5|t3             0.770    1.245    0.618    0.536
#>     O5|t4             1.341    1.340    1.001    0.317
#>     O5|t5             1.934    1.710    1.132    0.258
#>     C1|t1            -2.012    1.370   -1.468    0.142
#>     C1|t2            -1.432    1.242   -1.153    0.249
#>     C1|t3            -0.960    1.259   -0.763    0.446
#>     C1|t4            -0.268    1.390   -0.193    0.847
#>     C1|t5             0.736    1.773    0.415    0.678
#>     C2|t1            -1.880    1.430   -1.315    0.189
#>     C2|t2            -1.203    1.249   -0.963    0.335
#>     C2|t3            -0.784    1.246   -0.630    0.529
#>     C2|t4            -0.146    1.376   -0.106    0.915
#>     C2|t5             0.827    1.716    0.482    0.630
#>     C3|t1            -1.887    2.060   -0.916    0.360
#>     C3|t2            -1.191    1.499   -0.795    0.427
#>     C3|t3            -0.778    1.303   -0.597    0.550
#>     C3|t4            -0.022    1.253   -0.018    0.986
#>     C3|t5             0.921    1.577    0.584    0.559
#>     C4|t1            -0.550    2.051   -0.268    0.788
#>     C4|t2             0.208    1.521    0.137    0.891
#>     C4|t3             0.663    1.364    0.486    0.627
#>     C4|t4             1.289    1.259    1.024    0.306
#>     C4|t5             2.034    1.340    1.518    0.129
#>     C5|t1            -0.886    2.347   -0.377    0.706
#>     C5|t2            -0.254    1.622   -0.157    0.876
#>     C5|t3             0.047    1.397    0.034    0.973
#>     C5|t4             0.630    1.246    0.505    0.613
#>     C5|t5             1.302    1.380    0.944    0.345
#>     E1|t1            -0.728    1.283   -0.567    0.571
#>     E1|t2            -0.063    1.273   -0.049    0.961
#>     E1|t3             0.311    1.392    0.224    0.823
#>     E1|t4             0.789    1.682    0.469    0.639
#>     E1|t5             1.373    2.489    0.552    0.581
#>     E2|t1            -0.867    3.126   -0.277    0.782
#>     E2|t2            -0.153    1.969   -0.078    0.938
#>     E2|t3             0.168    1.644    0.102    0.918
#>     E2|t4             0.754    1.305    0.578    0.563
#>     E2|t5             1.320    1.289    1.024    0.306
#>     E3|t1            -1.625    2.404   -0.676    0.499
#>     E3|t2            -1.010    1.722   -0.587    0.557
#>     E3|t3            -0.518    1.487   -0.348    0.728
#>     E3|t4             0.267    1.272    0.210    0.834
#>     E3|t5             1.147    1.323    0.867    0.386
#>     E4|t1            -1.628    2.197   -0.741    0.459
#>     E4|t2            -1.061    1.636   -0.648    0.517
#>     E4|t3            -0.715    1.464   -0.488    0.625
#>     E4|t4            -0.262    1.285   -0.204    0.838
#>     E4|t5             0.640    1.252    0.511    0.609
#>     E5|t1            -1.813    2.753   -0.659    0.510
#>     E5|t2            -1.198    1.762   -0.680    0.496
#>     E5|t3            -0.798    1.470   -0.543    0.587
#>     E5|t4            -0.168    1.264   -0.133    0.894
#>     E5|t5             0.790    1.341    0.589    0.556
#>     A1|t1            -0.403    1.353   -0.298    0.766
#>     A1|t2             0.363    1.229    0.295    0.768
#>     A1|t3             0.780    1.249    0.624    0.533
#>     A1|t4             1.270    1.423    0.892    0.372
#>     A1|t5             1.889    1.867    1.012    0.312
#>     A2|t1            -2.157    1.593   -1.354    0.176
#>     A2|t2            -1.558    1.287   -1.211    0.226
#>     A2|t3            -1.222    1.230   -0.994    0.320
#>     A2|t4            -0.512    1.299   -0.394    0.693
#>     A2|t5             0.456    1.645    0.277    0.782
#>     A3|t1            -1.849    1.450   -1.275    0.202
#>     A3|t2            -1.332    1.240   -1.074    0.283
#>     A3|t3            -0.995    1.233   -0.807    0.420
#>     A3|t4            -0.349    1.342   -0.260    0.795
#>     A3|t5             0.593    1.744    0.340    0.734
#>     A4|t1            -1.716    1.477   -1.162    0.245
#>     A4|t2            -1.206    1.258   -0.959    0.338
#>     A4|t3            -0.927    1.240   -0.747    0.455
#>     A4|t4            -0.423    1.374   -0.308    0.758
#>     A4|t5             0.202    1.734    0.116    0.907
#>     A5|t1            -2.014    1.365   -1.476    0.140
#>     A5|t2            -1.363    1.243   -1.096    0.273
#>     A5|t3            -0.945    1.259   -0.751    0.453
#>     A5|t4            -0.277    1.404   -0.197    0.844
#>     A5|t5             0.652    1.755    0.371    0.710
#>     N1|t1            -0.710    1.168   -0.608    0.543
#>     N1|t2            -0.059    1.095   -0.054    0.957
#>     N1|t3             0.336    1.212    0.277    0.782
#>     N1|t4             0.892    1.131    0.789    0.430
#>     N1|t5             1.488    1.104    1.347    0.178
#>     N2|t1            -1.163    0.927   -1.255    0.210
#>     N2|t2            -0.472    0.931   -0.508    0.612
#>     N2|t3            -0.094    0.919   -0.102    0.919
#>     N2|t4             0.571    0.778    0.733    0.463
#>     N2|t5             1.259    0.730    1.726    0.084
#>     N3|t1            -0.923    0.853   -1.081    0.280
#>     N3|t2            -0.210    0.636   -0.330    0.741
#>     N3|t3             0.106    0.708    0.150    0.881
#>     N3|t4             0.685    0.656    1.044    0.296
#>     N3|t5             1.350    0.753    1.792    0.073
#>     N4|t1            -0.960    1.010   -0.950    0.342
#>     N4|t2            -0.234    0.822   -0.285    0.775
#>     N4|t3             0.155    0.699    0.222    0.825
#>     N4|t4             0.767    0.930    0.825    0.409
#>     N4|t5             1.339    0.676    1.981    0.048
#>     N5|t1            -0.714    0.689   -1.035    0.300
#>     N5|t2            -0.035    0.727   -0.048    0.962
#>     N5|t3             0.298    0.714    0.417    0.677
#>     N5|t4             0.824    0.704    1.171    0.242
#>     N5|t5             1.357    0.812    1.670    0.095
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .O1                1.000                           
#>    .O2                1.000                           
#>    .O3                1.000                           
#>    .O4                1.000                           
#>    .O5                1.000                           
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
#>    .A1                1.000                           
#>    .A2                1.000                           
#>    .A3                1.000                           
#>    .A4                1.000                           
#>    .A5                1.000                           
#>    .N1                1.000                           
#>    .N2                1.000                           
#>    .N3                1.000                           
#>    .N4                1.000                           
#>    .N5                1.000                           
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
#> 1 4.250
#> 2 2.178
#> 3 1.626
#> 4 1.273
#> 5 1.268
#> 6 1.306
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
