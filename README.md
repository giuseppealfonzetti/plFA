
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
models using the `lavaan` syntax. An example is given below.

## Fitting an ordinal CFA model

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

The `{lavaan}` model syntax to fit a confirmatory factor analysis model
with 5 latent factors, where the first 5 items load on the first factor
and so on, is:

``` r
mod <- "
  opn =~ O1 + O2 + O3 + O4 + O5  # Openness
  con =~ C1 + C2 + C3 + C4 + C5  # Conscientiousness
  ext =~ E1 + E2 + E3 + E4 + E5  # Extraversion
  agr =~ A1 + A2 + A3 + A4 + A5  # Agreeableness
  neu =~ N1 + N2 + N3 + N4 + N5  # Neuroticism
"
```

Users of `{lavaan}` will be familiar with this syntax. The call to the
`cfa()` function is made as usual:

``` r
# Fit the model (note, this cfa() function is from the plFA package)
fit <- cfa(model = mod, data = bfi, std.lv = TRUE)
summary(fit)
#> plFA 0.1.0 
#>   ⨉
#> lavaan 0.6-19 did NOT end normally after 107 iterations
#> ** WARNING ** Estimates below are most likely unreliable
#> 
#>   Estimator                                        PML
#>   Optimization method                           UCMINF
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
#>     O1                0.647    0.126    5.156    0.000
#>     O2               -0.444    0.264   -1.684    0.092
#>     O3                0.808    0.283    2.850    0.004
#>     O4                0.136    0.199    0.681    0.496
#>     O5               -0.486    0.307   -1.582    0.114
#>   con =~                                              
#>     C1                0.564    0.267    2.113    0.035
#>     C2                0.571    0.124    4.600    0.000
#>     C3                0.543    0.154    3.517    0.000
#>     C4               -0.760    0.222   -3.425    0.001
#>     C5               -0.706    0.233   -3.030    0.002
#>   ext =~                                              
#>     E1               -0.540    0.184   -2.931    0.003
#>     E2               -0.727    0.194   -3.753    0.000
#>     E3                0.690    0.138    5.005    0.000
#>     E4                0.742    0.122    6.093    0.000
#>     E5                0.620    0.225    2.760    0.006
#>   agr =~                                              
#>     A1               -0.335    0.115   -2.927    0.003
#>     A2                0.663    0.160    4.149    0.000
#>     A3                0.771    0.083    9.258    0.000
#>     A4                0.545    0.119    4.573    0.000
#>     A5                0.774    0.183    4.217    0.000
#>   neu =~                                              
#>     N1                0.837    0.264    3.167    0.002
#>     N2                0.801    0.147    5.454    0.000
#>     N3                0.755    0.145    5.223    0.000
#>     N4                0.666    0.278    2.400    0.016
#>     N5                0.570    0.213    2.680    0.007
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   opn ~~                                              
#>     con               0.322    0.092    3.507    0.000
#>     ext               0.466    0.071    6.593    0.000
#>     agr               0.269    0.090    2.978    0.003
#>     neu              -0.146    0.084   -1.738    0.082
#>   con ~~                                              
#>     ext               0.390    0.074    5.297    0.000
#>     agr               0.370    0.085    4.326    0.000
#>     neu              -0.301    0.089   -3.376    0.001
#>   ext ~~                                              
#>     agr               0.702    0.091    7.678    0.000
#>     neu              -0.267    0.091   -2.927    0.003
#>   agr ~~                                              
#>     neu              -0.232    0.057   -4.047    0.000
#> 
#> Thresholds:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>     O1|t1            -2.447    0.181  -13.507    0.000
#>     O1|t2            -1.695    0.186   -9.135    0.000
#>     O1|t3            -1.190    0.149   -8.008    0.000
#>     O1|t4            -0.430    0.120   -3.582    0.000
#>     O1|t5             0.455    0.217    2.093    0.036
#>     O2|t1            -0.565    0.194   -2.911    0.004
#>     O2|t2             0.131    0.224    0.583    0.560
#>     O2|t3             0.495    0.152    3.251    0.001
#>     O2|t4             1.011    0.221    4.570    0.000
#>     O2|t5             1.557    0.222    7.025    0.000
#>     O3|t1            -1.990    0.197  -10.098    0.000
#>     O3|t2            -1.451    0.270   -5.368    0.000
#>     O3|t3            -0.956    0.171   -5.586    0.000
#>     O3|t4            -0.122    0.318   -0.384    0.701
#>     O3|t5             0.839    0.298    2.818    0.005
#>     O4|t1            -2.187    0.356   -6.149    0.000
#>     O4|t2            -1.567    0.223   -7.035    0.000
#>     O4|t3            -1.246    0.264   -4.719    0.000
#>     O4|t4            -0.607    0.332   -1.830    0.067
#>     O4|t5             0.255    0.191    1.336    0.182
#>     O5|t1            -0.584    0.160   -3.661    0.000
#>     O5|t2             0.261    0.109    2.401    0.016
#>     O5|t3             0.771    0.109    7.080    0.000
#>     O5|t4             1.341    0.113   11.905    0.000
#>     O5|t5             1.931    0.181   10.682    0.000
#>     C1|t1            -2.011    0.361   -5.567    0.000
#>     C1|t2            -1.432    0.145   -9.857    0.000
#>     C1|t3            -0.960    0.218   -4.401    0.000
#>     C1|t4            -0.269    0.228   -1.178    0.239
#>     C1|t5             0.736    0.291    2.528    0.011
#>     C2|t1            -1.880    0.189   -9.943    0.000
#>     C2|t2            -1.203    0.110  -10.968    0.000
#>     C2|t3            -0.785    0.110   -7.104    0.000
#>     C2|t4            -0.146    0.114   -1.279    0.201
#>     C2|t5             0.828    0.184    4.490    0.000
#>     C3|t1            -1.887    0.219   -8.596    0.000
#>     C3|t2            -1.191    0.148   -8.033    0.000
#>     C3|t3            -0.778    0.156   -5.004    0.000
#>     C3|t4            -0.022    0.144   -0.156    0.876
#>     C3|t5             0.921    0.225    4.093    0.000
#>     C4|t1            -0.550    0.211   -2.614    0.009
#>     C4|t2             0.208    0.099    2.109    0.035
#>     C4|t3             0.664    0.115    5.750    0.000
#>     C4|t4             1.289    0.132    9.734    0.000
#>     C4|t5             2.033    0.145   14.022    0.000
#>     C5|t1            -0.886    0.220   -4.033    0.000
#>     C5|t2            -0.254    0.221   -1.151    0.250
#>     C5|t3             0.047    0.173    0.274    0.784
#>     C5|t4             0.630    0.182    3.459    0.001
#>     C5|t5             1.303    0.294    4.426    0.000
#>     E1|t1            -0.728    0.105   -6.932    0.000
#>     E1|t2            -0.063    0.123   -0.512    0.609
#>     E1|t3             0.311    0.125    2.487    0.013
#>     E1|t4             0.789    0.109    7.243    0.000
#>     E1|t5             1.373    0.190    7.242    0.000
#>     E2|t1            -0.867    0.145   -5.975    0.000
#>     E2|t2            -0.153    0.171   -0.892    0.373
#>     E2|t3             0.169    0.120    1.403    0.161
#>     E2|t4             0.754    0.141    5.351    0.000
#>     E2|t5             1.320    0.205    6.442    0.000
#>     E3|t1            -1.625    0.186   -8.737    0.000
#>     E3|t2            -1.010    0.115   -8.798    0.000
#>     E3|t3            -0.518    0.124   -4.179    0.000
#>     E3|t4             0.267    0.110    2.425    0.015
#>     E3|t5             1.148    0.121    9.463    0.000
#>     E4|t1            -1.628    0.183   -8.899    0.000
#>     E4|t2            -1.061    0.112   -9.444    0.000
#>     E4|t3            -0.716    0.103   -6.914    0.000
#>     E4|t4            -0.263    0.097   -2.697    0.007
#>     E4|t5             0.640    0.113    5.669    0.000
#>     E5|t1            -1.813    0.183   -9.913    0.000
#>     E5|t2            -1.199    0.217   -5.529    0.000
#>     E5|t3            -0.798    0.158   -5.051    0.000
#>     E5|t4            -0.168    0.148   -1.134    0.257
#>     E5|t5             0.790    0.228    3.462    0.001
#>     A1|t1            -0.403    0.306   -1.320    0.187
#>     A1|t2             0.363    0.129    2.815    0.005
#>     A1|t3             0.780    0.197    3.963    0.000
#>     A1|t4             1.270    0.265    4.790    0.000
#>     A1|t5             1.888    0.249    7.592    0.000
#>     A2|t1            -2.156    0.257   -8.382    0.000
#>     A2|t2            -1.558    0.150  -10.403    0.000
#>     A2|t3            -1.223    0.182   -6.711    0.000
#>     A2|t4            -0.513    0.149   -3.439    0.001
#>     A2|t5             0.456    0.225    2.022    0.043
#>     A3|t1            -1.849    0.250   -7.386    0.000
#>     A3|t2            -1.332    0.136   -9.771    0.000
#>     A3|t3            -0.995    0.160   -6.219    0.000
#>     A3|t4            -0.350    0.152   -2.297    0.022
#>     A3|t5             0.593    0.227    2.619    0.009
#>     A4|t1            -1.716    0.223   -7.687    0.000
#>     A4|t2            -1.206    0.141   -8.583    0.000
#>     A4|t3            -0.927    0.153   -6.051    0.000
#>     A4|t4            -0.423    0.164   -2.576    0.010
#>     A4|t5             0.202    0.222    0.907    0.365
#>     A5|t1            -2.015    0.240   -8.388    0.000
#>     A5|t2            -1.364    0.143   -9.505    0.000
#>     A5|t3            -0.946    0.154   -6.146    0.000
#>     A5|t4            -0.277    0.149   -1.862    0.063
#>     A5|t5             0.652    0.234    2.788    0.005
#>     N1|t1            -0.711    0.065  -10.933    0.000
#>     N1|t2            -0.060    0.068   -0.890    0.373
#>     N1|t3             0.335    0.057    5.854    0.000
#>     N1|t4             0.892    0.067   13.295    0.000
#>     N1|t5             1.489    0.066   22.629    0.000
#>     N2|t1            -1.163    0.079  -14.756    0.000
#>     N2|t2            -0.473    0.074   -6.395    0.000
#>     N2|t3            -0.094    0.080   -1.179    0.238
#>     N2|t4             0.570    0.070    8.128    0.000
#>     N2|t5             1.259    0.082   15.302    0.000
#>     N3|t1            -0.924    0.087  -10.628    0.000
#>     N3|t2            -0.210    0.091   -2.321    0.020
#>     N3|t3             0.106    0.086    1.235    0.217
#>     N3|t4             0.685    0.075    9.103    0.000
#>     N3|t5             1.351    0.080   16.910    0.000
#>     N4|t1            -0.961    0.077  -12.503    0.000
#>     N4|t2            -0.235    0.073   -3.214    0.001
#>     N4|t3             0.155    0.063    2.442    0.015
#>     N4|t4             0.768    0.076   10.070    0.000
#>     N4|t5             1.341    0.079   17.004    0.000
#>     N5|t1            -0.714    0.071  -10.109    0.000
#>     N5|t2            -0.035    0.074   -0.477    0.633
#>     N5|t3             0.297    0.076    3.925    0.000
#>     N5|t4             0.825    0.084    9.792    0.000
#>     N5|t5             1.358    0.073   18.612    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .O1                0.581                           
#>    .O2                0.803                           
#>    .O3                0.347                           
#>    .O4                0.982                           
#>    .O5                0.764                           
#>    .C1                0.681                           
#>    .C2                0.674                           
#>    .C3                0.705                           
#>    .C4                0.422                           
#>    .C5                0.502                           
#>    .E1                0.708                           
#>    .E2                0.471                           
#>    .E3                0.524                           
#>    .E4                0.450                           
#>    .E5                0.615                           
#>    .A1                0.888                           
#>    .A2                0.561                           
#>    .A3                0.406                           
#>    .A4                0.703                           
#>    .A5                0.401                           
#>    .N1                0.300                           
#>    .N2                0.359                           
#>    .N3                0.429                           
#>    .N4                0.556                           
#>    .N5                0.675                           
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
#>    0.647   -0.444    0.808    0.136   -0.486    0.564    0.571    0.543 
#>  con=~C4  con=~C5  ext=~E1  ext=~E2  ext=~E3  ext=~E4  ext=~E5  agr=~A1 
#>   -0.760   -0.706   -0.540   -0.727    0.690    0.742    0.620   -0.335 
#>  agr=~A2  agr=~A3  agr=~A4  agr=~A5  neu=~N1  neu=~N2  neu=~N3  neu=~N4 
#>    0.663    0.771    0.545    0.774    0.837    0.801    0.755    0.666 
#>  neu=~N5    O1|t1    O1|t2    O1|t3    O1|t4    O1|t5    O2|t1    O2|t2 
#>    0.570   -2.447   -1.695   -1.190   -0.430    0.455   -0.565    0.131 
#>    O2|t3    O2|t4    O2|t5    O3|t1    O3|t2    O3|t3    O3|t4    O3|t5 
#>    0.495    1.011    1.557   -1.990   -1.451   -0.956   -0.122    0.839 
#>    O4|t1    O4|t2    O4|t3    O4|t4    O4|t5    O5|t1    O5|t2    O5|t3 
#>   -2.187   -1.567   -1.246   -0.607    0.255   -0.584    0.261    0.771 
#>    O5|t4    O5|t5    C1|t1    C1|t2    C1|t3    C1|t4    C1|t5    C2|t1 
#>    1.341    1.931   -2.011   -1.432   -0.960   -0.269    0.736   -1.880 
#>    C2|t2    C2|t3    C2|t4    C2|t5    C3|t1    C3|t2    C3|t3    C3|t4 
#>   -1.203   -0.785   -0.146    0.828   -1.887   -1.191   -0.778   -0.022 
#>    C3|t5    C4|t1    C4|t2    C4|t3    C4|t4    C4|t5    C5|t1    C5|t2 
#>    0.921   -0.550    0.208    0.664    1.289    2.033   -0.886   -0.254 
#>    C5|t3    C5|t4    C5|t5    E1|t1    E1|t2    E1|t3    E1|t4    E1|t5 
#>    0.047    0.630    1.303   -0.728   -0.063    0.311    0.789    1.373 
#>    E2|t1    E2|t2    E2|t3    E2|t4    E2|t5    E3|t1    E3|t2    E3|t3 
#>   -0.867   -0.153    0.169    0.754    1.320   -1.625   -1.010   -0.518 
#>    E3|t4    E3|t5    E4|t1    E4|t2    E4|t3    E4|t4    E4|t5    E5|t1 
#>    0.267    1.148   -1.628   -1.061   -0.716   -0.263    0.640   -1.813 
#>    E5|t2    E5|t3    E5|t4    E5|t5    A1|t1    A1|t2    A1|t3    A1|t4 
#>   -1.199   -0.798   -0.168    0.790   -0.403    0.363    0.780    1.270 
#>    A1|t5    A2|t1    A2|t2    A2|t3    A2|t4    A2|t5    A3|t1    A3|t2 
#>    1.888   -2.156   -1.558   -1.223   -0.513    0.456   -1.849   -1.332 
#>    A3|t3    A3|t4    A3|t5    A4|t1    A4|t2    A4|t3    A4|t4    A4|t5 
#>   -0.995   -0.350    0.593   -1.716   -1.206   -0.927   -0.423    0.202 
#>    A5|t1    A5|t2    A5|t3    A5|t4    A5|t5    N1|t1    N1|t2    N1|t3 
#>   -2.015   -1.364   -0.946   -0.277    0.652   -0.711   -0.060    0.335 
#>    N1|t4    N1|t5    N2|t1    N2|t2    N2|t3    N2|t4    N2|t5    N3|t1 
#>    0.892    1.489   -1.163   -0.473   -0.094    0.570    1.259   -0.924 
#>    N3|t2    N3|t3    N3|t4    N3|t5    N4|t1    N4|t2    N4|t3    N4|t4 
#>   -0.210    0.106    0.685    1.351   -0.961   -0.235    0.155    0.768 
#>    N4|t5    N5|t1    N5|t2    N5|t3    N5|t4    N5|t5 opn~~con opn~~ext 
#>    1.341   -0.714   -0.035    0.297    0.825    1.358    0.322    0.466 
#> opn~~agr opn~~neu con~~ext con~~agr con~~neu ext~~agr ext~~neu agr~~neu 
#>    0.269   -0.146    0.390    0.370   -0.301    0.702   -0.267   -0.232
head(parameterEstimates(fit))
#>   lhs op rhs    est
#> 1 opn =~  O1  0.647
#> 2 opn =~  O2 -0.444
#> 3 opn =~  O3  0.808
#> 4 opn =~  O4  0.136
#> 5 opn =~  O5 -0.486
#> 6 con =~  C1  0.564
head(fitted(fit)$cov)  # model-implied covariance matrix for the UVs
#>             O1          O2         O3          O4          O5          C1
#> O1  1.00000000 -0.28766473  0.5229862  0.08776934 -0.31480932  0.11770367
#> O2 -0.28766473  1.00000000 -0.3589212 -0.06023538  0.21605108 -0.08077907
#> O3  0.52298624 -0.35892120  1.0000000  0.10951039 -0.39278969  0.14685966
#> O4  0.08776934 -0.06023538  0.1095104  1.00000000 -0.06591931  0.02464649
#> O5 -0.31480932  0.21605108 -0.3927897 -0.06591931  1.00000000 -0.08840154
#> C1  0.11770367 -0.08077907  0.1468597  0.02464649 -0.08840154  1.00000000
#>             C2          C3          C4          C5         E1         E2
#> O1  0.11902338  0.11323967 -0.15854985 -0.14714078 -0.1629614 -0.2193540
#> O2 -0.08168478 -0.07771546  0.10881147  0.10098152  0.1118391  0.1505409
#> O3  0.14850627  0.14128989 -0.19782371 -0.18358854 -0.2033280 -0.2736894
#> O4  0.02492283  0.02371175 -0.03319945 -0.03081046 -0.0341232 -0.0459315
#> O5 -0.08939271 -0.08504884  0.11907913  0.11051033  0.1223924  0.1647462
#> C1  0.32212968  0.30647641 -0.42910572 -0.39822776 -0.1189600 -0.1601260
#>             E3          E4          E5          A1          A2          A3
#> O1  0.20823024  0.22383148  0.18712833 -0.05848417  0.11559788  0.13439663
#> O2 -0.14290672 -0.15361373 -0.12842465  0.04013721 -0.07933389 -0.09223532
#> O3  0.25981026  0.27927604  0.23348126 -0.07297109  0.14423225  0.16768757
#> O4  0.04360225  0.04686906  0.03918362 -0.01224626  0.02420555  0.02814190
#> O5 -0.15639167 -0.16810901 -0.14054304  0.04392464 -0.08681999 -0.10093882
#> C1  0.15200579  0.16339452  0.13660162 -0.06997619  0.13831263  0.16080529
#>             A4          A5          N1          N2          N3          N4
#> O1  0.09506869  0.13494464 -0.07927433 -0.07588712 -0.07159264 -0.06315254
#> O2 -0.06524487 -0.09261141  0.05440533  0.05208071  0.04913345  0.04334107
#> O3  0.11861783  0.16837132 -0.09891111 -0.09468487 -0.08932661 -0.07879584
#> O4  0.01990685  0.02825665 -0.01659960 -0.01589034 -0.01499110 -0.01322379
#> O5 -0.07140150 -0.10135040  0.05953912  0.05699515  0.05376977  0.04743082
#> C1  0.11374949  0.16146098 -0.14224881 -0.13617084 -0.12846488 -0.11332008
#>             N5
#> O1 -0.05398299
#> O2  0.03704809
#> O3 -0.06735493
#> O4 -0.01130374
#> O5  0.04054401
#> C1 -0.09686636
```

## Session information

``` r
sessioninfo::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.4.1 (2024-06-14)
#>  os       macOS 15.3
#>  system   aarch64, darwin20
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       Asia/Brunei
#>  date     2025-02-05
#>  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package      * version    date (UTC) lib source
#>  cli            3.6.3      2024-06-21 [1] CRAN (R 4.4.0)
#>  colorspace     2.1-1      2024-07-26 [1] CRAN (R 4.4.0)
#>  digest         0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
#>  dplyr          1.1.4      2023-11-17 [1] CRAN (R 4.4.0)
#>  evaluate       1.0.3      2025-01-10 [1] CRAN (R 4.4.1)
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
#>  plFA         * 0.1.0      2025-02-05 [1] local
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
#>  xfun           0.50       2025-01-07 [1] CRAN (R 4.4.1)
#>  yaml           2.3.10     2024-07-26 [1] CRAN (R 4.4.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```

<!-- To fit a latent factor analysis model with 5 factors, where the first 5 items load on the first factor and so on, we can use the `fit_plFA()` function. -->
<!-- In preparation, we must build the following objects to be supplied as a list to the `CONSTR_LIST` argument: -->
<!-- - `CONSTRMAT`: A matrix of constraints on the factor loadings. -->
<!-- - `CORRFLAG`: A flag indicating whether the factor correlations are estimated. -->
<!-- - `CONSTRVAR`: A vector of constraints on the factor variances, if any. -->
<!-- - `STDLV`: A flag indicating whether the factor loadings are standardised. -->
<!-- - `LLC`: A list indicating linear constraints on factor loading estimations, if any. -->
<!-- The above ingredients are prepared using the code below. -->
<!-- ```{r} -->
<!-- p <- 25         # number of items -->
<!-- q <- 5          # number of factors -->
<!-- n <- nrow(bfi)  # number of observations -->
<!-- # Convert data frame to matrix, with each item taking values 0,1,2,...,5 -->
<!-- D <- as.matrix(data.frame(lapply(bfi[, 1:p], as.numeric))) - 1 -->
<!-- head(D) -->
<!-- # Constraint matrix: NA indicate free loadings -->
<!-- constrmat <- build_constrMat(P = p, Q = q, STRUCT = "simple") -->
<!-- head(constrmat, 12) -->
<!-- # Build constraint list -->
<!-- constr_list <- list( -->
<!--   CONSTRMAT = constrmat, -->
<!--   CORRFLAG = TRUE, -->
<!--   CONSTRVAR = rep(1, q),  # these options standardise the   -->
<!--   STDLV = TRUE,           # factor variances -->
<!--   LLC = NULL -->
<!-- ) -->
<!-- ``` -->
<!-- We can now fit the factor model: -->
<!-- ```{r} -->
<!-- fit <- fit_plFA( -->
<!--   DATA = D,  -->
<!--   CONSTR_LIST = constr_list, -->
<!--   METHOD = "ucminf", -->
<!--   VERBOSE = TRUE -->
<!-- ) -->
<!-- print(fit) -->
<!-- ``` -->
<!-- The default method is using the optimiser from the `{ucminf}` package, which works well for small data sizes. -->
<!-- Alternatively, the user may choose `METHOD = "SA"` to turn on the stochastic approximation algorithm. -->
<!-- Use the function `getPar()` to extract the parameter estimates: -->
<!-- ```{r} -->
<!-- getPar(fit) -->
<!-- ``` -->
<!-- Computing the standard errors is done separately using the `computeVar()` function: -->
<!-- ```{r} -->
<!-- var <- computeVar(OBJ = fit, DATA = D, VERBOSE = TRUE) -->
<!-- str(var) -->
<!-- # Standard errors -->
<!-- head(sqrt(diag(var$vcov) / n)) -->
<!-- ``` -->
<!-- ## `{lavaan}` wrapper -->
<!-- This package also provides a user-friendly `{lavaan}`-style interface to estimate factor analysis models using the `cfa()` function. -->
<!-- Let's rerun the above example, but this time using the stochastic approximation algorithm. -->

[^1]: Most :)–we’re working on it!
