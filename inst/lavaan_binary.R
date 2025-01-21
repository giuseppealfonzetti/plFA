library(plFA)
library(lavaan)
set.seed(123)

p <- 10  # number of items
q <- 2  # number of latent variables
n <- 2000  # number of observations

pop_model <- "
  # Factor loadings
  eta1 =~ 0.8 * y1 + 0.7 * y2 + 0.47 * y3 + 0.38 * y4 + 0.34 * y5
  eta2 =~ 0.8 * y6 + 0.7 * y7 + 0.47 * y8 + 0.38 * y9 + 0.34 * y10

  # Factor correlation
  eta1 ~~ 0.3 * eta2

  # Thresholds for ordinal items
  y1 | -1.43*t1
  y2 | -0.55*t1
  y3 | -0.13*t1
  y4 | -0.72*t1
  y5 | -1.13*t1
  y6 | -1.43*t1
  y7 | -0.55*t1
  y8 | -0.13*t1
  y9 | -0.72*t1
  y10 | -1.13*t1
"

DAT <- simulateData(pop_model, sample.nobs = n)
dat <- as.data.frame(lapply(DAT, ordered))
D <- as.matrix(DAT) - 1
mod <- "
  eta1 =~ y1 + y2 + y3 + y4 + y5
  eta2 =~ y6 + y7 + y8 + y9 + y10
"

tictoc::tic("DWLS Estimation")
fit_lav1 <- lavaan::cfa(mod, dat, std.lv = !TRUE)
coef(fit_lav1)
tictoc::toc()
tictoc::tic("Pairwise ML Estimation")
fit_lav2 <- lavaan::cfa(mod, dat, std.lv = TRUE, estimator = "PML")
coef(fit_lav2)
tictoc::toc()

A <- build_constrMat(P = p, Q = q, STRUCT = "simple")
fit_plFA <- fit_plFA(
  DATA = D,
  CONSTR_LIST = list("CONSTRMAT" = A, "CORRFLAG" = 1),
  METHOD = "hyper"
)

# These two error
getPar(fit_plFA, OPTION = "raw")
computeVar(OBJ = fit_plFA, DATA = D, NUMDERIV = TRUE)
