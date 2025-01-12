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
  y1 | -1.5*t1 + 0*t2 + 1.5*t3
  y2 | -1.5*t1 + 0*t2 + 1.5*t3
  y3 | -1.5*t1 + 0*t2 + 1.5*t3
  y4 | -1.5*t1 + 0*t2 + 1.5*t3
  y5 | -1.5*t1 + 0*t2 + 1.5*t3
  y6 | -1.5*t1 + 0*t2 + 1.5*t3
  y7 | -1.5*t1 + 0*t2 + 1.5*t3
  y8 | -1.5*t1 + 0*t2 + 1.5*t3
  y9 | -1.5*t1 + 0*t2 + 1.5*t3
  y10 | -1.5*t1 + 0*t2 + 1.5*t3
"

DAT <- simulateData(pop_model, sample.nobs = n)
dat <- as.data.frame(lapply(DAT, ordered))
D <- as.matrix(DAT) - 1
mod <- "
  eta1 =~ y1 + y2 + y3 + y4 + y5
  eta2 =~ y6 + y7 + y8 + y9 + y10
"

tictoc::tic("DWLS Estimation")
fit_lav1 <- cfa(mod, dat, std.lv = TRUE)
coef(fit_lav1)
tictoc::toc()
tictoc::tic("Pairwise ML Estimation")
fit_lav2 <- cfa(mod, dat, std.lv = TRUE, estimator = "PML")
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
