# test_that("check stoc_args() cleaning",{
#   p <- 50
#   args <- list()
#   ctrl <- stoc_args(ARGS = args, P = p)
#   expect_identical(ctrl,
#                    list(
#                      PAIRS_PER_ITERATION = 8,
#                      MAXT = round(50*p*(p-1)/ctrl$PAIRS_PER_ITERATION, 0),
#                      BURN = round(ctrl$MAXT/2, 0),
#                      ETA = .005,
#                      PAR1 = 1,
#                      PAR2 = .002,
#                      PAR3 = .75,
#                      STEPSIZEFLAG = 1,
#                      EACHCLOCK = round(ctrl$MAXT/10, 0),
#                      CHECKCONV = 0)
#
#                    )
# })
#
# test_that("check fit_plFA() input",{
#   set.seed(1)
#   p <- 8L; q <- 4L; n <- 100
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   Load <- gen_loadings(CONSTRMAT = A)
#   thr <- c(-1, 0, 1)
#   S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q)
#   D <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   cat <- apply(D, 2, max) + 1
#   theta <- get_theta(rep(thr, p), Load, S, cat, A)
#   f <- compute_frequencies(Y = D, C_VEC = cat)
#
#
#   lambda0_init <- c()
#   s <- 0
#
#   for (i in 1:length(cat)) {
#     vec <- 1:(cat[i]-1)
#     vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
#     lambda0_init[(s + 1):(s + cat[i] - 1)] <- vec
#     s <- s + cat[i] - 1
#   }
#   lambda_init = rep(0, sum(A))
#   transformed_rhos_init = rep(0, q*(q-1)/2)
#   #get_Sigma_u2(constrMat, transformed_rhos_init)
#
#   par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)
#
#
#   # tests for DATA
#   tst <- D; tst[1,1] <- NA
#   for (dt in list(NULL, NA, NaN, 3, 'hello', c(1,2,3), tst)) {
#     expect_error(fit_plFA(
#       DATA = dt,
#       CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'=1),
#       METHOD = 'ucminf'),
#       'DATA is not a numeric matrix.')
#   }
#
#   # tests for CONSTRMAT
#   tst <- A; tst[1,1] <- NA
#   for (cm in list(NULL, NA, NaN, 3, 'hello', c(1,2,3), tst)) {
#     expect_error(fit_plFA(
#       DATA = D,
#       CONSTR_LIST = list('CONSTRMAT' = cm, 'CORRFLAG'= 1),
#       METHOD = 'ucminf'),
#       'CONSTRMAT must be a binary matrix')
#   }
#
#   # tests for checking dimensions of DATA and CONSTRMAT
# tst <- list( c(8,7,4), c(7,8,4), c(3,3,4))
#   for (i in 1:length(tst)) {
#     expect_error(fit_plFA(
#       DATA = D[,1:tst[[i]][1]],
#       CONSTR_LIST = list('CONSTRMAT' = A[1:tst[[i]][2],1:tst[[i]][3]], 'CORRFLAG'= 1),
#       METHOD = 'ucminf'),
#       'CONSTRMAT dimensions not allowed. Check Items x Factors.')
#   }
#
# # tests for CORRFLAG
#   for (cf in c(NULL, NA, NaN, 3, 'hello')) {
#     expect_error(fit_plFA(
#       DATA = D,
#       CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= cf),
#       METHOD = 'ucminf'),
#       'CORRFLAG must be 0 or 1')
#   }
#
# # tests for METHOD
# for (lab in c(NULL, NA, NaN, 3, 'hello', 'ber', 'hyp')) {
#   expect_error(fit_plFA(
#     DATA = D,
#     CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#     METHOD = lab),
#     'Method not available.')
# }
#
# # test for NULL INIT
# expect_message(fit_plFA(
#   DATA = D,
#   CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#   METHOD = 'ucminf',
#   INIT = NULL),
#   '1. Initialising at default values')
#
# # test for invalid INIT
# tst <- par_init; tst[1] <- NA
# expect_message(fit_plFA(
#   DATA = D,
#   CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#   METHOD = 'ucminf',
#   INIT = tst),
#   '1. Initialising at default values')
#
# # test valid init
# expect_message(fit_plFA(
#   DATA = D,
#   CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#   METHOD = 'ucminf',
#   INIT = par_init),
#   '1. Initialising at INIT vector.')
#
# # test invalid INIT
# expect_error(fit_plFA(
#   DATA = D,
#   CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#   METHOD = 'ucminf',
#   INIT = par_init[1:4]),
#   paste0('init vector has length 4 instead of ', length(par_init) ,'.'))
#
# # test for METHOD
# for (lab in c('ucminf', 'bernoulli', 'hyper')) {
#   expect_message(fit_plFA(
#     DATA = D,
#     CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
#     METHOD = lab,
#     INIT = NULL),
#     paste0('3. Optimising with ', lab, '...'))
# }
# })

# p = number of items, q = number of latent variables, n = number of observations
p <- 6; q <- 2; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(0){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- check_cnstr_loadings(A, stdlv)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q)
  constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.5, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions

  d <- nthr + nload + ncorr + nvar
  d


  #
  theta <- get_theta(
    THRESHOLDS = rep(thr, p),
    LOADINGS = Load,
    LATENT_COV = S,
    CAT = cat,
    CONSTRMAT = A,
    CONSTRVAR = exp(constr_lsd)^2,
    CORRFLAG = corrflag,
    STDLV = stdlv
  )


  dat <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))

  numFit <- fit_plFA(
    DATA = dat,
    CONSTR_LIST = constr_list,
    METHOD = 'ucminf')

  numFit

  # extract estimated parameter vector
  numPar <- getPar(numFit, OPTION="raw")
  # extract list of parameter estimates
  numParList <- getPar(numFit, OPTION = 'list')

  # mean square error
  mean((numPar-theta)^2)

  modcov <- numParList$loadings%*%numParList$latent_correlations%*%t(numParList$loadings)
  abs(modcov[upper.tri(modcov)])>1
  #
  numVar <- computeVar(OBJ=numFit, DATA = dat)
  numVar$asymptotic_variance

}


