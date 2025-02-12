Sys.setenv("R_TESTS" = "")

# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

###########################
## Check METHOD="SA" init #
###########################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- NULL
  A <- check_cnstr_loadings(A, stdlv, LLC=llc)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)

  plFit <- fit_plFA(
    DATA = dat,
    CONSTR_LIST = constr_list,
    INIT_METHOD="standard",
    METHOD = 'SA',
    CPP_CONTROL_MAIN = list(BURNE=0, MAXE=0, PAIRS_PER_ITERATION=1, EPOCHS_PER_CHECK=1),
    VERBOSE = FALSE)

  plFit@stoFit@path$nll


  fit <- cpp_plSA2(
    FREQ = plFit@freq,
    VALFREQ = plFit@freq,
    N = plFit@dims@n,
    C_VEC = plFit@dims@cat,
    CONSTRMAT=plFit@cnstr@loadings,
    CONSTRLOGSD=plFit@cnstr@loglatsd,
    LLC = plFit@cnstr@llc,
    THETA_INIT = plFit@init,
    DIH = rep(1, length(plFit@init)),
    NTHR = plFit@dims@nthr,
    NLOAD= plFit@dims@nload,
    NCORR= plFit@dims@ncorr,
    NVAR= plFit@dims@nvar,
    PAIRS_PER_ITERATION=plFit@stoFit@control$PAIRS_PER_ITERATION,
    STEP0=plFit@stoFit@control$STEP0,
    STEP1=plFit@stoFit@control$STEP1,
    STEP2=plFit@stoFit@control$STEP2,
    STEP3=plFit@stoFit@control$STEP3,
    BURNE=plFit@stoFit@control$BURNE,
    MAXE=plFit@stoFit@control$MAXE,
    EPOCHS_PER_CHECK = plFit@stoFit@control$EPOCHS_PER_CHECK,
    TOL_END=plFit@stoFit@control$TOL_END,
    TOL_BURN = plFit@stoFit@control$TOL_BURN,
    CHECK_TOL=plFit@stoFit@control$CHECK_TOL,
    MAX_TOL_COUNTER = plFit@stoFit@control$MAX_TOL_COUNTER,
    SEED=plFit@stoFit@control$SEED,
    VERBOSE=FALSE,
    VERBOSE_ITER=FALSE
  )

  test_that("check initial nll from SA", {expect_equal(plFit@stoFit@path$nll[1],
                                                       fit$path_nll[1])})


}


###################################
## Check METHOD="SA" vs cpp_plSA ##
###################################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- NULL
  A <- check_cnstr_loadings(A, stdlv, LLC=llc)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)

  plFit <- fit_plFA(
    DATA = dat,
    CONSTR_LIST = constr_list,
    INIT_METHOD="standard",
    METHOD = 'SA',
    CPP_CONTROL_MAIN = list(BURNE=1, MAXE=3, PAIRS_PER_ITERATION=8, EPOCHS_PER_CHECK=1),
    VERBOSE = FALSE)

  plFit@stoFit@path$nll


  fit <- cpp_plSA2(
    FREQ = plFit@freq,
    VALFREQ = plFit@freq,
    N = plFit@dims@n,
    C_VEC = plFit@dims@cat,
    CONSTRMAT=plFit@cnstr@loadings,
    CONSTRLOGSD=plFit@cnstr@loglatsd,
    LLC = plFit@cnstr@llc,
    THETA_INIT = plFit@init,
    DIH = rep(1, length(plFit@init)),
    NTHR = plFit@dims@nthr,
    NLOAD= plFit@dims@nload,
    NCORR= plFit@dims@ncorr,
    NVAR= plFit@dims@nvar,
    PAIRS_PER_ITERATION=plFit@stoFit@control$PAIRS_PER_ITERATION,
    STEP0=plFit@stoFit@control$STEP0,
    STEP1=plFit@stoFit@control$STEP1,
    STEP2=plFit@stoFit@control$STEP2,
    STEP3=plFit@stoFit@control$STEP3,
    BURNE=plFit@stoFit@control$BURNE,
    MAXE=plFit@stoFit@control$MAXE,
    EPOCHS_PER_CHECK = plFit@stoFit@control$EPOCHS_PER_CHECK,
    TOL_END=plFit@stoFit@control$TOL_END,
    TOL_BURN = plFit@stoFit@control$TOL_BURN,
    CHECK_TOL=plFit@stoFit@control$CHECK_TOL,
    MAX_TOL_COUNTER = plFit@stoFit@control$MAX_TOL_COUNTER,
    SEED=plFit@stoFit@control$SEED,
    VERBOSE=FALSE,
    VERBOSE_ITER=FALSE
  )

  test_that("check nll against direct cpp call", {
    expect_true(all(plFit@stoFit@path$nll== fit$path_nll))})
}

#######################
## Check METHOD="SA" ##
#######################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- NULL
  A <- check_cnstr_loadings(A, stdlv, LLC=llc)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)

  plFit <- fit_plFA(
    DATA = dat,
    CONSTR_LIST = constr_list,
    INIT_METHOD="standard",
    METHOD = 'SA',
    CPP_CONTROL_MAIN = list(),
    VERBOSE = FALSE)

  test_that("finite nll", {expect_true(all(is.finite(plFit@stoFit@path$nll)))})
  test_that("finite par", {lapply(plFit@stoFit@path$avtheta, function(x){expect_true(all(is.finite(x)))})})
}

############################
## Check INIT_METHOD="SA" ##
############################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- NULL
  A <- check_cnstr_loadings(A, stdlv, LLC=llc)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)

  plFit <- fit_plFA(
    DATA = dat,
    CONSTR_LIST = constr_list,
    INIT_METHOD="SA",
    METHOD = 'SA',
    CPP_CONTROL_MAIN = list(),
    VERBOSE = FALSE)

  test_that("finite nll", {expect_true(all(is.finite(plFit@stoFit@path$nll)))})
  test_that("finite par", {lapply(plFit@stoFit@path$avtheta, function(x){expect_true(all(is.finite(x)))})})
}

# ###############
# ## MSE tests ##
# ###############
# #### (STDLV=FALSE, CORRFLAG=TRUE) ####
# #### free correlation matrix and latent variances ######
# if(0){
#   set.seed(123)
#   stdlv <- FALSE
#   corrflag <- TRUE
#
#   # Simple loading matrix constraints
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   llc <- NULL
#   A <- check_cnstr_loadings(A, stdlv, LLC=llc)
#   # Draw some random loadings according to constraints
#   Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
#   nload <- sum(is.na(A))
#
#   # Generate random latent correlation matrix
#   tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
#   ncorr <- if(corrflag)q*(q-1)/2 else 0
#   R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)
#
#
#   # Generate random latent variances
#   constr_var <- rep(NA, q)
#   constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
#   nvar  <- sum(is.na(constr_lsd))
#   tsdvec <- constr_lsd
#   tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
#   Dmat <- diag(exp(tsdvec),q,q)
#
#   S <- Dmat %*% R %*% Dmat
#
#   # Dimensions
#
#   d <- nthr + nload + ncorr + nvar
#   d
#
#
#   #
#   theta <- get_theta(
#     THRESHOLDS = rep(thr, p),
#     LOADINGS = Load,
#     LATENT_COV = S,
#     CAT = cat,
#     CONSTRMAT = A,
#     CONSTRVAR = exp(constr_lsd)^2,
#     CORRFLAG = corrflag,
#     STDLV = stdlv
#   )
#
#
#   dat <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)
#
#   numFit <- fit_plFA(
#     DATA = dat,
#     CONSTR_LIST = constr_list,
#     INIT_METHOD="standard",
#     METHOD = 'SA',
#     VERBOSE = TRUE)
#
#   numFit
#
#   # extract estimated parameter vector
#   numPar <- getPar(numFit, OPTION="raw")
#   # extract list of parameter estimates
#   numParList <- getPar(numFit, OPTION = 'list')
#
#   # mean square error
#   ## thresholds
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE thresholds"),{
#              expect_true(
#                mean((numParList$thresholds-rep(thr,p))^2) <1e-2
#              )
#            }
#   )
#
#   ## loadings
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE loadings"),{
#              expect_true(
#                mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
#                <1e-2
#              )
#            }
#   )
#
#   ## covariances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE covariances"),{
#              expect_true(
#                mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
#              )
#            }
#   )
#
#   ## variances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE sd"),{
#              expect_true(
#                mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
#              )
#            }
#   )
#
# }
#
# #### (STDLV=FALSE, CORRFLAG=TRUE) ####
# #### free correlation matrix and latent variances ######
# #### with linear constraints on loadings
# if(0){
#   set.seed(123)
#   stdlv <- FALSE
#   corrflag <- TRUE
#
#   # Simple loading matrix constraints
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   llc <- list(list(c(2,1), c(1,5,2)),
#               list(c(3,1), c(1,6,2)))
#   A <- check_cnstr_loadings(A, stdlv, LLC=llc)
#   # Draw some random loadings according to constraints
#   Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
#   nload <- sum(is.na(A))
#
#   # Generate random latent correlation matrix
#   tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
#   ncorr <- if(corrflag)q*(q-1)/2 else 0
#   R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)
#
#
#   # Generate random latent variances
#   constr_var <- rep(NA, q)
#   constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
#   nvar  <- sum(is.na(constr_lsd))
#   tsdvec <- constr_lsd
#   tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
#   Dmat <- diag(exp(tsdvec),q,q)
#
#   S <- Dmat %*% R %*% Dmat
#
#   # Dimensions
#
#   d <- nthr + nload + ncorr + nvar
#   d
#
#
#   #
#   theta <- get_theta(
#     THRESHOLDS = rep(thr, p),
#     LOADINGS = Load,
#     LATENT_COV = S,
#     CAT = cat,
#     CONSTRMAT = A,
#     CONSTRVAR = exp(constr_lsd)^2,
#     CORRFLAG = corrflag,
#     STDLV = stdlv
#   )
#
#
#   dat <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc)
#
#   numFit <- fit_plFA(
#     DATA = dat,
#     CONSTR_LIST = constr_list,
#     INIT_METHOD="standard",
#     METHOD = 'SA')
#
#   numFit
#
#   # extract estimated parameter vector
#   numPar <- getPar(numFit, OPTION="raw")
#   # extract list of parameter estimates
#   numParList <- getPar(numFit, OPTION = 'list')
#
#   # mean square error
#   ## thresholds
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE thresholds"),{
#              expect_true(
#                mean((numParList$thresholds-rep(thr,p))^2) <1e-2
#              )
#            }
#   )
#
#   ## loadings
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE loadings"),{
#              expect_true(
#                mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
#                <1e-2
#              )
#            }
#   )
#
#   ## covariances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE covariances"),{
#              expect_true(
#                mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
#              )
#            }
#   )
#
#   ## variances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE sd"),{
#              expect_true(
#                mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
#              )
#            }
#   )
#
# }
#
#
# #### (STDLV=TRUE, CORRFLAG=TRUE) ####
# #### free correlation matrix but all latent variances fixed ######
# if(0){
#   set.seed(123)
#   stdlv <- TRUE
#   corrflag <- TRUE
#
#   # Simple loading matrix constraints
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   A <- check_cnstr_loadings(A, stdlv)
#
#   # Draw some random loadings according to constraints
#   Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
#   nload <- sum(is.na(A))
#
#   # Generate random latent correlation matrix
#   tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
#   ncorr <- if(corrflag)q*(q-1)/2 else 0
#   R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)
#
#
#   # Generate random latent variances
#   constr_var <- rep(NA, q)
#   constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
#   nvar  <- sum(is.na(constr_lsd))
#   tsdvec <- constr_lsd
#   tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
#   Dmat <- diag(exp(tsdvec),q,q)
#
#   S <- Dmat %*% R %*% Dmat
#
#   # Dimensions
#
#   d <- nthr + nload + ncorr + nvar
#   d
#
#
#   #
#   theta <- get_theta(
#     THRESHOLDS = rep(thr, p),
#     LOADINGS = Load,
#     LATENT_COV = S,
#     CAT = cat,
#     CONSTRMAT = A,
#     CONSTRVAR = exp(constr_lsd)^2,
#     CORRFLAG = corrflag,
#     STDLV = stdlv
#   )
#
#
#   dat <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)
#
#   numFit <- fit_plFA(
#     DATA = dat,
#     CONSTR_LIST = constr_list,
#     INIT_METHOD='ucminf',
#     METHOD = 'SA')
#
#   numFit
#
#   # extract estimated parameter vector
#   numPar <- getPar(numFit, OPTION="raw")
#   # extract list of parameter estimates
#   numParList <- getPar(numFit, OPTION = 'list')
#
#   # mean square error
#   ## thresholds
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE thresholds"),{
#              expect_true(
#                mean((numParList$thresholds-rep(thr,p))^2) <1e-2
#              )
#            }
#   )
#
#   ## loadings
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE loadings"),{
#              expect_true(
#                mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
#                <1e-2
#              )
#            }
#   )
#
#   ## covariances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE covariances"),{
#              expect_true(
#                mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
#              )
#            }
#   )
#
#   ## variances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE sd"),{
#              expect_true(
#                mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
#              )
#            }
#   )
#
#   #
#   # numVar <- computeVar(OBJ=numFit, DATA = dat)
#   # numVar$asymptotic_variance
#
# }
#
# #### (STDLV=TRUE, CORRFLAG=TRUE) ####
# #### free correlation matrix but all latent variances fixed ######
# #### some variances specified by the user ####
# if(0){
#   set.seed(123)
#   stdlv <- TRUE
#   corrflag <- TRUE
#
#   # Simple loading matrix constraints
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   A <- check_cnstr_loadings(A, stdlv)
#
#   # Draw some random loadings according to constraints
#   Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
#   nload <- sum(is.na(A))
#
#   # Generate random latent correlation matrix
#   tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
#   ncorr <- if(corrflag)q*(q-1)/2 else 0
#   R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)
#
#
#   # Generate random latent variances
#   constr_var <- rep(NA, q); constr_var[1:2] <- .9
#   constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
#   nvar  <- sum(is.na(constr_lsd))
#   tsdvec <- constr_lsd
#   tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
#   Dmat <- diag(exp(tsdvec),q,q)
#
#   S <- Dmat %*% R %*% Dmat
#
#   # Dimensions
#
#   d <- nthr + nload + ncorr + nvar
#   d
#
#
#   #
#   theta <- get_theta(
#     THRESHOLDS = rep(thr, p),
#     LOADINGS = Load,
#     LATENT_COV = S,
#     CAT = cat,
#     CONSTRMAT = A,
#     CONSTRVAR = exp(constr_lsd)^2,
#     CORRFLAG = corrflag,
#     STDLV = stdlv
#   )
#
#
#   dat <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)
#
#   numFit <- fit_plFA(
#     DATA = dat,
#     CONSTR_LIST = constr_list,
#     METHOD = 'ucminf')
#
#   numFit
#
#   # extract estimated parameter vector
#   numPar <- getPar(numFit, OPTION="raw")
#   # extract list of parameter estimates
#   numParList <- getPar(numFit, OPTION = 'list')
#
#   # mean square error
#
#   ## thresholds
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE thresholds"),{
#              expect_true(
#                mean((numParList$thresholds-rep(thr,p))^2) <1e-2
#              )
#            }
#   )
#
#   ## loadings
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE loadings"),{
#              expect_true(
#                mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
#                <1e-2
#              )
#            }
#   )
#
#   ## covariances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE covariances"),{
#              expect_true(
#                mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
#              )
#            }
#   )
#
#   ## variances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE sd"),{
#              expect_true(
#                mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
#              )
#            }
#   )
#
# }
#
# #### (STDLV=FALSE, CORRFLAG=TRUE) ####
# #### no correlation matrix but all free latent variances ######
# if(0){
#   set.seed(123)
#   stdlv <- FALSE
#   corrflag <- FALSE
#
#   # Simple loading matrix constraints
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   A <- check_cnstr_loadings(A, stdlv)
#
#   # Draw some random loadings according to constraints
#   Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
#   nload <- sum(is.na(A))
#
#   # Generate random latent correlation matrix
#   tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
#   ncorr <- if(corrflag)q*(q-1)/2 else 0
#   R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)
#
#
#   # Generate random latent variances
#   constr_var <- rep(NA, q)
#   constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
#   nvar  <- sum(is.na(constr_lsd))
#   tsdvec <- constr_lsd
#   tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
#   Dmat <- diag(exp(tsdvec),q,q)
#
#   S <- Dmat %*% R %*% Dmat
#
#   # Dimensions
#
#   d <- nthr + nload + ncorr + nvar
#   d
#
#
#   #
#   theta <- get_theta(
#     THRESHOLDS = rep(thr, p),
#     LOADINGS = Load,
#     LATENT_COV = S,
#     CAT = cat,
#     CONSTRMAT = A,
#     CONSTRVAR = exp(constr_lsd)^2,
#     CORRFLAG = corrflag,
#     STDLV = stdlv
#   )
#
#
#   dat <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)
#
#   numFit <- fit_plFA(
#     DATA = dat,
#     CONSTR_LIST = constr_list,
#     METHOD = 'ucminf')
#
#   numFit
#
#   # extract estimated parameter vector
#   numPar <- getPar(numFit, OPTION="raw")
#   # extract list of parameter estimates
#   numParList <- getPar(numFit, OPTION = 'list')
#
#   # mean square error
#
#   ## thresholds
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE thresholds"),{
#              expect_true(
#                mean((numParList$thresholds-rep(thr,p))^2) <1e-2
#              )
#            }
#   )
#
#   ## loadings
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE loadings"),{
#              expect_true(
#                mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
#                <1e-2
#              )
#            }
#   )
#
#   ## covariances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE covariances"),{
#              expect_true(
#                mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
#              )
#            }
#   )
#
#   ## variances
#   test_that(
#     paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
#            ", MSE sd"),{
#              expect_true(
#                mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
#              )
#            }
#   )
#
#   #
#   # numVar <- computeVar(OBJ=numFit, DATA = dat)
#   # numVar$asymptotic_variance
#
# }
