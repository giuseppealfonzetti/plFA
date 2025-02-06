# p = number of items, q = number of latent variables, n = number of observations
p <- 16; q <- 4; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

###############
## MSE tests ##
###############
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
  ## thresholds
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE thresholds"),{
             expect_true(
               mean((numParList$thresholds-rep(thr,p))^2) <1e-2
             )
           }
  )

  ## loadings
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE loadings"),{
             expect_true(
               mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
               <1e-2
             )
           }
  )

  ## covariances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE covariances"),{
             expect_true(
               mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
             )
           }
  )

  ## variances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE sd"),{
             expect_true(
               mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
             )
           }
  )

}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
#### with linear constraints on loadings
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1,5,2)),
              list(c(3,1), c(1,6,2)))
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
  ## thresholds
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE thresholds"),{
             expect_true(
               mean((numParList$thresholds-rep(thr,p))^2) <1e-2
             )
           }
  )

  ## loadings
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE loadings"),{
             expect_true(
               mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
               <1e-2
             )
           }
  )

  ## covariances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE covariances"),{
             expect_true(
               mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
             )
           }
  )

  ## variances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE sd"),{
             expect_true(
               mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
             )
           }
  )

}


#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
if(1){
  set.seed(123)
  stdlv <- TRUE
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

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
  ## thresholds
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE thresholds"),{
             expect_true(
               mean((numParList$thresholds-rep(thr,p))^2) <1e-2
             )
           }
  )

  ## loadings
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE loadings"),{
             expect_true(
               mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
               <1e-2
             )
           }
  )

  ## covariances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE covariances"),{
             expect_true(
               mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
             )
           }
  )

  ## variances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE sd"),{
             expect_true(
               mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
             )
           }
  )

  #
  # numVar <- computeVar(OBJ=numFit, DATA = dat)
  # numVar$asymptotic_variance

}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
#### some variances specified by the user ####
if(1){
  set.seed(123)
  stdlv <- TRUE
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
  constr_var <- rep(NA, q); constr_var[1:2] <- .9
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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

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

  ## thresholds
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE thresholds"),{
             expect_true(
               mean((numParList$thresholds-rep(thr,p))^2) <1e-2
             )
           }
  )

  ## loadings
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE loadings"),{
             expect_true(
               mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
               <1e-2
             )
           }
  )

  ## covariances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE covariances"),{
             expect_true(
               mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
             )
           }
  )

  ## variances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE sd"),{
             expect_true(
               mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
             )
           }
  )

}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### no correlation matrix but all free latent variances ######
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- FALSE

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
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

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

  ## thresholds
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE thresholds"),{
             expect_true(
               mean((numParList$thresholds-rep(thr,p))^2) <1e-2
             )
           }
  )

  ## loadings
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE loadings"),{
             expect_true(
               mean((numParList$loadings[is.na(constr_list$CONSTRMAT)]-Load[is.na(constr_list$CONSTRMAT)])^2)
               <1e-2
             )
           }
  )

  ## covariances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE covariances"),{
             expect_true(
               mean((numParList$latent_correlations[lower.tri(S)]-S[lower.tri(S)])^2)<1e-2
             )
           }
  )

  ## variances
  test_that(
    paste0("CORRFLAG:", corrflag, ", STDLV:", stdlv,
           ", MSE sd"),{
             expect_true(
               mean((sqrt(diag(numParList$latent_correlations))-sqrt(diag(S)))^2)<1e-2
             )
           }
  )

  #
  # numVar <- computeVar(OBJ=numFit, DATA = dat)
  # numVar$asymptotic_variance

}
