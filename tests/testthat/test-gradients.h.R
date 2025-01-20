###################################################
# Test file for gradients.h based on wrappers for #
# pair contributions in pairs.h as exported by    #
# exportedFuns.h                                  #
###################################################

set.seed(123)

# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 100

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p
pair_nll <- function(PAR, K, L, OPTION=0){
  pair <- cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 0,
    OPTION      = OPTION
  )
  out <- pair$nll/n
  return(out)
}

pair_ngr <- function(PAR, K, L, OPTION=0){
  pair <- cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 1,
    OPTION      = OPTION
  )
  out <- pair$ngr/n
  return(out)
}

pair_nll_old <- function(PAR, K, L){
  pair <- cpp_compute_pair(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 0,
    OPTION      = 0
  )
  out <- pair$ll/n
  return(out)
}
pair_ngr_old <- function(PAR, K, L){
  pair <- cpp_compute_pair(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 1,
    OPTION      = 0
  )
  out <- pair$ngr/n
  return(out)
}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
{
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
{
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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.5, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
#### some variances specified by the user ####
{
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
  constr_var <- rep(NA, q); constr_var[1:2] <- .5
  constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.5, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}


{
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
  constr_var <- rep(NA, q); constr_var[2] <- .5
  constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -1, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}

#### (STDLV=FALSE, CORRFLAG=FALSE) ####
#### no correlation matrix but some free latent variances ######
{
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
  constr_var <- rep(NA, q); constr_var[2] <- .5
  constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -1, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### no correlation matrix but all free latent variances ######
{
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
  constr_var <- rep(NA, q);
  constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -1, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- check_data(dat)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- check_dims(dat, constr_list)

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })

      test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
        expect_equal(pair_ngr(PAR=theta, K=k, L=l, OPTION=1), numDeriv::grad(func=pair_nll, x=theta, K=k,L=l), tolerance = 1e-4)
      })
    }
  }
}









