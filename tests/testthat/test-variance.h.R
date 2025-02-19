# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

full_ngr <- function(par_vec){
  n <- dims$n
  mod <- cpp_multiThread_completePairwise(
    N          = dims$n,
    C_VEC      = dims$cat,
    CONSTRMAT  = constr_list$CONSTRMAT,
    CONSTRLOGSD= constr_list$CONSTRLOGSD,
    LLC        = constr_list$LLC,
    FREQ       = f,
    THETA      = par_vec,
    CORRFLAG   = constr_list$CORRFLAG,
    NTHR       = dims$nthr,
    NLOAD      = dims$nload,
    NCORR      = dims$ncorr,
    NVAR       = dims$nvar,
    GRFLAG     = 1,
    SILENTFLAG = 1
  )
  out <- mod$iter_ngr/dims$n
  return(out)
}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- NULL
  A <- check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)


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
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- check_dims(dat, constr_list)

  RcppParallel::setThreadOptions(numThreads = 1)
  HJ <- cpp_sample_estimators_HJ(
    THETA = theta,
    FREQ = f,
    DATA = dat,
    C_VEC = cat,
    CONSTRMAT = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    LLC = constr_list$LLC,
    N = dims$n,
    CORRFLAG = constr_list$CORRFLAG,
    NTHR = dims$nthr,
    NLOAD = dims$nload,
    NCORR= dims$ncorr,
    NVAR = dims$nvar
  )

  test_that("check total gradient", {
    expect_equal(
      -HJ$gradH/n,
      full_ngr(theta)
    )

    expect_equal(
      -HJ$gradJ/n,
      full_ngr(theta)
    )
  })

}
