test_that("check pair gradient ",{
  set.seed(1)
  p <- 8L; q <- 4L; n <- 100
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1, 0, 1)
  tcorrvec <- rnorm(q*(q-1)/2)
  S <- cpp_get_latvar_vec2mat(SVEC=tcorrvec, Q=q)
  D <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  cat <- apply(D, 2, max) + 1
  corrflag <- 1
  theta <- get_theta(rep(thr, p), Load, S, cat, A, CORRFLAG = corrflag)
  cat(length(par_init))
  f <- compute_frequencies(Y = D, C_VEC = cat)


  D <- check_data(D)
  constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag))
  dims <- check_dims(D, constr_list)

  lambda0_init <- init_thresholds(dims, constr_list)
  lambda_init <- init_loadings(dims, constr_list)
  transformed_rhos_init = init_transformed_latcorr(dims, constr_list)

  #get_Sigma_u2(constrMat, transformed_rhos_init)

  par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)
  cat(length(par_init))

  pair_nll <- function(par_vec, OPTION=0){
    pair <- cpp_compute_pair(
      A = A,
      C_VEC = cat,
      THETA = par_vec,
      CORRFLAG = 1,
      k = k,
      l = l,
      PAIRS_TABLE = f,
      SILENTFLAG = 1,
      GRADFLAG = 0,
      OPTION = OPTION
    )
    out <- pair$ll/n
    return(out)
  }

  pair_nll(par_init)
  pair_nll(theta)


  pair_ngr <- function(par_vec, OPTION=0){
    pair <- cpp_compute_pair(
      A = A,
      C_VEC = cat,
      THETA = par_vec,
      CORRFLAG = 1,
      k = k,
      l = l,
      PAIRS_TABLE = f,
      SILENTFLAG = 1,
      GRADFLAG = 1,
      OPTION = OPTION
    )
    out <- pair$ngradient/n
    return(out)
  }

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      expect_identical(sum(abs(pair_ngr(theta, OPTION=0) - numDeriv::grad(pair_nll, theta))>1e-4), 0L)
      expect_equal(sum(abs(pair_ngr(par_init, OPTION=0) - numDeriv::grad(pair_nll, par_init)) >1e-4), 0L)
    }
  }

  # for (l in 2:p) {
  #   for (k in 1:(l-1)) {
  #     expect_identical(sum(abs(pair_ngr(theta, OPTION=1) - numDeriv::grad(pair_nll, theta))>1e-4), 0L)
  #     expect_equal(sum(abs(pair_ngr(par_init, OPTION=1) - numDeriv::grad(pair_nll, par_init)) >1e-4), 0L)
  #   }
  # }


})

# l <- 2
# k <- 1
# microbenchmark::microbenchmark(
#   sample=pair_ngr(theta),
#   ntimespi=pair_ngr(theta, OPTION=1)
# )
