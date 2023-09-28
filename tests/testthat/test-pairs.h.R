test_that("check pair gradient ",{
  set.seed(1)
  p <- 8L; q <- 4L; n <- 100
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1, 0, 1)
  S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q)
  D <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  cat <- apply(D, 2, max) + 1
  theta <- get_theta(rep(thr, p), Load, S, cat, A)
  f <- compute_frequencies(Y = D, C_VEC = cat)


  lambda0_init <- c()
  s <- 0

  for (i in 1:length(cat)) {
    vec <- 1:(cat[i]-1)
    vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
    lambda0_init[(s + 1):(s + cat[i] - 1)] <- vec
    s <- s + cat[i] - 1
  }
  lambda_init = rep(0, sum(A))
  transformed_rhos_init = rep(0, q*(q-1)/2)
  #get_Sigma_u2(constrMat, transformed_rhos_init)

  par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)

  pair_nll <- function(par_vec){
    pair <- compute_pair(
      A = A,
      C_VEC = cat,
      THETA = par_vec,
      CORRFLAG = 1,
      k = k,
      l = l,
      PAIRS_TABLE = f,
      SILENTFLAG = 1,
      GRADFLAG = 0
    )
    out <- pair$ll/n
    return(out)
  }
  pair_ngr <- function(par_vec){
    pair <- compute_pair(
      A = A,
      C_VEC = cat,
      THETA = par_vec,
      CORRFLAG = 1,
      k = k,
      l = l,
      PAIRS_TABLE = f,
      SILENTFLAG = 1,
      GRADFLAG = 1
    )
    out <- pair$ngradient/n
    return(out)
  }

  for (l in 2:p) {
    for (k in 1:(l-1)) {
      expect_identical(sum(abs(pair_ngr(theta) - numDeriv::grad(pair_nll, theta))>1e-4), 0L)
      expect_equal(sum(abs(pair_ngr(par_init) - numDeriv::grad(pair_nll, par_init)) >1e-4), 0L)
    }
  }


})
