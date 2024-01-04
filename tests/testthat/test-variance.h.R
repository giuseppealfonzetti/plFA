test_that("check gradient from H",{
  set.seed(1)
  p <- 6L; q <- 2L; n <- 500
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

  # function for nll
  full_nll <- function(par_vec){
    mod <- multiThread_completePairwise(
      N = n,
      C_VEC = cat,
      CONSTRMAT = A,
      FREQ = f,
      THETA = par_vec,
      CORRFLAG = 1,
      SILENTFLAG = 1
    )
    out <- mod$iter_nll/n
    return(out)
  }


  # function for gradient
  full_ngr <- function(par_vec){
    mod <- multiThread_completePairwise(
      N = n,
      C_VEC = cat,
      CONSTRMAT = A,
      FREQ = f,
      THETA = par_vec,
      CORRFLAG = 1,
      SILENTFLAG = 1
    )

    out <- mod$iter_ngradient/n
    return(out)
  }

  Hgr <- function(par_vec){
    out <- estimate_H(
      C_VEC = cat,
      A = A,
      THETA = par_vec,
      FREQ = f,
      N = n,
      CORRFLAG = 1
    )$gradient
    out <- -out/n
    return(out)
  }

  H <- function(par_vec){
    out <- estimate_H(
      C_VEC = cat,
      A = A,
      THETA = par_vec,
      FREQ = f,
      N = n,
      CORRFLAG = 1
    )$est_H
    return(out)
  }
  expect_equal(full_ngr(theta), Hgr(theta))
  expect_equal(full_ngr(par_init), Hgr(par_init))
  # expect_equal(sum(abs(H(theta)+numDeriv::jacobian(Hgr, theta))>1e-3), 0)

})
