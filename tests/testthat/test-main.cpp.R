test_that("check full gradient ",{
  set.seed(1)
  p <- 8L; q <- 4L; n <- 100
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1, 0, 1)
  S <- get_S2(THETA = rnorm(q*(q-1)/2), Q = q)
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
    lambda0_ <- par_vec[1:length(lambda0_init)]
    lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
    mod <-multiThread_completePairwise(
      Y = D,
      C_VEC = cat,
      A = A,
      FREQ = f,
      TAU = lambda0_,
      LAMBDA = lambda_,
      TRANSFORMED_RHOS = transformed_rhos_,
      CORRFLAG = 1,
      GRFLAG = 0,
      SILENTFLAG = 1
    )
    out <- mod$iter_nll/n
    return(out)
  }


  # function for gradient
 full_ngr <- function(par_vec){
    lambda0_ <- par_vec[1:length(lambda0_init)]
    lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
    mod <-multiThread_completePairwise(
      Y = D,
      C_VEC = cat,
      A = A,
      FREQ = f,
      TAU = lambda0_,
      LAMBDA = lambda_,
      TRANSFORMED_RHOS = transformed_rhos_,
      CORRFLAG = 1,
      GRFLAG = 1,
      SILENTFLAG = 1
    )

    out <- mod$iter_ngradient/n
    return(out)
  }

 expect_identical(sum(abs(numDeriv::grad(full_nll, theta) - full_ngr(theta)) > 1e-4), 0L)
 expect_identical(sum(abs(numDeriv::grad(full_nll, par_init) - full_ngr(par_init)) > 1e-4), 0L)


})
