test_that("check fit_plFA input",{
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


  # tests for DATA
  tst <- D; tst[1,1] <- NA
  for (dt in list(NULL, NA, NaN, 3, 'hello', c(1,2,3), tst)) {
    expect_error(fit_plFA(
      DATA = dt,
      CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'=1),
      METHOD = 'ucminf'),
      'DATA is not a numeric matrix.')
  }

  # tests for CONSTRMAT
  tst <- A; tst[1,1] <- NA
  for (cm in list(NULL, NA, NaN, 3, 'hello', c(1,2,3), tst)) {
    expect_error(fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = cm, 'CORRFLAG'= 1),
      METHOD = 'ucminf'),
      'CONSTRMAT must be a binary matrix')
  }

  # tests for checking dimensions of DATA and CONSTRMAT
tst <- list( c(8,7,4), c(7,8,4), c(3,3,4))
  for (i in 1:length(tst)) {
    expect_error(fit_plFA(
      DATA = D[,1:tst[[i]][1]],
      CONSTR_LIST = list('CONSTRMAT' = A[1:tst[[i]][2],1:tst[[i]][3]], 'CORRFLAG'= 1),
      METHOD = 'ucminf'),
      'CONSTRMAT dimensions not allowed. Check Items x Factors.')
  }

# tests for CORRFLAG
  for (cf in c(NULL, NA, NaN, 3, 'hello')) {
    expect_error(fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= cf),
      METHOD = 'ucminf'),
      'CORRFLAG must be 0 or 1')
  }

# tests for METHOD
for (lab in c(NULL, NA, NaN, 3, 'hello', 'ber', 'hyp')) {
  expect_error(fit_plFA(
    DATA = D,
    CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
    METHOD = lab),
    'Method not available.')
}

# tests for INIT
expect_message(fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'ucminf',
  INIT = NULL),
  '1. Initialising at default values')

tst <- par_init; tst[1] <- NA
expect_message(fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'ucminf',
  INIT = tst),
  '1. Initialising at default values')

expect_message(fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'ucminf',
  INIT = par_init),
  '1. Initialising at INIT vector.')

expect_error(fit_plFA(
  DATA = D,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG'= 1),
  METHOD = 'ucminf',
  INIT = par_init[1:4]),
  paste0('init vector has length 4 instead of ', length(par_init) ,'.'))
})
