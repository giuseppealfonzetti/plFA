#### Test get_corr() ####
test_that("get_corr() input type", {
  expect_error(get_corr(THETA = rep(0, 10), Q = '10'), 'Q not numeric.')
  expect_error(get_corr(THETA = rep('0', 10), Q = 10), 'THETA not numeric.')
  expect_error(get_corr(THETA = c(rep(0, 10), NaN), Q = 10), 'THETA not numeric.')
  expect_error(get_corr(THETA = c(rep(0, 10), NA), Q = 10), 'THETA not numeric.')
  expect_error(get_corr(THETA = c(rep(0, 10), Inf), Q = 10), 'THETA not numeric.')
  expect_error(get_corr(THETA = rep(0, 10), Q = NA), 'Q not numeric.')
})

test_that("get_corr() input dimensions", {
  set.seed(123)
  expect_error(get_corr(THETA = rep(0, 10), Q = 10), 'Check THETA dimensions')
})

test_that("Extracted correlation vector length", {
  set.seed(123)
  q <- 10
  expect_equal(length(get_corr(THETA = rnorm(q*(q-1)/2), Q = q)), q*(q-1)/2)
})

test_that("Extracted correlation range", {
  set.seed(123)
  q <- 10; corr <- get_corr(THETA = rnorm(q*(q-1)/2), Q = q)
  expect_equal(sum(corr > -1 & corr < 1), q*(q-1)/2)
})

test_that("Extracted correlation positive definiteness", {
  set.seed(123)
  q <- 10; corr <- get_corr(THETA = rnorm(q*(q-1)/2), Q = q)
  R <- diag(1, q, q); R[upper.tri(R)] <- corr; R <- R + t(R)
  expect_identical(matrixcalc::is.positive.definite(R), T)
})




#### Test get_lambda() ####

test_that("get_lambda() input type", {
  expect_error(get_lambda(THETA = rep(0, 10), C = '30', P = 10, Q = 3), 'C not numeric.')
  expect_error(get_lambda(THETA = rep(0, 10), C = 30, P = NA, Q = 3), 'P not numeric.')
  expect_error(get_lambda(THETA = rep(0, 10), C = 30, P = 10, Q = Inf), 'Q not numeric.')
  expect_error(get_lambda(THETA = c(rep(0, 10), NaN), C = 30, P = 10, Q = 3), 'THETA not numeric.')
  expect_error(get_lambda(THETA = c(rep(0, 10), NA), C = 30, P = 10, Q = 3), 'THETA not numeric.')
})

test_that("get_lambda() input dimensions", {
  expect_error(get_lambda(THETA = rep(0, 10), C = 30, P = 10, Q = 3), 'Check THETA dimensions')
})

#### Test get_theta() ####
test_that("get_theta() input type", {
  set.seed(123)
  p <- 10; q <- 5; cat <- rep(2,p)
  tau <- rep(0, sum(cat)-p)
  A <- matrix(rbinom(p*q, 1, 1/q), p, q)
  load <- A*matrix(rnorm(p*q), p, q)
  latent_cov <- diag(1, q, q)

  expect_error(get_theta(TAU = rep(NA, sum(cat)-p),
                         LOADINGS = load, LATENT_COV = latent_cov,
                         CAT = cat, CONSTRMAT = A), 'TAU not numeric.')
  expect_error(get_theta(TAU = tau,
                         LOADINGS = NA, LATENT_COV = latent_cov,
                         CAT = cat, CONSTRMAT = A), 'LOADINGS not numeric.')
  expect_error(get_theta(TAU = tau,
                         LOADINGS = load, LATENT_COV = Inf,
                         CAT = cat, CONSTRMAT = A), 'LATENT_COV not numeric.')
  expect_error(get_theta(TAU = tau,
                         LOADINGS = load, LATENT_COV = latent_cov,
                         CAT = NaN, CONSTRMAT = A), 'CAT not numeric.')

  theta <- get_theta(TAU = tau,
                     LOADINGS = load, LATENT_COV = latent_cov,
                     CAT = cat, CONSTRMAT = A)

})

test_that("get_theta() basic output checks", {
  set.seed(123)
  p <- 10; q <- 5; cat <- rep(2,p)
  tau <- rep(0, sum(cat)-p)
  A <- matrix(rbinom(p*q, 1, 1/q), p, q)
  load <- A*matrix(rnorm(p*q), p, q)
  latent_cov <- diag(1, q, q)

  theta <- get_theta(TAU = tau,
                     LOADINGS = load, LATENT_COV = latent_cov,
                     CAT = cat, CONSTRMAT = A)
  expect_identical(sum(theta[1:(sum(cat)-p)]-tau), 0)
  expect_identical(sum(!(theta[(sum(cat)-p+1):(q*(q-1)/2)]%in%load)), 0L)
})

