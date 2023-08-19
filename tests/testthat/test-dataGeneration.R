#### test compute_frequencies() ####
test_that("compute_frequencies() output", {
  set.seed(123)
  p <- 10; q <- 2; n <- 10
  y <- matrix(as.numeric(rbinom(n*p, 1, .5)), n, p)
  f <- compute_frequencies(Y = y, C_VEC = rep(2, p))
  expect_equal(ncol(f), 4*p*(p-1)/2)
  expect_equal(nrow(f), 5)
  expect_identical(sum(f[5,]<0), 0L)
  expect_identical(sum(f[5,]>n), 0L)
  expect_identical(sum(!is.finite(f)), 0L)
})

test_that("compute_frequencies() errors", {
  set.seed(123)
  p <- 10; q <- 2; n <- 10
  y <- matrix(as.numeric(rbinom(n*p, 1, .5)), n, p)
  f <- compute_frequencies(Y = y, C_VEC = rep(2, p))
  y1 <- y; y1[3,5] <- NaN
  expect_error(compute_frequencies(Y = y1,  C_VEC = rep(2, p)), 'Y is not a numeric matrix.')
  expect_error(compute_frequencies(Y = rep(2,p),  C_VEC = rep(2, p)), 'Y is not a numeric matrix.')
  expect_error(compute_frequencies(Y = y,  C_VEC = c(rep(2, p),NA)), 'C_VEC is not a numeric vector.')
  expect_error(compute_frequencies(Y = y,  C_VEC = matrix(2, 3, 4)), 'C_VEC is not a numeric vector.')
  expect_error(compute_frequencies(Y = y,  C_VEC = rep(2, p-1)), 'Y and C_VEC dimensions do not match.')
})

#### test build_constrMat() ####
test_that("build_constrMat() output", {
  p <- 10; q <- 2
  A <- build_constrMat(P = p, Q = q, STRUCT = 'triangular')
  expect_equal(sum(A), (p*q)-(q*(q-1)/2))
  expect_equal(sum(A[lower.tri(A, diag = T)]), (p*q)-(q*(q-1)/2))
  expect_equal(sum(A[upper.tri(A, diag = F)]), 0)
  expect_equal(sum(A%in%c(0,1)), p*q)
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  expect_equal(sum(A), p)
  expect_equal(sum(rowSums(A)==1), p)
})

test_that("build_constrMat() errors", {
  p <- 10; q <- 2
  A <- build_constrMat(P = p, Q = q, STRUCT = 'triangular')
  expect_error(build_constrMat(P = 'p', Q = q, STRUCT = 'triangular'), 'P not numeric.')
  expect_error(build_constrMat(P = p, Q = 'q', STRUCT = 'triangular'), 'Q not numeric.')
  expect_error(build_constrMat(P = p, Q = q, STRUCT = 'tri'), 'STRUCT not available.')
})

#### test gen_loadings() ####
test_that("gen_loadings() output", {
  set.seed(123)
  p <- 10; q <- 2
  A <- build_constrMat(P = p, Q = q, STRUCT = 'triangular')
  Load <- gen_loadings(CONSTRMAT = A)
  expect_identical(ncol(A), ncol(Load))
  expect_identical(nrow(A), nrow(Load))
  expect_identical(A!=0, Load!=0)
  expect_identical(Load[1,1]==Load[2,1], FALSE)
  Load <- gen_loadings(FIXED = .5, CONSTRMAT = A)
  expect_identical(ncol(A), ncol(Load))
  expect_identical(nrow(A), nrow(Load))
  expect_identical(A!=0, Load!=0)
  expect_identical(Load[1,1]==Load[2,1], TRUE)
  expect_identical(Load[1,1]==.5, TRUE)
})
test_that("gen_loadings() errors", {
  p <- 10; q <- 2
  A <- build_constrMat(P = p, Q = q, STRUCT = 'triangular')
  expect_error(Load <- gen_loadings(CONSTRMAT = 2), 'CONSTRMAT must be a matrix')
  expect_error(Load <- gen_loadings(CONSTRMAT = NA), 'CONSTRMAT must be a matrix')

})

#### Test sim_data() ####
test_that("sim_data() errors", {
  set.seed(123)
  p <- 10; q <- 5; n <- 50
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1, 0, 1)
  S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q)
  expect_error(sim_data(
    SAMPLE_SIZE = 'n',
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'SAMPLE_SIZE must be a positive integer.')
  expect_error(sim_data(
    SAMPLE_SIZE = -1,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'SAMPLE_SIZE must be a positive integer.')
  expect_error(sim_data(
    SAMPLE_SIZE = NA,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'SAMPLE_SIZE must be a positive integer.')
  expect_error(sim_data(
    SAMPLE_SIZE = 3.5,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'SAMPLE_SIZE must be a positive integer.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = NA,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'LOADINGS is not a numeric matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = 3,
    THRESHOLDS = thr,
    LATENT_COV = S
  ), 'LOADINGS is not a numeric matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = NA,
    LATENT_COV = S
  ), 'THRESHOLDS is not a numeric vector.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = c(2,1),
    LATENT_COV = S
  ), 'THRESHOLDS must be sorted.')

  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = 1
  ), 'LATENT_COV is not a numeric matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = NA
  ), 'LATENT_COV is not a numeric matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = rep(1,2)
  ), 'LATENT_COV is not a numeric matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = matrix(1:9,3,3)
  ), 'LATENT_COV not symmetric.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = (matrix(1:9,3,3)+t(matrix(1:9,3,3)))
  ), 'LATENT_COV not positive definite.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = diag(3, q, q)
  ), 'LATENT_COV not a correlation matrix.')
  expect_error(sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = diag(1, 2, 2)
  ), 'LOADINGS and LATENT_COV dimensions not compatible.')
}
)

test_that("sim_data() output", {
  set.seed(123)
  p <- 10L; q <- 5L; n <- 50L
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1, 0, 1)
  S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q)
  D <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  expect_identical(ncol(D), p)
  expect_identical(nrow(D), n)
  expect_identical(round(D), D)
  expect_identical(D<=length(thr), matrix(TRUE, n, p))
  expect_identical(D>=0, matrix(TRUE, n, p))



}
)
