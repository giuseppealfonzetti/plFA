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

#### test build_constMat() ####
