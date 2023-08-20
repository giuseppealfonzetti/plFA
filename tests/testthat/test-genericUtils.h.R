test_that("get_S() check implementations", {
  set.seed(123)
  p <- 10; q <- 5; n <- 50
  par <- rnorm(q*(q-1)/2)

  S <- get_S(THETA = par, Q = q)
  S2 <- get_S2(THETA = par, Q = q)
  expect_equal(S, S2)
})
