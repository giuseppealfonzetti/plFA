test_that("[binary] cfa function works", {

  # Generate some data
  set.seed(456)
  mod <- "
  f1 =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5
  f2 =~ 0.7*y6 + 0.7*y7 + 0.7*y8 + 0.7*y9 + 0.7*y10
  f1 ~~ 0.3*f2
  f1 ~~ 1*f1
  f2 ~~ 1*f2

  y1 | 0.5*t1
  y2 | 0.6*t1
  y3 | 0.7*t1
  y4 | 0.8*t1
  y5 | 0.9*t1
  y6 | 0.5*t1
  y7 | 0.6*t1
  y8 | 0.7*t1
  y9 | 0.8*t1
  y10 | 0.9*t1
  "
  dat <- lavaan::simulateData(mod, sample.nobs = 2000)
  dat <- as.data.frame(lapply(dat, ordered))

  # lavaan syntax fit
  mod <- "
    eta1 =~ y1 + y2 + y3 + y4 + y5
    eta2 =~ y6 + y7 + y8 + y9 + y10
  "
  fit <- cfa(mod, dat, std.lv = TRUE)

  expect_s4_class(fit, "lavaan")
  expect_s4_class(fit, "plFAlavaan")
  expect_true(inherits(coef(fit), "lavaan.vector"))
})

test_that("[ordinal] cfa function works", {

  # Generate some data
  set.seed(456)
  mod <- "
  f1 =~ 0.7*y1 + 0.7*y2 + 0.7*y3
  f2 =~ 0.7*y4 + 0.7*y5 + 0.7*y6
  f3 =~ 0.7*y7 + 0.7*y8 + 0.7*y9
  f1 ~~ 0.3*f2
  f1 ~~ 0.2*f3
  f2 ~~ 0.4*f3
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f3 ~~ 1*f3

  y1 | -0.5*t1 + 0*t2 + 0.5*t3
  y2 | -0.5*t1 + 0*t2 + 0.5*t3
  y3 | -0.5*t1 + 0*t2 + 0.5*t3
  y4 | -0.5*t1 + 0*t2 + 0.5*t3
  y5 | -0.5*t1 + 0*t2 + 0.5*t3
  y6 | -0.5*t1 + 0*t2 + 0.5*t3
  y7 | -0.5*t1 + 0*t2 + 0.5*t3
  y8 | -0.5*t1 + 0*t2 + 0.5*t3
  y9 | -0.5*t1 + 0*t2 + 0.5*t3
  "
  dat <- lavaan::simulateData(mod, sample.nobs = 2000)
  dat <- as.data.frame(lapply(dat, ordered))

  # lavaan syntax fit
  mod <- "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
    eta3 =~ y7 + y8 + y9
  "
  fit <- cfa(mod, dat, std.lv = TRUE)

  expect_s4_class(fit, "lavaan")
  expect_s4_class(fit, "plFAlavaan")
  expect_true(inherits(coef(fit), "lavaan.vector"))
})





