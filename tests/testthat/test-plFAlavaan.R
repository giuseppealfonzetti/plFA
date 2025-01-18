test_that("[binary] cfa function works", {

  # Generate some data
  set.seed(456)
  p <- 10; q <- 2; n <- 1000
  A <- build_constrMat(P = p, Q = q, STRUCT = "simple")
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- 0.75
  S <- get_S(THETA = rnorm(q * (q - 1) / 2), Q = q, CORRFLAG = 1)
  D <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  )
  # Needs to be ordinal data for lavaan to recognise it
  dat <- as.data.frame(lapply(as.data.frame(D), ordered))
  names(dat) <- paste0("y", seq_len(p))

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
  p <- 10; q <- 2; n <- 250
  A <- build_constrMat(P = p, Q = q, STRUCT = "simple")
  Load <- gen_loadings(CONSTRMAT = A)
  thr <- c(-1.5, 0, 1.5)
  S <- get_S(THETA = rnorm(q * (q - 1) / 2), Q = q, CORRFLAG = 1)
  D <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S
  )
  # Needs to be ordinal data for lavaan to recognise it
  dat <- as.data.frame(lapply(as.data.frame(D), ordered))
  names(dat) <- paste0("y", seq_len(p))

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





