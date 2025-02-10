testthat::skip()

mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
fit1a <- plFA::cfa(mod, LSAT, std.lv = TRUE)
fit1b <- plFA::cfa(mod, LSAT, std.lv = TRUE, verbose = TRUE,
                   estimator.args = list(computevar_numderiv = TRUE))
fit2 <- lavaan::cfa(mod, LSAT, std.lv = TRUE, estimator = "PML")

# Compare coef
expect_equal(coef(fit1a), coef(fit2), tolerance = 1e-5)

# Compare optim value
expect_equal(fit1a@optim$fx, fit2@optim$fx)
# Not the same. I know in lavaan the pl function is divided by sample size for
# numerical reasons. But still very different
# G: I guess we differ only for some constant terms (non-parameter dependent). So it does not affect optimisation

# Compare Hinv matrix
Hinv1a <- with(fit1a@external, computeVar(plFA, D)$invH[idx_plFA2lav, idx_plFA2lav])
Hinv1b <- with(fit1b@external, computeVar(plFA, D, NUMDERIV = TRUE)$invH[idx_plFA2lav, idx_plFA2lav])
Hinv2 <- lavaan:::lav_model_information_observed(
  lavmodel =fit2@Model,
  lavsamplestats = fit2@SampleStats,
  lavdata = fit2@Data,
  lavoptions = fit2@Options,
  lavcache = fit2@Cache,
  inverted = TRUE
)
expect_equal(Hinv1a, Hinv1b, tolerance = 1e-5)  # ucminf vs numderiv
expect_equal(Hinv1a, Hinv2, tolerance = 1e-5)  # ucminf vs lavaan
expect_equal(Hinv1b, Hinv2, tolerance = 1e-5)  # numderiv vs lavaan

# Compare J matrix
J1 <- with(fit1a@external, computeVar(plFA, D)$J[idx_plFA2lav, idx_plFA2lav])
J2 <- lavaan:::lav_model_information_firstorder(
  lavmodel =fit2@Model,
  lavsamplestats = fit2@SampleStats,
  lavdata = fit2@Data,
  lavoptions = fit2@Options,
  lavcache = fit2@Cache
)
expect_equal(J1, J2, tolerance = 1e-5)
