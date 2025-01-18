is_ordinal_df <- function(D) {
  if (!is.data.frame(D)) return(FALSE)
  else all(sapply(D, is.ordered))
}

#' Fit Confirmatory Factor Analysis Models using plFA
#'
#' @param model A description of the user-specified model. Typically, the model
#'   is described using the lavaan model syntax. See [`lavaan::model.syntax`]
#'   for more information. Alternatively, a parameter table (eg. the output of
#'   the [`lavaan::lavaanify()`] function) is also accepted.
#' @param data A data frame containing the observed variables used in the model.
#'   Variables must be declared as ordered factors.
#' @param start Starting values to use.
#' @param estimator (For `lavaan` compatibility only). The estimator is `PML`.
#' @param estimator.args A list of arguments for [`plFA()`]. TBC.
#' @param information (For `lavaan` compatibility only). The information matrix is...
#' @param control A list of control parameters for the estimation algorithm. See [`plFA()`] for more information.
#' @param ... Additional arguments to be passed to [lavaan()].
#'
#' @returns A `plFAlavaan` object.
#' @export
cfa <- function(
    model,
    data,
    std.lv = TRUE,
    start = NULL,
    estimator = "PML",
    estimator.args = list(
      method = "ucminf",
      valdata = NULL,
      iterations_subset = NULL,
      ncores = 1
    ),
    information = "observed",
    control = list(),
    ...
) {

  lavargs <- list(...)

  # Validate arguments
  D <- sapply(data, as.numeric) - 1

  if (!is_ordinal_df(data)) {
    cli::cli_abort("Ordinal data required")
  }
  if (isFALSE(std.lv)) {
    cli::cli_alert_warning("Non-standardised latent variable estimation not implemented yet. Setting `std.lv = TRUE`.")
    std.lv <- TRUE
  }
  lavargs$std.lv <- std.lv
  lavargs$control <- control
  method <- estimator.args$method
  if (is.null(method)) {
    cli::cli_alert_warning("No estimation method specified. Setting `method = 'ucminf'`.")
    method <- "ucminf"
  }
  if (method != "ucminf") {
    cli::cli_alert_warning("Only 'ucminf' method implemented for now.")
    method <- "ucminf"
  }
  valdata <- estimator.args$valdata
  iterations_subset <- estimator.args$iterations_subset
  ncores <- estimator.args$ncores
  if ("verbose" %in% names(lavargs)) verbose <- lavargs$verbose
  else verbose <- FALSE

  # Options for computeVar()
  computevar_numderiv <- estimator.args$numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE
  computevar_option <- estimator.args$numderiv
  if (is.null(computevar_option)) computevar_option <- "transformed"

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE
  fit0 <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)

  # Fit plFA -------------------------------------------------------------------

  # Build A matrix
  FREE <- lavaan::inspect(fit0, what = "free")
  lambda <- lavaan::inspect(fit0, what = "est")$lambda
  p <- nrow(lambda)
  q <- ncol(lambda)
  # A <- build_constrMat(P = p, Q = q, STRUCT = "simple")
  A <- lambda
  class(A) <- "matrix"
  A[] <- NA
  A[FREE$lambda == 0] <- lambda[FREE$lambda == 0]

  # CORRFLAG
  psi <- FREE$psi
  LTpsi <- psi[lower.tri(psi)]
  if (length(LTpsi) > 0 & any(LTpsi > 0)) {
    corrflag <- 1
  } else {
    corrflag <- 0
  }

  fit1 <- fit_plFA(
    DATA = D,
    CONSTR_LIST = list(CONSTRMAT = A, CORRFLAG = corrflag),
    VALDATA = valdata,
    METHOD = method,
    CONTROL = control,
    INIT = start,
    ITERATIONS_SUBSET = iterations_subset,
    VERBOSEFLAG = as.numeric(verbose),
    NCORES = ncores
  )
  vars <- computeVar(
    OBJ = fit1,
    DATA = D,
    NUMDERIV = computevar_numderiv,
    OPTION = computevar_option
  )

  # list(fit0 = fit0, fit1 = fit1, vars = vars)
  out <- create_lav_from_fitplFA(fit0, fit1, vars, D)
  new("plFAlavaan", out)
}

create_lav_from_fitplFA <- function(fit0, fit1, vars, D) {

  # Get coefficients and standard errors
  parlist <- extract_par(
    THETA = fit1@theta,
    OPTION = "list",
    C = sum(fit1@dims@cat),
    P = fit1@dims@p,
    Q = fit1@dims@q,
    CONSTRMAT = fit1@cnstr@loadings,
    CORRFLAG = fit1@cnstr@corrflag
  )
  FREE <- lavaan::inspect(fit0, what = "free")
  lambda <- parlist$loadings[FREE$lambda > 0]
  tau <- parlist$thresholds[FREE$tau > 0]
  psi <- parlist$latent_correlations[FREE$psi > 0 & lower.tri(FREE$psi)]
  x <- c(lambda, tau, psi)

  SE <- sqrt(vars$asymptotic_variance)
  vcov <- vars$vcov

  # Change version slot
  fit0@version <- as.character(packageVersion("plFA"))

  # Change timing slot
  fit0@timing$optim <- fit0@timing$optim + fit1@RTime
  fit0@timing$total <- fit0@timing$total + fit1@RTime

  # Change Model and implied slots
  fit0@Model <- lavaan::lav_model_set_parameters(fit0@Model, x)
  fit0@implied <- lavaan::lav_model_implied(fit0@Model)

  # Change ParTable and pta slots
  pt <- lavaan::partable(fit0)
  pt$est[pt$free > 0] <- x
  pt$se <- 0
  pt$se[pt$free > 0] <- SE
  fit0@ParTable <- as.list(pt)
  fit0@pta$names <- names(pt)

  # Change Options slot
  fit0@Options$estimator <- "PML"
  fit0@Options$optim.method <- fit1@method
  # fit0@Options$estimator.args <- list()
  # fit0@Options$test <- "standard"
  fit0@Options$se <- "robust.huber.white"  # this is the sandwich
  fit0@Options$do.fit <- TRUE

  # Change Fit slot
  fit0@Fit@x <- x
  fit0@Fit@se <- SE
  fit0@Fit@iterations <- as.integer(fit1@numFit$info["neval"])
  fit0@Fit@converged <- fit1@numFit$convergence == 1L  # FIXME: Check!!

  # Change optim slot
  fit0@optim$x <- x
  # fit0@optim$dx <- 0
  fit0@optim$npar <- length(x)
  fit0@optim$fx <- fit0@Fit@fx
  fit0@optim$fx.group <- fit0@Fit@fx.group
  fit0@optim$iterations <- fit1@numFit$info["neval"]
  fit0@optim$converged <- fit1@numFit$convergence == 1L

  # Change loglik slot
  # fit0@loglik$estimator <-
  #   if (fit$estimator == "ML") "ML"
  #   else if (fit$estimator == "IBRM") "IMP-BR ML"
  #   else if (fit$estimator == "IBRMP") "IMP-BR ML"
  #   else if (fit$estimator == "EBRM") "EXP-BR ML"
  # Change vcov slot
  fit0@vcov$se <- "robust.sem"
  fit0@vcov$information <- "expected"
  fit0@vcov$vcov <- vcov

  # fit0@test <- fit_lav@test
  # fit0@baseline <- fit_lav@baseline

  # Include the entire output of fit_sem
  fit0@external <- list(plFA = fit1, D = D)

  fit0
}
