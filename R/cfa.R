is_ordinal_df <- function(D) {
  if (!is.data.frame(D)) return(FALSE)
  else all(sapply(D, is.ordered))
}

cfa <- function(
    model,
    data,
    start = NULL,
    estimator = "SPML",
    estimator.args = list(
      method = "hyper",
      valdata = NULL,
      iterations_subset = NULL,
      ncores = 1
    ),
    information = "observed",
    std.lv = TRUE,
    control = list(),
    ...
) {

  lavargs <- list(...)

  # Validate arguments
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
    cli::cli_alert_warning("No estimation method specified. Setting `method = 'hyper'.")
    method <- "hyper"
  }
  valdata <- estimator.args$valdata
  iterations_subset <- estimator.args$iterations_subset
  ncores <- estimator.args$ncores
  if ("verbose" %in% names(lavargs)) verbose <- lavargs$verbose
  else verbose <- FALSE

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE
  fit0 <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)

  # Fit plFA -------------------------------------------------------------------
  D <- as.matrix(DAT) - 1

  # Get A matrix
  FREE <- lavaan::inspect(fit0, what = "free")
  lambda <- FREE$lambda  # Loading matrix
  p <- nrow(lambda)
  q <- ncol(lambda)
  if (isTRUE(lavargs$std.lv)) {
    A <- build_constrMat(P = p, Q = q, STRUCT = "simple")
  } else {
    cli::cli_abort("Not implemented yet.")
  }

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
  return(fit1)

  # out <- create_lav_from_fitsem(fit, model, data, ...)
  # new("brlavaan", out)
}
