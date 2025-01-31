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
#' @param std.lv TBC
#' @param start Starting values to use.
#' @param estimator (For `lavaan` compatibility only). The estimator is `PML`.
#' @param estimator.args A list of arguments for [fit_plFA()]. TBC.
#' @param information (For `lavaan` compatibility only). The information matrix is...
#' @param control A list of control parameters for the estimation algorithm. See [fit_plFA()] for more information.
#' @param ... Additional arguments to be passed to [lavaan()].
#'
#' @returns A `plFAlavaan` object.
#' @export
cfa <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(
      method = c("ucminf", "SA"),
      init_method = c("SA", "custom", "standard"),
      cpp_control_init = NULL,
      ncores = 1,
      valdata = NULL,
      computevar_numderiv = FALSE
    ),
    start = NULL,
    information = "observed",
    control = list(),
    verbose = FALSE,
    ...
) {

  # Validate arguments ---------------------------------------------------------
  if (!is_ordinal_df(data)) {
    cli::cli_abort("Ordinal data required")
  }
  method <- estimator.args$method
  if (is.null(method)) method <- "ucminf"
  method <- rlang::arg_match(method, c("ucminf", "SA"))
  init_method <- estimator.args$init_method
  if (is.null(init_method)) init_method <- "SA"
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  valdata <- estimator.args$valdata
  if (method == "SA") {
    cpp_control_main <- control
    control <- list()
  }
  cpp_control_init <- estimator.args$cpp_control_init
  computevar_numderiv <- estimator.args$computevar_numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE

  # Initialise {lavaan} model object -------------------------------------------
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$std.lv <- std.lv
  lavargs$do.fit <- FALSE
  fit0 <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)

  if (fit0@Data@ngroups > 1L) {
    cli::cli_abort("Multigroup analysis is not currently supported.")
  }

  # Fit plFA -------------------------------------------------------------------
  D <- fit0@Data@X[[1]] - 1

  # Build the p x q loading constraints matrix
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
  psi <- lavaan::inspect(fit0, "est")$psi
  LTpsi <- FREE$psi[lower.tri(psi)]
  if (length(LTpsi) > 0 & any(LTpsi > 0)) {
    corrflag <- TRUE
  } else {
    corrflag <- FALSE
  }
  constrvar <- diag(psi)
  constrvar[diag(FREE$psi) > 0] <- NA

  # Linear loadings constraints
  llc <- NULL
  pt <- partable(fit0)
  ptlc <- pt[pt$op == "==", ]
  Lambdaid <- lavaan::inspect(fit0, "free")$lambda

  if (nrow(ptlc) > 0L) {
    llc <- list()
    coord1 <- ptlc$lhs
    coord2 <- strsplit(ptlc$rhs, "\\+")[[1]]
    for (i in seq_along(coord1)) {
      idx <- pt$id[pt$label == coord1[i]]
      llc[[i]] <- list(which(Lambdaid == idx, arr.ind = TRUE))

      for (j in seq_along(coord2)) {
        cz <- strsplit(coord2[j], "\\*")[[1]]
        idx <- pt$id[pt$label == cz[2]]
        tmp <- which(Lambdaid == idx, arr.ind = TRUE)
        llc[[i]] <- c(llc[[i]], list(as.numeric(c(cz[1], tmp))))
      }
    }
  }

  # Build constraint list
  constr_list <- list(
    CONSTRMAT = A,
    CONSTRVAR = constrvar,  # FIXME: Correct way to specify variance constraints?
    CORRFLAG = corrflag,
    STDLV = std.lv,
    LLC = llc
  )

  fit1 <- fit_plFA(
    DATA = D,
    CONSTR_LIST = constr_list,
    METHOD = method,
    INIT_METHOD = init_method,
    NCORES = ncores,
    VALDATA = valdata,
    CONTROL = control,
    CPP_CONTROL_MAIN = cpp_control_main,
    CPP_CONTROL_INIT = cpp_control_init,
    INIT = start,
    VERBOSE = verbose
  )
  vars <- computeVar(
    OBJ = fit1,
    DATA = D,
    NUMDERIV = computevar_numderiv,
    OPTION = "transformed"
  )

  # list(fit0 = fit0, fit1 = fit1, vars = vars)
  out <- create_lav_from_fitplFA(fit0, fit1, vars, D)
  new("plFAlavaan", out)
}

create_lav_from_fitplFA <- function(fit0, fit1, vars, D) {

  # Get coefficients and standard errors
  FREE <- lavaan::inspect(fit0, what = "free")

  parlist <- extract_par(
    THETA = fit1@theta,
    OPTION = "list",
    CONSTRMAT = fit1@cnstr@loadings,
    CONSTRLOGSD = fit1@cnstr@loglatsd,
    LLC = fit1@cnstr@llc,
    NTHR = fit1@dims@nthr,
    NLOAD = fit1@dims@nload,
    NCORR = fit1@dims@ncorr,
    NVAR = fit1@dims@nvar
  )
  lambda <- parlist$loadings[FREE$lambda > 0]
  tau <- parlist$thresholds[FREE$tau > 0]
  psi <- parlist$latent_correlations[FREE$psi > 0 & lower.tri(FREE$psi, diag = TRUE)]

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

