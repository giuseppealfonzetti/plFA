utils::globalVariables(c("clock"))

#' Fit factor models for ordinal data with pairwise likelihood methods
#'
#' @description
#'
#' fit_plFA() allows to fit models both numerically and stochastically
#'
#' @param DATA Integer data matrix of dimension \eqn{n*p}. Categories must be coded starting from zero.
#' @param CONSTR_LIST List of constraints. It must contain \tabular{ll}{
#'    \code{CONSTRMAT} \tab \eqn{p*q} binary matrix. Elements set to 1 refers to free loding paameters. Elements set to 0 refer to null loadings. \cr
#'    \tab \cr
#'    \code{CORRFLAG} \tab Binary indicator. Set it to 0 if the latent variables are assumed to be independent. Set it 1 otherwise. \cr
#' }
#' @param VALDATA Data matrix used as validation set to check for convergence.
#' @param METHOD Label for the method chosen. Possible values are: \tabular{ll}{
#'    \code{'ucminf'} \tab for estimation via the numerical optimiser from the \code{ucminf} package. \cr
#'    \tab \cr
#'    \code{'bernoulli'} \tab for stochastic estimation with weights following a Bernoulli distribution. \cr
#'     \tab \cr
#'    \code{'hyper'} \tab for stochastic estimation with weights following a multivariate hypergeometric distribution. \cr
#' }
#' @param NCORES Integer value setting the number of threads to carry out the estimation.
#' @param INIT Initialising vector. If not provided, optiisation starts at the zero vector.
#' @param ITERATIONS_SUBSET Vector with the indexes of iterations used to store trajectories.
#' @param CONTROL List of control options to pass to the optimiser: \tabular{ll}{
#'    \code{control} \tab  See \code{ucminf} documentation \cr
#'    \tab \cr
#'    \code{hessian} \tab  See \code{ucminf} documentation \cr
#'    \tab \cr
#'    \code{PAIRS_PER_ITERATION} \tab  (Average) number of pairs to draw stochastically at each iteration with method (\code{'bernoulli'}) \code{'hyper'}.   \cr
#'    \tab \cr
#'    \code{MAXT} \tab  Maximum number of iterations.\cr
#'    \tab \cr
#'    \code{BURN} \tab  Length of the burning period.\cr
#'    \tab \cr
#'    \code{ETA} \tab  Initial stepsize.\cr
#'    \tab \cr
#'    \code{PAR1} \tab Hyperparameter for stepsize scheduling by Xu (2011): Scaling. \cr
#'    \tab \cr
#'    \code{PAR2} \tab Hyperparameter for stepsize scheduling by Xu (2011): Smallest Hessian eigenvalue. \cr
#'    \tab \cr
#'    \code{PAR1} \tab Hyperparameter for stepsize scheduling by Xu (2011): Decay rate. \cr
#'    \tab \cr
#'    \code{STEPSIZEFLAG} \tab Choose stepsize scheduling: Set 0 for Polyak and Juditsky (1992), 1 for Xu (2011). \cr
#'    \tab \cr
#'    \code{EACHCLOCK} \tab How often (in terms of iterations) to measure single iteration computational times (using \code{RcppClock}. \cr
#'    \tab \cr
#'    \code{CHECKCONV} \tab Flag to check for convergence using complete pairwise likelihood on the validation set. \cr
#'    \tab \cr
#'    \code{EACHCHECK} \tab How often (in terms of iterations) to check for convergence. \cr
#'    \tab \cr
#'    \code{TOL} \tab Tolerance between consecutive evaluations of the validation composite likelihood. \cr
#'    \tab \cr
#'    \code{TOLCOUNT} \tab How many checks the algorithm is allowed to accept if the validation negative composite likelihood starts increasing. \cr
#'    \tab \cr
#' }
#' @param VERBOSEFLAG TRUE verbose output
#' @export
fit_plFA <- function(
    DATA,
    CONSTR_LIST,
    VALDATA = NULL,
    METHOD = c('ucminf', "SA"),
    INIT = NULL,
    INIT_METHOD = NULL,
    CONTROL = list(),
    CPP_CONTROL_MAIN = NULL,
    CPP_CONTROL_INIT = NULL,    VERBOSE = FALSE,
    NCORES = 1
){

  start_time <- Sys.time()

  dat <- check_data(DATA)
  constr_list <- check_cnstr(CONSTR_LIST)
  method <- match.arg(METHOD)


  # Identify model dimensions
  dims <- check_dims(dat, constr_list)

  tmp <- new('PlFaFit',
             cnstr = new('Constraints',
                         loadings = constr_list$CONSTRMAT,
                         corrflag = constr_list$CORRFLAG,
                         stdlv    = constr_list$STDLV,
                         loglatsd = constr_list$CONSTRLOGSD),
             dims = new('Dimensions',
                        n     = dims$n,
                        p     = dims$p,
                        q     = dims$q,
                        cat   = dims$cat,
                        pairs = dims$pairs,
                        nthr  = dims$nthr,
                        nload = dims$nload,
                        ncorr = dims$ncorr,
                        nvar  = dims$nvar,
                        npar  = dims$d),
             method = METHOD)



  # Set up multi-threads computations
  tmp@cores <- NCORES
  RcppParallel::setThreadOptions(numThreads = NCORES)


  if(VERBOSE) message('- Computing frequencies...')
  freq_start_time <- Sys.time()
  tmp@freq <- pairs_freq(DATA, dims$cat)
  freq_end_time <- Sys.time()
  tmp@freqTime <- as.numeric(difftime(freq_end_time, freq_start_time, units = 'secs')[1])

  # Check Initialisation
  tmp@init <- check_init(INIT, tmp@freq, dims, constr_list, INIT_METHOD, VERBOSE, CPP_ARGS = CPP_CONTROL_INIT)





  # Numerical optimisation
  if(METHOD == 'ucminf'){

    if(VERBOSE) message('- Optimising with ucminf...')
    tmp@method <- "ucminf"

    # Compute frequency table bivariate patterns

    Rwr_ncl <- function(par_vec){
      n <- dims$n
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = constr_list$CONSTRMAT,
        CONSTRLOGSD= constr_list$CONSTRLOGSD,
        FREQ       = tmp@freq,
        THETA      = par_vec,
        CORRFLAG   = constr_list$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        SILENTFLAG = 1
      )
      out <- mod$iter_nll/dims$n
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      n <- dims$n
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = constr_list$CONSTRMAT,
        CONSTRLOGSD= constr_list$CONSTRLOGSD,
        FREQ       = tmp@freq,
        THETA      = par_vec,
        CORRFLAG   = constr_list$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        GRFLAG     = 1,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient/dims$n
      return(out)
    }

    # list of ucminf args
    args <- list(
      'par' = tmp@init,
      'fn' = Rwr_ncl,
      'gr' = Rwr_ngr,
      'control' = ifelse(is.null(CONTROL$ctrl), list(), CONTROL$ctrl),
      'hessian' = ifelse(is.null(CONTROL$hessian), 2, CONTROL$hessian))

    # optimisation
    opt <- do.call(ucminf::ucminf, args)


    end_time <- Sys.time()
    tmp@numFit <- opt
    tmp@theta <- opt$par
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])


    if(VERBOSE) message('Done! (', round(tmp@RTime,2),' secs)')
    RcppParallel::setThreadOptions(numThreads = 1)
    return(tmp)
  }

  if(METHOD=="SA"){
    if(!is.null(VALDATA)){
      if(VERBOSE) message('- Computing frequencies on validation data...')
      tmp@valfreq <- pairs_freq(VALDATA, dims$cat)
    }

    if(VERBOSE) message('- Optimising via pairwise stochastic approximation...')
    stoFit <- new('StoFit')
    tmp@method <- "SA"


    sa_args <- check_plSA_args(FREQ        = tmp@freq,
                               VALFREQ     = if(is.null(VALDATA)) NULL else tmp@valfreq,
                               DIMS        = dims,
                               CONSTR_LIST = constr_list,
                               LIST        = CPP_CONTROL_MAIN,
                               SETTING     = "main")

    sa_args <-  c(list(N           = dims$n,
                       C_VEC       = dims$cat,
                       CONSTRMAT   = constr_list$CONSTRMAT,
                       CONSTRLOGSD = constr_list$CONSTRLOGSD,
                       THETA_INIT  = tmp@init,
                       NTHR        = dims$nthr,
                       NLOAD       = dims$nload,
                       NCORR       = dims$ncorr,
                       NVAR        = dims$nvar
    ), sa_args)

    opt <- do.call(cpp_plSA, sa_args)
    end_time <- Sys.time()
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

    if(!opt$convergence_burn) warning(paste0("Burn-in period did not reach tolerance level ", sa_args$TOL_NLL*10, ". Stopped at ", sa_args$BURN, " iterations. Try increasing BURN or STEP0."))
    if(opt$proj_after_burn)   warning("Projections performed after the burn-in. Trajectories might be unstable.")
    if(opt$convergence==0)    warning(paste0("Burn-in period did not reach tolerance level ", sa_args$TOL_NLL, ". Stopped at ", sa_args$MAXT, " iterations. Try increasing MAXT or STEP0.."))
    if(opt$convergence==-1)   warning(paste0("Divergent trajectories. Decrease STEP0."))
    if(opt$neg_pdiff)         warning("Possible divergent trajectories detected. Try decreasing STEP0.")


    # stoFit@path_iters   <- opt$path_iters
    # stoFit@path_theta   <- opt$path_theta
    # stoFit@path_avtheta <- opt$path_avtheta
    # stoFit@path_nll     <- opt$path_nll
    stoFit@nll    <- opt$nll
    stoFit@last_iter    <- opt$last_iter
    # stoFit@convergence  <- opt$convergence
    stoFit@cppTime <- summary(clock, units = 's')
    tmp@stoFit <- stoFit
    tmp@theta <- opt$avtheta


    if(VERBOSE) message('Done! (', round(tmp@RTime,2),' secs)')
    RcppParallel::setThreadOptions(numThreads = 1)
    return(tmp)







  }
}

































