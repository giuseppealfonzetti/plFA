utils::globalVariables(c("clock"))
#' Internal function to check arguments in fit_plFA
#' @param ARGS List of arguments
#' @param P Number of pairs
stoc_args <- function(ARGS, P){

  pairs <- P*(P-1)/2
  tmp <- ARGS

  tmp$PAIRS_PER_ITERATION <- ifelse(!is.numeric(ARGS$PAIRS_PER_ITERATION), 8, round(ARGS$PAIRS_PER_ITERATION,0))
  tmp$MAXT <- ifelse(!is.numeric(ARGS$MAXT), round(100*pairs/tmp$PAIRS_PER_ITERATION,0), round(ARGS$MAXT,0))
  tmp$BURN <- ifelse(!is.numeric(ARGS$BURN), round(tmp$MAXT/2,0), round(ARGS$BURN,0))
  tmp$ETA <- ifelse(!is.numeric(ARGS$ETA), .005, ARGS$ETA)
  tmp$PAR1 <- ifelse(!is.numeric(ARGS$PAR1), 1, ARGS$PAR1)
  tmp$PAR2 <- ifelse(!is.numeric(ARGS$PAR2), .002, ARGS$PAR2)
  tmp$PAR3 <- ifelse(!is.numeric(ARGS$PAR3), 3/4, ARGS$PAR3)
  tmp$STEPSIZEFLAG <- ifelse(!is.numeric(ARGS$STEPSIZEFLAG), 1, max(1, round(abs(ARGS$STEPSIZEFLAG),0)))
  tmp$EACHCLOCK <- ifelse(!is.numeric(ARGS$EACHCLOCK), round(tmp$MAXT/10,0), round(ARGS$EACHCLOCK,0))
  tmp$CHECKCONV <- ifelse(!is.numeric(ARGS$CHECKCONV), 0, ARGS$CHECKCONV)


  out <- tmp

  return(out)
}

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
    METHOD = 'ucminf',
    CONTROL = list(),
    INIT = NULL,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0,
    NCORES = 1
){

  start_time <- Sys.time()

  if(sum(!is.finite(DATA))>0 | !is.matrix(DATA)) stop('DATA is not a numeric matrix.')
  if(sum(!is.finite(VALDATA))>0 | !is.matrix(VALDATA)) VALDATA <- DATA
  if(sum(!is.finite(CONSTR_LIST$CONSTRMAT))>0 | !is.matrix(CONSTR_LIST$CONSTRMAT) | sum(!(CONSTR_LIST$CONSTRMAT %in% c(0,1)))!=0) stop('CONSTRMAT must be a binary matrix')
  if(nrow(CONSTR_LIST$CONSTRMAT)!=ncol(DATA) | ncol(CONSTR_LIST$CONSTRMAT)>=nrow(CONSTR_LIST$CONSTRMAT)) stop('CONSTRMAT dimensions not allowed. Check Items x Factors.')
  if(is.null(CONSTR_LIST$CORRFLAG) | sum(is.finite(CONSTR_LIST$CORRFLAG))==0 | !(CONSTR_LIST$CORRFLAG %in% c(0,1))) stop('CORRFLAG must be 0 or 1')
  if(!(METHOD %in% c('ucminf','bernoulli', 'hyper'))) stop('Method not available.')
  # Identify model dimensions
  p <- ncol(DATA)
  n <- nrow(DATA)
  q <- ncol(CONSTR_LIST$CONSTRMAT)
  categories <- apply(DATA, 2, max, na.rm = T) + 1
  d = sum(categories)-p + sum(CONSTR_LIST$CONSTRMAT) + q*(q-1)/2

  tmp <- new('PlFaFit',
             cnstr = new('Constraints', loadings = CONSTR_LIST$CONSTRMAT, corrflag = CONSTR_LIST$CORRFLAG),
             dims = new('Dimensions', n = n, p = p, q = q, cat = categories, pairs = p*(p-1)/2),
             method = METHOD)
  # Check Initialisation

    if(is.vector(INIT)){
      if(length(INIT)!=d)
        stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
      else if((sum(!is.finite(INIT))==0)){
        message('1. Initialising at INIT vector.')
        tmp@init <-  INIT
      }else{
        message('1. Initialising at default values')
        lambda0_init <- c()
        s <- 0

        for (i in 1:length(categories)) {
          vec <- 1:(categories[i]-1)
          vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
          lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
          s <- s + categories[i] - 1
        }
        lambda_init = rep(0.5, sum(CONSTR_LIST$CONSTRMAT))
        transformed_rhos_init = rep(0, q*(q-1)/2)
        tmp@init <-  c(lambda0_init, lambda_init, transformed_rhos_init)
      }
    }

    if(!is.vector(INIT)){
      message('1. Initialising at default values')
      lambda0_init <- c()
      s <- 0

      for (i in 1:length(categories)) {
        vec <- 1:(categories[i]-1)
        vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
        lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
        s <- s + categories[i] - 1
      }
      lambda_init = rep(0.5, sum(CONSTR_LIST$CONSTRMAT))
      transformed_rhos_init = rep(0, q*(q-1)/2)
      tmp@init <-  c(lambda0_init, lambda_init, transformed_rhos_init)
    }

 message('2. Computing frequencies...')
  freq_start_time <- Sys.time()
  tmp@freq <- pairs_freq(DATA, categories)
  tmp@valfreq <- pairs_freq(VALDATA, categories)
  freq_end_time <- Sys.time()

  tmp@freqTime <- as.numeric(difftime(freq_end_time, freq_start_time, units = 'secs')[1])

  tmp@cores <- NCORES
  RcppParallel::setThreadOptions(numThreads = NCORES)


  # Numerical optimisation
  if(METHOD == 'ucminf'){

    message('3. Optimising with ucminf...')

    # Compute frequency table bivariate patterns

    Rwr_ncl <- function(par_vec){
      mod <- multiThread_completePairwise(
        N = n,
        C_VEC = categories,
        CONSTRMAT = CONSTR_LIST$CONSTRMAT,
        FREQ = tmp@freq,
        THETA = par_vec,
        CORRFLAG = CONSTR_LIST$CORRFLAG,
        SILENTFLAG = 1
      )
      out <- mod$iter_nll/n
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <- multiThread_completePairwise(
        N = n,
        C_VEC = categories,
        CONSTRMAT = CONSTR_LIST$CONSTRMAT,
        FREQ = tmp@freq,
        THETA = par_vec,
        CORRFLAG = CONSTR_LIST$CORRFLAG,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient/n
      return(out)
    }

    # list of ucminf args
    args <- list(
      'par' = tmp@init,
      'fn' = Rwr_ncl,
      'gr' = Rwr_ngr,
      'control' = ifelse(is.null(CONTROL$ctrl), list(), CONTROL$ctrl),
      'hessian' = ifelse(is.null(CONTROL$hessian), 0, CONTROL$hessian))

    # optimisation
    opt <- do.call(ucminf::ucminf, args)


    end_time <- Sys.time()
    tmp@numFit <- opt
    tmp@theta <- opt$par
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])


    message('4. Done! (', round(tmp@RTime,2),' secs)')
    return(tmp)
  }

  # Stochastic approximation
  if(METHOD == 'bernoulli' | METHOD == 'hyper'){

    stoFit <- new('StoFit')
    message(paste0('3. Optimising with ', METHOD, '...'))

    # Check stochastic control parameters
    cpp_ctrl <- stoc_args(CONTROL, P = p)

    # Check iterations selected
    if(!is.null(ITERATIONS_SUBSET)){
      stoFit@trajSubset <- unique(c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT))
    }else{
      stoFit@trajSubset <- 0:cpp_ctrl$MAXT
    }

    # Guarantee reproducibility stochastic optimisation
    # Note: R set.seed() is only needed by Bernoulli sampling.
    # For Hypergeometric sampling the seed is directly passed to cpp
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(
      list( 'FREQ' = tmp@freq,
            'VALFREQ' = tmp@valfreq,
            'N' = tmp@dims@n,
            'THETA_INIT' = tmp@init,
            'C_VEC' = tmp@dims@cat
            ),
      c( CONSTR_LIST, cpp_ctrl) )

    args$METHODFLAG <- ifelse(METHOD == 'hyper', 0, 1)
    fit <- do.call(plFA, args)

    end_time <- Sys.time()
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

    stoFit@projIters <- fit$post_index
    stoFit@trajSubset <- c(stoFit@trajSubset[stoFit@trajSubset<=fit$last_iter], fit$last_iter)
    stoFit@pathTheta <- fit$path_theta[stoFit@trajSubset + 1,]
    stoFit@pathAvTheta <- fit$path_av_theta[stoFit@trajSubset + 1,]
    stoFit@pathGrad <- fit$path_grad[stoFit@trajSubset,]
    stoFit@control <- cpp_ctrl
    stoFit@lastIter <- fit$last_iter
    stoFit@pathValNll <- cbind(fit$check_val_iter, fit$check_val_nll)
    stoFit@convergence <- fit$convergence
    tmp@theta <- stoFit@pathAvTheta[nrow(stoFit@pathAvTheta),]
    stoFit@cppTime <- summary(clock, units = 's')
    tmp@stoFit <- stoFit

    message('4. Done! (',round(tmp@RTime,2),' secs)')


    return(tmp)
  }
}
