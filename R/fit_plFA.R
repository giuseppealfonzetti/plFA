#' Internal function to check arguments in fit_plFA
#' @export
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

  out <- tmp

  return(out)
}

#' Fit factor models for ordinal data with pairwise likelihood methods
#'
#' fit_plFA() allows to fit models both numerically and stochastically
#'
#'
#' @export
fit_plFA2 <- function(
    DATA_LIST = list('DATA', 'CONSTRMAT', 'CORRFLAG'),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
    INIT = NULL,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0,
    NCORES = 1
){

  out <- list()
  start_time <- Sys.time()
  # Identify model dimensions
  p <- ncol(DATA_LIST$DATA)
  n <- nrow(DATA_LIST$DATA)
  q <- ncol(DATA_LIST$CONSTRMAT)
  categories <- apply(DATA_LIST$DATA, 2, max, na.rm = T) + 1 ### number of categories in each item
  d = sum(categories)-p + sum(DATA_LIST$CONSTRMAT) + q*(q-1)/2

  out$categories <- categories
  # Check Constraints
  if(is.null(DATA_LIST$CONSTRMAT)){
    stop('CONSTRMAT not declared')
  }else if(nrow(DATA_LIST$CONSTRMAT)!=ncol(DATA_LIST$DATA) || ncol(DATA_LIST$CONSTRMAT)>=nrow(DATA_LIST$CONSTRMAT)){
    stop('CONSTRMAT dimensions not acceptable. Check Items x Factors.')
  }
  out$constraints <- DATA_LIST$CONSTRMAT

  # Check Initialisation
  if(is.vector(INIT)){
    if(length(INIT)!=d)
      stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
    else
      message('1. Initialising at init vector.')
    out$theta_init <-  INIT
  }else{
    if(is.null(INIT))
      message('1. Initialising at default values')
    lambda0_init <- c()
    s <- 0

    for (i in 1:length(categories)) {
      vec <- 1:(categories[i]-1)
      vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
      lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
      s <- s + categories[i] - 1
    }
    lambda_init = rep(0.5, sum(DATA_LIST$CONSTRMAT))
    transformed_rhos_init = rep(0, q*(q-1)/2)
    out$theta_init <-  c(lambda0_init, lambda_init, transformed_rhos_init)
  }

  # Check if method entry is correct
  if(!(METHOD %in% c('ucminf', 'hyper'))) stop('Method not available.')
  out$method <- METHOD
  RcppParallel::setThreadOptions(numThreads = NCORES)


  # Numerical optimisation
  if(METHOD == 'ucminf'){

    message('2. Optimising with ucminf...')

    # Compute frequency table bivariate patterns
    freq_tab <- pairs_freq(DATA_LIST$DATA, categories)

    Rwr_ncl <- function(par_vec){
      lambda0_ <- par_vec[1:length(lambda0_init)]
      lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
      transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
      mod <-multiThread_completePairwise(
        Y = DATA_LIST$DATA,
        C_VEC = categories,
        A = DATA_LIST$CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = DATA_LIST$CORRFLAG,
        GRFLAG = 0,
        SILENTFLAG = 1
      )
      out <- mod$iter_nll
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      lambda0_ <- par_vec[1:length(lambda0_init)]
      lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
      transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
      mod <-multiThread_completePairwise(
        Y = DATA_LIST$DATA,
        C_VEC = categories,
        A = DATA_LIST$CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = DATA_LIST$CORRFLAG,
        GRFLAG = 1,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient
      return(out)
    }

    # list of ucminf args
    args <- list(
      'par' = out$theta_init,
      'fn' = Rwr_ncl,
      'gr' = Rwr_ngr,
      'control' = UCMINF_CONTROL$ctrl,
      'hessian' = UCMINF_CONTROL$hessian)

    # optimisation
    start_opt <- Sys.time()
    opt <- do.call(ucminf::ucminf, args)
    end_opt <- Sys.time()
    out$num_time <- as.numeric(difftime(end_opt, start_opt, units = 'secs')[1])

    out$fit <- opt

    out$control <- UCMINF_CONTROL
    out$theta   <- opt$par

    end_time <- Sys.time()
    out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
    message('3. Done! (', round(out$time,2),' secs)')
    return(out)
  }

  # Stochastic approximation of numerical optimiser
  if(METHOD == 'bernoulli' | METHOD == 'hyper' | METHOD == 'recycle_hyper' ){

    message(paste0('2. Optimising with ', METHOD, '...'))

    # Check stochastic control parameters
    # cpp_ctrl <- check_args(CPP_CONTROL, N = n)

    cpp_ctrl <- CPP_CONTROL
    # Check iterations selected
    if(!is.null(ITERATIONS_SUBSET)){
      out$iterations_subset <- unique(c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT))
    }else{
      out$iterations_subset <- 0:cpp_ctrl$MAXT
    }


    # Guarantee reproducibility stochastic optimisation
    # Note: R set.seed() is only needed by Bernoulli sampling.
    # For Hypergeometric sampling the seed is directly passed to cpp
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    # args <- append(
    #   list( 'TAU' = out$theta_init[1:(sum(out$categories)-p)],
    #         'LAMBDA' = get_lambda(out$theta_init, sum(out$categories), p, q),
    #         'TRANSFORMED_RHOS' = get_corr(out$theta_init, q),
    #         'C_VEC' = out$categories),
    #   c( DATA_LIST, cpp_ctrl) )
    args <- append(
      list( 'THETA_INIT' = out$theta_init,
            'C_VEC' = out$categories),
      c( DATA_LIST, cpp_ctrl) )

    args$METHODFLAG <- dplyr::if_else(METHOD == 'hyper', 0, 1)

    if(METHOD == 'hyper'){args$METHODFLAG <- 0} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 1}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 2}
    fit <- do.call(plFA, args)
    message('\n3. Rearranging output...')

    end_time <- Sys.time()
    out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

    out$post_index <- fit$post_index
    out$iterations_subset <- c(out$iterations_subset, fit$last_iter)#[c(out$iterations_subset, fit$las_iter-1)<=(fit$las_iter-1)]
    out$iterations_subset <- out$iterations_subset[out$iterations_subset<=fit$last_iter]
    fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,]
    fit$path_av_theta <- fit$path_av_theta[out$iterations_subset + 1,]
    fit$path_grad     <- fit$path_grad[out$iterations_subset,]
    # fit$path_theta    <- purrr::reduce(fit$path_theta, rbind)
    # fit$path_av_theta <- purrr::reduce(fit$path_av_theta, rbind)
    # fit$path_grad     <- purrr::reduce(fit$path_grad, rbind)
    out$control <- cpp_ctrl
    out$fit <- fit
    out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]

    if('RcppClock'%in% (.packages())) out$clock <- summary(clock, units = 's')
    message('4. Done! (', round(out$time,2),' secs)')


    return(out)
  }

}



#' Fit factor models for ordinal data with pairwise likelihood methods
#'
#' fit_plFA() allows to fit models both numerically and stochastically
#'
#' @param DATA Integer data matrix of dimension \eqn{n*p}. Categories must be coded starting from zero.
#' @param CONSTR_LIST List of constraints. It must contain \tabular{ll}{
#'    \code{CONSTRMAT} \tab \eqn{p*q} binary matrix. Elements set to 1 refers to free loding paameters. Elements set to 0 refer to null loadings. \cr
#'    \tab \cr
#'    \code{CORRFLAG} \tab Binary indicator. Set it to 0 if the latent variables are assumed to be independent. Set it 1 otherwise. \cr
#' }
#' @param METHOD Label for the method chosen. Possible values are: \tabular{ll}{
#'    \code{ucminf} \tab for estimation via the numerical optimiser from the \code{ucminf} package. \cr
#'    \tab \cr
#'    \code{bernoulli} \tab for stochastic estimation with weights following a Bernoulli distribution. \cr
#'     \tab \cr
#'    \code{hyper} \tab for stochastic estimation with weights following a multivariate hypergeometric distribution. \cr
#' }
#' @param NCORES Integer value setting the number of threads to carry out the estimation.
#'
#'
#' @export
fit_plFA <- function(
    DATA,
    CONSTR_LIST,
    METHOD = 'ucminf',
    CONTROL = list(),
    INIT = NULL,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0,
    NCORES = 1
){

  start_time <- Sys.time()

  if(sum(!is.finite(DATA))>0 | !is.matrix(DATA)) stop('DATA is not a numeric matrix.')
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





  RcppParallel::setThreadOptions(numThreads = NCORES)


  # Numerical optimisation
  if(METHOD == 'ucminf'){

    message('2. Optimising with ucminf...')

    # Compute frequency table bivariate patterns
    freq_tab <- pairs_freq(DATA, categories)

    Rwr_ncl <- function(par_vec){
      lambda0_ <- par_vec[1:(sum(categories)-p)]
      lambda_ <- par_vec[(sum(categories)-p+1):(sum(categories)-p+sum(CONSTR_LIST$CONSTRMAT))]
      transformed_rhos_ <- par_vec[(sum(categories)-p+1+sum(CONSTR_LIST$CONSTRMAT)):length(par_vec)]
      mod <-multiThread_completePairwise(
        Y = DATA,
        C_VEC = categories,
        A = CONSTR_LIST$CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = CONSTR_LIST$CORRFLAG,
        GRFLAG = 0,
        SILENTFLAG = 1
      )
      out <- mod$iter_nll
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      lambda0_ <- par_vec[1:(sum(categories)-p)]
      lambda_ <- par_vec[(sum(categories)-p+1):(sum(categories)-p+sum(CONSTR_LIST$CONSTRMAT))]
      transformed_rhos_ <- par_vec[(sum(categories)-p+1+sum(CONSTR_LIST$CONSTRMAT)):length(par_vec)]
      mod <-multiThread_completePairwise(
        Y = DATA,
        C_VEC = categories,
        A = CONSTR_LIST$CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = CONSTR_LIST$CORRFLAG,
        GRFLAG = 1,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient
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
    # start_opt <- Sys.time()
    opt <- do.call(ucminf::ucminf, args)


    end_time <- Sys.time()
    tmp@numFit <- opt
    tmp@theta <- opt$par
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])


    message('3. Done! (', round(tmp@RTime,2),' secs)')
    return(tmp)
  }

  # Stochastic approximation
  if(METHOD == 'bernoulli' | METHOD == 'hyper'){

    stoFit <- new('StoFit')
    message(paste0('2. Optimising with ', METHOD, '...'))

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
      list( 'DATA' = DATA,
            'THETA_INIT' = tmp@init,
            'C_VEC' = tmp@dims@cat
            ),
      c( CONSTR_LIST, cpp_ctrl) )

    args$METHODFLAG <- ifelse(METHOD == 'hyper', 0, 1)
    fit <- do.call(plFA, args)
    message('\n3. Rearranging output...')


    end_time <- Sys.time()
    tmp@RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
    stoFit@projIters <- fit$post_index
    stoFit@trajSubset <- c(stoFit@trajSubset[stoFit@trajSubset<=fit$last_iter], fit$last_iter)
    stoFit@pathTheta <- fit$path_theta[stoFit@trajSubset + 1,]
    stoFit@pathAvTheta <- fit$path_av_theta[stoFit@trajSubset + 1,]
    stoFit@pathGrad <- fit$path_grad[stoFit@trajSubset,]
    stoFit@control <- cpp_ctrl
    tmp@theta <- stoFit@pathAvTheta[nrow(stoFit@pathAvTheta),]
    stoFit@cppTime <- summary(clock, units = 's')
    tmp@stoFit <- stoFit

    message('4. Done! (',round(tmp@RTime,2),' secs)')


    return(tmp)
  }
}
