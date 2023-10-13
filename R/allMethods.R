#' Show method for the output of fit_plFA()
setMethod('show', 'PlFaFit', function(object){
  if(object@method == 'ucminf'){
    cat(
      "Numerical estimate obtained via ucminf.\n"
    )
  }else{
    cat(
      "- Stochastic estimate\n\n",
      "    Sampling scheme: ", object@method, "\n",
      "    Pairs per iteration: ", object@stoFit@control$PAIRS_PER_ITERATION, " out of ", object@dims@pairs, "\n",
      "    Iterations: ", object@stoFit@lastIter, "\n",
      "\n- Use getPar() to extract parameter estimates."
    )
  }
})


#' Extract parameters trajectory along the optimisation
#'
#' @param OBJ Object of class StoFit
#' @param LAB Can take values \code{pathTheta} for raw trajectories or \code{pathAvTheta} for averaged ones. Set by default at \code{pathAvTheta}.
#' @param OPTION Parameterisation option.

#' @export
setGeneric('getThetaPath', function(OBJ, LAB = 'pathAvTheta', OPTION = 'transformed') standardGeneric('getThetaPath'))
setMethod('getThetaPath', 'PlFaFit', function(OBJ, LAB, OPTION) extract_theta_path(OBJ, LAB, OPTION))

#' Extract parameters trajectory along the optimisation
#'
#' @param OBJ Object of class FitPlFa
#' @param LAB Can take values \code{pathTheta} for raw trajectories or \code{pathAvTheta} for averaged ones. Set by default at \code{pathAvTheta}.
#' @param OPTION Parameterisation option.
extract_theta_path <- function(OBJ, LAB = 'pathAvTheta', OPTION = 'transformed'){
  if(!(LAB%in%c('pathTheta', 'pathAvTheta'))) stop('The extraction of this path is not implemented.')
  if(!(OPTION%in%c('raw', 'transformed'))) stop('OPTION choice not implemented.')

  # out <- list()
  iters <- OBJ@stoFit@trajSubset
  path  <- slot(OBJ@stoFit, LAB)

  out <- dplyr::tibble(iter = iters)
  parList <- split(t(path), rep(1:nrow(path), each = ncol(path)))
  newParList <- lapply(parList, function(x) extract_par(THETA= x,
                                            OPTION = OPTION,
                                            C = sum(OBJ@dims@cat),
                                            P = OBJ@dims@p,
                                            Q = OBJ@dims@q,
                                            CONSTRMAT = OBJ@cnstr@loadings))

  out$par <- newParList
  return(out)
}


#' Extract parameters
#'
#' @param OBJ Object of class StoFit, PlFaFit. Also usable on raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param OPTION Can take values \tabular{ll}{
#'    \code{'raw'} \tab Raw vector of parameters. Latent correlations are reparametrised, thus not directly interpretable. \cr
#'    \tab \cr
#'    \code{'transformed'} \tab Transformed vector of parameters, with interpretable latent correlations. \cr
#'     \tab \cr
#'    \code{'list'} \tab List containing the ordered thresholds vector, the loading matrix and the latent correlation matrix. \cr
#' }
#' @export
setGeneric('getPar', function(OBJ, ...) standardGeneric('getPar'))
setMethod('getPar', 'PlFaFit', function(OBJ, ...) extract_par(THETA = OBJ@theta, C = sum(OBJ@dims@cat), P = OBJ@dims@p, Q = OBJ@dims@q, CONSTRMAT = OBJ@cnstr@loadings, ...))
setMethod('getPar', 'vector', function(OBJ, ...) extract_par(THETA = OBJ, ...))

#' Extract parameters
#'
#' @param THETA Raw theta vector.
#' @param OPTION Can take values \tabular{ll}{
#'    \code{'raw'} \tab Returns \code{THETA}. \cr
#'    \tab \cr
#'    \code{'transformed'} \tab Transformed vector of parameters, with interpretable latent correlations. \cr
#'     \tab \cr
#'    \code{'list'} \tab List containing the ordered thresholds vector, the loading matrix and the latent correlation matrix. \cr
#' }
#' @param C Sum of the number of categories for each item.
#' @param P Number if items
#' @param Q Number of latent variables
#' @param CONSTRMAT Matrix of dimension \eqn{p * q} with binary loading constraints. 1 for free loadings, 0 otherwise.
extract_par <- function(THETA, OPTION = 'transformed', C, P, Q, CONSTRMAT){
  if(!(OPTION %in% c('raw', 'transformed', 'list'))) stop('The OPTION chosen is not implemented')


  if(OPTION == 'raw') return(THETA)
  if(OPTION == 'transformed') return(c(THETA[1:(length(THETA)-Q*(Q-1)/2)], get_corr(THETA = THETA, Q = Q)))
  if(OPTION == 'list') {
    out <- list()
    thr = THETA[1:(C-P)]
    ld <- CONSTRMAT
    ld_vec <- get_lambda(THETA = THETA, C = C, P = P, Q = Q)
    s <- 1
    for (j in 1:ncol(ld)) {
      for (i in 1:nrow(ld)) {
        if(CONSTRMAT[i,j]!=0){
          ld[i,j] <- ld_vec[s]
          s = s+1
        }
      }
    }
    S <- get_S(THETA = THETA, Q = Q)

    out <- list(thresholds = thr, loadings = ld, latent_correlations = S)
    return(out)
  }
}

#' Compute estimates' variance
#' @export
setGeneric('computeVar', function(OBJ, DATA, NUMDERIV = F) standardGeneric('computeVar'), signature = 'OBJ')
setMethod('computeVar', 'PlFaFit',
          function(OBJ, DATA, NUMDERIV){ if(OBJ@method == 'ucminf'){
            compute_var(THETA     = OBJ@theta,
                       C_VEC     = OBJ@dims@cat,
                       N         = OBJ@dims@n,
                       IT        = NULL,
                       PAIRS     = NULL,
                       PPI       = NULL,
                       CONSTRMAT = OBJ@cnstr@loadings,
                       CORRFLAG  = OBJ@cnstr@corrflag,
                       FREQ      = OBJ@freq,
                       DATA      = DATA,
                       METHOD    = OBJ@method,
                       NUMDERIV  = NUMDERIV)
          }else{
            compute_var(THETA     = OBJ@theta,
                       C_VEC     = OBJ@dims@cat,
                       N         = OBJ@dims@n,
                       IT        = OBJ@stoFit@lastIter - OBJ@stoFit@control$BURN,
                       PAIRS     = OBJ@dims@pairs,
                       PPI       = OBJ@stoFit@control$PAIRS_PER_ITERATION,
                       CONSTRMAT = OBJ@cnstr@loadings,
                       CORRFLAG  = OBJ@cnstr@corrflag,
                       FREQ      = OBJ@freq,
                       DATA      = DATA,
                       METHOD    = OBJ@method,
                       NUMDERIV  = NUMDERIV)
          }}
)

#setMethod('computeVar', 'vector', function(OBJ, DATA, ...) compute_var(OBJ, DATA, ...))
#' @export
compute_var <- function(THETA, C_VEC, N, IT = NULL, PAIRS = NULL, PPI = NULL,  CONSTRMAT, CORRFLAG, FREQ, DATA, METHOD, NUMDERIV = F){

  opt_noise <- NA
  Hhat <- NA
  if(NUMDERIV){
    message('1. Computing H numerically...')
    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <-multiThread_completePairwise(
        N = N,
        C_VEC = C_VEC,
        CONSTRMAT = CONSTRMAT,
        FREQ = FREQ,
        THETA = par_vec,
        CORRFLAG = CORRFLAG,
        GRFLAG = 1,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient/N
      return(out)
    }

    Hhat <- numDeriv::jacobian(Rwr_ngr, THETA)
  }else{
    message('1. Estimating H...')
    Hhat <- estimate_H(
      C_VEC = C_VEC,
      A = CONSTRMAT,
      THETA = THETA,
      FREQ = FREQ,
      N = N,
      CORRFLAG = CORRFLAG
    )$est_H
  }

  message('2. Estimating J...')
  Jhat <- estimate_J(
    Y = DATA,
    C_VEC = C_VEC,
    A = CONSTRMAT,
    THETA = THETA,
    CORRFLAG = CORRFLAG
  )$est_J

  if(!(det(Hhat)>0)) stop('H not invertible.')
  message('3. Inverting H...')
  invH <- solve(Hhat)
  sandwich <- invH %*% Jhat %*% invH

  message('3. Computing the variances...')
  if(METHOD =='ucminf'){
    asy_var <- diag(sandwich/N)
  }else{
    asy_var <- diag(sandwich/N)
    a1 <- PAIRS*(PAIRS - PPI)/(PPI*(PAIRS-1))
    a2 <- (PAIRS-PPI)/(PPI*(PAIRS-1))
    opt_noise <- diag(invH%*%(a1*Hhat - a2*Jhat)%*%invH/(N*IT))

  }
  message('Done!')

  return(
    list(
      H = Hhat,
      J = Jhat,
      asymptotic_variance = asy_var,
      optimisation_noise = opt_noise

    )
  )

}

