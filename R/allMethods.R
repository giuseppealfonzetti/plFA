#' Show method for the output of fit_plFA()
#' @param object Object of class StoFit
setMethod('show', 'PlFaFit', function(object){
  if(object@method == 'ucminf'){
    cat(
      "- Numerical estimate obtained via ucminf.\n\n",
      "   Total time:", round(object@RTime, 2), "s (Data reduction:", round(object@freqTime,2), "s)\n",
      "   Cores used:", object@cores, "\n",

      "\n- Use getPar() to extract parameter estimates."

    )
  }else{
    cat(
      "- Stochastic estimate\n\n",
      "    Sampling scheme:", object@method, "\n",
      "    Pairs per iteration:", object@stoFit@control$PAIRS_PER_ITERATION, " out of ", object@dims@pairs, "\n",
      "    Iterations:", object@stoFit@lastIter, "\n",
     # "    Convergence:", object@stoFit@convergence, "\n",
      "    Total time:", round(object@RTime, 2), "s (Data reduction:", round(object@freqTime,2), "s)\n",
      "    Cores used:", object@cores, "\n",


      "\n- Use getPar() to extract parameter estimates."
    )
  }
})


#' Extract parameters trajectory along the optimisation
#'
#' @param OBJ Object of class StoFit
#' @param LAB Can take values \code{pathTheta} for raw trajectories or \code{pathAvTheta} for averaged ones. Set by default at \code{pathAvTheta}.
#' @param OPTION Parameterisation option.
#'
#' @export
setGeneric('getThetaPath', function(OBJ, LAB = 'pathAvTheta', OPTION = 'transformed') standardGeneric('getThetaPath'))

#' Extract parameters trajectory along the optimisation
#'
#' @param OBJ Object of class StoFit
#' @param LAB Can take values \code{pathTheta} for raw trajectories or \code{pathAvTheta} for averaged ones. Set by default at \code{pathAvTheta}.
#' @param OPTION Parameterisation option.
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
                                            CONSTRMAT = OBJ@cnstr@loadings,
                                            CORRFLAG = OBJ@cnstr@corrflag))

  out$par <- newParList
  return(out)
}


#' Extract parameters
#'
#' @param OBJ Object of class StoFit, PlFaFit. Also usable on raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
#' @export
setGeneric('getPar', function(OBJ, ...) standardGeneric('getPar'))
#' Extract parameters
#'
#' @param OBJ Object of class PlFaFit.
#' @param ... Additional arguments
setMethod('getPar', 'PlFaFit', function(OBJ, ...) extract_par(THETA = OBJ@theta, C = sum(OBJ@dims@cat), P = OBJ@dims@p, Q = OBJ@dims@q, CONSTRMAT = OBJ@cnstr@loadings, CORRFLAG = OBJ@cnstr@corrflag, ...))
#' Extract parameters
#'
#' @param OBJ Raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
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
extract_par <- function(THETA, OPTION = c('list', 'raw', 'transformed'), C, P, Q, CONSTRMAT, CORRFLAG){
  OPTION <- match.arg(OPTION)

  if(OPTION == 'raw'){
    return(THETA)
    }else if(OPTION == 'transformed'){
      ncorr <- 0; if(CORRFLAG) ncorr <- Q*(Q-1)/2
      if(CORRFLAG){
        return(c(THETA[1:(C-P+sum(is.na(CONSTRMAT)))], get_corr(THETA = THETA, Q = Q, CORRFLAG=CORRFLAG)))
      }else{
        return(c(THETA[1:(C-P+sum(is.na(CONSTRMAT)))]))
      }

    }else if(OPTION == 'list') {
    out <- list()
    thr = THETA[1:(C-P)]
    ld <- cpp_get_loadings_theta2mat(
      THETA = THETA,
      CONSTRMAT = CONSTRMAT,
      P = P,
      C = C,
      NLOAD = sum(is.na(CONSTRMAT))
    )
    S <- cpp_get_latvar_theta2mat(
      THETA = THETA,
      Q=Q,
      D=length(THETA),
      CORRFLAG = CORRFLAG
    )

    out <- list(thresholds = thr, loadings = ld, latent_correlations = S)
    return(out)
  }
}

#' Compute variance of the estimates
#' @param OBJ Object of class plFaFit
#' @param DATA Original data
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param OPTION `transformed` if correlations are of interest. `raw` for inference on the reparametrised level.
#' @export
setGeneric('computeVar', function(OBJ, DATA, NUMDERIV = F, OPTION = 'transformed') standardGeneric('computeVar'), signature = 'OBJ')
#' Compute variance of the estimates
#' @param OBJ Object of class plFaFit
#' @param DATA Original data
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param OPTION `transformed` if correlations are of interest. `raw` for inference on the reparametrised level.
setMethod('computeVar', 'PlFaFit',
          function(OBJ, DATA, NUMDERIV, OPTION){ if(OBJ@method == 'ucminf'){
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
                       NUMDERIV  = NUMDERIV,
                       OPTION    = OPTION,
                       INVHAPPRX = OBJ@numFit$invhessian)
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
                       NUMDERIV  = NUMDERIV,
                       OPTION    = OPTION,
                       INVHAPPRX = OBJ@numFit$invhessian)
          }}
)

#setMethod('computeVar', 'vector', function(OBJ, DATA, ...) compute_var(OBJ, DATA, ...))
compute_var <- function(THETA, C_VEC, N, IT = NULL, PAIRS = NULL, PPI = NULL,  CONSTRMAT, CORRFLAG, FREQ, DATA, METHOD, NUMDERIV = F, OPTION = 'transformed', INVHAPPRX=NULL){

  Rwr_getPar <- function(par_vec){
    getPar(par_vec, OPTION = OPTION, C = sum(C_VEC),
           P = nrow(CONSTRMAT),
           Q = ncol(CONSTRMAT),
           CONSTRMAT = CONSTRMAT,
           CORRFLAG = CORRFLAG)
  }
  trJacob <- numDeriv::jacobian(Rwr_getPar, THETA)

  opt_noise <- NA
  Hhat <- NA
  if(NUMDERIV){
    message('1. Computing H numerically...')
    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <-cpp_multiThread_completePairwise(
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
    if(is.null(INVHAPPRX)){
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

  }

  message('2. Estimating J...')
  Jhat <- estimate_J(
    Y = DATA,
    C_VEC = C_VEC,
    A = CONSTRMAT,
    THETA = THETA,
    CORRFLAG = CORRFLAG
  )$est_J

  invH <- INVHAPPRX
  if(is.null(INVHAPPRX)){
    message('3. Inverting H...')
    if(!(det(Hhat)>0)) stop('H not invertible.')
    invH <- solve(Hhat)
  }

  sandwich <- invH %*% Jhat %*% invH

  message('3. Computing the variances...')
  vcov <- (trJacob %*% sandwich %*% t(trJacob)) / N
  if(METHOD =='ucminf'){
    asy_var <- diag(vcov)
  } else {
    asy_var <- diag(vcov)
    a1 <- PAIRS*(PAIRS - PPI)/(PPI*(PAIRS-1))
    a2 <- (PAIRS-PPI)/(PPI*(PAIRS-1))
    opt_noise <- diag( trJacob %*% (invH%*%(a1*Hhat - a2*Jhat)%*%invH/(N*IT) ) %*% t(trJacob) )
  }
  message('Done!')

  return(
    list(
      trJacob = trJacob,
      H = Hhat,
      J = Jhat,
      vcov = vcov,
      asymptotic_variance = asy_var,
      optimisation_noise = opt_noise
    )
  )

}

