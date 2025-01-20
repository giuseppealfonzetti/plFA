#' Show method for the output of fit_plFA()
#' @param object Object of class StoFit
setMethod('show', 'PlFaFit', function(object){

  cat( "- Dimensions:\n",
       "  - Sample size: ", object@dims@n, "\n",
       "  - Items: ", object@dims@p,  " (", object@dims@pairs, " pairs)\n",
       "  - Latent traits: ", object@dims@q,"\n\n",
       "- Free parameters:\n",
       "  - Thresholds: ", object@dims@nthr, "\n",
       "  - Loadings: ", object@dims@nload, "\n",
       "  - Latent correlations: ", object@dims@ncorr, "\n",
       "  - Latent variances: ", object@dims@nvar, "\n",
       "  - Total: ", object@dims@npar, "\n\n")

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
setMethod('getPar', 'PlFaFit', function(OBJ, ...) extract_par(THETA       = OBJ@theta,
                                                              CONSTRMAT   = OBJ@cnstr@loadings,
                                                              CONSTRLOGSD = OBJ@cnstr@loglatsd,
                                                              NTHR        = OBJ@dims@nthr,
                                                              NLOAD       = OBJ@dims@nload,
                                                              NCORR       = OBJ@dims@ncorr,
                                                              NVAR        = OBJ@dims@nvar,
                                                              ...))
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
extract_par <- function(THETA, OPTION = c('list', 'raw', 'transformed'), CONSTRMAT, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR){
  OPTION <- match.arg(OPTION)

  q <- ncol(CONSTRMAT)
  if(OPTION == 'raw'){
    return(THETA)
    }else if(OPTION == 'transformed'){

      vec <- THETA[1:(NTHR+NLOAD)]
      if(NCORR>0){
        vec <- c(vec, get_corr(THETA = THETA,
                               NTHR  = NTHR,
                               NLOAD = NLOAD,
                               NCORR = NCORR,
                               NVAR  = NVAR,
                               Q     = q))
      }
      if(NVAR>0){
        vec <- c(vec, exp(THETA[(NTHR+NLOAD+NCORR+1):(NTHR+NLOAD+NCORR+NVAR)]))
      }

      return(vec)


    }else if(OPTION == 'list') {
    out <- list()
    thr = THETA[1:NTHR]
    ld <- cpp_loadings_theta2mat(
      THETA     = THETA,
      CONSTRMAT = CONSTRMAT,
      NTHR      = NTHR,
      NLOAD     = NLOAD)

    S <- cpp_latvar_theta2mat(
      THETA       = THETA,
      CONSTRLOGSD = CONSTRLOGSD,
      NTHR        = NTHR,
      NLOAD       = NLOAD,
      NCORR       = NCORR,
      NVAR        = NVAR,
      Q           = q
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
            compute_var(THETA      = OBJ@theta,
                       C_VEC       = OBJ@dims@cat,
                       N           = OBJ@dims@n,
                       IT          = NULL,
                       PAIRS       = NULL,
                       PPI         = NULL,
                       CONSTRMAT   = OBJ@cnstr@loadings,
                       CONSTRLOGSD = OBJ@cnstr@loglatsd,
                       NTHR        = OBJ@dims@nthr,
                       NLOAD       = OBJ@dims@nload,
                       NCORR       = OBJ@dims@ncorr,
                       NVAR        = OBJ@dims@nvar,
                       FREQ        = OBJ@freq,
                       DATA        = DATA,
                       METHOD      = OBJ@method,
                       NUMDERIV    = NUMDERIV,
                       INVHAPPRX   = OBJ@numFit$invhessian)
          }else{
            compute_var(THETA    = OBJ@theta,
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
compute_var <- function(THETA, C_VEC, N, IT = NULL, PAIRS = NULL, PPI = NULL,
                        CONSTRMAT, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR,
                        FREQ, DATA, METHOD, NUMDERIV = F, INVHAPPRX=NULL){

  Rwr_getPar <- function(x){
    getPar(OBJ         = x,
           OPTION      = "transformed",
           CONSTRMAT   = CONSTRMAT,
           CONSTRLOGSD = CONSTRLOGSD,
           NTHR        = NTHR,
           NLOAD       = NLOAD,
           NCORR       = NCORR,
           NVAR        = NVAR)
  }
  trJacob <- numDeriv::jacobian(Rwr_getPar, x=THETA)

  opt_noise <- NA
  Hhat <- NA
  if(NUMDERIV){
    message('- Computing H numerically...')
    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = N,
        C_VEC      = C_VEC,
        CONSTRMAT  = CONSTRMAT,
        CONSTRLOGSD= CONSTRLOGSD,
        FREQ       = FREQ,
        THETA      = par_vec,
        CORRFLAG   = (NCORR>0),
        NTHR       = NTHR,
        NLOAD      = NLOAD,
        NCORR      = NCORR,
        NVAR       = NVAR,
        GRFLAG     = 1,
        SILENTFLAG = 1
      )}

    Hhat <- numDeriv::jacobian(Rwr_ngr, THETA)
  }else{
    if(is.null(INVHAPPRX)){
      message('- Estimating H...')
      # Hhat <- estimate_H(
      #   C_VEC = C_VEC,
      #   A = CONSTRMAT,
      #   THETA = THETA,
      #   FREQ = FREQ,
      #   N = N,
      #   CORRFLAG = CORRFLAG
      # )$est_H
    }

  }

  message('- Estimating J...')
  Jhat <- estimate_J(
    Y           = DATA,
    C_VEC       = C_VEC,
    A           = CONSTRMAT,
    CONSTRLOGSD = CONSTRLOGSD,
    THETA       = THETA,
    CORRFLAG    = as.numeric((NCORR>0)),
    NTHR        = NTHR,
    NLOAD       = NLOAD,
    NCORR       = NCORR,
    NVAR        = NVAR
  )$est_J

  invH <- INVHAPPRX
  if(is.null(INVHAPPRX)){
    message('3. Inverting H...')
    if(!(det(Hhat)>0)) stop('H not invertible.')
    invH <- solve(Hhat)
  }

  sandwich <- invH %*% Jhat %*% invH

  message('- Computing the variances...')
  if(METHOD =='ucminf'){
    asy_var <- diag((trJacob%*% sandwich %*% t(trJacob))/N)
  }else{
    asy_var <- diag((trJacob%*% sandwich %*% t(trJacob))/N)
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
      invH = invH,
      asymptotic_variance = asy_var,
      optimisation_noise = opt_noise
    )
  )

}

