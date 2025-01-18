#' @export
get_S <- function(THETA, Q, CORRFLAG){
  THETA <- check_theta(THETA)

  return(cpp_get_latvar_theta2mat(THETA, Q, length(THETA), CORRFLAG))
}



#' Extract vector of correlations from parameter vector
#'
#' get_corr() extracts the latent correlations from the
#' parameter vector. The vector has length \eqn{q*(q-1)/2},
#' where \eqn{q} is the number of latent variables.
#'
#' @param Q Number of latent variables.
#' @param THETA Parameter vector
#' @export
get_corr <- function(THETA, Q, CORRFLAG){
  Q <- check_q(Q)
  THETA <- check_theta(THETA)
  # if(length(THETA)<(Q*(Q-1)/2))stop('Check THETA dimensions.')

  s <- get_S(THETA = THETA, Q = Q, CORRFLAG = CORRFLAG)
  out <- s[upper.tri(s)]

  return(out)
}

#' Extract free loadings from parameter vector
#'
#' get_lambda() extracts the vector of free loadings
#' from the parameter vector
#'
#' @param THETA Raw vector of parameters.
#' @param C Sum of the number of categories for each item.
#' @param P Number if items
#' @param Q Number of latent variables
#'
#'@export
get_lambda <- function(THETA, C, P, Q){
  THETA <- check_theta(THETA)
  P     <- check_p(P)
  Q     <- check_q(Q)
  if(!is.finite(C))stop('C not numeric.')


  d <- length(THETA)
  ncorr <- Q*(Q-1)/2

  if(length(THETA)<=(C-P+ncorr))stop('Check THETA dimensions.')
  if(C<(2*P))stop('Check THETA dimensions.')

  lambda <- THETA[(C-P+1):(d-ncorr)]
  lambda
}

#' Get full parameter vector
#'
#' Assemble the parameter vector starting from the threshold vector,
#' the latent covariance matrix and the loading matrix
#'
#' @param TAU Numerical vector of thresholds parameters.
#' If all the \eqn{p} items have \eqn{c} possible categories, it contains
#' \eqn{p(c-1)} elements.
#' @param LOADINGS Numerical matrix of dimension \eqn{p*q},
#' where \eqn{q} is the number of latent variables.
#' @param LATENT_COV Correlation matrix of dimension \eqn{q}.
#' @param CAT Integer vector storing the possible number of categories
#' for each item. Values must be ordered following items order in the dataset.
#' @param CONSTRMAT Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
#' that the corresponding element in the loading matrix is free to be estimated.
#' A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
#' to \eqn{0}.
#' @param CORRFLAG Logical flag indicating whether the latent covariance matrix
#'
#' @export
get_theta <- function(TAU, LOADINGS, LATENT_COV=NULL, CAT, CONSTRMAT, CORRFLAG){

  TAU <- check_thresholds(TAU, CAT)
  LOADINGS <- check_loadings(LOADINGS)
  CONSTRMAT <- check_cnstr_loadings(CONSTRMAT)

  if(CORRFLAG)LATENT_COV <- check_latcov(LATENT_COV)

  stopifnot(dim(LOADINGS)==dim(CONSTRMAT))
  stopifnot(ncol(LOADINGS)==ncol(LATENT_COV))
  stopifnot(length(TAU)==(sum(CAT)-nrow(CONSTRMAT)))



  load_vec <- cpp_get_loadings_mat2vec(LOADINGS, CONSTRMAT, sum(is.na(CONSTRMAT)))

  corr_vec <- NULL
  if(CORRFLAG)  corr_vec <- cpp_get_latvar_mat2vec(LATENT_COV)
  theta <- c(TAU, load_vec, corr_vec)
  return(theta)
}


#' Extract estimation trajectories
#'
#' Extract trajectories of interest along the stochastic optimisation
#' from a fitted object
#'
#' @param MOD_OBJ Output of [fit_plFA].
#' @param PATH_LAB Label for the trajectory of interest. Use "path_av_theta"
#' for the averaged parameters trajectories.
#'
#'@export
get_tidy_path <- function(MOD_OBJ, PATH_LAB){
  iters <- MOD_OBJ$iterations_subset
  path  <- MOD_OBJ$fit[[PATH_LAB]]

  out <- dplyr::mutate(dplyr::tibble(iter = iters),
      path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
    )

  colnames(out) <- c('iter', PATH_LAB)

  return(out)
}

