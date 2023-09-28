#' Extract vector of correlations from parameter vector
#'
#' get_corr() extracts the latent correlations from the
#' parameter vector. The vector has length \eqn{q*(q-1)/2},
#' where \eqn{q} is the number of latent variables.
#'
#' @param CONSTRMAT Matrix of dimension \eqn{p * q},
#' with \eqn{p} the number of items.
#' @param Q Number of latent variables.
#'
#' @export
get_corr <- function(THETA, Q){
  if(sum(!is.finite(THETA))>0)stop('THETA not numeric.')
  if(!is.finite(Q))stop('Q not numeric.')
  if(length(THETA)<(Q*(Q-1)/2))stop('Check THETA dimensions.')

  s <- get_S(THETA = THETA, Q = Q)
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
  if(sum(!is.finite(THETA))>0) stop('THETA not numeric.')
  if(!is.finite(C))stop('C not numeric.')
  if(!is.finite(P))stop('P not numeric.')
  if(!is.finite(Q))stop('Q not numeric.')

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
#' @param A Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
#' that the corresponding element in the loading matrix is free to be estimated.
#' A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
#' to \eqn{0}.
#'
#' @export
get_theta <- function(TAU, LOADINGS, LATENT_COV, CAT, CONSTRMAT){

  if(sum(!is.finite(TAU))>0)stop('TAU not numeric.')
  if(sum(!is.finite(LOADINGS))>0)stop('LOADINGS not numeric.')
  if(sum(!is.finite(LATENT_COV))>0)stop('LATENT_COV not numeric.')
  if(sum(!is.finite(CAT))>0)stop('CAT not numeric.')
  if(!matrixcalc::is.positive.definite(LATENT_COV))stop('Latent correlation matrix not positive definite!')

  if(sum(dim(LOADINGS)==dim(CONSTRMAT))<2)stop('LOADINGS and A dimensions do not match.')
  if(ncol(LOADINGS)!=ncol(LATENT_COV))stop('LOADINGS and LATENT_COV dimensions do not match.')
  if(length(TAU)!=(sum(CAT)-nrow(CONSTRMAT)))stop('TAU and CAT dimensions do not match.')

  load_vec <- c()
  s <- 1
  for (j in 1:ncol(LOADINGS)) {
    for (i in 1:nrow(LOADINGS)) {
      if(CONSTRMAT[i,j]!=0){
        load_vec[s] <- LOADINGS[i,j]
        s = s+1
      }
    }
  }

  corr_vec <- get_par_from_S(LATENT_COV)
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

  out <- dplyr::tibble(iter = iters) %>%
    dplyr::mutate(
      path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
    )

  colnames(out) <- c('iter', PATH_LAB)

  return(out)
}

# Custom call to pvnorm to pass to rcpp functions
pmvnorm_custom <- function(lower, upper, sigma, seed = 123){
  set.seed(seed)
  out <- pmvnorm(lower = lower, upper = upper, sigma = sigma, keepAttr = F)
  return(out)
}
