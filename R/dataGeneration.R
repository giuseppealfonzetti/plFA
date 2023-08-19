#' Compute pairwise frequencies
#'
#' It returns a 5-rows matrix, with one column for each combination of items and
#' categories.
#' Row0: item k, Row1: item l, Row2; category item k, Row3: category item l, Row4: freq
#'
#' @param Y Integer matrix of dimension \eqn{n*p}, where \eqn{n} is the sample size
#' and \eqn{p} is the number of items considered. Categories must be coded starting from zero.
#' For example, an item with three categories can only accept values contained in
#' \eqn{{0, 1, 2}}.
#' @param C_VEC Integer vector indicating how many possible categories are associated to
#' each item in 'Y'.
#'
#' @export
compute_frequencies <- function(Y, C_VEC){
  if(sum(!is.finite(Y))>0 | !is.matrix(Y))stop('Y is not a numeric matrix.')
  if(sum(!is.finite(C_VEC))>0 | is.matrix(C_VEC))stop('C_VEC is not a numeric vector.')
  if(length(C_VEC) != ncol(Y))stop('Y and C_VEC dimensions do not match.')

  pairs_freq(Y = Y, C_VEC = C_VEC)
}

#' Automatic construction of simple constraints matrices
#'
#' build_constrMat() creates a binary matrix with entries equal to 1
#' when referring to free loadings, and 0 entries otherwise.
#'
#' @param P Number of items.
#' @param Q Number of latent vaiables.
#' @param STRUCT Label for the constraint structure chosen. Only two possible values:
#' 'triangular' and 'simple'.
#'
#' @export
build_constrMat <- function(P, Q, STRUCT = 'triangular'){
  if(!is.finite(P))stop('P not numeric.')
  if(!is.finite(Q))stop('Q not numeric.')
  if(!(STRUCT %in% c('simple', 'triangular')))stop('STRUCT not available.')
  if(STRUCT == 'triangular'){
    # Build lower trinagular matrix
    constrMat <- matrix(1, nrow = P, ncol = Q)
    for(j in 1:P){
      for (h in 1:Q) {
        if(h > j) constrMat[j,h] <- 0
      }
    }

  } else if ( STRUCT == 'simple'){
    # Each item to one factor
    loadings_per_factor <- P %/% Q
    remaining_loadings <- P %% Q

    constrMat <- matrix(0, nrow = P, ncol = Q)
    for (h in 1:(Q-1)) {
      constrMat[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor), h] <- rep(1,loadings_per_factor)
    }
    h <- Q
    constrMat[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor + remaining_loadings), h] <- rep(1,loadings_per_factor + remaining_loadings)

  }
  return(constrMat)
}


#' Construct loading matrix
#'
#' gen_loadings() construct a (possibly) random loading matrix following the
#' constraints passed.
#'
#' @param FIXED Fixed value to assign to all free loadings. If 'NULL' it draws
#' them randomly from Unif(0,1)
#' @param CONSTRMAT Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
#' that the corresponding element in the loading matrix is free to be estimated.
#' A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
#' to \eqn{0}.
#' @param SEED Random seed.
#'
#' @export
# generate the matrix of loadings [pxq] drawing from normal rv
gen_loadings <- function(FIXED = NULL, CONSTRMAT, SEED = 123){
  if(!is.matrix(CONSTRMAT))stop('CONSTRMAT must be a matrix')

  set.seed(SEED)
  p <- nrow(CONSTRMAT);
  q <- ncol(CONSTRMAT)

  out <-  CONSTRMAT
  for (j in 1:p) {
    for (h in 1:q) {
      if(out[j,h] != 0) out[j,h] <- runif(1, 0,1)
    }
  }

  if(is.null(FIXED)==FALSE){
    for (j in 1:p) {
      for (h in 1:q) {
        if(out[j,h] != 0) out[j,h] <- FIXED
      }
    }
  }

  return(out)
}

#' Simulate new data
#'
#' Simulate data given the values of the parameters for thresholds, loadings and latent correlations
#'
#' @param SAMPLE_SIZE Number of observations to simulate.
#' @param LOADINGS Matrix of loadings, with dimensions \eqn{p*q}.
#' @param THRESHOLDS Vector of thresholds fo a single item.
#' It is assumed to be the same for all the items.
#' @param LATENT_COV Latent correlation matrix. Dimension \eqn{q*q}.
#' @param SEED Numeric seed for random generation.
#'
#'@export
sim_data <- function(SAMPLE_SIZE, LOADINGS, THRESHOLDS, LATENT_COV, SEED = 123){
  if(!is.finite(SAMPLE_SIZE) | length(SAMPLE_SIZE)!=1 )stop('SAMPLE_SIZE must be a positive integer.')
  if(SAMPLE_SIZE<=0 | round(SAMPLE_SIZE)!=SAMPLE_SIZE )stop('SAMPLE_SIZE must be a positive integer.')
  if(sum(!is.finite(LOADINGS))!=0 | !is.matrix(LOADINGS))stop('LOADINGS is not a numeric matrix.')
  if(sum(!is.finite(THRESHOLDS))!=0 | !is.vector(THRESHOLDS))stop('THRESHOLDS is not a numeric vector.')
  if(sum(!is.finite(LATENT_COV))!=0 | !is.matrix(LATENT_COV))stop('LATENT_COV is not a numeric matrix.')

  if(!isSymmetric.matrix(LATENT_COV))stop('LATENT_COV not symmetric.')
  if(!matrixcalc::is.positive.definite(LATENT_COV))stop('LATENT_COV not positive definite.')
  if(is.unsorted(THRESHOLDS))stop('THRESHOLDS must be sorted.')
  if(!(sum(round(diag(LATENT_COV))==1)==ncol(LATENT_COV)))stop('LATENT_COV not a correlation matrix.')
  if(ncol(LOADINGS)!=ncol(LATENT_COV))stop('LOADINGS and LATENT_COV dimensions not compatible.')

  set.seed(SEED)

  p <- nrow(LOADINGS)
  error_variance <- diag(1, p, p) - diag(diag(LOADINGS%*%LATENT_COV%*%t(LOADINGS)),p,p)

  errors <- mvtnorm::rmvnorm(n = SAMPLE_SIZE, sigma = error_variance)
  factors <- mvtnorm::rmvnorm(n = SAMPLE_SIZE, sigma = LATENT_COV)

  URV <- factors %*% t(LOADINGS) + errors

  out <- matrix(0, SAMPLE_SIZE, p)
  for(i in 1:length(THRESHOLDS)){
    out[URV>THRESHOLDS[i]] <- i
  }

  return(out)
}
