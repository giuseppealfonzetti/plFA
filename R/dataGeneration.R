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
#'@export
build_constrMat <- function(P, Q, STRUCT = 'triangular'){
  if(STRUCT == 'triangular'){
    # Build lower trinagular matrix
    constrMat <- matrix(1, nrow = P, ncol = Q)
    for(j in 1:p){
      for (h in 1:q) {
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

#'@export
# generate normal random latent scores [nxq]
gen_scores <- function(N, Q, RHO = 0, SEED = 123){
  set.seed(SEED)
  cov_mat <- matrix(RHO, nrow = Q, ncol = Q)
  diag(cov_mat) <- 1

  out <-  mvtnorm::rmvnorm(N, mean = rep(0, nrow(cov_mat)), sigma = cov_mat)
  return(out)
}

#'@export
# generate the matrix of loadings [pxq] drawing from normal rv
gen_loadings <- function(FIXED = NULL, CONSTRAINT_MAT = NULL, SEED = 123){
  set.seed(SEED)
  p <- nrow(CONSTRAINT_MAT);
  q <- ncol(CONSTRAINT_MAT)

  out <-  CONSTRAINT_MAT
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



#'@export
gen_URV_data <- function(SAMPLE_SIZE, LOADINGS, THRESHOLDS, LATENT_COV, SEED = 123){

  set.seed(SEED)

  p <- nrow(LOADINGS)
  error_variance <- diag(1, p, p) - diag(diag(LOADINGS%*%LATENT_COV%*%t(LOADINGS)),p,p)

  errors <- mvtnorm::rmvnorm(n = SAMPLE_SIZE, sigma = error_variance)
  factors <- mvtnorm::rmvnorm(n = SAMPLE_SIZE, sigma = LATENT_COV)

  URV <- factors %*% t(LOADINGS) + errors

  out <- matrix(0, SAMPLE_SIZE, p)
  for(i in 1:length(THRESHOLDS)){
    # out <- purrr::modify2(out, URV, ~ifelse(.y > THRESHOLDS[i], .x + 1, .x))
    out[URV>THRESHOLDS[i]] <- i
  }

  return(out)
}
