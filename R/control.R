check_data <- function(DATA){
  stopifnot(is.matrix(DATA))
  stopifnot(is.finite(DATA))
  stopifnot(is.numeric(DATA))
  stopifnot(min(DATA)==0)
  return(round(DATA,0))
}

check_cnstr_loadings <- function(MAT){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(!is.infinite(sum(MAT, na.rm = TRUE)))
  return(MAT)
}

check_cnstr <- function(LIST){
  ls <- list()
  ls$CONSTRMAT <- check_cnstr_loadings(LIST$CONSTRMAT)
  if(is.null(LIST$CORRFLAG)) ls$CORRFLAG <- 1
  stopifnot(LIST$CORRFLAG%in%c(0,1))
  ls$CORRFLAG <- LIST$CORRFLAG
  return(ls)
}

check_p <- function(P){
  stopifnot(is.numeric(P))
  stopifnot(is.finite(P))
  stopifnot(P>1)
  return(as.integer(P))
}

check_q <- function(Q){
  stopifnot(is.numeric(Q))
  stopifnot(is.finite(Q))
  stopifnot(Q>0)
  return(as.integer(Q))
}

check_n <- function(N){
  stopifnot(is.numeric(N))
  stopifnot(is.finite(N))
  stopifnot(N>0)
  return(as.integer(N))
}

check_loadings <- function(MAT){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(is.finite(MAT))
  stopifnot(nrow(MAT)>ncol(MAT))
  return(MAT)
}

check_thresholds_j <- function(TH){
  stopifnot(is.vector(TH))
  stopifnot(is.numeric(TH))
  stopifnot(!is.unsorted(TH))
  return(TH)
}

check_thresholds <- function(THF, CAT){
  s <- 0
  for(j in 1:length(CAT)){
    check_thresholds_j(THF[(s + 1):(s + CAT[j] - 1)])
    s <- s + CAT[j] - 1
  }
  return(THF)
}

check_latcov <- function(MAT){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(is.finite(MAT))
  stopifnot(nrow(MAT)==ncol(MAT))
  stopifnot(isSymmetric.matrix(MAT))
  #check also positive definiteness
  return(MAT)
}

check_dims <- function(DATA, CONSTR_LIST){
  p <- ncol(DATA)
  n <- nrow(DATA)
  q <- ncol(CONSTR_LIST$CONSTRMAT)
  categories <- apply(DATA, 2, max, na.rm = T) + 1
  d <- sum(categories)-p + sum(is.na(CONSTR_LIST$CONSTRMAT))
  if(CONSTR_LIST$CORRFLAG==1) d <- d + q*(q-1)/2
  stopifnot(q<p)
  stopifnot(categories>1)

  return(list(
    n=n,
    p=p,
    q=q,
    d=d,
    cat=categories,
    pairs=p*(p-1)/2
  ))
}

check_theta <- function(THETA){

  stopifnot(is.vector(THETA))
  stopifnot(is.numeric(THETA))
  stopifnot(is.finite(THETA))

  return(THETA)
}


check_init <- function(INIT, DIMS, CONSTR_LIST){
  if(is.null(INIT)){
    message('1. Initialising at default values')

    # initialise thresholds
    lambda0_init <- init_thresholds(DIMS, CONSTR_LIST)
    lambda_init  <- init_loadings(DIMS, CONSTR_LIST)
    transformed_rhos_init  <- init_transformed_latcorr(DIMS, CONSTR_LIST)
    return(c(lambda0_init, lambda_init, transformed_rhos_init))

  }else{
    stopifnot(is.vector(INIT))
    stopifnot(is.numeric(INIT))
    stopifnot(is.finite(INIT))
    stopifnot(length(INIT)==DIMS$d)
    message('1. Initialising at INIT vector.')
    return(INIT)
  }

}
