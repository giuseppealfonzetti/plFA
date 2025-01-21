check_data <- function(DATA){
  stopifnot(is.matrix(DATA))
  stopifnot(is.finite(DATA))
  stopifnot(is.numeric(DATA))
  stopifnot(min(DATA)==0)
  return(round(DATA,0))
}
check_cnstr_corrflag <- function(CORRFLAG){
  if(is.null(CORRFLAG)) CORRFLAG <- 1
  stopifnot(CORRFLAG%in%c(0,1))
  return(CORRFLAG)
}
check_cnstr_stdlv <- function(STDLV){
  if(is.null(STDLV)) STDLV <- TRUE
  stopifnot(is.logical(STDLV))
  return(STDLV)
}
check_cnstr_loadings <- function(MAT, STDLV){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(!is.infinite(sum(MAT, na.rm = TRUE)))

  if(!STDLV){
    for(j in 1:ncol(MAT)){
      wna <- which(is.na(MAT[,j]))
      fna <- wna[1]
      if(is.null(fna)) break
      if(fna==1){MAT[fna,j] <- 1}else{
        if(all(MAT[1:(fna-1),j]==0))MAT[fna,j] <- 1
      }
    }
  }

  return(MAT)
}
check_cnstr_latvar <- function(VEC, Q, STDLV){
  if(is.null(VEC) & STDLV) VEC <- rep(1, Q)
  if(is.null(VEC) & !STDLV) VEC <- rep(NA, Q)

  stopifnot(is.vector(VEC))
  stopifnot(!is.infinite(sum(VEC, na.rm = TRUE)))
  stopifnot(VEC[!is.na(VEC)]>0)

  if(STDLV){
    VEC[is.na(VEC)] <- 1
  }

  return(log(sqrt(VEC)))
}
check_cnstr <- function(LIST){
  ls <- list()
  ls$CORRFLAG <- check_cnstr_corrflag(LIST$CORRFLAG)
  ls$STDLV    <- check_cnstr_stdlv(LIST$STDLV)
  ls$CONSTRMAT <- check_cnstr_loadings(LIST$CONSTRMAT, STDLV = ls$STDLV)
  q <- ncol(ls$CONSTRMAT)
  ls$CONSTRLOGSD <- check_cnstr_latvar(LIST$CONSTRVAR, Q=q, STDLV = ls$STDLV)

  stopifnot(q==length(ls$CONSTRLOGSD))
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

  nthr <- sum(categories)-p
  nload <- sum(is.na(CONSTR_LIST$CONSTRMAT))
  ncorr <- 0; if(CONSTR_LIST$CORRFLAG==1) ncorr <- q*(q-1)/2
  nvar  <- sum(is.na(CONSTR_LIST$CONSTRLOGSD))
  d <- nthr + nload + ncorr + nvar
  stopifnot(q<p)
  stopifnot(categories>1)




  return(list(
    n=n,
    p=p,
    q=q,
    nthr=nthr,
    nload=nload,
    ncorr=ncorr,
    nvar=nvar,
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
    transformed_latsd_init <- init_transformed_latsd(DIMS, CONSTR_LIST)
    INIT <- c(lambda0_init, lambda_init, transformed_rhos_init, transformed_latsd_init)
   }else{
    stopifnot(is.vector(INIT))
    stopifnot(is.numeric(INIT))
    stopifnot(is.finite(INIT))
    stopifnot(length(INIT)==DIMS$d)
    message('1. Initialising at INIT vector.')
   }

  stopifnot(length(INIT)==DIMS$d)
  return(INIT)




}
