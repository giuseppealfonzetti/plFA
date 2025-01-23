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

check_init_par <- function(PAR, DIMS, CONSTR_LIST){
  stopifnot(is.vector(PAR))
  stopifnot(is.numeric(PAR))
  stopifnot(is.finite(PAR))
  stopifnot(length(PAR)==DIMS$d)

  L <- cpp_loadings_theta2mat(THETA     = PAR,
                              CONSTRMAT = CONSTR_LIST$CONSTRMAT,
                              NTHR      = DIMS$nthr,
                              NLOAD     = DIMS$nload)
  check_loadings(L)
  S <- cpp_latvar_theta2mat(THETA       = PAR,
                            CONSTRLOGSD = CONSTR_LIST$CONSTRLOGSD,
                            NTHR        = DIMS$nthr,
                            NLOAD       = DIMS$nload,
                            NCORR       = DIMS$ncorr,
                            NVAR        = DIMS$nvar,
                            Q           = DIMS$q)
  check_latcov(S)
  Sy <- L%*%S%*%t(L)
  # stopifnot(diag(Sy)<=1)
  # if(!all(diag(Sy)<=1))warning("Negative error variances implied by init parameters")
}

check_init <- function(INIT, FREQ, DIMS, CONSTR_LIST, INIT_METHOD=c("plSA", "custom", "standard"), VERBOSE=FALSE, CPP_ARGS=NULL){

  if(!is.null(INIT)) INIT_METHOD <- "custom"
  INIT_METHOD <- match.arg(INIT_METHOD)

  if(INIT_METHOD=="custom"){
    if(VERBOSE) cat('- Initialised at INIT vector.\n')
  }else if(INIT_METHOD=="standard"){
    if(VERBOSE) cat('- Initialising at default values ...')
    lambda0_init <- init_thresholds(DIMS, CONSTR_LIST)
    lambda_init  <- init_loadings(DIMS, CONSTR_LIST)
    transformed_rhos_init  <- init_transformed_latcorr(DIMS, CONSTR_LIST)
    transformed_latsd_init <- init_transformed_latsd(DIMS, CONSTR_LIST)
    INIT <- c(lambda0_init, lambda_init, transformed_rhos_init, transformed_latsd_init)
    INIT <- cpp_sa_proj(THETA=INIT,
                        CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                        CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
                        C_VEC=DIMS$cat,
                        CORRFLAG=CONSTR_LIST$CORRFLAG,
                        NTHR=DIMS$nthr,
                        NLOAD=DIMS$nload,
                        NCORR=DIMS$ncorr,
                        NVAR=DIMS$nvar)
    if(VERBOSE) cat('- Done.\n')
  }else if(INIT_METHOD=="plSA"){
    if(VERBOSE) cat('- Initialising with plSA ...')
    CPP_ARGS <- check_plSA_args(FREQ=FREQ, DIMS=DIMS, CONSTR_LIST=CONSTR_LIST, LIST=CPP_ARGS)
    INIT <- init_par_with_plSA(DIMS, CONSTR_LIST, CPP_ARGS)
    if(VERBOSE) cat('- Done.\n')

  }

  check_init_par(INIT, DIMS, CONSTR_LIST)
  return(INIT)
}





check_plSA_args <- function(FREQ, DIMS, CONSTR_LIST, LIST=NULL, SETTING=c("init", "main")){
  out <- list()
  out$FREQ <- out$VALFREQ <- FREQ

  if(is.null(LIST$SAMPLER)) LIST$SAMPLER <- 0
  stopifnot(LIST$SAMPLER%in%c(0,1))
  out$SAMPLER <- LIST$SAMPLER

  if(is.null(LIST$PAIRS_PER_ITERATION)) LIST$PAIRS_PER_ITERATION <- 16
  stopifnot(is.numeric(LIST$PAIRS_PER_ITERATION))
  stopifnot(LIST$PAIRS_PER_ITERATION>0)
  out$PAIRS_PER_ITERATION <- as.integer(LIST$PAIRS_PER_ITERATION)

  if(is.null(LIST$SCHEDULE)) LIST$SCHEDULE <- 1
  stopifnot(LIST$SCHEDULE%in%c(0,1))
  out$SCHEDULE <- LIST$SCHEDULE

  if(is.null(LIST$STEP0)) LIST$STEP0 <- sqrt(1/DIMS$d)
  stopifnot(is.numeric(LIST$STEP0))
  stopifnot(LIST$STEP0>0)
  out$STEP0 <- LIST$STEP0

  if(is.null(LIST$STEP1)) LIST$STEP1 <- 1
  stopifnot(is.numeric(LIST$STEP1))
  stopifnot(LIST$STEP1>0)
  out$STEP1 <- LIST$STEP1

  if(is.null(LIST$STEP2)) LIST$STEP2 <- 1e-5
  stopifnot(is.numeric(LIST$STEP2))
  stopifnot(LIST$STEP2>0)
  out$STEP2 <- LIST$STEP2

  if(is.null(LIST$STEP3)) LIST$STEP3 <- .75
  stopifnot(is.numeric(LIST$STEP3))
  stopifnot(LIST$STEP3>0)
  out$STEP3 <- LIST$STEP3

  if(is.null(LIST$BURN)) LIST$BURN <- 150
  stopifnot(is.numeric(LIST$BURN))
  stopifnot(LIST$BURN>0)
  out$BURN <- as.integer(LIST$BURN)

  if(is.null(LIST$MAXT)) LIST$MAXT <- 200
  stopifnot(is.numeric(LIST$MAXT))
  stopifnot(LIST$MAXT>0)
  out$MAXT <- as.integer(LIST$MAXT)

  if(is.null(LIST$TOL_WINDOW)) LIST$TOL_WINDOW <- 1
  stopifnot(is.numeric(LIST$TOL_WINDOW))
  stopifnot(LIST$TOL_WINDOW>0)
  out$TOL_WINDOW <- as.integer(LIST$TOL_WINDOW)

  if(is.null(LIST$TOL)) LIST$TOL <- 1e-8
  stopifnot(is.numeric(LIST$TOL))
  stopifnot(LIST$TOL>0)
  out$TOL <- LIST$TOL

  if(is.null(LIST$CHECK_TOL)) LIST$CHECK_TOL <- FALSE
  stopifnot(is.logical(LIST$CHECK_TOL))
  out$CHECK_TOL <- LIST$CHECK_TOL


  if(is.null(LIST$CHECK_WINDOW)) LIST$CHECK_WINDOW <- 250
  stopifnot(is.numeric(LIST$CHECK_WINDOW))
  stopifnot(LIST$CHECK_WINDOW>0)
  out$CHECK_WINDOW <- as.integer(LIST$CHECK_WINDOW)

  if(is.null(LIST$PATH_WINDOW)) LIST$PATH_WINDOW <- out$CHECK_WINDOW
  stopifnot(is.numeric(LIST$PATH_WINDOW))
  stopifnot(LIST$PATH_WINDOW>0)
  out$PATH_WINDOW <- as.integer(LIST$PATH_WINDOW)

  if(is.null(LIST$CLOCK_WINDOW)) LIST$CLOCK_WINDOW <- out$CHECK_WINDOW
  stopifnot(is.numeric(LIST$CLOCK_WINDOW))
  stopifnot(LIST$CLOCK_WINDOW>0)
  out$CLOCK_WINDOW <- as.integer(LIST$CLOCK_WINDOW)

  if(is.null(LIST$SEED)) LIST$SEED <- 123
  stopifnot(is.numeric(LIST$SEED))
  stopifnot(LIST$SEED>0)
  out$SEED <- as.integer(LIST$SEED)

  if(is.null(LIST$VERBOSE)) LIST$VERBOSE <- FALSE
  stopifnot(is.logical(LIST$VERBOSE))
  out$VERBOSE <- LIST$VERBOSE


  return(out)

}



































































