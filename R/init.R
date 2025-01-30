init_thresholds <- function(DIMS, CONSTR_LIST){
  out <- c()
  s <- 0
  for (i in 1:length(DIMS$cat)) {
    if(DIMS$cat[i]==2){
      out[(s + 1):(s + DIMS$cat[i] - 1)] <- 0
    }else{
      vec <- 1:(DIMS$cat[i]-1)
      vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
      out[(s + 1):(s + DIMS$cat[i] - 1)] <- vec
    }
    s <- s + DIMS$cat[i] - 1
  }
  return(out)
}

init_loadings <- function(DIMS, CONSTR_LIST, FXD=0.25){
  out <-  rep(FXD, sum(is.na(CONSTR_LIST$CONSTRMAT)))
  return(out)
}

init_transformed_latcorr <- function(DIMS, CONSTR_LIST){

  out <- NULL
  if(CONSTR_LIST$CORRFLAG==1) out <- rep(0, DIMS$q*(DIMS$q-1)/2)
  return(out)
}

init_transformed_latsd <- function(DIMS, CONSTR_LIST){
  rep(0, sum(is.na(CONSTR_LIST$CONSTRLOGSD)))
}

init_par_with_plSA <- function(DIMS, CONSTR_LIST, CPP_ARGS=NULL){
  lambda0_init <- init_thresholds(DIMS, CONSTR_LIST)
  lambda_init  <- init_loadings(DIMS, CONSTR_LIST)
  transformed_rhos_init  <- init_transformed_latcorr(DIMS, CONSTR_LIST)
  transformed_latsd_init <- init_transformed_latsd(DIMS, CONSTR_LIST)
  start_par <- c(lambda0_init, lambda_init, transformed_rhos_init, transformed_latsd_init)
  start_par <- cpp_sa_proj(THETA=start_par,
                          CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                          CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
                          LLC = CONSTR_LIST$LLC,
                          C_VEC=DIMS$cat,
                          CORRFLAG=CONSTR_LIST$CORRFLAG,
                          NTHR=DIMS$nthr,
                          NLOAD=DIMS$nload,
                          NCORR=DIMS$ncorr,
                          NVAR=DIMS$nvar)

  cpp_args <- c(list(N           = DIMS$n,
                     C_VEC       = DIMS$cat,
                     CONSTRMAT   = CONSTR_LIST$CONSTRMAT,
                     CONSTRLOGSD = CONSTR_LIST$CONSTRLOGSD,
                     LLC         = CONSTR_LIST$LLC,
                     THETA_INIT  = start_par,
                     NTHR        = DIMS$nthr,
                     NLOAD       = DIMS$nload,
                     NCORR       = DIMS$ncorr,
                     NVAR        = DIMS$nvar
  ), CPP_ARGS)

  sto_init <- do.call(cpp_plSA, cpp_args)


  return(sto_init$avtheta)
}
