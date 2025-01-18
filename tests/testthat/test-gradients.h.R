set.seed(1)

# p = number of items, q = number of latent variables, n = number of observations
p <- 6; q <- 2; n <- 100

stdlv <- FALSE
corrflag <- TRUE

# Thresholds vector for each item
thr <- c(-.5, .5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

# Simple loading matrix constraints
A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
A <- check_cnstr_loadings(A, stdlv)

# Draw some random loadings according to constraints
Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
nload <- sum(is.na(A))

# Generate random latent correlation matrix
tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
ncorr <- if(corrflag)q*(q-1)/2 else 0
R <- cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


# Generate random latent variances
constr_var <- rep(NA, q)
constr_lsd <- check_cnstr_latvar(constr_var, q, stdlv)
nvar  <- sum(is.na(constr_lsd))
tsdvec <- constr_lsd
tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.5, .1)
Dmat <- diag(exp(tsdvec),q,q)

S <- Dmat %*% R %*% Dmat

# Dimensions
d <- nthr + nload + ncorr + nvar
d
#
theta <- get_theta(
  THRESHOLDS = rep(thr, p),
  LOADINGS = Load,
  LATENT_COV = S,
  CAT = cat,
  CONSTRMAT = A,
  CONSTRVAR = exp(constr_lsd)^2,
  CORRFLAG = corrflag,
  STDLV=stdlv
)
length(theta)
cpp_latvar_mat2vec(S=S,
                   CONSTRLOGSD = constr_lsd,
                   NCORR = ncorr,
                   NVAR = nvar)

dat <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)
f <- compute_frequencies(Y = dat, C_VEC = cat)

dat <- check_data(dat)
constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
dims <- check_dims(dat, constr_list)

lambda0_init <- init_thresholds(dims, constr_list)
lambda_init  <- init_loadings(dims, constr_list, FXD=0.5)
transformed_rhos_init <- init_transformed_latcorr(dims, constr_list)
transformed_sd_init <- init_transformed_latsd(dims, constr_list)

par_init <- c(lambda0_init, lambda_init, transformed_rhos_init, transformed_sd_init)
length(par_init)


pair_nll <- function(PAR, K, L){
  pair <- cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 0,
    OPTION      = 0
  )
  out <- pair$nll/n
  return(out)
}

pair_ngr <- function(PAR, K, L){
  pair <- cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 1,
    OPTION      = 0
  )
  out <- pair$ngr/n
  return(out)
}

pair_nll_old <- function(PAR, K, L){
  pair <- cpp_compute_pair(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 0,
    OPTION      = 0
  )
  out <- pair$ll/n
  return(out)
}
pair_ngr_old <- function(PAR, K, L){
  pair <- cpp_compute_pair(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    C_VEC       = cat,
    THETA       = PAR,
    CORRFLAG    = corrflag,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    SILENTFLAG  = 1,
    GRADFLAG    = 1,
    OPTION      = 0
  )
  out <- pair$ngr/n
  return(out)
}
k <- 2; l <- 1
pair_nll(theta, K=k, L=l)
pair_nll_old(theta, K=k, L=l)

numDeriv::grad(func = pair_nll, x=theta, K=k,L=l)
pair_ngr(PAR=theta, K=k, L=l)

numDeriv::grad(func = pair_nll_old, x=theta, K=k,L=l)
pair_ngr_old(PAR=theta, K=k, L=l)





















