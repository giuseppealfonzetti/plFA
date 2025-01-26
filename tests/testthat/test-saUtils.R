# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 2; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5,-0.5,.5, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

set.seed(1)
stdlv <- FALSE
corrflag <- TRUE

# Simple loading matrix constraints
A <- build_constrMat(P = p, Q = q, STRUCT = 'cross')
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
tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
Dmat <- diag(exp(tsdvec),q,q)

S <- Dmat %*% R %*% Dmat
Sy <- Load %*% S %*% t(Load)
diag(Sy)

which(diag(Sy)>1)
# Dimensions

d <- nthr + nload + ncorr + nvar

theta <- get_theta(
  THRESHOLDS = rep(thr, p),
  LOADINGS = Load,
  LATENT_COV = S,
  CAT = cat,
  CONSTRMAT = A,
  CONSTRVAR = exp(constr_lsd)^2,
  CORRFLAG = corrflag,
  STDLV = stdlv
)
constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))


proj_theta <- cpp_sa_proj(
  THETA = theta,
  CONSTRMAT = constr_list$CONSTRMAT,
  CONSTRLOGSD = constr_list$CONSTRLOGSD,
  C_VEC = cat,
  CORRFLAG =constr_list$CORRFLAG,
  NTHR  = nthr,
  NLOAD = nload,
  NCORR = ncorr,
  NVAR = nvar
)


pLoad <- cpp_loadings_theta2mat(
  THETA = proj_theta,
  CONSTRMAT = constr_list$CONSTRMAT,
  NTHR = nthr,
  NLOAD = nload
)

pS <- cpp_latvar_theta2mat(
  THETA = proj_theta,
  CONSTRLOGSD = constr_list$CONSTRLOGSD,
  NTHR = nthr,
  NLOAD=nload,
  NCORR = ncorr,
  NVAR = nvar,
  Q=q
)


Load
pLoad

S
pS
