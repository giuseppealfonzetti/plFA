# p = number of items, q = number of latent variables, n = number of observations
p <- 6; q <- 2; n <- 5000

# Thresholds vector for each item
thr <- c(-1.5, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

set.seed(1)
stdlv <- TRUE
corrflag <- TRUE

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
tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
Dmat <- diag(exp(tsdvec),q,q)

S <- Dmat %*% R %*% Dmat

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


dat <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)
constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
f <- compute_frequencies(Y = dat, C_VEC = cat)
dims <- check_dims(dat, constr_list)
init <- check_init(NULL, f, dims, constr_list, INIT_METHOD = "standard", FALSE)
init

test_that("init par",{
  expect_true(all(is.finite(init)))
  expect_true(all(is.finite(init)))}
)


Rwr_ncl <- function(par_vec){
  n <- dims$n
  mod <- cpp_multiThread_completePairwise(
    N          = dims$n,
    C_VEC      = dims$cat,
    CONSTRMAT  = constr_list$CONSTRMAT,
    CONSTRLOGSD= constr_list$CONSTRLOGSD,
    LLC        = constr_list$LLC,
    FREQ       = f,
    THETA      = par_vec,
    CORRFLAG   = constr_list$CORRFLAG,
    NTHR       = dims$nthr,
    NLOAD      = dims$nload,
    NCORR      = dims$ncorr,
    NVAR       = dims$nvar,
    SILENTFLAG = 1
  )
  out <- mod$iter_nll/n
  return(out)
}
init_nll <- Rwr_ncl(init)
test_that("init nll",{
  expect_true(is.finite(init_nll))
  expect_true(init_nll<1e3)
})



Rwr_ngr <- function(par_vec){
  n <- dims$n
  mod <- cpp_multiThread_completePairwise(
    N          = dims$n,
    C_VEC      = dims$cat,
    CONSTRMAT  = constr_list$CONSTRMAT,
    CONSTRLOGSD= constr_list$CONSTRLOGSD,
    LLC        = constr_list$LLC,
    FREQ       = f,
    THETA      = par_vec,
    CORRFLAG   = constr_list$CORRFLAG,
    NTHR       = dims$nthr,
    NLOAD      = dims$nload,
    NCORR      = dims$ncorr,
    NVAR       = dims$nvar,
    SILENTFLAG = 1
  )
  out <- mod$iter_ngradient
  return(out)
}
init_ngr <- Rwr_ngr(init)
test_that("init ngr",{
  expect_true(all(is.finite(init_ngr)))
})
