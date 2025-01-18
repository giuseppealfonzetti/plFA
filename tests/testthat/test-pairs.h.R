#### CORRFLAG ON ####
set.seed(1)
p <- 8L; q <- 4L; n <- 100
A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
Load <- gen_loadings(CONSTRMAT = A)
thr <- c(-1,1)
corrflag <- 1
S <- diag(1,q,q)
if(corrflag) S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q, CORRFLAG = corrflag)
D <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)
cat <- apply(D, 2, max) + 1
theta <- get_theta(rep(thr, p), Load, S, cat, A,CORRFLAG = corrflag)
f <- compute_frequencies(Y = D, C_VEC = cat)


D <- check_data(D)
constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag))
dims <- check_dims(D, constr_list)
dims$d

lambda0_init <- init_thresholds(dims, constr_list)
lambda_init  <- init_loadings(dims, constr_list, FXD=0)
transformed_rhos_init <- init_transformed_latcorr(dims, constr_list)

par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)
length(par_init)

pair_nll <- function(par_vec, OPTION){
  pair <- cpp_compute_pair(
    A = A,
    C_VEC = cat,
    THETA = par_vec,
    CORRFLAG = corrflag,
    k = k,
    l = l,
    PAIRS_TABLE = f,
    SILENTFLAG = 1,
    GRADFLAG = 0,
    OPTION = OPTION
  )
  out <- -pair$ll/n
  return(out)
}
pair_ngr <- function(par_vec, OPTION){
  pair <- cpp_compute_pair(
    A = A,
    C_VEC = cat,
    THETA = par_vec,
    CORRFLAG = corrflag,
    k = k,
    l = l,
    PAIRS_TABLE = f,
    SILENTFLAG = 1,
    GRADFLAG = 1,
    OPTION=OPTION
  )
  out <- -pair$ngradient/n
  return(out)
}

for (l in 2:p) {
  for (k in 1:(l-1)) {
    test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
      expect_equal(pair_ngr(theta, OPTION=0), numDeriv::grad(pair_nll, theta, OPTION=0), tolerance = 1e-4)
    })

    test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
      expect_equal(pair_ngr(theta, OPTION=1), numDeriv::grad(pair_nll, theta, OPTION=1), tolerance = 1e-4)
    })
  }
}


#### CORRFLAG OFF ####
set.seed(1)
p <- 8L; q <- 4L; n <- 100
A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
Load <- gen_loadings(CONSTRMAT = A)
thr <- c(-1,1)
corrflag <- 0
S <- diag(1,q,q)
if(corrflag) S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q, CORRFLAG = corrflag)
D <- sim_data(
  SAMPLE_SIZE = n,
  LOADINGS = Load,
  THRESHOLDS = thr,
  LATENT_COV = S)
cat <- apply(D, 2, max) + 1
theta <- get_theta(rep(thr, p), Load, S, cat, A,CORRFLAG = corrflag)
f <- compute_frequencies(Y = D, C_VEC = cat)


D <- check_data(D)
constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag))
dims <- check_dims(D, constr_list)
dims$d

lambda0_init <- init_thresholds(dims, constr_list)
lambda_init  <- init_loadings(dims, constr_list, FXD=0)
transformed_rhos_init <- init_transformed_latcorr(dims, constr_list)

par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)
length(par_init)

pair_nll <- function(par_vec, OPTION){
  pair <- cpp_compute_pair(
    A = A,
    C_VEC = cat,
    THETA = par_vec,
    CORRFLAG = corrflag,
    k = k,
    l = l,
    PAIRS_TABLE = f,
    SILENTFLAG = 1,
    GRADFLAG = 0,
    OPTION = OPTION
  )
  out <- -pair$ll/n
  return(out)
}
pair_ngr <- function(par_vec, OPTION){
  pair <- cpp_compute_pair(
    A = A,
    C_VEC = cat,
    THETA = par_vec,
    CORRFLAG = corrflag,
    k = k,
    l = l,
    PAIRS_TABLE = f,
    SILENTFLAG = 1,
    GRADFLAG = 1,
    OPTION=OPTION
  )
  out <- -pair$ngradient/n
  return(out)
}

for (l in 2:p) {
  for (k in 1:(l-1)) {
    test_that(paste0("check gradient pair (", l, ",", k, "):") ,{
      expect_equal(pair_ngr(theta, OPTION=0), numDeriv::grad(pair_nll, theta, OPTION=0), tolerance = 1e-4)
    })

    test_that(paste0("check gradient (pi version) pair (", l, ",", k, "):") ,{
      expect_equal(pair_ngr(theta, OPTION=1), numDeriv::grad(pair_nll, theta, OPTION=1), tolerance = 1e-4)
    })
  }
}

