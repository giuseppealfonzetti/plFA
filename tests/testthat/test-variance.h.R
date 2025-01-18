# test_that("check gradient from H",{
#   set.seed(1)
#   p <- 8L; q <- 4L; n <- 100
#   A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
#   Load <- gen_loadings(CONSTRMAT = A)
#   thr <- c(-1,1)
#   corrflag <- 1
#   S <- diag(1,q,q)
#   if(corrflag) S <- get_S(THETA = rnorm(q*(q-1)/2), Q = q, CORRFLAG = corrflag)
#   D <- sim_data(
#     SAMPLE_SIZE = n,
#     LOADINGS = Load,
#     THRESHOLDS = thr,
#     LATENT_COV = S)
#   cat <- apply(D, 2, max) + 1
#   theta <- get_theta(rep(thr, p), Load, S, cat, A,CORRFLAG = corrflag)
#   f <- compute_frequencies(Y = D, C_VEC = cat)
#
#
#
#   D <- check_data(D)
#   constr_list <- check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag))
#   dims <- check_dims(D, constr_list)
#   dims$d
#
#   lambda0_init <- init_thresholds(dims, constr_list)
#   lambda_init  <- init_loadings(dims, constr_list, FXD=0)
#   transformed_rhos_init <- init_transformed_latcorr(dims, constr_list)
#
#   par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)
#   length(par_init)
#
#
#   # function for nll
#   full_nll <- function(par_vec){
#     mod <- cpp_multiThread_completePairwise(
#       N = n,
#       C_VEC = cat,
#       CONSTRMAT = A,
#       FREQ = f,
#       THETA = par_vec,
#       CORRFLAG = 1,
#       SILENTFLAG = 1
#     )
#     out <- mod$iter_nll/n
#     return(out)
#   }
#
#
#   # function for gradient
#   full_ngr <- function(par_vec){
#     mod <- cpp_multiThread_completePairwise(
#       N = n,
#       C_VEC = cat,
#       CONSTRMAT = A,
#       FREQ = f,
#       THETA = par_vec,
#       CORRFLAG = 1,
#       SILENTFLAG = 1
#     )
#
#     out <- mod$iter_ngradient/n
#     return(out)
#   }
#
#   Hgr <- function(par_vec){
#     out <- estimate_H(
#       C_VEC = cat,
#       A = A,
#       THETA = par_vec,
#       FREQ = f,
#       N = n,
#       CORRFLAG = corrflag
#     )$gradient
#     out <- -out/n
#     return(out)
#   }
#
#   H <- function(par_vec){
#     out <- estimate_H(
#       C_VEC = cat,
#       A = A,
#       THETA = par_vec,
#       FREQ = f,
#       N = n,
#       CORRFLAG = corrflag
#     )$est_H
#     return(out)
#   }
#   expect_equal(full_ngr(theta), Hgr(theta))
#   expect_equal(full_ngr(par_init), Hgr(par_init))
#   # expect_equal(sum(abs(H(theta)+numDeriv::jacobian(Hgr, theta))>1e-3), 0)
#
# })
