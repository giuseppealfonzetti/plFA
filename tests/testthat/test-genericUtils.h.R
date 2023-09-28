# test_that("get_S() check implementations", {
#   set.seed(123)
#   p <- 10; q <- 5;
#   par <- rnorm(q*(q-1)/2)
#
#   S <- get_S(THETA = par, Q = q)
#   S2 <- get_S2(THETA = par, Q = q)
#   expect_equal(S, S2)
#
#   R_S <- function(TRHO){
#     vec <- par
#     vec[idx+1] <- TRHO
#     S <- get_S2(THETA = vec, Q = q)
#     S
#   }
#
#   for (idx in 0:(length(par)-1)) {
#     numgr <- numDeriv::jacobian(R_S,par[idx+1])
#     numgr <- matrix(as.vector(numgr),q,q)
#     gr <- grad_S(A = matrix(1,p,q), THETA = par, IDX = idx)
#     gr2 <- grad_S2(A = matrix(1,p,q), THETA = par, IDX = idx)
#
#     expect_equal(numgr, gr)
#     expect_equal(gr, gr2)
#   }
# })
