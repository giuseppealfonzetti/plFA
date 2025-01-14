#ifndef exportedFuns_H
#define exportedFuns_H

#include "pairs.h"


//' Get loading matrix from theta
//'
//' get_Lam() constructs the loading matrix from the parameter vector
//'
//' @param THETA Numerical vector of parameters.
//' @param A Binary matrix of dimension \eqn{p*q} where \eqn{p} is the number
//' of items and \eqn{q} the number of latent variables. Entries equal to
//' \eqn{1} refer to free loadings, while entries equal to \eqn{0} indicate
//' loadings constrained to be null.
//' @param C Sum of the number of categories for each item.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_Lam(Eigen::Map<Eigen::MatrixXd> A,
                        const unsigned int C,
                        Eigen::Map<Eigen::VectorXd> THETA
){
  return params::get_Lam(A, C, THETA);
}

//' Get latent correlation matrix from theta
//'
//' get_S() extracts the latent correlation matrix from theta assuming
//' theta elements to be reparametrised following the
//' Lewandowski-Kurowicka-Joe (2009) transform.
//'
//' @param THETA Numerical vector of parameters.
//' @param Q Number of latent variables.
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_S(Eigen::Map<Eigen::VectorXd> THETA,
                      const unsigned int Q
){
  return params::get_S(THETA, Q);
}

//' Get transformed parameters from latent correlation matrix
//'
//' Use Lewandowski-Kurowicka-Joe (2009) transformation on latent correlations
//' with Choleski decomposition.
//'
//' @param S Latent correlation matrix.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_par_from_S(Eigen::Map<Eigen::MatrixXd> S){
  return params::get_par_from_S(S);
}

//' Single pair contribution
//'
//' @description
//' Wrapper of pair_contribution() used for unit tests
//'
//' @param A Constraint matrix. Loadings free to be estimated are identified by a 1.
//' @param C_VEC Vector containing the number of categories for each item
//' @param THETA Parameter vector
//' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
//' @param k first index identifying the pair
//' @param l second index identifying the pair
//' @param PAIRS_TABLE output from [pairs_freq()]
//' @param SILENTFLAG optional for verbose output
//' @param GRADFLAG 1 to compute gradient
// [[Rcpp::export]]
Rcpp::List cpp_compute_pair(
     Eigen::Map<Eigen::MatrixXd> A,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int CORRFLAG,
     const unsigned int k,
     const unsigned int l,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int SILENTFLAG,
     const unsigned int GRADFLAG
 ){
   const unsigned int p = A.rows();
   const unsigned int q = A.cols();
   const unsigned int d = THETA.size();
   const unsigned int c = C_VEC.sum();
   const unsigned int nthr = c-p;
   const unsigned int ncorr = q*(q-1)/2;
   const unsigned int nload = d-nthr-ncorr;

   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);

   pairs::pair_contribution(A, C_VEC, THETA, CORRFLAG, k, l, PAIRS_TABLE, SILENTFLAG, GRADFLAG, ll, gradient);

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("ll") = -ll,
       Rcpp::Named("ngradient") = -gradient
     );
   return(output);

 }

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_thresholds_theta2vec(Eigen::Map<Eigen::VectorXd> THETA,
                                             const unsigned int P,
                                             const unsigned int C){
  return params::get_thresholds_theta2vec(THETA, P, C);
}
//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_loadings_mat2vec(Eigen::Map<Eigen::MatrixXd> LOADINGS,
                                         Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
                                         const int  NLOAD){
  return params::get_loadings_mat2vec(LOADINGS, CONSTRMAT, NLOAD);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_loadings_vec2mat(Eigen::Map<Eigen::VectorXd> LOADINGS,
                                         Eigen::Map<Eigen::MatrixXd> CONSTRMAT){
  return params::get_loadings_vec2mat(LOADINGS, CONSTRMAT);
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S){
  return params::get_latvar_mat2vec(S);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_latvar_vec2mat(Eigen::Map<Eigen::VectorXd> SVEC,
                                       const unsigned int Q){
  return params::get_latvar_vec2mat(SVEC, Q);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int P,
    const int Q,
    const int D,
    const int C
){
  return params::get_loadings_theta2mat(THETA, CONSTRMAT, P, Q, D, C);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_latvar_theta2mat(Eigen::Map<Eigen::VectorXd> THETA,
                                         const int Q,
                                         const int D
){
  return params::get_latvar_theta2mat(THETA, Q, D);
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_latvar_theta2vec(Eigen::Map<Eigen::VectorXd> THETA,
                                         const unsigned int NTHR,
                                         const unsigned int NLOAD,
                                         const unsigned int NCORR){
  return params::get_latvar_theta2vec(THETA, NTHR, NLOAD, NCORR);
}
#endif
