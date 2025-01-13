#ifndef exportedFuns_H
#define exportedFuns_H

#include "pairs.h"

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
#endif
