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
Rcpp::List cpp_compute_pair_ext(
     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
     Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int CORRFLAG,
     const int NTHR,
     const int NLOAD,
     const int NCORR,
     const int NVAR,
     const unsigned int K,
     const unsigned int L,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int SILENTFLAG,
     const unsigned int GRADFLAG,
     const int OPTION = 0
 ){
   const int p = CONSTRMAT.rows();
   const int q = CONSTRMAT.cols();
   const int c = C_VEC.sum();
   const int d = NTHR+NLOAD+NCORR+NVAR;

   if(THETA.size()!=(NTHR+NLOAD+NCORR+NVAR))Rcpp::stop("Check theta size");


   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
   pairs::pair_contribution_extended(CONSTRMAT, CONSTRLOGSD, C_VEC, THETA, CORRFLAG,
                                     NTHR, NLOAD, NCORR, NVAR, K, L,
                                     PAIRS_TABLE, SILENTFLAG, GRADFLAG,
                                     ll, gradient);

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("nll") = -ll,
       Rcpp::Named("ngradient") = -gradient
     );
   return(output);

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
     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
     Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int CORRFLAG,
     const int NTHR,
     const int NLOAD,
     const int NCORR,
     const int NVAR,
     const unsigned int K,
     const unsigned int L,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int SILENTFLAG,
     const unsigned int GRADFLAG,
     const int OPTION = 0
 ){
   const unsigned int p = CONSTRMAT.rows();
   const unsigned int q = CONSTRMAT.cols();
   const unsigned int d = THETA.size();
   const unsigned int c = C_VEC.sum();
   const unsigned int nthr = c-p;
   unsigned int ncorr = 0; if(CORRFLAG==1) ncorr = q*(q-1)/2;
   const unsigned int nload = d-nthr-ncorr;

   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);

   if(OPTION==0){
     pairs::pair_contribution(CONSTRMAT, CONSTRLOGSD, C_VEC, THETA, CORRFLAG,
                              NTHR, NLOAD, NCORR, NVAR, K, L,
                              PAIRS_TABLE, SILENTFLAG, GRADFLAG, ll, gradient);
   }else if(OPTION==1){
     pairs::pair_contribution2(CONSTRMAT, C_VEC, THETA, CORRFLAG, K, L, PAIRS_TABLE, SILENTFLAG, GRADFLAG, ll, gradient);
   }

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
                                       const int Q){
  return params::get_latvar_vec2mat(SVEC, Q);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int P,
    const int C,
    const int NLOAD
){
  return params::get_loadings_theta2mat(THETA, CONSTRMAT, P, C, NLOAD);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_get_latvar_theta2mat(Eigen::Map<Eigen::VectorXd> THETA,
                                         const int Q,
                                         const int D,
                                         const int CORRFLAG
){
  return params::get_latvar_theta2mat(THETA, Q, D, CORRFLAG);
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_get_latvar_theta2vec(Eigen::Map<Eigen::VectorXd> THETA,
                                         const int NTHR,
                                         const int NLOAD,
                                         const int NCORR,
                                         const int CORRFLAG){
  return params::get_latvar_theta2vec(THETA, NTHR, NLOAD, NCORR, CORRFLAG);
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_loadings_theta2vec(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2vec(THETA, NTHR, NLOAD);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2mat(THETA, CONSTRMAT, NTHR, NLOAD);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_mat2vec(
    Eigen::Map<Eigen::MatrixXd> LOADINGS,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int  NLOAD){
  return params::loadings::mat2vec(LOADINGS, CONSTRMAT, NLOAD);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_vec2cmat(Eigen::Map<Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q){
  return params::latvar::vec2cmat(VEC, NCORR, Q);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_vec2dmat(Eigen::Map<Eigen::VectorXd> VEC,
                             Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                             const int NCORR,
                             const int NVAR,
                             const int Q){
  return params::latvar::vec2dmat(VEC, CONSTRLOGSD, NCORR, NVAR, Q);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2cmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q){

  Eigen::VectorXd vec = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR);
  return params::latvar::vec2cmat(vec, NCORR, Q);
}


//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q){
  return params::latvar::theta2mat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q);
}

//' @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S,
                                   Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                                   const int NCORR,
                                   const int NVAR){
  return params::latvar::mat2vec(S, CONSTRLOGSD, NCORR, NVAR);
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_mat2cmat(Eigen::Map<Eigen::MatrixXd> S){
  return params::latvar::mat2cmat(S);
}
#endif
