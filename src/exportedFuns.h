#ifndef exportedFuns_H
#define exportedFuns_H

#include "pairs.h"
#include "optimisationUtils.h"


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
   if(OPTION==0){
     pairs::pair_contribution_extended(CONSTRMAT, CONSTRLOGSD, C_VEC, THETA, CORRFLAG,
                                       NTHR, NLOAD, NCORR, NVAR, K, L,
                                       PAIRS_TABLE, SILENTFLAG, GRADFLAG,
                                       ll, gradient);
     }else if(OPTION==1){
       pairs::pair_contribution2(CONSTRMAT, CONSTRLOGSD, C_VEC, THETA, CORRFLAG,
                                         NTHR, NLOAD, NCORR, NVAR, K, L,
                                         PAIRS_TABLE, SILENTFLAG, GRADFLAG,
                                         ll, gradient);
       }else if(OPTION==2){
         // pairs::pair_contribution( CONSTRMAT, CONSTRLOGSD, C_VEC, THETA, CORRFLAG,
         //                           NTHR, NLOAD, NCORR, NVAR, K, L,
         //                           PAIRS_TABLE, SILENTFLAG, GRADFLAG,
         //                           ll, gradient);
       }


   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("nll") = -ll,
       Rcpp::Named("ngradient") = -gradient
     );
   return(output);

 }

// [[Rcpp::export]]
Eigen::VectorXd cpp_loadings_theta2vec(
    Eigen::Map<Eigen::VectorXd> THETA,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2vec(THETA, NTHR, NLOAD);
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_theta2mat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int NTHR,
    const int NLOAD
){
  return params::loadings::theta2mat(THETA, CONSTRMAT, NTHR, NLOAD);
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_loadings_mat2vec(
    Eigen::Map<Eigen::MatrixXd> LOADINGS,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const int  NLOAD){
  return params::loadings::mat2vec(LOADINGS, CONSTRMAT, NLOAD);
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_vec2cmat(Eigen::Map<Eigen::VectorXd> VEC,
                             const int NCORR,
                             const int Q){
  return params::latvar::vec2cmat(VEC, NCORR, Q);
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_vec2dmat(Eigen::Map<Eigen::VectorXd> VEC,
                             Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                             const int NCORR,
                             const int NVAR,
                             const int Q){
  return params::latvar::vec2dmat(VEC, CONSTRLOGSD, NCORR, NVAR, Q);
}

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

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_theta2dmat(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int Q){


  return params::latvar::theta2dmat(THETA,CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q);
}

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

// [[Rcpp::export]]
Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S,
                                   Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
                                   const int NCORR,
                                   const int NVAR){
  return params::latvar::mat2vec(S, CONSTRLOGSD, NCORR, NVAR);
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_latvar_mat2cmat(Eigen::Map<Eigen::MatrixXd> S){
  return params::latvar::mat2cmat(S);
}


// [[Rcpp::export]]
Eigen::VectorXd cpp_sa_proj(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
){
  bool flag = false;
  Eigen::VectorXd theta = THETA;
  sa::proj2(CONSTRMAT,
           CONSTRLOGSD,
           C_VEC,
           CORRFLAG,
           NTHR,
           NLOAD,
           NCORR,
           NVAR,
           theta,
           flag);
  return theta;
}




#endif
