#ifndef pairs_H
#define pairs_H
#include <RcppEigen.h>

#include "genericUtils.h"
#include "bivariateProbs.h"
#include "gradients.h"

namespace pairs{

  // Single pair contribution to
  // 1. Log-likelihood
  // 2. Gradient
  void pair_contribution_extended(
      // Parameters
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> LAM,
      const Eigen::Ref<const Eigen::MatrixXd> RU,
      const Eigen::Ref<const Eigen::MatrixXd> DU,
      const Eigen::Ref<const Eigen::MatrixXd> SIGMAU,
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int K,
      const int L,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const int SILENTFLAG,
      const int GRADFLAG,
      double &LL,
      Eigen::VectorXd &GRADIENT);
  // ){
  //   const int p = A.rows();
  //   const int q = A.cols();
  //
  //   // Identifies quantities related to pair (k,l)
  //   const int ck          = C_VEC(K);
  //   const int cl          = C_VEC(L);
  //   const Eigen::VectorXd lambdak  = LAM.row(K);
  //   const Eigen::VectorXd lambdal  = LAM.row(L);
  //   const double rho_kl            = lambdak.transpose() * SIGMAU * lambdal;
  //
  //   // Initialize pairs table
  //   Eigen::MatrixXd pairs_tab     = FREQ;
  //   pairs_tab.conservativeResize(FREQ.rows() + 1, Eigen::NoChange_t() );
  //
  //   // identify column index in freq table
  //   // i1: starting index item k
  //   int i1 = 0;
  //   if(K > 1){
  //     for(int u = 1; u < K; u++){
  //       const int cu = C_VEC(u);
  //       i1 += cu * C_VEC.segment(0,u).sum();
  //     }
  //   }
  //
  //   // i2 starting index from i1 for item l
  //   int i2 = 0;
  //   if(L > 0){
  //     i2 = C_VEC.segment(0,L).sum() * C_VEC(K);
  //   }
  //
  //   ////////////////////////////
  //   /* LIKELIHOOD COMPUTATION */
  //   ////////////////////////////
  //   for(int sk = 0; sk < ck; sk ++){
  //
  //     // i3: starting index from i2 for cat sk
  //     const int i3 = sk * cl;
  //
  //     for(int sl = 0; sl < cl; sl ++){
  //
  //       // final column index for pairs_tab. Print to check
  //       const int r = i1 + i2 + i3 + sl;
  //
  //       // read frequency
  //       const int n_sksl = FREQ(4, r);
  //
  //       // identify thresholds
  //       const Eigen::VectorXd pi_thresholds = params::extract_thresholds(TAU, C_VEC, K, L, sk, sl);
  //
  //       // compute pi
  //       const double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, K, L, sk, sl);
  //       pairs_tab(5,r) = pi_sksl;
  //
  //       // update ll
  //       LL += n_sksl * log(pi_sksl+1e-8);
  //     }
  //   }
  //
  //   //////////////////////////
  //   /* GRADIENT COMPUTATION */
  //   /////////////////////////
  //
  //   if(GRADFLAG == 1){
  //
  //     int iter = 0;
  //
  //     /////////////////////////////////////////////////////
  //     // (k,l)-pair likelihood derivative wrt thresholds //
  //     /////////////////////////////////////////////////////
  //     grads::thresholds(GRADIENT, pairs_tab, C_VEC, TAU, rho_kl, K, L, ck, cl, i1, i2, iter);
  //
  //     ///////////////////////////////////////////////////////////
  //     // (k,l)-pair likelihood derivative wrt URV correlation: //
  //     // intermediate step for derivatives wrt                 //
  //     // loadings and factor correlations                      //
  //     ///////////////////////////////////////////////////////////
  //     double tmp_kl = grads::rho_urv(pairs_tab, C_VEC, TAU, rho_kl, K, L, ck, cl, i1, i2);
  //
  //
  //     ///////////////////////////////////////////////////
  //     // (k,l)-pair likelihood derivative wrt loadings //
  //     ///////////////////////////////////////////////////
  //     grads::loadings(GRADIENT, A, LLC, SIGMAU, lambdak, lambdal, tmp_kl, p, q, K, L, iter);
  //
  //     /////////////////////////////////////////////////////////////////////////////
  //     // (k,l)-pair likelihood derivative wrt reparametrised latent correlations //
  //     /////////////////////////////////////////////////////////////////////////////
  //     if(CORRFLAG == 1){
  //       Eigen::VectorXd transformed_rhos = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR).segment(0, NCORR);
  //       grads::lat_corr(GRADIENT, A, lambdak.transpose()*DU, DU*lambdal, transformed_rhos, tmp_kl, q, NCORR, iter);
  //     }
  //
  //     ///////////////////////////////////////////////////////////////////////////
  //     // (k,l)-pair likelihood derivative wrt reparameterised latent variances //
  //     ///////////////////////////////////////////////////////////////////////////
  //     if(NVAR >=1){
  //       for(int j=0; j<q; j++){
  //         if(!std::isfinite(CONSTRLOGSD(j))){
  //           // Eigen::VectorXd ej = Eigen::VectorXd::Zero(q); ej(j) = 1;
  //           Eigen::MatrixXd dD = Eigen::MatrixXd::Zero(q,q); dD(j,j)=DU(j,j);
  //           GRADIENT(iter) = tmp_kl*lambdak.transpose()*(dD*RU*DU+DU*RU*dD)*lambdal;
  //           // Rcpp::Rcout<<"idx:"<<iter<<", gr:"<< GRADIENT(iter)<<"\n";
  //           iter++;
  //         }
  //       }
  //
  //     }
  //
  //   }
  //
  // }
    // Single pair contribution to
  // 1. Log-likelihood
  // 2. Gradient
  void pair_contribution2(
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> LAM,
      const Eigen::Ref<const Eigen::MatrixXd> RU,
      const Eigen::Ref<const Eigen::MatrixXd> DU,
      const Eigen::Ref<const Eigen::MatrixXd> SIGMAU,
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int K,
      const int L,
      const Eigen::MatrixXd &PAIRS_TABLE,
      const int SILENTFLAG,
      const int GRADFLAG,
      double &LL,
      Eigen::VectorXd &GRADIENT);
  // ){
  //   const int p = A.rows();
  //   const int q = A.cols();
  //   const int d = NTHR+NLOAD+NCORR+NVAR;
  //
  //   // // rearrange parameters
  //   // Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, LLC, NTHR, NLOAD);
  //   // Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
  //   // Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
  //   // Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  //   // Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);
  //
  //   // Identifies quantities related to pair (k,l)
  //   const int ck          = C_VEC(K);
  //   const int cl          = C_VEC(L);
  //   const Eigen::VectorXd lambdak  = LAM.row(K);
  //   const Eigen::VectorXd lambdal  = LAM.row(L);
  //   const double rho_kl            = lambdak.transpose() * SIGMAU * lambdal;
  //
  //   // Initialize pairs table
  //   Eigen::MatrixXd pairs_tab     = PAIRS_TABLE;
  //   pairs_tab.conservativeResize(PAIRS_TABLE.rows() + 1, Eigen::NoChange_t() );
  //
  //   // identify column index in freq table
  //   // i1: starting index item k
  //   int i1 = 0;
  //   if(K > 1){
  //     for(int u = 1; u < K; u++){
  //       const int cu = C_VEC(u);
  //       i1 += cu * C_VEC.segment(0,u).sum();
  //     }
  //   }
  //
  //   // i2 starting index from i1 for item l
  //   int i2 = 0;
  //   if(L > 0){
  //     i2 = C_VEC.segment(0,L).sum() * C_VEC(K);
  //   }
  //
  //
  //   ////////////////////////////
  //   /* LIKELIHOOD AND GRADIENT COMPUTATION */
  //   ////////////////////////////
  //   for(int sk = 0; sk < ck; sk ++){
  //
  //     // i3: starting index from i2 for cat sk
  //     const int i3 = sk * cl;
  //
  //     for(int sl = 0; sl < cl; sl ++){
  //
  //       // final column index for pairs_tab. Print to check
  //       const int r = i1 + i2 + i3 + sl;
  //
  //       // read frequency
  //       const int n_sksl = PAIRS_TABLE(4, r);
  //
  //       // identify thresholds
  //       const Eigen::VectorXd pi_thresholds = params::extract_thresholds(TAU, C_VEC, K, L, sk, sl);
  //
  //       // compute pi
  //       const double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, K, L, sk, sl);
  //       if(SILENTFLAG == 0)Rcpp::Rcout << "("<<K<<","<<L<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
  //       pairs_tab(5,r) = pi_sksl;
  //
  //
  //       // update ll
  //       LL += n_sksl * log(pi_sksl+1e-8);
  //
  //       if(GRADFLAG == 1){
  //
  //         GRADIENT += (n_sksl/(pi_sksl+1e-8))* grads::pi(THETA,
  //                      A,
  //                      CONSTRLOGSD,
  //                      LLC,
  //                      C_VEC,
  //                      pi_thresholds,
  //                      SIGMAU,
  //                      DU,
  //                      RU,
  //                      lambdak,
  //                      lambdal,
  //                      rho_kl,
  //                      d,
  //                      p,
  //                      q,
  //                      K,
  //                      L,
  //                      sk,
  //                      sl,
  //                      CORRFLAG,
  //                      NTHR,
  //                      NLOAD,
  //                      NCORR,
  //                      NVAR);
  //
  //       }
  //     }
  //   }
  //
  // }

}





#endif
