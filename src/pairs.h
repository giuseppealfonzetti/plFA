#ifndef pairs_H
#define pairs_H
#include <RcppParallel.h>

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
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const unsigned int K,
      const unsigned int L,
      const Eigen::MatrixXd &PAIRS_TABLE,
      const unsigned int SILENTFLAG,
      const unsigned int GRADFLAG,
      double &LL,
      Eigen::VectorXd &GRADIENT
  ){
    const unsigned int p = A.rows();
    const unsigned int q = A.cols();

    // rearrange parameters
    Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, NTHR, NLOAD);
    Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
    Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
    Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
    Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);

    // Identifies quantities related to pair (k,l)
    const unsigned int ck          = C_VEC(K);
    const unsigned int cl          = C_VEC(L);
    const Eigen::VectorXd lambdak  = Lam.row(K);
    const Eigen::VectorXd lambdal  = Lam.row(L);
    const double rho_kl            = lambdak.transpose() * Sigma_u * lambdal;

    // Rcpp::Rcout<<"rho_kl:"<<Ru<<"\n";
    // Initialize pairs table
    Eigen::MatrixXd pairs_tab     = PAIRS_TABLE;
    pairs_tab.conservativeResize(PAIRS_TABLE.rows() + 1, Eigen::NoChange_t() );

    // identify column index in freq table
    // i1: starting index item k
    unsigned int i1 = 0;
    if(K > 1){
      for(unsigned int u = 1; u < K; u++){
        const unsigned int cu = C_VEC(u);
        i1 += cu * C_VEC.segment(0,u).sum();
      }
    }

    // i2 starting index from i1 for item l
    unsigned int i2 = 0;
    if(L > 0){
      i2 = C_VEC.segment(0,L).sum() * C_VEC(K);
    }

    ////////////////////////////
    /* LIKELIHOOD COMPUTATION */
    ////////////////////////////
    for(unsigned int sk = 0; sk < ck; sk ++){

      // i3: starting index from i2 for cat sk
      const unsigned int i3 = sk * cl;

      for(unsigned int sl = 0; sl < cl; sl ++){

        // final column index for pairs_tab. Print to check
        const unsigned int r = i1 + i2 + i3 + sl;

        // read frequency
        const unsigned int n_sksl = PAIRS_TABLE(4, r);

        // identify thresholds
        const Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, C_VEC, K, L, sk, sl);

        // compute pi
        const double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, K, L, sk, sl);
        pairs_tab(5,r) = pi_sksl;

        // update ll
        LL += n_sksl * log(pi_sksl+1e-8);
      }
    }

    //////////////////////////
    /* GRADIENT COMPUTATION */
    /////////////////////////

    if(GRADFLAG == 1){

      unsigned int iter = 0;

      /////////////////////////////////////////////////////
      // (k,l)-pair likelihood derivative wrt thresholds //
      /////////////////////////////////////////////////////
      grads::thresholds(GRADIENT, pairs_tab, C_VEC, tau, rho_kl, K, L, ck, cl, i1, i2, iter);

      ///////////////////////////////////////////////////////////
      // (k,l)-pair likelihood derivative wrt URV correlation: //
      // intermediate step for derivatives wrt                 //
      // loadings and factor correlations                      //
      ///////////////////////////////////////////////////////////
      double tmp_kl = grads::rho_urv(pairs_tab, C_VEC, tau, rho_kl, K, L, ck, cl, i1, i2);


      ///////////////////////////////////////////////////
      // (k,l)-pair likelihood derivative wrt loadings //
      ///////////////////////////////////////////////////
      grads::loadings(GRADIENT, A, Sigma_u, lambdak, lambdal, tmp_kl, p, q, K, L, iter);

      /////////////////////////////////////////////////////////////////////////////
      // (k,l)-pair likelihood derivative wrt reparametrised latent correlations //
      /////////////////////////////////////////////////////////////////////////////
      if(CORRFLAG == 1){
        Eigen::VectorXd transformed_rhos = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR).segment(0, NCORR);
        grads::lat_corr(GRADIENT, A, lambdak.transpose()*Du, Du*lambdal, transformed_rhos, tmp_kl, q, NCORR, iter);
      }

      ///////////////////////////////////////////////////////////////////////////
      // (k,l)-pair likelihood derivative wrt reparameterised latent variances //
      ///////////////////////////////////////////////////////////////////////////
      if(NVAR >=1){
        for(int j=0; j<q; j++){
          if(!std::isfinite(CONSTRLOGSD(j))){
            // Eigen::VectorXd ej = Eigen::VectorXd::Zero(q); ej(j) = 1;
            Eigen::MatrixXd dD = Eigen::MatrixXd::Zero(q,q); dD(j,j)=Du(j,j);
            GRADIENT(iter) = tmp_kl*lambdak.transpose()*(dD*Ru*Du+Du*Ru*dD)*lambdal;
            // Rcpp::Rcout<<"idx:"<<iter<<", gr:"<< GRADIENT(iter)<<"\n";
            iter++;
          }
        }

      }

    }

  }
    // Single pair contribution to
  // 1. Log-likelihood
  // 2. Gradient
  void pair_contribution2(
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const unsigned int K,
      const unsigned int L,
      const Eigen::MatrixXd &PAIRS_TABLE,
      const unsigned int SILENTFLAG,
      const unsigned int GRADFLAG,
      double &LL,
      Eigen::VectorXd &GRADIENT
  ){
    const unsigned int p = A.rows();
    const unsigned int q = A.cols();
    const unsigned int d = NTHR+NLOAD+NCORR+NVAR;

    // rearrange parameters
    Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, NTHR, NLOAD);
    Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
    Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
    Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
    Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);

    // Identifies quantities related to pair (k,l)
    const unsigned int ck          = C_VEC(K);
    const unsigned int cl          = C_VEC(L);
    const Eigen::VectorXd lambdak  = Lam.row(K);
    const Eigen::VectorXd lambdal  = Lam.row(L);
    const double rho_kl            = lambdak.transpose() * Sigma_u * lambdal;

    // Initialize pairs table
    Eigen::MatrixXd pairs_tab     = PAIRS_TABLE;
    pairs_tab.conservativeResize(PAIRS_TABLE.rows() + 1, Eigen::NoChange_t() );

    // identify column index in freq table
    // i1: starting index item k
    unsigned int i1 = 0;
    if(K > 1){
      for(unsigned int u = 1; u < K; u++){
        const unsigned int cu = C_VEC(u);
        i1 += cu * C_VEC.segment(0,u).sum();
      }
    }

    // i2 starting index from i1 for item l
    unsigned int i2 = 0;
    if(L > 0){
      i2 = C_VEC.segment(0,L).sum() * C_VEC(K);
    }


    ////////////////////////////
    /* LIKELIHOOD AND GRADIENT COMPUTATION */
    ////////////////////////////
    for(unsigned int sk = 0; sk < ck; sk ++){

      // i3: starting index from i2 for cat sk
      const unsigned int i3 = sk * cl;

      for(unsigned int sl = 0; sl < cl; sl ++){

        // final column index for pairs_tab. Print to check
        const unsigned int r = i1 + i2 + i3 + sl;

        // read frequency
        const unsigned int n_sksl = PAIRS_TABLE(4, r);

        // identify thresholds
        const Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, C_VEC, K, L, sk, sl);

        // compute pi
        const double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, K, L, sk, sl);
        if(SILENTFLAG == 0)Rcpp::Rcout << "("<<K<<","<<L<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
        pairs_tab(5,r) = pi_sksl;


        // update ll
        LL += n_sksl * log(pi_sksl+1e-8);

        if(GRADFLAG == 1){

          GRADIENT += (n_sksl/(pi_sksl+1e-8))* grads::pi(THETA,
                       A,
                       CONSTRLOGSD,
                       C_VEC,
                       pi_thresholds,
                       Sigma_u,
                       Du,
                       Ru,
                       lambdak,
                       lambdal,
                       rho_kl,
                       d,
                       p,
                       q,
                       K,
                       L,
                       ck,
                       cl,
                       sk,
                       sl,
                       i1,
                       i2,
                       CORRFLAG,
                       NTHR,
                       NLOAD,
                       NCORR,
                       NVAR);

        }
      }
    }

  }

  // RcppParallel Worker to compute pair_contribution
  // on multiple pairs in parallel
  struct SubsetWorker : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Map<Eigen::MatrixXd> constrmat;
    const Eigen::Map<Eigen::VectorXd> constrsd;
    const Eigen::Map<Eigen::VectorXd> c_vec;
    const Eigen::MatrixXd &pairs_table;
    const Eigen::MatrixXd &items_pairs;
    const unsigned int corrFLAG;
    const int nthr;
    const int nload;
    const int ncorr;
    const int nvar;
    const unsigned int silentFLAG;
    const unsigned int gradFLAG;

    //// Iteration:
    const Eigen::VectorXd &theta;
    const std::vector<int> &index_vector;

    // Output quantities:
    double subset_ll;
    Eigen::VectorXd subset_gradient, subset_gradient2;

    // Constructor 1:
    SubsetWorker(
      const Eigen::Map<Eigen::MatrixXd> CONSTRMAT_,
      const Eigen::Map<Eigen::VectorXd> CONSTRSD_,
      const Eigen::Map<Eigen::VectorXd> C_VEC_,
      const Eigen::MatrixXd &PAIRS_TABLE_,
      const Eigen::MatrixXd &ITEMS_PAIRS_,
      const unsigned int CORRFLAG_,
      const int NTHR_,
      const int NLOAD_,
      const int NCORR_,
      const int NVAR_,
      const unsigned int SILENTFLAG_,
      const unsigned int GRADFLAG_,
      const Eigen::VectorXd &THETA_,
      const std::vector<int> &INDEX_VECTOR_
    ):
      constrmat(CONSTRMAT_), constrsd(CONSTRSD_), c_vec(C_VEC_), pairs_table(PAIRS_TABLE_), items_pairs(ITEMS_PAIRS_),
      corrFLAG(CORRFLAG_), nthr(NTHR_), nload(NLOAD_), ncorr(NCORR_), nvar(NVAR_), silentFLAG(SILENTFLAG_), gradFLAG(GRADFLAG_),
      theta(THETA_), index_vector(INDEX_VECTOR_), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(THETA_.size())){}

    // Constructor 2:
    SubsetWorker(const SubsetWorker &OBJ_, RcppParallel::Split):
      constrmat(OBJ_.constrmat), constrsd(OBJ_.constrsd), c_vec(OBJ_.c_vec), pairs_table(OBJ_.pairs_table), items_pairs(OBJ_.items_pairs),
      corrFLAG(OBJ_.corrFLAG), nthr(OBJ_.nthr), nload(OBJ_.nload), ncorr(OBJ_.ncorr), nvar(OBJ_.nvar),
      silentFLAG(OBJ_.silentFLAG), gradFLAG(OBJ_.gradFLAG),
      theta(OBJ_.theta), index_vector(OBJ_.index_vector), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const SubsetWorker &RHS);

  };

  void SubsetWorker::operator()(std::size_t BEGIN, std::size_t END){
    // Run along the pairs identified by indexes in index_vector
    // computing their contribution to nll, gradient and eventually Hessian.
    unsigned int d = theta.size();
    //subset_gradient.resize(d); subset_gradient.setZero();
    for (unsigned int h = BEGIN; h < END; h++){

      // identify corresponding column in pairs_table
      unsigned int col = index_vector[h];

      // identify the pair
      unsigned int k = items_pairs(0, col);
      unsigned int l = items_pairs(1, col);

      // initialize empty pair-output
      double pair_ll = 0;
      Eigen::VectorXd pair_gradient = Eigen::VectorXd::Zero(d);

      // computation of log-likelihood, gradient
      pair_contribution_extended(constrmat, constrsd, c_vec, theta, corrFLAG, nthr, nload, ncorr, nvar, k, l, pairs_table, silentFLAG, gradFLAG, pair_ll, pair_gradient);

      // update
      {
        subset_ll += pair_ll;
        subset_gradient += pair_gradient;
      }
    }
  }

  void SubsetWorker::join(const SubsetWorker &RHS){
    subset_ll += RHS.subset_ll;
    for(unsigned i = 0; i < subset_gradient.size(); i++){
      subset_gradient(i) += RHS.subset_gradient(i);

    }
  }
}





#endif
