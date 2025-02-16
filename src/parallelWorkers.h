#ifndef paralellWorkers_H
#define paralellWorkers_H
#include <RcppParallel.h>
#include "pairs.h"

namespace parallelWorkers{

  // RcppParallel Worker to compute pair_contribution
  // on multiple pairs in parallel
  struct pairsGradient : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> constrmat;
    const Eigen::Ref<const Eigen::VectorXd> constrsd;
    const std::vector<std::vector<std::vector<double>>> llc;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> pairs_table;
    const Eigen::Ref<const Eigen::MatrixXd> items_pairs;
    const int corrFLAG;
    const int nthr;
    const int nload;
    const int ncorr;
    const int nvar;
    const int silentFLAG;
    const int gradFLAG;

    //// Iteration:
    const Eigen::VectorXd &theta;
    const std::vector<int> &index_vector;

    // Output quantities:
    double subset_ll;
    Eigen::VectorXd subset_gradient, subset_gradient2;

    // Constructor 1:
    pairsGradient(
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT_,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRSD_,
      const std::vector<std::vector<std::vector<double>>> LLC_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> PAIRS_TABLE_,
      const Eigen::Ref<const Eigen::MatrixXd> ITEMS_PAIRS_,
      const int CORRFLAG_,
      const int NTHR_,
      const int NLOAD_,
      const int NCORR_,
      const int NVAR_,
      const int SILENTFLAG_,
      const int GRADFLAG_,
      const Eigen::VectorXd &THETA_,
      const std::vector<int> &INDEX_VECTOR_
    ):
      constrmat(CONSTRMAT_), constrsd(CONSTRSD_), llc(LLC_), c_vec(C_VEC_), pairs_table(PAIRS_TABLE_), items_pairs(ITEMS_PAIRS_),
      corrFLAG(CORRFLAG_), nthr(NTHR_), nload(NLOAD_), ncorr(NCORR_), nvar(NVAR_), silentFLAG(SILENTFLAG_), gradFLAG(GRADFLAG_),
      theta(THETA_), index_vector(INDEX_VECTOR_), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(THETA_.size())){}

    // Constructor 2:
    pairsGradient(const pairsGradient &OBJ_, RcppParallel::Split):
      constrmat(OBJ_.constrmat), constrsd(OBJ_.constrsd), llc(OBJ_.llc), c_vec(OBJ_.c_vec), pairs_table(OBJ_.pairs_table), items_pairs(OBJ_.items_pairs),
      corrFLAG(OBJ_.corrFLAG), nthr(OBJ_.nthr), nload(OBJ_.nload), ncorr(OBJ_.ncorr), nvar(OBJ_.nvar),
      silentFLAG(OBJ_.silentFLAG), gradFLAG(OBJ_.gradFLAG),
      theta(OBJ_.theta), index_vector(OBJ_.index_vector), subset_ll(0.0),
      subset_gradient(Eigen::VectorXd::Zero(theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const pairsGradient &RHS);

  };

  void pairsGradient::operator()(std::size_t BEGIN, std::size_t END){
    // Run along the pairs identified by indexes in index_vector
    // computing their contribution to nll, gradient and eventually Hessian.
    int d = theta.size();
    //subset_gradient.resize(d); subset_gradient.setZero();
    for (int h = BEGIN; h < END; h++){

      // identify corresponding column in pairs_table
      const int col = index_vector[h];

      // identify the pair
      const int k = items_pairs(0, col);
      const int l = items_pairs(1, col);

      // initialize empty pair-output
      double pair_ll = 0;
      Eigen::VectorXd pair_gradient = Eigen::VectorXd::Zero(d);

      // computation of log-likelihood, gradient
      // rearrange parameters
      const int q = constrmat.cols();
      const Eigen::MatrixXd Lam             = params::loadings::theta2mat(theta, constrmat, llc, nthr, nload);
      const Eigen::MatrixXd Ru              = params::latvar::theta2cmat(theta, nthr, nload, ncorr, nvar, q);
      const Eigen::MatrixXd Du              = params::latvar::theta2dmat(theta, constrsd, nthr, nload, ncorr, nvar, q);
      const Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
      const Eigen::VectorXd tau             = params::thresholds::theta2vec(theta, nthr);
      pairs::pair_contribution_extended(constrmat, constrsd, llc, c_vec, theta,
                                 Lam, Ru, Du, Sigma_u, tau,
                                 corrFLAG, nthr, nload, ncorr, nvar, k, l, pairs_table, silentFLAG, gradFLAG, pair_ll, pair_gradient);

                                  // update
                                 {
                                   subset_ll += pair_ll;
                                   subset_gradient += pair_gradient;
                                 }

    }
  }

  void pairsGradient::join(const pairsGradient &RHS){
    subset_ll += RHS.subset_ll;
    for(int i = 0; i < subset_gradient.size(); i++){
      subset_gradient(i) += RHS.subset_gradient(i);

    }
  }

  // RcppParallel Worker to compute gradient outer product of
  // bivariate response patterns. Used to evaluate the
  // sample estimator of H
  struct patternwiseOuterProd : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> constrmat;
    const Eigen::Ref<const Eigen::VectorXd> constrsd;
    const std::vector<std::vector<std::vector<double>>> llc;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;
    const Eigen::Ref<const Eigen::MatrixXd> freq;
    const int corrFLAG;
    const int nthr;
    const int nload;
    const int ncorr;
    const int nvar;

    //// Iteration:
    const Eigen::VectorXd &theta;
    Eigen::MatrixXd &gradmat;

    // Output quantities:
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;


    // Constructor 1:
    patternwiseOuterProd(
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT_,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRSD_,
      const std::vector<std::vector<std::vector<double>>> LLC_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ_,
      const int CORRFLAG_,
      const int NTHR_,
      const int NLOAD_,
      const int NCORR_,
      const int NVAR_,
      const Eigen::VectorXd &THETA_,
      Eigen::MatrixXd &GRADMAT_

    ):
      constrmat(CONSTRMAT_), constrsd(CONSTRSD_), llc(LLC_), c_vec(C_VEC_), freq(FREQ_),
      corrFLAG(CORRFLAG_), nthr(NTHR_), nload(NLOAD_), ncorr(NCORR_), nvar(NVAR_),
      theta(THETA_),
      gradmat(GRADMAT_),
      gradient(Eigen::VectorXd::Zero(THETA_.size())),
      hessian(Eigen::MatrixXd::Zero(THETA_.size(), THETA_.size()))
      {}

    // Constructor 2:
    patternwiseOuterProd(const patternwiseOuterProd &OBJ_, RcppParallel::Split):
      constrmat(OBJ_.constrmat), constrsd(OBJ_.constrsd), llc(OBJ_.llc), c_vec(OBJ_.c_vec), freq(OBJ_.freq),
      corrFLAG(OBJ_.corrFLAG), nthr(OBJ_.nthr), nload(OBJ_.nload), ncorr(OBJ_.ncorr), nvar(OBJ_.nvar),
      theta(OBJ_.theta),
      gradmat(OBJ_.gradmat),
      gradient(Eigen::VectorXd::Zero(theta.size())),
      hessian(Eigen::MatrixXd::Zero(theta.size(), theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const patternwiseOuterProd &RHS);

  };

  void patternwiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){
    // Run along the pairs identified by indexes in index_vector
    // computing their contribution to nll, gradient and eventually Hessian.
    const int d = theta.size();
    const int q = constrmat.cols();
    const int p = constrmat.rows();

    // rearrange parameters
    const Eigen::MatrixXd Lam             = params::loadings::theta2mat(theta, constrmat, llc, nthr, nload);
    const Eigen::MatrixXd Ru              = params::latvar::theta2cmat(theta, nthr, nload, ncorr, nvar, q);
    const Eigen::MatrixXd Du              = params::latvar::theta2dmat(theta, constrsd, nthr, nload, ncorr, nvar, q);
    const Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
    const Eigen::VectorXd tau             = params::thresholds::theta2vec(theta, nthr);

    //subset_gradient.resize(d); subset_gradient.setZero();
    for (int h = BEGIN; h < END; h++){


      // identify the pair
      const int k      = freq(0, h);
      const int l      = freq(1, h);
      const int sk     = freq(2, h);
      const int sl     = freq(3, h);
      const int n_sksl = freq(4, h);

      // identify thresholds
      Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, c_vec, k, l, sk, sl);

      // implied model correlation
      double rho_kl = Sigma_u(k,l);

      // compute pi
      double pi_sksl = biprobs::compute_pi(c_vec, pi_thresholds, rho_kl, k, l, sk, sl);

      // initialize empty pair-output
      double pair_ll = 0;
      Eigen::VectorXd pair_gradient = Eigen::VectorXd::Zero(d);


      Eigen::VectorXd pi_grad = grads::pi(theta,
                                          constrmat,
                                          constrsd,
                                          llc,
                                          c_vec,
                                          pi_thresholds,
                                          Sigma_u,
                                          Du,
                                          Ru,
                                          Lam.row(k),
                                          Lam.row(l),
                                          rho_kl,
                                          d,
                                          p,
                                          q,
                                          k,
                                          l,
                                          sk,
                                          sl,
                                          ncorr>0,
                                          nthr,
                                          nload,
                                          ncorr,
                                          nvar);

      gradmat.col(h) = (1/(pi_sksl+1e-8))*pi_grad;

      // update
      {
       gradient += (n_sksl/(pi_sksl+1e-8))*pi_grad;
       hessian +=  (n_sksl/(pow(pi_sksl,2)+1e-8))*pi_grad*pi_grad.transpose();
      }


    }
  }

  void patternwiseOuterProd::join(const patternwiseOuterProd &RHS){
    gradient += RHS.gradient;
    hessian += RHS.hessian;
  }


  // RcppParallel Worker to estimate H and J
  struct casewiseOuterProd : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Ref<const Eigen::MatrixXd> data;
    const Eigen::Ref<const Eigen::MatrixXd> gradmat;
    const Eigen::Ref<const Eigen::VectorXd> c_vec;

    // Output quantities:
    Eigen::VectorXd gradient;
    Eigen::MatrixXd variability;



    // Constructor 1:
    casewiseOuterProd(
      const Eigen::Ref<const Eigen::MatrixXd> DATA_,
      const Eigen::Ref<const Eigen::MatrixXd> GRADMAT_,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC_
    ):
      data(DATA_), gradmat(GRADMAT_), c_vec(C_VEC_){}

    // Constructor 2:
    casewiseOuterProd(const casewiseOuterProd &OBJ_, RcppParallel::Split):
      data(OBJ_.data), gradmat(OBJ_.gradmat), c_vec(OBJ_.c_vec){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const casewiseOuterProd &RHS);

  };

  void casewiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){

    const int d = gradmat.rows();
    const int p = data.cols();
    //subset_gradient.resize(d); subset_gradient.setZero();
    for (int i = BEGIN; i < END; i++){

      Eigen::VectorXd gradienti = Eigen::VectorXd::Zero(d);

      for(int k = 1; k < p; k++){
        const int ck = c_vec(k);

        // identify column index in freq table
        // i1: starting index item k
        int i1 = 0;
        if(k > 1){
          for(unsigned int u = 1; u < k; u++){
            const int cu = c_vec(u);
            i1 += cu * c_vec.segment(0,u).sum();
          }
        }



        const int sk = data(i, k);

        for(int l = 0; l < k; l++){
          const int cl = c_vec(l);

          // i2 starting index from i1 dor item l
          int i2 = 0;
          if(l > 0){
            i2 = c_vec.segment(0,l).sum() * c_vec(k);
          }

          // i3: starting index from i2 for cat sk
          unsigned int i3 = sk * cl;


          const int sl = data(i, l);

          // final column index for pairs_tab. Print to check
          const int r = i1 + i2 + i3 + sl;

          gradienti += gradmat.col(r);






    }
  }

      gradient += gradienti;
      variability += gradienti * gradienti.transpose();

    }
  }

  void casewiseOuterProd::join(const casewiseOuterProd &RHS){
    gradient += RHS.gradient;
    variability += RHS.variability;
  }
}

#endif
