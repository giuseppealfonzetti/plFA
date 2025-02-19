#include "parallelWorkers.h"

void parallelWorkers::pairsGradient::operator()(std::size_t BEGIN, std::size_t END){
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

void parallelWorkers::pairsGradient::join(const pairsGradient &RHS){
  subset_ll += RHS.subset_ll;
  for(int i = 0; i < subset_gradient.size(); i++){
    subset_gradient(i) += RHS.subset_gradient(i);

  }
}




void parallelWorkers::patternwiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){
  // Run along the pairs identified by indexes in index_vector
  // computing their contribution to nll, gradient and eventually Hessian.
  const int d = theta.size();
  const int q = constrmat.cols();
  const int p = constrmat.rows();

  // Rcpp::Rcout << "d:"<<d<<",q:"<<q<<", p:"<<p<<"\n";
  // rearrange parameters
  const Eigen::MatrixXd Lam             = params::loadings::theta2mat(theta, constrmat, llc, nthr, nload);
  const Eigen::MatrixXd Ru              = params::latvar::theta2cmat(theta, nthr, nload, ncorr, nvar, q);
  const Eigen::MatrixXd Du              = params::latvar::theta2dmat(theta, constrsd, nthr, nload, ncorr, nvar, q);
  const Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  const Eigen::VectorXd tau             = params::thresholds::theta2vec(theta, nthr);
  const Eigen::MatrixXd Sigma_y         = Lam*Sigma_u*Lam.transpose();

  // Rcpp::Rcout << "theta:\n";
  // Rcpp::Rcout << theta.transpose()<<"\n";
  // Rcpp::Rcout << "Lam:\n";
  // Rcpp::Rcout << Lam<<"\n";
  // Rcpp::Rcout << "Sigma_u:\n";
  // Rcpp::Rcout << Sigma_u<<"\n";
  // Rcpp::Rcout << "Sigma_y:\n";
  // Rcpp::Rcout << Sigma_y<<"\n";
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
    double rho_kl = Sigma_y(k,l);

    // compute pi
    double pi_sksl = biprobs::compute_pi(c_vec, pi_thresholds, rho_kl, k, l, sk, sl);

    // Rcpp::Rcout << "k:"<<k<<", l:"<<l<<", sk:"<<sk<<", sl:"<< sl <<", n:"<<n_sksl<<", rho:"<<rho_kl<<", pi:"<<pi_sksl<<"\n";

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
      // gradient +=pi_grad;

      hessian +=  (n_sksl/(pow(pi_sksl,2)+1e-8))*pi_grad*pi_grad.transpose();
      // hessian +=  pi_grad*pi_grad.transpose();

    }


  }
}

void parallelWorkers::patternwiseOuterProd::join(const patternwiseOuterProd &RHS){
  gradient += RHS.gradient;
  hessian += RHS.hessian;
}


void parallelWorkers::casewiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){

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

void parallelWorkers::casewiseOuterProd::join(const casewiseOuterProd &RHS){
  gradient += RHS.gradient;
  variability += RHS.variability;
}
