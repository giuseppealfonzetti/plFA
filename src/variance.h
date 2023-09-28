#ifndef variance_H
#define variance_H

#include "bivariateProbs.h"
//' Estimate of H
//'
//' Compute a sample estimate of the expected negative Hessian by taking
//' advantage of the second Bartlett's identity at the single pair level
//'
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_H(
    Eigen::Map<Eigen::VectorXd> C_VEC,                // Vector containing the number of categories for each item
    Eigen::Map<Eigen::MatrixXd> A,                    // Constraint matrix. Loadings free to be estimated are identified by a 1
    Eigen::VectorXd &THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    int N,
    int CORRFLAG
){
  unsigned int p = A.rows();                                             // number of items
  unsigned int q = A.cols();                                             // number of latents
  unsigned int d = THETA.size();                                         // number of parameters
  unsigned int ncorr = q*(q-1)/2;                                        // number of correlations
  unsigned int c = C_VEC.sum();                                          // total number of categories
  unsigned int nthr = c-p;                                               // number of thresholds
  unsigned int nload = d - ncorr - nthr;                                 // number of loadings
  unsigned int R = p*(p-1)/2;                                            // number of pairs of items


  // Copy frequencies, and build pair dictionary
  Eigen::MatrixXd pairs_table = FREQ;
  pairs_table.conservativeResize(pairs_table.rows() + 1, Eigen::NoChange_t() );


  // rearrange parameters
  Eigen::MatrixXd Lam            = get_Lam(A, c, THETA);
  Eigen::MatrixXd Sigma_u        = get_S(THETA, q);
  Eigen::VectorXd tau            = THETA.segment(0,c-p);

  double ll = 0;
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
  Eigen::MatrixXd est_H = Eigen::MatrixXd::Zero(d,d);
  for(unsigned int k = 1; k < p; k++){
    unsigned int ck = C_VEC(k);
    Eigen::VectorXd lambdak = Lam.row(k);

    // identify column index in freq table
    // i1: starting index item k
    int i1 = 0;
    if(k > 1){
      for(int u = 1; u < k; u++){
        int cu = C_VEC(u);
        //if(silentFLAG == 0)Rcpp::Rcout << "u: " << u << ", cu: "<< cu << "\n";
        i1 += cu * C_VEC.segment(0,u).sum();
      }
    }

    for(unsigned int l = 0; l < k; l++){
      unsigned int cl = C_VEC(l);
      Eigen::VectorXd lambdal = Lam.row(l);

      double rho_kl = lambdak.transpose() * Sigma_u * lambdal;


      // i2 starting index from i1 dor item l
      int i2 = 0;
      if(l > 0){
        i2 = C_VEC.segment(0,l).sum() * C_VEC(k);
      }
      for(unsigned int sk = 0; sk < ck; sk ++){

        // i3: starting index from i2 for cat sk
        unsigned int i3 = sk * cl;

        for(unsigned int sl = 0; sl < cl; sl ++){

          // final column index for pairs_tab. Print to check
          unsigned int r = i1 + i2 + i3 + sl;

          // read frequency
          unsigned int n_sksl = pairs_table(4, r);

          // identify thresholds
          Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);

          // compute pi
          double pi_sksl = compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
          //if(SILENTFLAG == 0)Rcpp::Rcout << "("<<k<<","<<l<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
          pairs_table(5,r) = pi_sksl;
          Eigen::VectorXd pi_grad = compute_pi_grad(A, C_VEC, pi_thresholds, Sigma_u, Lam, THETA, rho_kl, k, l, sk, sl, CORRFLAG);

          gradient += (n_sksl/(pi_sksl+1e-8))*pi_grad;
          // Eigen::MatrixXd tempH = pi_grad*pi_grad.transpose();
          est_H +=  (n_sksl/(pow(pi_sksl,2)+1e-8))*pi_grad*pi_grad.transpose();

          // update ll
          ll += n_sksl * log(pi_sksl+1e-8);
        }
      }
    }
  }


  est_H /=N;
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("ll") = ll,
      Rcpp::Named("gradient") = gradient,
      Rcpp::Named("est_H") = est_H
    );

  return(output);
}

//' Estimate of J
//'
//' Compute a sample estimate of the variability matrix via the sample average outer product
//' of the composite score
//'
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_J(
    Eigen::Map<Eigen::MatrixXd> Y,                    // Manifest data
    Eigen::Map<Eigen::VectorXd> C_VEC,                // Vector containing the number of categories for each item
    Eigen::Map<Eigen::MatrixXd> A,                    // Constraint matrix. Loadings free to be estimated are identified by a 1
    Eigen::VectorXd &THETA,
    int CORRFLAG
){
  unsigned int p = A.rows();                                             // number of items
  unsigned int q = A.cols();                                             // number of latents
  unsigned int d = THETA.size();                                         // number of parameters
  unsigned int ncorr = q*(q-1)/2;                                        // number of correlations
  unsigned int c = C_VEC.sum();                                          // total number of categories
  unsigned int nthr = c-p;                                               // number of thresholds
  unsigned int nload = d - ncorr - nthr;                                 // number of loadings
  unsigned int R = p*(p-1)/2;                                            // number of pairs of items
  unsigned int n = Y.rows();


  // rearrange parameters
  Eigen::MatrixXd Lam            = get_Lam(A, c, THETA);
  Eigen::MatrixXd Sigma_u        = get_S(THETA, q);
  Eigen::VectorXd tau            = THETA.segment(0,c-p);

  double ll = 0;
  Eigen::MatrixXd est_J = Eigen::MatrixXd::Zero(d,d);
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);

  for(unsigned int i = 0; i < n; i++){
    Eigen::VectorXd gradienti = Eigen::VectorXd::Zero(d);

    for(unsigned int k = 1; k < p; k++){
      const unsigned int ck = C_VEC(k);
      const Eigen::VectorXd lambdak = Lam.row(k);

      // identify column index in freq table
      // i1: starting index item k
      int i1 = 0;
      if(k > 1){
        for(unsigned int u = 1; u < k; u++){
          const unsigned int cu = C_VEC(u);
          //if(silentFLAG == 0)Rcpp::Rcout << "u: " << u << ", cu: "<< cu << "\n";
          i1 += cu * C_VEC.segment(0,u).sum();
        }
      }

      unsigned int sk = Y(i, k);

      for(unsigned int l = 0; l < k; l++){
        const unsigned int cl = C_VEC(l);
        const Eigen::VectorXd lambdal = Lam.row(l);

        const double rho_kl = lambdak.transpose() * Sigma_u * lambdal;


        // i2 starting index from i1 dor item l
        int i2 = 0;
        if(l > 0){
          i2 = C_VEC.segment(0,l).sum() * C_VEC(k);
        }

        unsigned int sl = Y(i, l);


        // i3: starting index from i2 for cat sk
        const unsigned int i3 = sk * cl;

        // final column index for pairs_tab. Print to check
        const unsigned int r = i1 + i2 + i3 + sl;

        // read frequency
        // const unsigned int n_sksl = pairs_table(4, r);

        // identify thresholds
        const Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);

        // compute pi
        const double pi_sksl = compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);

        // compute pi gradient contribution
        Eigen::VectorXd pi_grad = compute_pi_grad(A, C_VEC, pi_thresholds, Sigma_u, Lam, THETA, rho_kl, k, l, sk, sl, CORRFLAG);

        // Update gradient
        gradienti += (1/(pi_sksl+1e-8))*pi_grad;
      }
    }
    gradient += gradienti;
    est_J += gradienti * gradienti.transpose();
  }


  est_J /=n;
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gradient") = gradient,
      Rcpp::Named("est_J") = est_J
    );

  return(output);
}

#endif
