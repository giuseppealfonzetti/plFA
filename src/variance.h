#ifndef variance_H
#define variance_H

#include "bivariateProbs.h"
#include "genericUtils.h"
#include "gradients.h"

//' Estimate of H
//'
//' @description
//' Compute a sample estimate of the expected negative Hessian by taking
//' advantage of the second Bartlett's identity at the single pair level
//'
//' @param A Constraint matrix. Loadings free to be estimated are identified by a 1.
//' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
//' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
//' @param C_VEC Vector containing the number of categories for each item
//' @param THETA Parameter vector
//' @param FREQ output from [pairs_freq()]
//' @param N sample size
//' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
//' @param NTHR Number of thresholds parameters.
//' @param NLOAD Number of free loadings parameters
//' @param NCORR Number of free latent correlations parameters.
//' @param NVAR Number of free latent variance parameters.
//'
// [[Rcpp::export]]
Rcpp::List estimate_H(
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> A,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    int N,
    int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
){

  // Copy frequencies, and build pair dictionary
  Eigen::MatrixXd pairs_table = FREQ;
  pairs_table.conservativeResize(pairs_table.rows() + 1, Eigen::NoChange_t() );


  const int p = A.rows();
  const int q = A.cols();
  const int d = NTHR+NLOAD+NCORR+NVAR;

  // rearrange parameters
  Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, LLC, NTHR, NLOAD);
  Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);

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
          Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, C_VEC, k, l, sk, sl);

          // compute pi
          double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
          pairs_table(5,r) = pi_sksl;
          Eigen::VectorXd pi_grad = grads::pi(THETA,
                                              A,
                                              CONSTRLOGSD,
                                              LLC,
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
                                              k,
                                              l,
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
          gradient += (n_sksl/(pi_sksl+1e-8))*pi_grad;
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
//' @description
//' Compute a sample estimate of the variability matrix via the sample average outer product
//' of the composite score
//'
//' @param A Constraint matrix. Loadings free to be estimated are identified by a 1.
//' @param C_VEC Vector containing the number of categories for each item
//' @param THETA Parameter vector
//' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
//' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
//' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
//' @param Y data matrix
//' @param NTHR Number of thresholds parameters.
//' @param NLOAD Number of free loadings parameters
//' @param NCORR Number of free latent correlations parameters.
//' @param NVAR Number of free latent variance parameters.
//'
// [[Rcpp::export]]
Rcpp::List estimate_J(
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> A,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
    Eigen::VectorXd &THETA,
    int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
){


  const int p = A.rows();
  const int q = A.cols();
  const int d = NTHR+NLOAD+NCORR+NVAR;
  const int n = Y.rows();

  // rearrange parameters
  Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, LLC, NTHR, NLOAD);
  Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);

  Eigen::MatrixXd est_J = Eigen::MatrixXd::Zero(d,d);
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);

  for(unsigned int i = 0; i < n; i++){
    Eigen::VectorXd gradienti = Eigen::VectorXd::Zero(d);
    Rcpp::checkUserInterrupt();

    for(unsigned int k = 1; k < p; k++){
      const unsigned int ck = C_VEC(k);
      const Eigen::VectorXd lambdak = Lam.row(k);

      // identify column index in freq table
      // i1: starting index item k
      int i1 = 0;
      if(k > 1){
        for(unsigned int u = 1; u < k; u++){
          const unsigned int cu = C_VEC(u);
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

        // read frequency
        // const unsigned int n_sksl = pairs_table(4, r);

        // identify thresholds
        const Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, C_VEC, k, l, sk, sl);

        // compute pi
        const double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);


        // compute pi gradient contribution
        Eigen::VectorXd pi_grad = grads::pi(THETA,
                                            A,
                                            CONSTRLOGSD,
                                            LLC,
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
                                            k,
                                            l,
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


//' Estimate Diagonal of H
//'
//' @description
//' Compute a sample estimate of the expected negative Hessian by taking
//' advantage of the second Bartlett's identity at the single pair level
//'
//' @param A Constraint matrix. Loadings free to be estimated are identified by a 1.
//' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
//' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
//' @param C_VEC Vector containing the number of categories for each item
//' @param THETA Parameter vector
//' @param FREQ output from [pairs_freq()]
//' @param N sample size
//' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
//' @param NTHR Number of thresholds parameters.
//' @param NLOAD Number of free loadings parameters
//' @param NCORR Number of free latent correlations parameters.
//' @param NVAR Number of free latent variance parameters.
//'
// [[Rcpp::export]]
Eigen::VectorXd cpp_DH(
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> A,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    int N,
    int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR
){

  // Copy frequencies, and build pair dictionary
  Eigen::MatrixXd pairs_table = FREQ;
  pairs_table.conservativeResize(pairs_table.rows() + 1, Eigen::NoChange_t() );


  const int p = A.rows();
  const int q = A.cols();
  const int d = NTHR+NLOAD+NCORR+NVAR;

  // rearrange parameters
  Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, A, LLC, NTHR, NLOAD);
  Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
  Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
  Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);

  double ll = 0;
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
  Eigen::VectorXd diagH = Eigen::VectorXd::Zero(d);

  for(unsigned int k = 1; k < p; k++){
    unsigned int ck = C_VEC(k);
    Eigen::VectorXd lambdak = Lam.row(k);

    // identify column index in freq table
    // i1: starting index item k
    int i1 = 0;
    if(k > 1){
      for(int u = 1; u < k; u++){
        int cu = C_VEC(u);
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
          Eigen::VectorXd pi_thresholds = params::extract_thresholds(tau, C_VEC, k, l, sk, sl);

          // compute pi
          double pi_sksl = biprobs::compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
          pairs_table(5,r) = pi_sksl;
          Eigen::VectorXd pi_grad = grads::pi(THETA,
                                              A,
                                              CONSTRLOGSD,
                                              LLC,
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
                                              k,
                                              l,
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


          // gradient += (n_sksl/(pi_sksl+1e-8))*pi_grad;
          diagH    +=  (n_sksl/(pow(pi_sksl,2)+1e-8))*Eigen::VectorXd(pi_grad.array().square());

          // update ll
        }
      }
    }
  }


  diagH /=N;
  return(diagH);
}
#endif
