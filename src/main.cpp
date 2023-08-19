#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <RcppClock.h>
#include <random>
#include <math.h>

// [[Rcpp::depends(RcppEigen, RcppParallel, RcppClock)]]
#include "genericUtils.h"
#include "frequencies.h"
#include "pairs.h"
#include "optimisationUtils.h"

//' Complete pairiwse iteration with multithreading option//'
//' Used by external optimisers
// [[Rcpp::export]]
Rcpp::List multiThread_completePairwise(
    Eigen::Map<Eigen::MatrixXd> Y,                    // Manifest data
    Eigen::Map<Eigen::VectorXd> C_VEC,                // Vector containing the number of categories for each item
    Eigen::Map<Eigen::MatrixXd> A,                    // Constraint matrix. Loadings free to be estimated are identified by a 1
    Eigen::Map<Eigen::VectorXd> TAU,                  // Initial values for thresholds parameters
    Eigen::Map<Eigen::VectorXd> LAMBDA,               // Initial values for loadings parameters
    Eigen::Map<Eigen::VectorXd> TRANSFORMED_RHOS,     // Initial values for latent correlations reparameterized trough Fisher's transformation
    Eigen::Map<Eigen::MatrixXd> FREQ,
    int CORRFLAG,
    int GRFLAG,
    int SILENTFLAG
){


  Eigen::VectorXd rep_tau = TAU;
  unsigned int n = Y.rows();                                             // number of units
  unsigned int p = A.rows();                                             // number of items
  unsigned int q = A.cols();                                             // number of latents
  unsigned int ncorr = TRANSFORMED_RHOS.size();                          // number of correlations
  unsigned int nthr = TAU.size();                                        // number of thresholds
  unsigned int nload = LAMBDA.size();                                    // number of loadings
  unsigned int d = nload + nthr + ncorr;                                 // number of parameters
  unsigned int c = C_VEC.sum();                                          // total number of categories
  unsigned int R = p*(p-1)/2;                                            // number of pairs of items
  unsigned int DFLAG, gradFLAG = 0;

  if(GRFLAG == 1){
    gradFLAG = 1; DFLAG = 0;
  }else if(GRFLAG==2){
    gradFLAG = 1; DFLAG = 1;
  }

  // Copy frequencies, and build pair dictionary
  Eigen::MatrixXd pairs_table = FREQ;
  Eigen::MatrixXd items_pairs(2,R);
  unsigned int r = 0;
  for(unsigned int k = 1; k < p; k ++){
    for(unsigned int l = 0; l < k; l++){
      items_pairs(0, r) = k;
      items_pairs(1, r) = l;
      r++;
    }
  }

  // Rearrange parameters
  Eigen::VectorXd theta(d);
  theta << rep_tau, LAMBDA, TRANSFORMED_RHOS;                                 // Complete parameters vector

  // Initialize vector of indeces for entries in pairs_table
  std::vector<int> vector_pairs(items_pairs.cols()) ;
  std::iota (std::begin(vector_pairs), std::end(vector_pairs), 0);

  double iter_ll = 0;
  Eigen::VectorXd iter_gradient(d); iter_gradient.fill(0.0);
  Eigen::VectorXd iter_gradient2(d); iter_gradient2.fill(0.0);

  SubsetWorker iteration_subset(A, C_VEC, pairs_table, items_pairs, CORRFLAG, SILENTFLAG, gradFLAG, theta, vector_pairs);
  RcppParallel::parallelReduce(0, R, iteration_subset);
  iter_ll = iteration_subset.subset_ll;
  iter_gradient = iteration_subset.subset_gradient;
  iter_gradient2 = iteration_subset.subset_gradient2;
  Eigen::MatrixXd H_approx = iter_gradient2.asDiagonal();
  // output list
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("iter_nll") = -iter_ll,
      Rcpp::Named("iter_ngradient") = -iter_gradient,
      Rcpp::Named("H_approx") = -H_approx
    );
  return(output);
}

