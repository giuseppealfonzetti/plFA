#ifndef fullPairwise_H
#define fullPairwise_H
#include "pairs.h"

namespace fullPairwise{

  Rcpp::List multiThreadContribution(
      const int N,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> FREQ,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR,
      const int GRFLAG,
      const int SILENTFLAG
  ){

    Rcpp::checkUserInterrupt();
    const int d = THETA.size();
    const int p = CONSTRMAT.rows();                                             // number of items
    int R = p*(p-1)/2;                                            // number of pairs of items

    // Copy frequencies, and build pair dictionary
    Eigen::MatrixXd items_pairs(2,R);
    int r = 0;
    for(int k = 1; k < p; k ++){
      for(int l = 0; l < k; l++){
        items_pairs(0, r) = k;
        items_pairs(1, r) = l;
        r++;
      }
    }

    // Rearrange parameters
    Eigen::VectorXd theta(d);
    theta << THETA;                                 // Complete parameters vector

    // Initialize vector of indices for entries in pairs_table
    std::vector<int> vector_pairs(items_pairs.cols()) ;
    std::iota (std::begin(vector_pairs), std::end(vector_pairs), 0);

    double iter_ll = 0;
    Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

    pairs::SubsetWorker iteration_subset(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, FREQ, items_pairs, CORRFLAG, NTHR, NLOAD, NCORR, NVAR,
                                         SILENTFLAG, GRFLAG, theta, vector_pairs);
    RcppParallel::parallelReduce(0, R, iteration_subset);
    iter_ll = iteration_subset.subset_ll;
    iter_gradient = iteration_subset.subset_gradient;

    // output list
    Rcpp::List output =
      Rcpp::List::create(
        Rcpp::Named("iter_nll") = -iter_ll,
        Rcpp::Named("iter_ngradient") = -iter_gradient
      );
    return(output);
  }
}

#endif
