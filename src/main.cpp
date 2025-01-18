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
#include "variance.h"
#include "exportedFuns.h"

//' Full pairwise iteration
//'
//' @description
//' Evaluate negative loglikelihood or gradient of the complete pool of pairs.
//' Used by external optimisers. Multithreading options via `RcppParallel`.
//'
//' @param N Number of observations
//' @param C_VEC Vector containing the number of categories for each item
//' @param CONSTRMAT Constraint matrix. Loadings free to be estimated are identified by a NA.
//' @param THETA Parameter vector
//' @param FREQ Frequency table
//' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
//' @param GRFLAG 0 to only compute the likelihood. 1 to also compute the gradient.
//' @param SILENTFLAG optional for verbose output
//'
// [[Rcpp::export]]
Rcpp::List cpp_multiThread_completePairwise(
    const unsigned int N,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRSD,
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::MatrixXd> FREQ,
    const int CORRFLAG,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int GRFLAG = 0,
    const int SILENTFLAG = 1
){

  const unsigned int d = THETA.size();
  const unsigned int n = N;                                             // number of units
  const unsigned int p = CONSTRMAT.rows();                                             // number of items
  const unsigned int q = CONSTRMAT.cols();                                             // number of latent variables
  const unsigned int c = C_VEC.sum();                                          // total number of categories
  const unsigned int nthr = c-p;                                        // number of thresholds
  unsigned int ncorr = 0; if(CORRFLAG==1) ncorr = q*(q-1)/2;                  // number of correlations
  const unsigned int nload = d - nthr - ncorr;                                    // number of loadings

  unsigned int R = p*(p-1)/2;                                            // number of pairs of items

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
  theta << THETA;                                 // Complete parameters vector

  // Initialize vector of indices for entries in pairs_table
  std::vector<int> vector_pairs(items_pairs.cols()) ;
  std::iota (std::begin(vector_pairs), std::end(vector_pairs), 0);

  double iter_ll = 0;
  Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

  pairs::SubsetWorker iteration_subset(CONSTRMAT, CONSTRSD, C_VEC, pairs_table, items_pairs, CORRFLAG, NTHR, NLOAD, NCORR, NVAR,
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

// //' Full pairwise likelihood
// //'
// //' @description
// //' Evaluate negative loglikelihood. Same structure of `multiThread_completePairwise`.
// //' Used to monitor validation log-likelihood. Multithreading options via `RcppParallel`.
// //'
// //' @param N Number of observations
// //' @param C_VEC Vector containing the number of categories for each item
// //' @param CONSTRMAT Constraint matrix. Loadings free to be estimated are identified by a 1.
// //' @param THETA Parameter vector
// //' @param FREQ Frequency table
// //' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
// //' @param SILENTFLAG optional for verbose output
// //'
// // [[Rcpp::export]]
// double multiThread_completePairwise_nll(
//     const unsigned int N,
//     Eigen::Map<Eigen::VectorXd> C_VEC,                // Vector containing the number of categories for each item
//     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,                    // Constraint matrix. Loadings free to be estimated are identified by a 1
//     Eigen::VectorXd THETA,
//     Eigen::Map<Eigen::MatrixXd> FREQ,
//     const unsigned int CORRFLAG,
//     const unsigned int SILENTFLAG
// ){
//
//   const unsigned int d = THETA.size();
//   const unsigned int n = N;                                             // number of units
//   const unsigned int p = CONSTRMAT.rows();                                             // number of items
//   const unsigned int q = CONSTRMAT.cols();                                             // number of latent variables
//   const unsigned int c = C_VEC.sum();                                          // total number of categories
//   const unsigned int nthr = c-p;                                        // number of thresholds
//   const unsigned int ncorr = q*(q-1)/2;                          // number of correlations
//   const unsigned int nload = d - nthr - ncorr;                                    // number of loadings
//
//   unsigned int R = p*(p-1)/2;                                            // number of pairs of items
//   unsigned int DFLAG, gradFLAG = 0;
//
//   // Copy frequencies, and build pair dictionary
//   Eigen::MatrixXd pairs_table = FREQ;
//   Eigen::MatrixXd items_pairs(2,R);
//   unsigned int r = 0;
//   for(unsigned int k = 1; k < p; k ++){
//     for(unsigned int l = 0; l < k; l++){
//       items_pairs(0, r) = k;
//       items_pairs(1, r) = l;
//       r++;
//     }
//   }
//
//   // Rearrange parameters
//   Eigen::VectorXd theta(d);
//   theta << THETA;                                 // Complete parameters vector
//
//   // Initialize vector of indeces for entries in pairs_table
//   std::vector<int> vector_pairs(items_pairs.cols()) ;
//   std::iota (std::begin(vector_pairs), std::end(vector_pairs), 0);
//
//   double iter_ll = 0;
//   Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);
//
//   pairs::SubsetWorker iteration_subset(CONSTRMAT, C_VEC, pairs_table, items_pairs, CORRFLAG, SILENTFLAG, 0, theta, vector_pairs);
//   RcppParallel::parallelReduce(0, R, iteration_subset);
//   iter_ll = -iteration_subset.subset_ll;
//
//   return(iter_ll/n);
// }

// //' Stochastic optimiser
// //'
// //' @description
// //' Core function to evaluate stochastic estimators
// //'
// //' @param FREQ Frequency table.
// //' @param VALFREQ Frequency table validation data.
// //' @param N Number of observations.
// //' @param C_VEC Vector containing the number of categories for each item.
// //' @param CONSTRMAT Constraint matrix. Loadings free to be estimated are identified by a 1.
// //' @param THETA_INIT Initial parameter vector
// //' @param CORRFLAG 1 to estimate latent correlations. 0 for orthogonal latent factors.
// //' @param METHODFLAG 0 for hypergeometric, 1 for Bernoulli.
// //' @param PAIRS_PER_ITERATION Number of pairs to draw per iteration.
// //' @param ETA Initial stepsize.
// //' @param BURN Initial burn-in period.
// //' @param MAXT Maximum number of iterations.
// //' @param TOLCOUNT Tolerance count for convergence with validation nll.
// //' @param SILENTFLAG Silent output.
// //' @param TOL Tolerance level.
// //' @param EACHCLOCK How often (in terms of iterations) to measure single iteration computational times (using \code{RcppClock}..
// //' @param PAR1 Hyperparameter for stepsize scheduling by Xu (2011): Scaling.
// //' @param PAR2 Hyperparameter for stepsize scheduling by Xu (2011): Smallest Hessian eigenvalue.
// //' @param PAR3 Hyperparameter for stepsize scheduling by Xu (2011): Decay rate.
// //' @param STEPSIZEFLAG Choose stepsize scheduling: Set 0 for Polyak and Juditsky (1992), 1 for Xu (2011).
// //' @param CHECKCONV Flag to check for convergence using complete pairwise likelihood on the validation set.
// //' @param EACHCHECK  How often (in terms of iterations) to check for convergence.
// //' @param SEED Randomising seed.
// //'
// // [[Rcpp::export]]
// Rcpp::List plFA(
//     Eigen::Map<Eigen::MatrixXd> FREQ,
//     Eigen::Map<Eigen::MatrixXd> VALFREQ,
//     const unsigned int N,
//     Eigen::Map<Eigen::VectorXd> C_VEC,
//     Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
//     Eigen::Map<Eigen::VectorXd> THETA_INIT,
//     const unsigned int CORRFLAG,
//     unsigned int METHODFLAG,
//     const unsigned int PAIRS_PER_ITERATION,                                          // Pairs drawn Per Iteration
//     double ETA,                                       // Proposed stepsize
//     const unsigned int BURN,
//     const unsigned int MAXT,                                      // Maximum number of iterations during the stochastic optimization
//     const unsigned int TOLCOUNT = 50,                                 // How may consecutive iterations need to satisfy the convergence check in order to declare convergence
//     const unsigned int SILENTFLAG = 1,
//     const double TOL = 1e-6,
//     const int EACHCLOCK= 100,
//     const double PAR1 = 1,
//     const double PAR2 = 1,
//     const double PAR3 = .75,
//     const unsigned int STEPSIZEFLAG = 0,
//     const unsigned int CHECKCONV = 0,
//     const unsigned int EACHCHECK = 100,
//     const unsigned int SEED = 123                                 // Random seed for sampling reproducibility
// ){
//
//   // Set up clock monitor to export to R session trough RcppClock
//   Rcpp::Clock clock;
//   clock.tick("Main");
//
//   if(SILENTFLAG == 0)Rcpp::Rcout << "Started!\n";
//
//   const unsigned int d = THETA_INIT.size();
//   const unsigned int n = N;                                             // number of units
//   const unsigned int p = CONSTRMAT.rows();                                             // number of items
//   const unsigned int q = CONSTRMAT.cols();                                             // number of latent variables
//   const unsigned int c = C_VEC.sum();                                          // total number of categories
//   const unsigned int nthr = c-p;                                        // number of thresholds
//   const unsigned int ncorr = q*(q-1)/2;                          // number of correlations
//   const unsigned int nload = d - nthr - ncorr;                                    // number of loadings
//   const unsigned int P = p*(p-1)/2;                                            // number of pairs of items
//
//   const double prob = static_cast<double>(PAIRS_PER_ITERATION)/static_cast<double>(P);
//   unsigned int convergence = 0;                                          // convergence check results
//   unsigned int tolGrad_counter = 0;
//   unsigned int tolPar_counter = 0;
//   unsigned int tolObj_counter = 0;
//
//   // Compute frequencies, items_pairs and items_pools
//   // Eigen::MatrixXd pairs_table(5,1);
//   Eigen::MatrixXd items_pairs(2,P);
//   // std::vector<std::vector<int>> items_pools(p);
//   // pairs_freq_cpp(DATA, C_VEC, pairs_table, items_pairs, items_pools);
//
//   unsigned int idx = 0;
//   for(unsigned int k = 1; k < p; k++){
//     const unsigned int ck = C_VEC(k);
//     for(unsigned int l = 0; l < k; l ++){
//
//       items_pairs(0, idx) = k;
//       items_pairs(1, idx) = l;
//       idx++;
//     }
//   }
//
//   Eigen::MatrixXd pairs_table = FREQ;
//
//   // Rearrange parameters
//   Eigen::VectorXd theta = THETA_INIT;
//   // theta << rep_tau, LAMBDA, TRANSFORMED_RHOS;                                 // Complete parameters vector
//   Eigen::MatrixXd Lam = params::get_Lam(CONSTRMAT, c, theta);                             // Loading matrix
//   Eigen::MatrixXd Sigma_u = params::get_S(theta, q);                            // Latent variable covariance matrix
//
//   // Initialize vector of indexes for entries in pairs_table
//   std::vector<int> outloop_pool;
//   std::vector<int> full_pool(P) ;
//   std::iota (std::begin(full_pool), std::end(full_pool), 0);
//   std::vector<int> chosen_pairs;                                              // Vector recording pairs chosen at each iteration
//
//   // Initialize storage for iterations quantities
//   Eigen::MatrixXd path_theta    = Eigen::MatrixXd::Zero(MAXT + 1, d); path_theta.row(0)    = theta;
//   Eigen::MatrixXd path_av_theta = Eigen::MatrixXd::Zero(MAXT + 1, d); path_av_theta.row(0) = theta;
//   Eigen::MatrixXd path_grad     = Eigen::MatrixXd::Zero(MAXT,     d);
//   Eigen::VectorXd path_nll      = Eigen::VectorXd::Zero(MAXT);
//   std::vector<int> post_index;                                                              // Track which iterations required post-update projection
//
//   // Initialize single iteration quantities
//   double nll = 0;
//   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
//   unsigned int last_iter = MAXT;
//
//   // std::vector<double> checkGrad;
//   std::vector<double> checkPar;
//   std::vector<double> check_val_nll;
//   std::vector<double> check_val_iter;
//
//
//   // Compute scaling constant
//   double scale;
//   switch(METHODFLAG){
//   case 0:
//     scale = static_cast<double>(P)/static_cast<double>(n*PAIRS_PER_ITERATION) ;
//     break;
//   case 1:
//     scale = prob/static_cast<double>(n);
//     break;
//   }
//
//   ///////////////////////
//   /* OPTIMISATION LOOP */
//   ///////////////////////
//   unsigned int conv_incr_counter = 0;
//   unsigned int conv_check_iterator = 1;
//   unsigned int sampling_window_iterator = 0;
//   for(unsigned int iter = 1; iter <= MAXT; iter++){
//     // check user interruption
//     Rcpp::checkUserInterrupt();
//     if(SILENTFLAG == 0) Rcpp::Rcout << "\rIteration:" << iter << " ";
//     if(iter % EACHCLOCK == 0) clock.tick("Iteration");
//
//     // Initialize empty contributions for iteration iter
//     double iter_nll = 0;
//     Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);
//
//     /////////////////////
//     // SAMPLING SCHEME //
//     /////////////////////
//     if(iter % EACHCLOCK == 0) clock.tick("Sampling_step");
//     std::vector<int> iter_chosen_pairs;
//     if(PAIRS_PER_ITERATION>=P){
//       // Complete pairwise likelihood
//       iter_chosen_pairs = full_pool;
//     }else{
//       // Stochastic sampling
//       iter_chosen_pairs = sampling_step(full_pool, METHODFLAG, prob, PAIRS_PER_ITERATION, p, SEED, SILENTFLAG, iter);
//     }
//
//     if(iter % EACHCLOCK == 0) clock.tock("Sampling_step");
//
//     ///////////////////////////
//     /* GRADIENT COMPUTATION  */
//     ///////////////////////////
//     if(iter % EACHCLOCK == 0) clock.tick("Stochastic_gradient");
//     pairs::SubsetWorker iteration_subset(CONSTRMAT, C_VEC, pairs_table, items_pairs, CORRFLAG, SILENTFLAG, 1, theta, iter_chosen_pairs);
//     RcppParallel::parallelReduce(0, iter_chosen_pairs.size(), iteration_subset);
//
//     iter_nll = -iteration_subset.subset_ll;
//     iter_gradient = -iteration_subset.subset_gradient;
//     iter_nll *= scale; iter_gradient *= scale;
//     if(iter % EACHCLOCK == 0) clock.tock("Stochastic_gradient");
//
//     ///////////////////////////
//     /*   PARAMETERS UPDATE  */
//     ///////////////////////////
//     if(iter % EACHCLOCK == 0) clock.tick("Update");
//     // double stepsize = ETA*pow(iter, -.5-1e-2);
//     // double stepsize = ETA * pow(1 + 1*pow(ETA,2)*iter, -.75);
//     // double stepsize = ETA*(1+.0002*ETA*pow(iter, -.75));
//
//     double stepsize = ETA;
//     switch(STEPSIZEFLAG){
//     case 0:
//       stepsize *= pow(iter, -PAR3);
//       break;
//     case 1:
//       stepsize *= PAR1 * pow(1 + PAR2*stepsize*iter, -PAR3);
//       break;
//     }
//     theta -= stepsize * iter_gradient;
//     if(iter % EACHCLOCK == 0) clock.tock("Update");
//
//     // 4. Post update loadings projection, if needed
//     bool checkevent = false;
//     stabilise_loadings(CONSTRMAT, C_VEC, theta, checkevent);
//     if(checkevent)post_index.push_back(iter);
//
//     /////////////////////////////////
//     /* STORE ITERATION QUANTITIES  */
//     /////////////////////////////////
//     path_theta.row(iter) = theta;
//     path_grad.row(iter-1) = iter_gradient;
//     path_nll(iter) = iter_nll;
//
//     /////////////////////////////////
//     /*    ONLINE AVERAGE UPDATE    */
//     /////////////////////////////////
//     if(iter <= BURN){
//       path_av_theta.row(iter) = path_theta.row(iter);
//     }else{
//       path_av_theta.row(iter) = ( (iter - BURN - 1) * path_av_theta.row(iter - 1) + path_theta.row(iter) ) / (iter - BURN);
//     }
//
//     if(iter % EACHCLOCK == 0) clock.tock("Iteration");
//     if(CHECKCONV == 1){
//       if(iter == BURN){
//         Eigen::VectorXd check_theta = path_av_theta.row(iter);
//         double val_nll = multiThread_completePairwise_nll(n, C_VEC, CONSTRMAT, check_theta, VALFREQ, CORRFLAG, SILENTFLAG);
//         check_val_nll.push_back(val_nll);
//         check_val_iter.push_back(iter);
//       }else if(iter > BURN){
//         if(conv_check_iterator == EACHCHECK){
//           Eigen::VectorXd check_theta = path_av_theta.row(iter);
//
//           double prev_check_nll = check_val_nll.back();
//           int prev_check_iter = check_val_iter.back();
//
//           double val_nll = multiThread_completePairwise_nll(n, C_VEC, CONSTRMAT, check_theta, VALFREQ, CORRFLAG, SILENTFLAG);
//           check_val_nll.push_back(val_nll);
//           check_val_iter.push_back(iter);
//
//           // If validation nll keeps decreasing then check tolerance
//           if(val_nll < prev_check_nll){
//             if((prev_check_nll - val_nll) < TOL){
//               convergence = 1;
//               last_iter = iter;
//               theta = path_av_theta.row(last_iter);
//               convergence = 1;
//               break;
//             }
//           }else{
//             // If validation nll increases
//             conv_incr_counter ++;
//             if(conv_incr_counter == TOLCOUNT){
//               auto minit = std::distance(std::begin(check_val_nll), std::min_element(std::begin(check_val_nll), std::end(check_val_nll)));
//               last_iter = check_val_iter.at(minit);
//               theta = path_av_theta.row(last_iter);
//               convergence = 2;
//               break;
//             }
//
//           }
//
//
//
//           conv_check_iterator = 1;
//
//         }else{
//           conv_check_iterator ++;
//         }
//       }
//
//     }
//   }
//
//
//   // output list
//   path_theta.conservativeResize(last_iter+1, d);
//   path_grad.conservativeResize(last_iter, d);
//   path_nll.conservativeResize(last_iter);
//
//   clock.tock("Main");
//   clock.stop("clock");
//
//   // output list
//   Rcpp::List output =
//     Rcpp::List::create(
//       Rcpp::Named("scale") = scale,
//       Rcpp::Named("path_theta") = path_theta,
//       Rcpp::Named("path_av_theta") = path_av_theta,
//       Rcpp::Named("path_grad") = path_grad,
//       Rcpp::Named("check_val_nll") = check_val_nll,
//       Rcpp::Named("check_val_iter") = check_val_iter,
//       Rcpp::Named("path_nll") = path_nll,
//       Rcpp::Named("post_index") = post_index,
//       Rcpp::Named("last_iter") = last_iter,
//       Rcpp::Named("theta") = theta,
//       Rcpp::Named("convergence") = convergence
//     );
//   return(output);
// }




