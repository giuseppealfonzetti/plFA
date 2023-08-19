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
    const unsigned int CORRFLAG,
    const unsigned int GRFLAG,
    const unsigned int SILENTFLAG
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

/* MAIN FUNCTION  for STOCHASTIC OPTIMIZATION*/
//' @export
// [[Rcpp::export]]
Rcpp::List plFA(
    Eigen::Map<Eigen::MatrixXd> DATA,                    // Manifest data
    Eigen::Map<Eigen::VectorXd> C_VEC,                // Vector containing the number of categories for each item
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,                    // Constraint matrix. Loadings free to be estimated are identified by a 1
    Eigen::Map<Eigen::VectorXd> THETA_INIT,                  // Initial values for thresholds parameters
    const unsigned int CORRFLAG,
    unsigned int METHODFLAG,
    const unsigned int PAIRS_PER_ITERATION = 1,                                          // Pairs drawn Per Iteration
    double ETA = 1e-5,                                       // Proposed stepsize
    const unsigned int BURN = 200,
    const unsigned int MAXT = 1000,                                      // Maximum number of iterations during the stochastic optimization
    const unsigned int TOLCOUNT = 10,                                 // How may consecutive iterations need to satisfy the convergence check in order to declare convergence
    const unsigned int SEED = 123 ,                                  // Random seed for sampling reproducibility
    const unsigned int SILENTFLAG = 1,
    const bool CHECKCONVERGENCE = false,
    const double TOL = 1e-6,
    const int TOLCOUNTER = 10,
    const int EACHCLOCK= 100,
    const double PAR1 = 1,
    const double PAR2 = 1,
    const double PAR3 = .75,
    const unsigned int SAMPLING_WINDOW = 1,
    const unsigned int STEPSIZEFLAG = 0
){

  // Set up clock monitor to export to R session trough RcppClock
  Rcpp::Clock clock;
  clock.tick("Main");

  if(SILENTFLAG == 0)Rcpp::Rcout << "Started!\n";

  const unsigned int d = THETA_INIT.size();
  const unsigned int n = DATA.rows();                                             // number of units
  const unsigned int p = CONSTRMAT.rows();                                             // number of items
  const unsigned int q = CONSTRMAT.cols();                                             // number of latent variables
  const unsigned int c = C_VEC.sum();                                          // total number of categories
  const unsigned int nthr = c-p;                                        // number of thresholds
  const unsigned int ncorr = q*(q-1)/2;                          // number of correlations
  const unsigned int nload = d - nthr - ncorr;                                    // number of loadings
  const unsigned int P = p*(p-1)/2;                                            // number of pairs of items

  const double prob = static_cast<double>(PAIRS_PER_ITERATION)/static_cast<double>(P);
  unsigned int convergence = 0;                                          // convergence check results
  unsigned int tolGrad_counter = 0;
  unsigned int tolPar_counter = 0;
  unsigned int tolObj_counter = 0;

  // Compute frequencies, items_pairs and items_pools
  Eigen::MatrixXd pairs_table(5,1);
  Eigen::MatrixXd items_pairs(2,P);
  std::vector<std::vector<int>> items_pools(p);
  pairs_freq_cpp(DATA, C_VEC, pairs_table, items_pairs, items_pools);

  // Rearrange parameters
  Eigen::VectorXd theta = THETA_INIT;
  // theta << rep_tau, LAMBDA, TRANSFORMED_RHOS;                                 // Complete parameters vector
  Eigen::MatrixXd Lam = get_Lam(CONSTRMAT, c, theta);                             // Loading matrix
  Eigen::MatrixXd Sigma_u = get_S(theta, q);                            // Latent variable covariance matrix

  // Initialize vector of indexes for entries in pairs_table
  std::vector<int> outloop_pool;
  std::vector<int> full_pool(items_pairs.cols()) ;
  std::iota (std::begin(full_pool), std::end(full_pool), 0);
  std::vector<int> chosen_pairs;                                              // Vector recording pairs chosen at each iteration

  // Initialize storage for iterations quantities
  Eigen::MatrixXd path_theta    = Eigen::MatrixXd::Zero(MAXT + 1, d); path_theta.row(0)    = theta;
  Eigen::MatrixXd path_av_theta = Eigen::MatrixXd::Zero(MAXT + 1, d); path_av_theta.row(0) = theta;
  Eigen::MatrixXd path_grad     = Eigen::MatrixXd::Zero(MAXT,     d);
  Eigen::VectorXd path_nll      = Eigen::VectorXd::Zero(MAXT);
  std::vector<int> post_index;                                                              // Track which iterations required post-update projection

  // Initialize single iteration quantities
  double nll = 0;
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
  unsigned int last_iter = MAXT;

  // std::vector<double> checkGrad;
  std::vector<double> checkPar;
  // std::vector<double> checkObj;

  // Compute scaling constant
  double scale;
  switch(METHODFLAG){
  case 0:
    scale = static_cast<double>(P)/static_cast<double>(n*PAIRS_PER_ITERATION) ;
    break;
  case 1:
    scale = prob/static_cast<double>(n);
    break;
  }

  ///////////////////////
  /* OPTIMISATION LOOP */
  ///////////////////////
  unsigned int sampling_window_iterator = 0;
  for(unsigned int iter = 1; iter <= MAXT; iter++){
    // check user interruption
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "\rIteration:" << iter << " ";
    if(iter % EACHCLOCK == 0) clock.tick("Iteration");

    // Initialize empty contributions for iteration iter
    double iter_nll = 0;
    Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

    /////////////////////
    // SAMPLING SCHEME //
    /////////////////////
    if(iter % EACHCLOCK == 0) clock.tick("Sampling_step");
    std::vector<int> iter_chosen_pairs;
    // if(PAIRS_PER_ITERATION>=P){
    //   // Complete pairwise likelihood
    //   iter_chosen_pairs = full_pool;
    // }else{
    //   // Stochastic sampling
    //   iter_chosen_pairs = sampling_step(full_pool, METHODFLAG, prob, PAIRS_PER_ITERATION, p, SEED, SILENTFLAG, iter);
    // }

    if(PAIRS_PER_ITERATION>=P) METHODFLAG = 99;
    switch(METHODFLAG){
      case 99:
        {// Complete pairwise likelihood
        iter_chosen_pairs = full_pool;
        break;}
      case 0:
        {std::vector<int> tmp = hyper_sampling(P, SEED + iter);
         iter_chosen_pairs = {tmp.begin(), tmp.begin()+int(PAIRS_PER_ITERATION)};
         break;}
      case 1:
        {iter_chosen_pairs = bernoulli_sampling(P, prob);
         break;}
      case 2:
        {if(sampling_window_iterator == 0){
          outloop_pool = hyper_sampling(P, SEED + iter);
        }

        int tmp_ind = sampling_window_iterator*PAIRS_PER_ITERATION;
        iter_chosen_pairs = {outloop_pool.begin() + tmp_ind, outloop_pool.begin() + tmp_ind + int(PAIRS_PER_ITERATION)};
        sampling_window_iterator++;
        if(sampling_window_iterator==SAMPLING_WINDOW) sampling_window_iterator = 0;
        break;}


    }

    if(iter % EACHCLOCK == 0) clock.tock("Sampling_step");

    ///////////////////////////
    /* GRADIENT COMPUTATION  */
    ///////////////////////////
    if(iter % EACHCLOCK == 0) clock.tick("Stochastic_gradient");
    SubsetWorker iteration_subset(CONSTRMAT, C_VEC, pairs_table, items_pairs, CORRFLAG, SILENTFLAG, 1, theta, iter_chosen_pairs);
    RcppParallel::parallelReduce(0, iter_chosen_pairs.size(), iteration_subset);

    iter_nll = -iteration_subset.subset_ll;
    iter_gradient = -iteration_subset.subset_gradient;
    iter_nll *= scale; iter_gradient *= scale;
    if(iter % EACHCLOCK == 0) clock.tock("Stochastic_gradient");

    if(SILENTFLAG == 0){Rcpp::Rcout <<" |_ iter_nll:"<<iter_nll<<"\n";}

    ///////////////////////////
    /*   PARAMETERS UPDATE  */
    ///////////////////////////
    if(iter % EACHCLOCK == 0) clock.tick("Update");
    // double stepsize = ETA*pow(iter, -.5-1e-2);
    // double stepsize = ETA * pow(1 + 1*pow(ETA,2)*iter, -.75);
    // double stepsize = ETA*(1+.0002*ETA*pow(iter, -.75));

    double stepsize = ETA;
    switch(STEPSIZEFLAG){
    case 0:
      stepsize *= pow(iter, -PAR3);
      break;
    case 1:
      stepsize *= PAR1 * pow(1 + PAR2*stepsize*iter, -PAR3);
      break;
    }
    theta -= stepsize * iter_gradient;
    if(iter % EACHCLOCK == 0) clock.tock("Update");

    // 4. Post update loadings normalisation, if needed
    if(SILENTFLAG == 0)Rcpp::Rcout <<" |_ post-";
    bool checkevent = false;
    stabilise_loadings(CONSTRMAT, C_VEC, theta, checkevent);
    if(checkevent)post_index.push_back(iter);

    /////////////////////////////////
    /* STORE ITERATION QUANTITIES  */
    /////////////////////////////////
    path_theta.row(iter) = theta;
    path_grad.row(iter-1) = iter_gradient;
    path_nll(iter) = iter_nll;

    /////////////////////////////////
    /*    ONLINE AVERAGE UPDATE    */
    /////////////////////////////////
    if(iter <= BURN){
      path_av_theta.row(iter) = path_theta.row(iter);
    }else{
      path_av_theta.row(iter) = ( (iter - BURN - 1) * path_av_theta.row(iter - 1) + path_theta.row(iter) ) / (iter - BURN);
    }
    /////////////////////////////////
    /*      CHECK CONVERGENCE      */
    /////////////////////////////////
    if((iter > BURN) & CHECKCONVERGENCE){

      // check norm diff parameter vector
      {
        const double parNorm = (path_av_theta.row(iter)-path_av_theta.row(iter-1)).norm()/(path_av_theta.row(iter-1)).norm();
        if(parNorm < TOL ){tolPar_counter ++; }else{tolPar_counter = 0;}
        checkPar.push_back(parNorm);
        if(tolPar_counter == TOLCOUNTER){last_iter = iter; convergence = 1; break;}
      }




    }
    if(iter % EACHCLOCK == 0) clock.tock("Iteration");

  }


  // output list
  path_theta.conservativeResize(last_iter+1, d);
  path_grad.conservativeResize(last_iter, d);
  path_nll.conservativeResize(last_iter);

  clock.tock("Main");
  clock.stop("clock");

  // output list
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("scale") = scale,
      Rcpp::Named("path_theta") = path_theta,
      Rcpp::Named("path_av_theta") = path_av_theta,
      Rcpp::Named("path_grad") = path_grad,
      Rcpp::Named("checkPar") = checkPar,
      Rcpp::Named("path_nll") = path_nll,
      Rcpp::Named("post_index") = post_index,
      Rcpp::Named("last_iter") = last_iter,
      Rcpp::Named("theta") = theta,
      Rcpp::Named("convergence") = convergence
    );
  return(output);
}
