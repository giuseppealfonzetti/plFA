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
//' Evaluate negative log pairwise likelihood or gradient of the complete pool of pairs.
//' Used by external optimisers. Multithreading options via `RcppParallel`.
//'
//' @param N Number of observations
//' @param C_VEC Vector containing the number of categories for each item
//' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA` refers to free loading parameters. Elements set to numerical values denote fixed values constraints.
//' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
//' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
//' @param THETA Parameter vector
//' @param FREQ Frequency table
//' @param CORRFLAG TRUE to estimate latent correlations. 0 for orthogonal latent factors.
//' @param NTHR Number of thresholds parameters.
//' @param NLOAD Number of free loadings parameters
//' @param NCORR Number of free latent correlations parameters.
//' @param NVAR Number of free latent variance parameters.
//' @param GRFLAG 0 to only compute the likelihood. 1 to also compute the gradient.
//' @param SILENTFLAG optional for verbose output
//'
// [[Rcpp::export]]
Rcpp::List cpp_multiThread_completePairwise(
    const unsigned int N,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
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

  Rcpp::checkUserInterrupt();
  const unsigned int d = THETA.size();
  const unsigned int p = CONSTRMAT.rows();                                             // number of items
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

  pairs::SubsetWorker iteration_subset(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, pairs_table, items_pairs, CORRFLAG, NTHR, NLOAD, NCORR, NVAR,
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






// [[Rcpp::export]]
Rcpp::List cpp_plSA(
    Eigen::Map<Eigen::MatrixXd> FREQ,
    Eigen::Map<Eigen::MatrixXd> VALFREQ,
    const int N,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
    const std::vector<std::vector<std::vector<double>>> LLC,
    Eigen::Map<Eigen::VectorXd> THETA_INIT,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int SAMPLER,
    const int PAIRS_PER_ITERATION,
    const int SCHEDULE,
    const double STEP0,
    const double STEP1,
    const double STEP2,
    const double STEP3,
    const int BURN,
    const int MAXT,
    const int TOL_WINDOW,
    const double TOL_NLL,
    const int CHECK_TOL,
    const int CHECK_WINDOW,
    const int PATH_WINDOW,
    const int CLOCK_WINDOW,
    const int SEED,
    const bool VERBOSE
){

    // Set up clock monitor to export to R session trough RcppClock
    Rcpp::Clock clock;
    clock.tick("Main");

    if(VERBOSE) Rcpp::Rcout << "Started!\n";


    const int d = THETA_INIT.size();
    const int n = N;
    const int p = CONSTRMAT.rows();
    // const int q = CONSTRMAT.cols();
    // const int c = C_VEC.sum();
    const int pairs = p*(p-1)/2;
    const int corrflag = NCORR>0;
    if(d!=NTHR+NLOAD+NCORR+NVAR) Rcpp::stop("check theta dimensions");
    if(PAIRS_PER_ITERATION>pairs)Rcpp::stop("too many pairs per iterations");

    int maxt               = MAXT;
    int burn               = BURN;
    const double prob      = static_cast<double>(PAIRS_PER_ITERATION)/static_cast<double>(pairs);
    int tol_counter        = 0;
    int tol_counter_burn   = 0;

    double stepsize   = STEP0;

    // Read frequencies, and initialise items_pairs
    // Eigen::MatrixXd pairs_table = FREQ;
    Eigen::MatrixXd items_pairs(2,pairs);
    unsigned int idx = 0;
    for(unsigned int k = 1; k < p; k++){
      for(unsigned int l = 0; l < k; l ++){
        items_pairs(0, idx) = k;
        items_pairs(1, idx) = l;
        idx++;
      }
    }

    // Initialize vector of indexes for entries in pairs_table
    std::vector<int> outloop_pool;
    std::vector<int> full_pool(pairs) ;
    std::iota (std::begin(full_pool), std::end(full_pool), 0);
    std::vector<int> chosen_pairs;

    // Initialize storage for iterations quantities
    Eigen::VectorXd theta   = THETA_INIT;
    Eigen::VectorXd avtheta = THETA_INIT;

    double          nll   = std::numeric_limits<double>::infinity();
    std::vector<Eigen::VectorXd> path_theta;
    std::vector<Eigen::VectorXd> path_avtheta;
    std::vector<Eigen::VectorXd> path_grad;
    std::vector<double>          path_nll;
    std::vector<int>             path_iters;
    std::vector<int>             post_index;


    // Initialize single iteration quantities
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
    int last_iter = maxt;

    // Compute scaling constant
    double scale = prob/static_cast<double>(n);


    // Diagnostics
    bool neg_pdiff        = false;
    bool proj_after_burn  = false;
    bool convergence_burn = true;
    int  convergence_full = 0;



    for(int t = 0; t <= maxt; t++){
      Rcpp::checkUserInterrupt();

      if(t % CLOCK_WINDOW == 0) clock.tick("Iteration");



      /////////////////////////
      // STORE PREVIOUS ITER //
      /////////////////////////
      if(((t) % PATH_WINDOW == 0) || (CHECK_TOL && tol_counter == CHECK_WINDOW) || (t == maxt)){
        path_iters.push_back(t);
        path_theta.push_back(theta);
        path_avtheta.push_back(avtheta);
        path_nll.push_back(nll);
      }

      // Initialize empty contributions for iteration
      Eigen::VectorXd iter_gradient = Eigen::VectorXd::Zero(d);

      ///////////////////////
      // COMPUTE FULL NLL //
      ///////////////////////
      if((CHECK_TOL && (t % CHECK_WINDOW == 0)) || (t == maxt)){
        if(t % CLOCK_WINDOW == 0) clock.tick("Obj eval");

        double prev_nll = nll;
        pairs::SubsetWorker fullpool_worker(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, VALFREQ, items_pairs, corrflag, NTHR, NLOAD, NCORR, NVAR,
                                             1, 0, theta, full_pool);
        RcppParallel::parallelReduce(0, pairs, fullpool_worker);
        nll = -fullpool_worker.subset_ll/n;


        // Check convergence
        double pdiff = - (nll - prev_nll)/prev_nll;
        if(pdiff<TOL_NLL){tol_counter++;}else{tol_counter=0;}
        if(pdiff<(TOL_NLL*10)){tol_counter_burn++;}else{tol_counter_burn=0;}
        if(pdiff<0 & t<=burn) neg_pdiff=true;



        // Verbose
        if(VERBOSE){
          Rcpp::Rcout << "Iter: " << t << " | step: " << stepsize << " | full npll: "<< nll << " | pdiff: " << pdiff << " | tol_counter:" << tol_counter<< "\n";
        }

        if((t>0) && !std::isfinite(nll)){
          last_iter       = t;
          convergence_full= -1;
          break;
        }
        if(t % CLOCK_WINDOW == 0) clock.tock("Obj eval");
      }

      /////////////////
      // AUTO BURNIN //
      /////////////////
      if(t==burn){
        if(VERBOSE) Rcpp::Rcout << "Burn-in ended at iter " << t << "\n";
        convergence_burn=false;
      }else if((CHECK_TOL && (tol_counter_burn >= TOL_WINDOW) && (t <= burn))){
        if(VERBOSE) Rcpp::Rcout << "Burn-in stopped at iter " << t << " | full npll: "<< nll << "\n";
        burn = t;
        tol_counter_burn=0;
        convergence_burn=true;
      }




      /////////////
      // EXITING //
      /////////////
      if(CHECK_TOL && (tol_counter>=TOL_WINDOW) && ( t>=(burn+CHECK_WINDOW))){
        if(VERBOSE) Rcpp::Rcout << "Converged at iter " << t << " | full npll: "<< nll << "\n";
        convergence_full = 1;
        last_iter = t;
        break;
      }else if(t==maxt){
        if(VERBOSE) Rcpp::Rcout << "Reached iter " << t << " | full npll: "<< nll << "\n";
        last_iter = t;
      }

      /////////////////////
      // SAMPLING SCHEME //
      /////////////////////
      if(t % CLOCK_WINDOW == 0) clock.tick("Sampling_step");
      std::vector<int> iter_chosen_pairs;
      iter_chosen_pairs = sa::sampling_step(full_pool, SAMPLER, prob, PAIRS_PER_ITERATION, p, SEED, false, t);
      if(t % CLOCK_WINDOW == 0) clock.tock("Sampling_step");



      ///////////////////////////
      // GRADIENT COMPUTATION  //
      ///////////////////////////
      if(t % CLOCK_WINDOW == 0) clock.tick("Stochastic_gradient");
      pairs::SubsetWorker iteration_subset(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, FREQ, items_pairs, corrflag,
                                           NTHR, NLOAD, NCORR, NVAR,
                                           1, 1, theta, iter_chosen_pairs);
      RcppParallel::parallelReduce(0, iter_chosen_pairs.size(), iteration_subset);

      iter_gradient = -iteration_subset.subset_gradient;
      iter_gradient *= scale;
      if(t % CLOCK_WINDOW == 0) clock.tock("Stochastic_gradient");


      //////////////////////////
      //   PARAMETERS UPDATE  //
      //////////////////////////
      if(t % CLOCK_WINDOW == 0) clock.tick("Update");
      switch(SCHEDULE){
      case 0:
        stepsize *= pow(t+1, -STEP3);
        break;
      case 1:
        stepsize *= STEP1 * pow(1 + STEP2*stepsize*(t+1), -STEP3);
        break;
      }


      // iter_gradient.segment(0, NTHR) /=10;
      iter_gradient.segment(NTHR, NLOAD) *= 5;
      if(NCORR>0) iter_gradient.segment(NTHR+NLOAD, NCORR) *= 100;
      if(NVAR>0)  iter_gradient.segment(NTHR+NLOAD+NCORR, NVAR) *= 5;


      theta -= stepsize * iter_gradient;
      if(t % CLOCK_WINDOW == 0) clock.tock("Update");

      ////////////////////////////////////////////////
      // Post update loadings projection, if needed //
      ///////////////////////////////////////////////

      bool checkevent;
      sa::proj2(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, corrflag, NTHR, NLOAD, NCORR, NVAR, theta, checkevent);
      if(checkevent) post_index.push_back(t);
      if(VERBOSE && checkevent && t>burn) proj_after_burn=true;



      if(t <= burn){
        avtheta = theta;
      }else{
        avtheta = ( (t - burn) * avtheta + theta ) / (t - burn + 1);
      }

      if(t % CLOCK_WINDOW == 0) clock.tock("Iteration");


    }















    clock.tock("Main");
    clock.stop("clock");


    Rcpp::List output =
      Rcpp::List::create(
        Rcpp::Named("scale") = scale,
        Rcpp::Named("nll")   = nll,
        Rcpp::Named("path_theta") = path_theta,
        Rcpp::Named("path_avtheta") = path_avtheta,
        Rcpp::Named("path_iters") = path_iters,
        Rcpp::Named("path_nll") = path_nll,
        Rcpp::Named("post_index") = post_index,
        Rcpp::Named("last_iter") = last_iter,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("avtheta") = avtheta,
        Rcpp::Named("convergence") = convergence_full,
        Rcpp::Named("convergence_burn") = convergence_burn,
        Rcpp::Named("proj_after_burn") = proj_after_burn,
        Rcpp::Named("neg_pdiff") = neg_pdiff
      );
    return(output);



}

