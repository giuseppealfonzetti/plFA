#ifndef optimisationUtils_H
#define optimisationUtils_H
#include <random>

#include "genericUtils.h"

namespace sa{
  /* SAMPLING PAIRS */
  std::vector<int> sampling_step(
      std::vector<int> &FULL_POOL,
      const unsigned int METHOD,                    // 0 hypergeometric, 1 Bernoulli
      const double PROB,
      const int PAIRS_PER_ITERATION,
      const unsigned int N_ITEMS,
      const unsigned int SEED,
      const bool VERBOSE,
      const unsigned int ITER
  ){
    // Set-up the randomizer
    std::mt19937 randomizer(SEED + ITER);

    // empty vector to collect chosen indices
    std::vector<int> chosen_index;

    if(METHOD == 0){
      std::shuffle(FULL_POOL.begin(), FULL_POOL.end(), randomizer);

      for(unsigned int draw = 0; draw < PAIRS_PER_ITERATION; draw++){
        if(VERBOSE)Rcpp::Rcout << "Draw "<< draw + 1 <<":\n |_ Shuffling the pairs...\n";

        // Collect the index drawn
        unsigned int pair_index = FULL_POOL[draw];
        chosen_index.push_back(pair_index);
        if(VERBOSE)Rcpp::Rcout << " |_ Position drawn:" << pair_index << "\n";

      }

    }else if(METHOD == 1){
      for(unsigned int pair = 0; pair < FULL_POOL.size(); pair++){
        if(R::runif(0,1) < PROB ) chosen_index.push_back(pair);

      }

    }

    return chosen_index;
  }

  std::vector<int> hyper_sampling(const unsigned int K, const unsigned int SEED){
    std::mt19937 randomizer(SEED);
    std::vector<int> pool(K);
    std::iota (std::begin(pool), std::end(pool), 0);
    std::shuffle(pool.begin(), pool.end(), randomizer);
    return(pool);
  }

  std::vector<int> bernoulli_sampling(const unsigned int K, const double PROB){
    std::vector<int> pool;
    for( int iterator = 0; iterator < K; iterator++){
      if(R::runif(0,1) < PROB ) pool.push_back(iterator);
    }
    return(pool);
  }

  void proj2(
        Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
        Eigen::Map<Eigen::VectorXd> CONSTRLOGSD,
        const std::vector<std::vector<std::vector<double>>> LLC,
        Eigen::Map<Eigen::VectorXd> C_VEC,
        const int CORRFLAG,
        const int NTHR,
        const int NLOAD,
        const int NCORR,
        const int NVAR,
        Eigen::VectorXd &THETA,
        bool &CHECKEVENT){
      //dimensions
      const int p = CONSTRMAT.rows(); // number of items
      const int q = CONSTRMAT.cols(); // number of latents

      // rearrange parameters
      Eigen::MatrixXd Lam             = params::loadings::theta2mat(THETA, CONSTRMAT, LLC, NTHR, NLOAD);
      Eigen::MatrixXd Ru              = params::latvar::theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, q);
      Eigen::MatrixXd Du              = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q);
      Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
      Eigen::VectorXd tau             = params::thresholds::theta2vec(THETA, NTHR);


      // Normalisation if needed
      for(unsigned int k = 0; k < p; k++){
        const Eigen::VectorXd lambdak = Lam.row(k);
        const double vk = lambdak.transpose()*Sigma_u*lambdak;
        if(vk>1){
          CHECKEVENT = true;

          for(int h=0; h<q; h++){

            if(!std::isfinite(CONSTRMAT(k,h))){
              // Free loading


              if(!std::isfinite(CONSTRLOGSD(h))){
                // Free loading and free variance
                Lam(k,h)/=sqrt(sqrt(vk));
                Du(h,h)/=sqrt(sqrt(vk));

              }else{

                // Free loading but fixed variance
                Lam(k,h)/= sqrt(vk);

              }
            }else{
              // Fixed loading

              if(CONSTRMAT(k,h)!=0.0){
                // Non-zero loading
                Du(h,h)/=sqrt(vk);
              }

            }
          }
        }
      }



      if(CHECKEVENT){

        int idx = 0;
        if(NVAR>0){
          Eigen::VectorXd norm_logsd(NVAR);
          for(unsigned int h = 0; h < q; h++){
            if(!std::isfinite(CONSTRLOGSD(h))){
              norm_logsd(idx) = log(Du(h,h));
              idx++;
            }
          }
          THETA.segment(NTHR+NLOAD+NCORR, NVAR) = norm_logsd;
        }

        idx=0;
        Eigen::VectorXd norm_lambda(NLOAD);
        for(unsigned int h = 0; h < q; h++){
          for(unsigned int j = 0; j < p; j++){
            if(!std::isfinite(CONSTRMAT(j,h)))
            {
              norm_lambda(idx) = Lam(j, h);
              idx ++;
            };
          };
        }
        THETA.segment(NTHR, NLOAD) = norm_lambda;

      }
    }




}

































#endif


