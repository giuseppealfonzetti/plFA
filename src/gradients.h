#ifndef gradients_H
#define gradients_H
#include "genericUtils.h"

namespace grads{
  void thresholds(
      Eigen::VectorXd &GRADIENT,
      const Eigen::Ref<const Eigen::MatrixXd> PAIRS_TAB,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const double RHO_KL,
      const unsigned int K,
      const unsigned int L,
      const unsigned int CK,
      const unsigned int CL,
      const unsigned int I1,
      const unsigned int I2,
      unsigned int &IDX
  ){
    int SILENTFLAG=1;
    for(unsigned int s = 0; s < TAU.size(); s++){
      double grs = 0; // temporary location for gradient related to s-th element of tau
      if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ gradient("<< s<< ")\n";

      // List three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
      if(s >= (C_VEC.segment(0, K).sum()) - (K) & s < C_VEC.segment(0, K + 1).sum() - (K + 1)){
        // [CASE 1]: threshold related to item k

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item k:\n";
        const unsigned int sk = s - (C_VEC.segment(0, K).sum()) + (K);

        // i3: starting index from i2 for cat sk and sk+1
        const unsigned int i3 = sk * CL;
        const unsigned int i3suc = (sk+1) * CL;
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ sk: " << sk << ". Summing over categories item l: ";

        // iterate over categories of item l
        for(unsigned int sl = 0; sl < CL; sl ++){
          if(SILENTFLAG == 0)Rcpp::Rcout << " ... cat" << sl ;

          // identify pairs_tab column for (sk,sl) and (sk+1, sl)
          const unsigned int r = I1 + I2 + i3 + sl;
          const unsigned int rsuc = I1 + I2 + i3suc + sl;

          // read frequences
          const unsigned int n_sksl = PAIRS_TAB(4, r);
          const unsigned int n_sksucsl = PAIRS_TAB(4, rsuc);

          // read probabilities
          const double pi_sksl = PAIRS_TAB(5, r);
          const double pi_sksucsl = PAIRS_TAB(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = extract_thresholds(TAU, C_VEC, K, L, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);

          // compute gradient
          const double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_sksucsl/(pi_sksucsl+1e-8)));
          const double tmp2 = R::dnorm(t_sk, 0, 1, 0);
          const double tmp3 = R::pnorm((t_sl-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          const double tmp4 = R::pnorm((t_sl_prev-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          grs += tmp1 * tmp2 * (tmp3 - tmp4);
        }
        if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

      }else if(s >= (C_VEC.segment(0, L).sum())-(L) & s<C_VEC.segment(0, L + 1).sum()-(L + 1)){
        // [CASE 2]: threshold related to item l

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item l\n";
        const unsigned int sl = s - (C_VEC.segment(0, L).sum()) + (L);

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  sl: " << sl << ". Summing over categories item k: ";

        // iterate over categories item k
        for(unsigned int sk = 0; sk < CK; sk ++){

          // i3: starting index from i2 for cat sk
          const unsigned int i3 = sk * CL;

          // identify pairs_tab column for (sk,sl) and (sk, sl + 1)
          const unsigned int r = I1 + I2 + i3 + sl;
          const unsigned int rsuc = I1 + I2 + i3 + sl + 1;

          // read frequences
          const unsigned int n_sksl = PAIRS_TAB(4, r);
          const unsigned int n_skslsuc = PAIRS_TAB(4, rsuc);

          // read probabilities
          const double pi_sksl = PAIRS_TAB(5, r);
          const double pi_skslsuc = PAIRS_TAB(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = extract_thresholds(TAU, C_VEC, K, L, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);


          if(SILENTFLAG == 0)Rcpp::Rcout<<"\n  |    |   |_ sk:"<< sk << ", r: "<< r<<", n_sksl:"
                                        << n_sksl<< ", n_sksl+1:" << n_skslsuc << ", pi_sksl:"
                                        << pi_sksl << ", pi_sksl+1:"<< pi_skslsuc << ", t_sk:"
                                        << t_sk<< ", t_sl:" << t_sl << "t_sk-1:"<< t_sk_prev;
          // compute gradient
          const double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_skslsuc/(pi_skslsuc+1e-8)));
          const double tmp2 = R::dnorm(t_sl, 0, 1, 0);
          const double tmp3 = R::pnorm((t_sk-RHO_KL*t_sl)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          const double tmp4 = R::pnorm((t_sk_prev-RHO_KL*t_sl)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          if(SILENTFLAG == 0)Rcpp::Rcout<<" => out" << sk << ":" << tmp1 * tmp2 * (tmp3 - tmp4);
          grs += tmp1 * tmp2 * (tmp3 - tmp4);
        }
        if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

      }else{
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  tau of other item\n";
      }

      GRADIENT(IDX) += grs;
      IDX ++;
    }

  }
}
#endif
