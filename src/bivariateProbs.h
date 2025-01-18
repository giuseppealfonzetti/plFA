#ifndef bivariateProbs_H
#define bivariateProbs_H

#include "bivariateNormal.h"
#include "gradients.h"

namespace biprobs{
  /*  PI (the probability of a bivariate response pattern) */
  // Compute specific pi_sksl
  double compute_pi(
      const Eigen::VectorXd c_vec,
      const Eigen::VectorXd pi_thresholds,
      double rho_kl,

      const unsigned int k,
      const unsigned int l,
      const unsigned int sk,
      const unsigned int sl
  ){
    unsigned int ck = c_vec(k);
    unsigned int cl = c_vec(l);

    // read pi related thresholds
    double t_sk = pi_thresholds(0);
    double t_sl = pi_thresholds(1);
    double t_sk_prev = pi_thresholds(2);
    double t_sl_prev = pi_thresholds(3);

    // Phi(t_sk, t_sl; rho_kl)
    double cum1;
    if(sk == (ck-1) & sl == (cl-1)){
      cum1 = 1;
    } else if(sk == (ck-1)){
      cum1 = R::pnorm(t_sl, 0, 1, 1, 0);
    } else if(sl == (cl-1)){
      cum1 = R::pnorm(t_sk, 0, 1, 1, 0);
    } else {
      cum1 = binorm::pbvnorm( t_sk, t_sl, rho_kl);
    }

    // Phi(t_sk, t_sl-1; rho_kl)
    double cum2;
    if(sl == 0){
      cum2 = 0;
    } else {
      cum2 = binorm::pbvnorm( t_sk, t_sl_prev, rho_kl);
    }
    // Phi(t_sk-1, t_sl; rho_kl)
    double cum3;
    if(sk == 0){
      cum3 = 0;
    } else{
      cum3 = binorm::pbvnorm( t_sk_prev, t_sl, rho_kl);
    }
    // Phi(t_sk-1, t_sl-1; rho_kl)
    double cum4;
    if(sl == 0 | sk == 0){
      cum4 = 0;
    }else{
      cum4 = binorm::pbvnorm( t_sk_prev, t_sl_prev, rho_kl);
    }

    double pi_sksl = cum1 - cum2 - cum3 + cum4;

    return pi_sksl;
  }

  double pi(
      const Eigen::VectorXd C_VEC,
      const Eigen::VectorXd PI_THRESHOLDS,
      double RHO_KL,
      const unsigned int K,
      const unsigned int L,
      const unsigned int SK,
      const unsigned int SL
  ){
    unsigned int ck = C_VEC(K);
    unsigned int cl = C_VEC(L);

    // read pi related thresholds
    double t_sk = PI_THRESHOLDS(0);
    double t_sl = PI_THRESHOLDS(1);
    double t_sk_prev = PI_THRESHOLDS(2);
    double t_sl_prev = PI_THRESHOLDS(3);

    // Phi(t_sk, t_sl; rho_kl)
    double cum1;
    if(SK == (ck-1) & SL == (cl-1)){
      cum1 = 1;
    } else if(SK == (ck-1)){
      cum1 = R::pnorm(t_sl, 0, 1, 1, 0);
    } else if(SL == (cl-1)){
      cum1 = R::pnorm(t_sk, 0, 1, 1, 0);
    } else {
      cum1 = binorm::pbvnorm( t_sk, t_sl, RHO_KL);
    }

    // Phi(t_sk, t_sl-1; rho_kl)
    double cum2;
    if(SL == 0){
      cum2 = 0;
    } else {
      cum2 = binorm::pbvnorm( t_sk, t_sl_prev, RHO_KL);
    }
    // Phi(t_sk-1, t_sl; rho_kl)
    double cum3;
    if(SK == 0){
      cum3 = 0;
    } else{
      cum3 = binorm::pbvnorm( t_sk_prev, t_sl, RHO_KL);
    }
    // Phi(t_sk-1, t_sl-1; rho_kl)
    double cum4;
    if(SL == 0 | SK == 0){
      cum4 = 0;
    }else{
      cum4 = binorm::pbvnorm( t_sk_prev, t_sl_prev, RHO_KL);
    }

    double pi_sksl = cum1 - cum2 - cum3 + cum4;

    return pi_sksl;
  }
  /* GRADIENT OF PI */
  // Compute gradient of specific pi_sksl
  Eigen::VectorXd compute_pi_grad(Eigen::Map<Eigen::MatrixXd> A,
                                  const Eigen::VectorXd c_vec,
                                  const Eigen::VectorXd &pi_thresholds,
                                  const Eigen::MatrixXd &Sigma_u,
                                  const Eigen::MatrixXd &Lam,
                                  const Eigen::VectorXd &theta,
                                  const double &rho_kl,
                                  const unsigned int k,
                                  const unsigned int l,
                                  const unsigned int sk,
                                  const unsigned int sl,
                                  const unsigned int corrFLAG){
    unsigned int d = theta.size();
    unsigned int p = A.rows();
    unsigned int q = A.cols();
    unsigned int ck = c_vec(k);
    unsigned int cl = c_vec(l);
    unsigned int ncorr = q*(q-1)/2;
    Eigen::VectorXd transformed_rhos = theta.segment(d-ncorr, ncorr);
    Eigen::VectorXd lambdak = Lam.row(k);
    Eigen::VectorXd lambdal = Lam.row(l);


    // read pi related thresholds
    double t_sk = pi_thresholds(0);
    double t_sl = pi_thresholds(1);
    double t_sk_prev = pi_thresholds(2);
    double t_sl_prev = pi_thresholds(3);

    //Intermediate derivative pi wrt to kl correlation
    // phi(t_sk, t_sl; rho_kl)
    double d1 = binorm::dbvnorm( t_sk, t_sl, rho_kl, 0);

    // phi(t_sk, t_sl-1; rho_kl)
    double d2 = binorm::dbvnorm( t_sk, t_sl_prev, rho_kl, 0);

    // phi(t_sk-1, t_sl; rho_kl)
    double d3 = binorm::dbvnorm( t_sk_prev, t_sl, rho_kl, 0);

    // phi(t_sk-1, t_sl-1; rho_kl)
    double d4 = binorm::dbvnorm( t_sk_prev, t_sl_prev, rho_kl, 0);

    double tmp_sksl = d1 - d2 - d3 + d4;

    //Rcpp::Rcout << "  |    |_ d1:"<<d1<<", d2:"<<d2<<", d3:"<<d3<< ", d4:"<<d4<<", tmp_sksl: "<< tmp_sksl<< "\n";
    // for each (sk, sl) compute the gradient of pi wrt the parameter vector
    Eigen::VectorXd dpi(d); dpi.fill(0.0);

    // iterator over parameter vector
    unsigned int iter_th = 0;

    ////////////////////////////////////////////////
    // gradient pi_sksl wrt to thresholds
    ///////////////////////////////////////////////

    // loop: iterate over elements of thresholds vector
    for(unsigned int s = 0; s < c_vec.sum()-p; s++){

      // Elicit three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
      if(s >= (c_vec.segment(0, k).sum()) - (k) & s < c_vec.segment(0, k + 1).sum() - (k + 1)){
        // [CASE 1]: threshold related to item k
        unsigned int sk_a = s - (c_vec.segment(0, k).sum()) + (k);
        //Rcpp::Rcout << "  |    |_ tau item k. sk_a:"<< sk_a <<"\n";
        if(sk_a == sk){
          double tmp1 = R::dnorm(t_sk, 0, 1, 0);
          double tmp2 = R::pnorm((t_sl-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sl_prev-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);;
          dpi(iter_th) = tmp1*tmp2-tmp1*tmp3;
        }else if(sk_a == (sk-1)){
          double tmp1 = R::dnorm(t_sk_prev, 0, 1, 0);
          double tmp2 = R::pnorm((t_sl-rho_kl*t_sk_prev)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sl_prev-rho_kl*t_sk_prev)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);;
          dpi(iter_th) = -tmp1*tmp2+tmp1*tmp3;
        }else{
          dpi(iter_th) = 0;
        }

      }else if(s >= (c_vec.segment(0, l).sum())-(l) & s<c_vec.segment(0, l + 1).sum()-(l + 1)){
        // [CASE 2]: threshold related to item l
        unsigned int sl_a = s - (c_vec.segment(0, l).sum()) + (l);
        //Rcpp::Rcout << "  |    |_ tau item l. sl_a:"<< sl_a<<"\n";
        if(sl_a == sl){
          double tmp1 = R::dnorm(t_sl, 0, 1, 0);
          double tmp2 = R::pnorm((t_sk-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sk_prev-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);;
          dpi(iter_th) = tmp1*tmp2-tmp1*tmp3;
        }else if(sl_a == (sl-1)){
          double tmp1 = R::dnorm(t_sl_prev, 0, 1, 0);
          double tmp2 = R::pnorm((t_sk-rho_kl*t_sl_prev)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sk_prev-rho_kl*t_sl_prev)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);;
          dpi(iter_th) = -tmp1*tmp2+tmp1*tmp3;
        }else{
          dpi(iter_th) = 0;
        }
      } else {
        // [CASE 3]: threshold non related to (k,l)
        //Rcpp::Rcout << "  |    |_ tau non related to (k,l)\n";
        dpi(iter_th) = 0;
      }
      iter_th ++;
    }

    //////////////////////////////////////////////
    // gradient pi_sksl wrt to loadings
    //////////////////////////////////////

    // double loop: iterate over elements of loadings matrix
    for(unsigned int j = 0; j < p; j++){
      for(unsigned int v = 0; v < q; v++){
        //if(silentFLAG == 0)Rcpp::Rcout << "  |_ visiting lambda_"<< j << v <<":\n";

        // elicit three cases: 1. free loading item k, 2. free loading l, 3. other
        if(j == k){
          if(A(j,v)!=0 ){
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ item k, free loading:\n";
            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
            double d_rho_kl = ev.transpose() * Sigma_u * lambdal;
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            dpi(iter_th) = tmp_sksl * d_rho_kl;
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter_th - tau.size() << "/"<< lambda.size() <<". Tot:: "<< iter_th << "/"<< d << "\n";

            iter_th ++;
          }
        }else if (j == l){
          if(A(j,v)!=0 ){
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ item l, free loading:\n";
            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
            double d_rho_kl = lambdak.transpose() * Sigma_u * ev;
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            //if(silentFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter_th - tau.size() << "/"<< lambda.size() <<". Tot:: "<< iter_th << "/"<< d << "\n";
            dpi(iter_th) = tmp_sksl * d_rho_kl;

            iter_th ++;
          }
        }else if(A(j,v)!=0){

          //if(silentFLAG == 0)Rcpp::Rcout << "  |  |_ Loadings:: " << iter_th - tau.size() << "/"<< lambda.size() <<". Tot:: "<< iter_th << "/"<< d << " [not included]\n";
          iter_th ++;
        }
      }
    }

    /////////////////////////////////////////
    // gradient pi_sksl wrt correlations
    ///////////////////////////////////////
    if(corrFLAG == 1){
      for(unsigned int thro_idx = 0; thro_idx < ncorr; thro_idx ++){

        Eigen::MatrixXd dSigma = grads::S(A, transformed_rhos, q, thro_idx);
        // Rcpp::Rcout << "dsigma wrt corr idx " << thro_idx << ":\n"<< dSigma << "\n";

        double d_rho_kl = lambdak.transpose() * dSigma * lambdal;
        dpi(iter_th) += tmp_sksl * d_rho_kl;
        iter_th++;
      }
    }

    return dpi;
  }

}




#endif
