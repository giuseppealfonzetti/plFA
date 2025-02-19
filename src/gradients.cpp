#include "gradients.h"

namespace grads{
  void thresholds(
      Eigen::VectorXd &GRADIENT,
      const Eigen::Ref<const Eigen::MatrixXd> PAIRS_TAB,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const double RHO_KL,
      const int K,
      const int L,
      const int CK,
      const int CL,
      const int I1,
      const int I2,
      int &IDX
  ){
    int SILENTFLAG=1;
    for(int s = 0; s < TAU.size(); s++){
      double grs = 0; // temporary location for gradient related to s-th element of tau
      if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ gradient("<< s<< ")\n";

      // List three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
      if(s >= (C_VEC.segment(0, K).sum()) - (K) && s < C_VEC.segment(0, K + 1).sum() - (K + 1)){
        // [CASE 1]: threshold related to item k

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item k:\n";
        const int sk = s - (C_VEC.segment(0, K).sum()) + (K);

        // i3: starting index from i2 for cat sk and sk+1
        const int i3 = sk * CL;
        const int i3suc = (sk+1) * CL;
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ sk: " << sk << ". Summing over categories item l: ";

        // iterate over categories of item l
        for(int sl = 0; sl < CL; sl ++){
          if(SILENTFLAG == 0)Rcpp::Rcout << " ... cat" << sl ;

          // identify pairs_tab column for (sk,sl) and (sk+1, sl)
          const int r = I1 + I2 + i3 + sl;
          const int rsuc = I1 + I2 + i3suc + sl;

          // read frequences
          const int n_sksl = PAIRS_TAB(4, r);
          const int n_sksucsl = PAIRS_TAB(4, rsuc);

          // read probabilities
          const double pi_sksl = PAIRS_TAB(5, r);
          const double pi_sksucsl = PAIRS_TAB(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = params::extract_thresholds(TAU, C_VEC, K, L, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1);  const double t_sl_prev = pi_thresholds(3);

          // compute gradient
          const double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_sksucsl/(pi_sksucsl+1e-8)));
          const double tmp2 = R::dnorm(t_sk, 0, 1, 0);
          const double tmp3 = R::pnorm((t_sl-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          const double tmp4 = R::pnorm((t_sl_prev-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          grs += tmp1 * tmp2 * (tmp3 - tmp4);
        }
        if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

      }else if(s >= (C_VEC.segment(0, L).sum())-(L) && s<C_VEC.segment(0, L + 1).sum()-(L + 1)){
        // [CASE 2]: threshold related to item l

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item l\n";
        const int sl = s - (C_VEC.segment(0, L).sum()) + (L);

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  sl: " << sl << ". Summing over categories item k: ";

        // iterate over categories item k
        for(int sk = 0; sk < CK; sk ++){

          // i3: starting index from i2 for cat sk
          const int i3 = sk * CL;

          // identify pairs_tab column for (sk,sl) and (sk, sl + 1)
          const int r = I1 + I2 + i3 + sl;
          const int rsuc = I1 + I2 + i3 + sl + 1;

          // read frequences
          const int n_sksl = PAIRS_TAB(4, r);
          const int n_skslsuc = PAIRS_TAB(4, rsuc);

          // read probabilities
          const double pi_sksl = PAIRS_TAB(5, r);
          const double pi_skslsuc = PAIRS_TAB(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = params::extract_thresholds(TAU, C_VEC, K, L, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2);


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

  double rho_urv(
      const Eigen::Ref<const Eigen::MatrixXd> PAIRS_TAB,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> TAU,
      const double RHO_KL,
      const int K,
      const int L,
      const int CK,
      const int CL,
      const int I1,
      const int I2
  ){
    double out = 0;
    int SILENTFLAG=1;

    for(int sk = 0; sk < CK; sk ++){
      for(int sl = 0; sl < CL; sl ++){
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ sk: "<< sk << ", sl: " << sl << ": \n";

        // identify pairs_tab column for (sk,sl)
        const int i3 = sk * CL;
        const int r = I1 + I2 + i3 + sl;

        // read freq
        const int n_sksl = PAIRS_TAB(4, r);

        // read prob
        const double pi_sksl = PAIRS_TAB(5, r);

        // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
        const Eigen::VectorXd pi_thresholds = params::extract_thresholds(TAU, C_VEC, K, L, sk, sl);
        const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);

        // phi(t_sk, t_sl; rho_kl)
        const double d1 = binorm::dbvnorm( t_sk, t_sl, RHO_KL, 0);

        // phi(t_sk, t_sl-1; rho_kl)
        const double d2 = binorm::dbvnorm( t_sk, t_sl_prev, RHO_KL, 0);

        // phi(t_sk-1, t_sl; rho_kl)
        const double d3 = binorm::dbvnorm( t_sk_prev, t_sl, RHO_KL, 0);

        // phi(t_sk-1, t_sl-1; rho_kl)
        const double d4 = binorm::dbvnorm( t_sk_prev, t_sl_prev, RHO_KL, 0);

        out += (n_sksl/(pi_sksl+1e-8)) * ( d1 - d2 - d3 + d4);
      }
    }

    return out;
  }

  void loadings(
      Eigen::VectorXd &GRADIENT,
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::MatrixXd> SIGMA_U,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAK,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAL,
      const double PRHO_URV,
      const int P,
      const int Q,
      const int K,
      const int L,
      int &IDX
  ){
    int SILENTFLAG=1;

    Eigen::MatrixXd dict=Eigen::MatrixXd::Zero(P, Q);
    for(int j = 0; j < P; j++){
      for(int v = 0; v < Q; v++){
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ visiting lambda_"<< j << v <<":\n";

        // elicit three cases: 1. free loading item k, 2. free loading l, 3. other
        if(j == K){
          if(!std::isfinite(A(j,v))){
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item k, free loading:\n";
            Eigen::VectorXd ev(Q); ev.fill(0.0); ev(v) = 1;
            const double d_rho_kl = ev.transpose() * SIGMA_U * LAMBDAL;
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            GRADIENT(IDX) += PRHO_URV * d_rho_kl;
            dict(j,v) = IDX;

            IDX ++;
          }
        }else if (j == L){
          if(!std::isfinite(A(j,v))){
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item l, free loading:\n";
            Eigen::VectorXd ev(Q); ev.fill(0.0); ev(v) = 1;
            const double d_rho_kl = LAMBDAK.transpose() * SIGMA_U * ev;
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            GRADIENT(IDX) += PRHO_URV * d_rho_kl;
            dict(j,v) = IDX;

            IDX ++;
          }
        }else if(!std::isfinite(A(j,v))){
          dict(j,v) = IDX;
          IDX ++;
        }
      }
    }

    for(int idx_lc = 0; idx_lc < LLC.size(); idx_lc++){
      const int target_row = LLC.at(idx_lc).at(0).at(0)-1;
      const int target_col = LLC.at(idx_lc).at(0).at(1)-1;
      const int n_loads_involved = LLC.at(idx_lc).size();
      // Rcpp::Rcout << "Constrain number " << idx_lc+1 << ", L("<<target_row+1<<","<<target_col+1
      //             <<"), Number of params involved: "<< n_loads_involved-1<<"\n";

      if(target_row==K){

        // Rcpp::Rcout<< "Target row == K \n";
        for(int idx_var = 1; idx_var<n_loads_involved; idx_var++){
          Eigen::VectorXd ev=Eigen::VectorXd::Zero(Q);
          const double coeff = LLC.at(idx_lc).at(idx_var).at(0);
          const int target_row_var = LLC.at(idx_lc).at(idx_var).at(1)-1;
          const int target_col_var = LLC.at(idx_lc).at(idx_var).at(2)-1;
          const int target_idx_var = dict(target_row_var, target_col_var);

          // Rcpp::Rcout << " + " << coeff << "*" "L("<<target_row_var+1<<","<<target_col_var+1<<"), idx:" << target_idx_var<< "\n";

          ev(target_col) = coeff;
          const double d_rho_kl = ev.transpose() * SIGMA_U * LAMBDAL;

          GRADIENT(target_idx_var) +=  PRHO_URV * d_rho_kl;
          }
      } else if(target_row==L){

        // Rcpp::Rcout<< "Target row ==L \n";

        for(int idx_var = 1; idx_var<n_loads_involved; idx_var++){
          Eigen::VectorXd ev=Eigen::VectorXd::Zero(Q);
          const double coeff = LLC.at(idx_lc).at(idx_var).at(0);
          const int target_row_var = LLC.at(idx_lc).at(idx_var).at(1)-1;
          const int target_col_var = LLC.at(idx_lc).at(idx_var).at(2)-1;
          const int target_idx_var = dict(target_row_var, target_col_var);
          ev(target_col) = coeff;
          const double d_rho_kl = LAMBDAK.transpose() * SIGMA_U * ev;

          GRADIENT(target_idx_var) +=  PRHO_URV * d_rho_kl;
        }
      }



    }
  }

  Eigen::MatrixXd S(const Eigen::Ref<const Eigen::MatrixXd> A,
                    const Eigen::Ref<const Eigen::VectorXd> TRANSFORMED_RHOS,
                    const int Q,
                    const int IDX
  ){

    // inverse of Fisher's transformation
    Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*TRANSFORMED_RHOS)).array().exp() - 1)/
      ((Eigen::VectorXd(2*TRANSFORMED_RHOS)).array().exp() + 1);

    // place tanh_entries in upper triangular matrix
    // with unit diagonal
    Eigen::MatrixXd Z(Q,Q); Z.setIdentity();
    int iterator = 0; int i,j;
    for(int col = 1; col < Q; col ++){
      for(int row = 0; row < col; row++){
        if(iterator==IDX){i = row; j = col;}
        Z(row, col) = tanh_entries(iterator);
        iterator ++;
      }
    }

    // construct upper Cholesky from Z
    Eigen::MatrixXd U(Q,Q);
    for( int col = 0; col < Q; col++){
      for( int row = 0; row < Q; row++){
        if(row > col) U(row, col) = 0;
        else if ((row == col) && (row == 0)) U(row, col) = 1;
        else if ((row > 0) && (row == col)) {
          if(Z(row-1,col)==0){
            double prod = 1;
            for(int row1 = 0; row1 < row; row1++){
              prod *= pow(1-pow(Z(row1, col), 2),.5);
            }
            U(row, col) = prod;
          }else{
            U(row, col) = U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5)/Z(row-1, col);
          }
        }
        else if ((row == 0) && (row < col)) U(row, col) = Z(row, col);
        else if ((row > 0) && (row < col)) {
          if(Z(row-1,col)==0){
            double prod = Z(row, col);
            for( int row1 = 0; row1 < row; row1++){
              prod *= pow(1-pow(Z(row1, col), 2),.5);
            }
            U(row, col) = prod;
          }else{
            U(row, col) = Z(row, col)/Z(row-1, col)*U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5);
          }
        }
      }
    }

    Eigen::MatrixXd dU(Q,Q); dU.setZero();

    for( int row = 0; row < Q; row++){
      if(row > j) dU(row, j) = 0;
      else if ((row == j) && (row == 0)) dU(row, j) = 0;
      else if ((row == 0) && (row < j) && (row == i)) dU(row, j) = 1;
      else if ((row > 0) && (row < j) && (row == i)) {
        if(Z(row-1, j)==0){
          double prod = 1;
          for( int row1 = 0; row1 < row; row1++){
            prod *= pow(1-pow(Z(row1, j), 2),.5);
          }
          dU(row, j) = prod;
        }else{
          dU(row, j) = U(row-1,j)*pow(1-pow(Z(row-1, j),2), .5)/Z(row-1, j);
        }}
      else if ((row > 0) && (row < j) && (i < row)) {
        if(Z(row-1, j)==0){
          double prod = Z(row, j);
          for( int row1 = 0; row1 < row; row1++){
            prod *= pow(1-pow(Z(row1, j), 2),.5);
          }
          dU(row, j) = - prod * Z(i, j) / (1-pow(Z(i,j), 2));
        }else{
          if(row-i==1){
            dU(row, j) = -Z(row, j)*pow(Z(row-1, j),-2)*U(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) + Z(row, j)/Z(row-1, j)*(dU(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) - U(row-1, j)* pow(1-pow(Z(row-1, j),2), -.5)*Z(row-1, j));
          }else{
            dU(row, j) = Z(row,j)*pow(Z(row-1,j),-1)*dU(row-1,j)*pow(1-pow(Z(row-1, j),2), .5);
          }
        }

      }
      else if ((row > 0) && (row == j) && (row == i)) dU(row, j) = 0;
      else if ((row > 0) && (row == j) && (i < row)){
        dU(row, j) = - U(row, j) * Z(i, j)/ (1-pow(Z(i,j), 2));
      }
    }

    double dz = pow(cosh(TRANSFORMED_RHOS(IDX)),-2);
    dU = dU*dz;

    Eigen::MatrixXd dS = dU.transpose()*U + U.transpose()*dU;

    return dS;
  }

  void lat_corr(
      Eigen::VectorXd &GRADIENT,
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAK,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAL,
      const Eigen::Ref<const Eigen::VectorXd> TRANSFORMED_RHOS,
      const double PRHO_URV,
      const int Q,
      const int NCORR,
      int &IDX
  ){
    for(int thro_idx = 0; thro_idx < NCORR; thro_idx ++){
      const Eigen::MatrixXd dSigma = grads::S(A, TRANSFORMED_RHOS, Q, thro_idx);
      const double d_rho_kl = LAMBDAK.transpose() * dSigma * LAMBDAL;
      GRADIENT(IDX) += PRHO_URV * d_rho_kl;
      IDX++;
    }
  }

  Eigen::VectorXd pi(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> A,
      const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
      const std::vector<std::vector<std::vector<double>>> LLC,
      const Eigen::Ref<const Eigen::VectorXd> C_VEC,
      const Eigen::Ref<const Eigen::VectorXd> PI_THRESHOLDS,
      const Eigen::Ref<const Eigen::MatrixXd> SIGMA_U,
      const Eigen::Ref<const Eigen::MatrixXd> D_U,
      const Eigen::Ref<const Eigen::MatrixXd> R_U,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAK,
      const Eigen::Ref<const Eigen::VectorXd> LAMBDAL,
      const double RHO_KL,
      const int D,
      const int P,
      const int Q,
      const int K,
      const int L,
      const int SK,
      const int SL,
      const int CORRFLAG,
      const int NTHR,
      const int NLOAD,
      const int NCORR,
      const int NVAR
  ){
    // read pi related thresholds
    const double t_sk = PI_THRESHOLDS(0);
    const double t_sl = PI_THRESHOLDS(1);
    const double t_sk_prev = PI_THRESHOLDS(2);
    const double t_sl_prev = PI_THRESHOLDS(3);

    //Intermediate derivative pi wrt to kl correlation
    // phi(t_sk, t_sl; rho_kl)
    const double d1 = binorm::dbvnorm( t_sk, t_sl, RHO_KL, 0);

    // phi(t_sk, t_sl-1; rho_kl)
    const double d2 = binorm::dbvnorm( t_sk, t_sl_prev, RHO_KL, 0);

    // phi(t_sk-1, t_sl; rho_kl)
    const double d3 = binorm::dbvnorm( t_sk_prev, t_sl, RHO_KL, 0);

    // phi(t_sk-1, t_sl-1; rho_kl)
    const double d4 = binorm::dbvnorm( t_sk_prev, t_sl_prev, RHO_KL, 0);

    const double tmp_sksl = d1 - d2 - d3 + d4;

    Eigen::VectorXd grad=Eigen::VectorXd::Zero(D);


    // iterator over parameter vector
    int idx = 0;

    ////////////////////////////////////////////////
    // gradient pi_sksl wrt to thresholds
    ///////////////////////////////////////////////

    // loop: iterate over elements of thresholds vector
    for(int s = 0; s < C_VEC.sum()-P; s++){

      // Elicit three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
      if(s >= (C_VEC.segment(0, K).sum()) - (K) && s < C_VEC.segment(0, K + 1).sum() - (K + 1)){
        // [CASE 1]: threshold related to item k
        int sk_a = s - (C_VEC.segment(0, K).sum()) + (K);
        //Rcpp::Rcout << "  |    |_ tau item k. sk_a:"<< sk_a <<"\n";
        if(sk_a == SK){
          double tmp1 = R::dnorm(t_sk, 0, 1, 0);
          double tmp2 = R::pnorm((t_sl-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sl_prev-RHO_KL*t_sk)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);;
          grad(idx) = tmp1*tmp2-tmp1*tmp3;
        }else if(sk_a == (SK-1)){
          double tmp1 = R::dnorm(t_sk_prev, 0, 1, 0);
          double tmp2 = R::pnorm((t_sl-RHO_KL*t_sk_prev)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sl_prev-RHO_KL*t_sk_prev)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);;
          grad(idx) = -tmp1*tmp2+tmp1*tmp3;
        }else{
          grad(idx) = 0;
        }

      }else if(s >= (C_VEC.segment(0, L).sum())-(L) && s<C_VEC.segment(0, L + 1).sum()-(L + 1)){
        // [CASE 2]: threshold related to item l
        int sl_a = s - (C_VEC.segment(0, L).sum()) + (L);
        //Rcpp::Rcout << "  |    |_ tau item l. sl_a:"<< sl_a<<"\n";
        if(sl_a == SL){
          double tmp1 = R::dnorm(t_sl, 0, 1, 0);
          double tmp2 = R::pnorm((t_sk-RHO_KL*t_sl)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sk_prev-RHO_KL*t_sl)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);;
          grad(idx) = tmp1*tmp2-tmp1*tmp3;
        }else if(sl_a == (SL-1)){
          double tmp1 = R::dnorm(t_sl_prev, 0, 1, 0);
          double tmp2 = R::pnorm((t_sk-RHO_KL*t_sl_prev)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);
          double tmp3 = R::pnorm((t_sk_prev-RHO_KL*t_sl_prev)/(pow(1-pow(RHO_KL,2), .5)), 0, 1, 1, 0);;
          grad(idx) = -tmp1*tmp2+tmp1*tmp3;
        }else{
          grad(idx) = 0;
        }
      } else {
        // [CASE 3]: threshold non related to (k,l)
        //Rcpp::Rcout << "  |    |_ tau non related to (k,l)\n";
        grad(idx) = 0;
      }
      idx ++;
    }

    //////////////////////////////////////////////
    // gradient pi_sksl wrt to loadings
    //////////////////////////////////////

    // double loop: iterate over elements of loadings matrix
    // for(int j = 0; j < P; j++){
    //   for(int v = 0; v < Q; v++){
    //
    //     // elicit three cases: 1. free loading item k, 2. free loading l, 3. other
    //     if(j == K){
    //       if(!std::isfinite(A(j,v))){
    //         Eigen::VectorXd ev(Q); ev.fill(0.0); ev(v) = 1;
    //         double d_rho_kl = ev.transpose() * SIGMA_U * LAMBDAL;
    //         grad(idx) = tmp_sksl * d_rho_kl;
    //         idx ++;
    //       }
    //     }else if (j == L){
    //       if(!std::isfinite(A(j,v))){
    //         Eigen::VectorXd ev(Q); ev.fill(0.0); ev(v) = 1;
    //         double d_rho_kl = LAMBDAK.transpose() * SIGMA_U * ev;
    //         grad(idx) = tmp_sksl * d_rho_kl;
    //         idx ++;
    //       }
    //     }else if(!std::isfinite(A(j,v))){
    //       idx ++;
    //     }
    //   }
    // }

    grads::loadings(grad, A, LLC, SIGMA_U, LAMBDAK, LAMBDAL, tmp_sksl, P, Q, K, L, idx);


    /////////////////////////////////////////
    // gradient pi_sksl wrt correlations
    ///////////////////////////////////////
    if(CORRFLAG == 1){
      Eigen::VectorXd transformed_rhos = params::latvar::theta2vec(THETA, NTHR, NLOAD, NCORR, NVAR).segment(0, NCORR);

      for(int thro_idx = 0; thro_idx < NCORR; thro_idx ++){

        Eigen::MatrixXd dR = grads::S(A, transformed_rhos, Q, thro_idx);
        // Rcpp::Rcout << "dsigma wrt corr idx " << thro_idx << ":\n"<< dSigma << "\n";

        double d_rho_kl = LAMBDAK.transpose()*D_U * dR * D_U*LAMBDAL;
        grad(idx) += tmp_sksl * d_rho_kl;
        idx++;
      }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (k,l)-pair likelihood derivative wrt reparameterised latent variances //
    ///////////////////////////////////////////////////////////////////////////
    if(NVAR > 0){
      for(int j=0; j<Q; j++){
        if(!std::isfinite(CONSTRLOGSD(j))){
          // Eigen::VectorXd ej = Eigen::VectorXd::Zero(q); ej(j) = 1;
          Eigen::MatrixXd dD = Eigen::MatrixXd::Zero(Q,Q); dD(j,j)=D_U(j,j);
          grad(idx) = tmp_sksl*LAMBDAK.transpose()*(dD*R_U*D_U+D_U*R_U*dD)*LAMBDAL;
          // Rcpp::Rcout<<"idx:"<<iter<<", gr:"<< GRADIENT(iter)<<"\n";
          idx++;
        }
      }

    }




















    return grad;


  }

}
