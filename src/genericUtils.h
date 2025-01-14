#ifndef utils_H
#define utils_H

namespace params{
  //' Get loading matrix from theta
  //'
  //' get_Lam() constructs the loading matrix from the parameter vector
  //'
  //' @param THETA Numerical vector of parameters.
  //' @param A Binary matrix of dimension \eqn{p*q} where \eqn{p} is the number
  //' of items and \eqn{q} the number of latent variables. Entries equal to
  //' \eqn{1} refer to free loadings, while entries equal to \eqn{0} indicate
  //' loadings constrained to be null.
  //' @param C Sum of the number of categories for each item.
  Eigen::MatrixXd get_Lam(const Eigen::Ref<const Eigen::MatrixXd> A,
                          const unsigned int C,
                          const Eigen::Ref<const Eigen::VectorXd> THETA
  ){
    const unsigned int p = A.rows();      // number of items
    const unsigned int q = A.cols();      // number of latent variables
    const unsigned int d = THETA.size();  // number of parameters

    const Eigen::VectorXd lambda = THETA.segment(C-p, d-C+p-q*(q-1)/2);
    Eigen::MatrixXd Lam = A;
    unsigned int iter = 0;
    for(unsigned int h = 0; h < q; h++){
      for(unsigned int j = 0; j < p; j++){
        if(A(j, h) != 0.0)
        {
          Lam(j, h) = lambda(iter) ;
          iter ++;
        };
      };
    }

    return Lam;
  }

  //' Get latent correlation matrix from theta
  //'
  //' get_S() extracts the latent correlation matrix from theta assuming
  //' theta elements to be reparametrised following the
  //' Lewandowski-Kurowicka-Joe (2009) transform.
  //'
  //' @param THETA Numerical vector of parameters.
  //' @param Q Number of latent variables.
  Eigen::MatrixXd get_S(const Eigen::Ref<const Eigen::VectorXd> THETA,
                        const unsigned int Q
  ){
    unsigned int d = THETA.size(); // number of parameters

    // extract uncostrained parameters related to latent correlations
    Eigen::VectorXd transformed_rhos = THETA.segment(d-Q*(Q-1)/2, Q*(Q-1)/2);

    // inverse of Fisher's transformation
    Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*transformed_rhos)).array().exp() - 1)/
      ((Eigen::VectorXd(2*transformed_rhos)).array().exp() + 1);

    // place tanh_entries in upper triangular matrix
    // with unit diagonal
    Eigen::MatrixXd Z(Q,Q); Z.setIdentity();
    unsigned int iterator = 0;
    for(unsigned int col = 1; col < Q; col ++){
      for(unsigned int row = 0; row < col; row++){
        Z(row, col) = tanh_entries(iterator);
        iterator ++;
      }
    }

    // construct upper Cholesky from Z
    Eigen::MatrixXd U(Q,Q);

    for( int col = 0; col < Q; col++){
      for( int row = 0; row < Q; row++){
        if(row > col) U(row, col) = 0;
        else if(row == col & row == 0) U(row, col) = 1;
        else if(row > 0 & row == col) {
          if(Z(row-1,col)==0){
            double prod = 1;
            for(unsigned int row1 = 0; row1 < row; row1++){
              prod *= pow(1-pow(Z(row1, col), 2),.5);
            }
            U(row, col) = prod;
          }else{
            U(row, col) = U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5)/Z(row-1, col);
          }
        }
        else if(row == 0 & row < col) U(row, col) = Z(row, col);
        else if(row > 0 & row < col) {
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


    // latent variable covariance matrix
    Eigen::MatrixXd S = U.transpose()*U;

    return S;
  }

  Eigen::VectorXd get_tau(const Eigen::Ref<const Eigen::VectorXd> THETA,
                          const unsigned int C,
                          const unsigned int P){
    return THETA.segment(0,C-P);
  }

  Eigen::VectorXd get_thros(const Eigen::Ref<const Eigen::VectorXd> THETA,
                            const unsigned int NTHR,
                            const unsigned int NLOAD,
                            const unsigned int NCORR){
    return THETA.segment(NTHR+NLOAD, NCORR);
  }



   //' Get transformed parameters from latent correlation matrix
   //'
   //' Use Lewandowski-Kurowicka-Joe (2009) transformation on latent correlations
   //' with Choleski decomposition.
   //'
   //' @param S Latent correlation matrix.
   Eigen::VectorXd get_par_from_S(const Eigen::Ref<const Eigen::MatrixXd> S){
     const unsigned int q = S.cols(); // number of latent variables

     Eigen::LLT<Eigen::MatrixXd> llt(S);
     Eigen::MatrixXd U = llt.matrixU();

     Eigen::VectorXd y(q*(q-1)/2);
     Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(q,q);
     Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(q,q);
     unsigned int iterator = 0;
     for(unsigned int col = 1; col < q; col ++){
       for(unsigned int row = 0; row < col; row++){
         if(row == 0) Z(row, col) = U(row, col);
         else if(row > 0){
           double prod = 1;
           for(unsigned int i = 0; i < row; i++){
             prod *= pow(1-pow(Z(i, col), 2),-.5);
             Z(row, col) = U(row, col)*prod;
           }
         }
         Y(row, col) = atanh(Z(row, col));
         y(iterator) = atanh(Z(row, col));
         iterator++;
       }
     }

     return y;
   }


/////////////////////////////////////////////////////////////////////////////////////////

  /* EXTRACT PI THRESHOLDS */
  // Extract from tau the trhesholds related to pi_sksl
  Eigen::VectorXd extract_thresholds(const Eigen::Ref<const Eigen::VectorXd> TAU,
                                     const Eigen::Ref<const Eigen::VectorXd> C_VEC,
                                     const unsigned int K,
                                     const unsigned int L,
                                     const unsigned int SK,
                                     const unsigned int SL){

    unsigned int ck = C_VEC(K);
    unsigned int cl = C_VEC(L);

    // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
    unsigned int sk_tau_index = C_VEC.segment(0, K).sum() - (K) + SK;   // index tau_sk in tau vector
    unsigned int sl_tau_index = C_VEC.segment(0, L).sum() - (L) + SL;   // index tau_sl in tau vector
    double t_sk;                                               // tau_sk
    if(SK == ck-1){
      t_sk = 100;
    } else {
      t_sk = TAU(sk_tau_index);
    }
    double t_sk_prev;                                           // tau_sk-1
    if(SK == 0){
      t_sk_prev = -100;
    }else{
      t_sk_prev = TAU(sk_tau_index-1);
    }

    double t_sl;                                               // tau_sl
    if(SL == cl-1){
      t_sl = 100;
    } else {
      t_sl = TAU(sl_tau_index);
    }
    double t_sl_prev;                                           // tau_sl-1
    if(SL == 0){
      t_sl_prev = -100;
    }else{
      t_sl_prev = TAU(sl_tau_index-1);
    }

    Eigen::VectorXd pi_thresholds(4);
    pi_thresholds << t_sk, t_sl, t_sk_prev, t_sl_prev;
    return pi_thresholds;
  }

  Eigen::VectorXd get_loadings_mat2vec(const Eigen::Ref<const Eigen::MatrixXd> LOADINGS,
                                       const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
                                       const int  NLOAD){

    Eigen::VectorXd out(NLOAD);
    const int p = LOADINGS.rows();
    const int q = LOADINGS.cols();
    int idx = 0;

    // construct vector with parameters free to be estimated
    for(int h = 0; h < q; h++){
      for(int j = 0; j < p; j++){
        if(!std::isfinite(CONSTRMAT(j,h))){
          out(idx)=LOADINGS(j,h);
          idx ++;
        };
      };
    }

    return out;


  }

  Eigen::MatrixXd get_loadings_vec2mat(const Eigen::Ref<const Eigen::VectorXd> LOADINGS,
                                       const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT){

    const int p = CONSTRMAT.rows();
    const int q = CONSTRMAT.cols();
    Eigen::MatrixXd out(p, q);

    int idx = 0;

    // construct vector with parameters free to be estimated
    for(int h = 0; h < q; h++){
      for(int j = 0; j < p; j++){
        if(!std::isfinite(CONSTRMAT(j,h))){
          out(j, h)=LOADINGS(idx);
          idx ++;
        }else{
          out(j, h)= CONSTRMAT(j, h);
        };
      };
    }

    return out;


  }

  Eigen::VectorXd get_latvar_mat2vec(const Eigen::Ref<const Eigen::MatrixXd> S){
    const unsigned int q = S.cols(); // number of latent variables

    Eigen::LLT<Eigen::MatrixXd> llt(S);
    Eigen::MatrixXd U = llt.matrixU();

    Eigen::VectorXd y(q*(q-1)/2);
    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(q,q);
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(q,q);
    unsigned int iterator = 0;
    for(unsigned int col = 1; col < q; col ++){
      for(unsigned int row = 0; row < col; row++){
        if(row == 0) Z(row, col) = U(row, col);
        else if(row > 0){
          double prod = 1;
          for(unsigned int i = 0; i < row; i++){
            prod *= pow(1-pow(Z(i, col), 2),-.5);
            Z(row, col) = U(row, col)*prod;
          }
        }
        Y(row, col) = atanh(Z(row, col));
        y(iterator) = atanh(Z(row, col));
        iterator++;
      }
    }

    return y;
  }

  Eigen::MatrixXd get_latvar_vec2mat(const Eigen::Ref<const Eigen::VectorXd> SVEC,
                                     const unsigned int Q){

     // extract uncostrained parameters related to latent correlations
     Eigen::VectorXd transformed_rhos = SVEC;

     // inverse of Fisher's transformation
     Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*transformed_rhos)).array().exp() - 1)/
       ((Eigen::VectorXd(2*transformed_rhos)).array().exp() + 1);

     // place tanh_entries in upper triangular matrix
     // with unit diagonal
     Eigen::MatrixXd Z(Q,Q); Z.setIdentity();
     unsigned int iterator = 0;
     for(unsigned int col = 1; col < Q; col ++){
       for(unsigned int row = 0; row < col; row++){
         Z(row, col) = tanh_entries(iterator);
         iterator ++;
       }
     }

     // construct upper Cholesky from Z
     Eigen::MatrixXd U(Q,Q);

     for( int col = 0; col < Q; col++){
       for( int row = 0; row < Q; row++){
         if(row > col) U(row, col) = 0;
         else if(row == col & row == 0) U(row, col) = 1;
         else if(row > 0 & row == col) {
           if(Z(row-1,col)==0){
             double prod = 1;
             for(unsigned int row1 = 0; row1 < row; row1++){
               prod *= pow(1-pow(Z(row1, col), 2),.5);
             }
             U(row, col) = prod;
           }else{
             U(row, col) = U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5)/Z(row-1, col);
           }
         }
         else if(row == 0 & row < col) U(row, col) = Z(row, col);
         else if(row > 0 & row < col) {
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


     // latent variable covariance matrix
     Eigen::MatrixXd S = U.transpose()*U;

     return S;
   }

  //' Get loading matrix from theta
  //'
  //' get_Lam() constructs the loading matrix from the parameter vector
  //'
  //' @param THETA Numerical vector of parameters.
  //' @param A Binary matrix of dimension \eqn{p*q} where \eqn{p} is the number
  //' of items and \eqn{q} the number of latent variables. Entries equal to
  //' \eqn{1} refer to free loadings, while entries equal to \eqn{0} indicate
  //' loadings constrained to be null.
  //' @param C Sum of the number of categories for each item.
  Eigen::MatrixXd get_loadings_theta2mat(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
      const int P,
      const int Q,
      const int D,
      const int C
  ){

    const Eigen::VectorXd lambda = THETA.segment(C-P, D-C+P-Q*(Q-1)/2);


    return get_loadings_vec2mat(lambda, CONSTRMAT);
  }


  Eigen::MatrixXd get_latvar_theta2mat(const Eigen::Ref<const Eigen::VectorXd> THETA,
                                       const int Q,
                                       const int D
  ){
    Eigen::VectorXd transformed_rhos = THETA.segment(D-Q*(Q-1)/2, Q*(Q-1)/2);

    return get_latvar_vec2mat(transformed_rhos, Q);
  }

  Eigen::VectorXd get_thresholds_theta2vec(const Eigen::Ref<const Eigen::VectorXd> THETA,
                          const unsigned int P,
                          const unsigned int C){
    return THETA.segment(0,C-P);
  }

  Eigen::VectorXd get_latvar_theta2vec(const Eigen::Ref<const Eigen::VectorXd> THETA,
                                   const unsigned int NTHR,
                                   const unsigned int NLOAD,
                                   const unsigned int NCORR){
    return THETA.segment(NTHR+NLOAD, NCORR);
  }


}


#endif
