#ifndef utils_H
#define utils_H

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
Eigen::MatrixXd get_Lam(Eigen::Map<Eigen::MatrixXd> A,
                        const unsigned int C,
                        const Eigen::VectorXd &THETA
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

/* EXTRACT PI THRESHOLDS */
// Extract from tau the trhesholds related to pi_sksl
Eigen::VectorXd extract_thresholds(const Eigen::VectorXd tau,
                                   const Eigen::VectorXd c_vec,
                                   const unsigned int k,
                                   const unsigned int l,
                                   const unsigned int sk,
                                   const unsigned int sl){

  unsigned int ck = c_vec(k);
  unsigned int cl = c_vec(l);

  // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
  unsigned int sk_tau_index = c_vec.segment(0, k).sum() - (k) + sk;   // index tau_sk in tau vector
  unsigned int sl_tau_index = c_vec.segment(0, l).sum() - (l) + sl;   // index tau_sl in tau vector
  double t_sk;                                               // tau_sk
  if(sk == ck-1){
    t_sk = 100;
  } else {
    t_sk = tau(sk_tau_index);
  }
  double t_sk_prev;                                           // tau_sk-1
  if(sk == 0){
    t_sk_prev = -100;
  }else{
    t_sk_prev = tau(sk_tau_index-1);
  }

  double t_sl;                                               // tau_sl
  if(sl == cl-1){
    t_sl = 100;
  } else {
    t_sl = tau(sl_tau_index);
  }
  double t_sl_prev;                                           // tau_sl-1
  if(sl == 0){
    t_sl_prev = -100;
  }else{
    t_sl_prev = tau(sl_tau_index-1);
  }

  Eigen::VectorXd pi_thresholds(4);
  pi_thresholds << t_sk, t_sl, t_sk_prev, t_sl_prev;
  return pi_thresholds;
}
//
// //' Get latent correlation matrix from theta
// //'
// //' get_S() extracts the latent correlation matrix from theta assuming
// //' theta elements to be reparametrised following the
// //' Lewandowski-Kurowicka-Joe (2009) transform.
// //'
// //' @param THETA Numerical vector of parameters.
// //' @param Q Number of latent variables.
// //'
// //' @export
//  // [[Rcpp::export]]
//  Eigen::MatrixXd get_S(const Eigen::VectorXd &THETA,
//                              const unsigned int Q
//  ){
//    unsigned int d = THETA.size(); // number of parameters
//
//    // extract uncostrained parameters related to latent correlations
//    Eigen::VectorXd transformed_rhos = THETA.segment(d-Q*(Q-1)/2, Q*(Q-1)/2);
//
//    // inverse of Fisher's transformation
//    Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*transformed_rhos)).array().exp() - 1)/
//      ((Eigen::VectorXd(2*transformed_rhos)).array().exp() + 1);
//
//    // place tanh_entries in upper triangular matrix
//    // with unit diagonal
//    Eigen::MatrixXd Z(Q,Q); Z.setIdentity();
//    unsigned int iterator = 0;
//    for(unsigned int col = 1; col < Q; col ++){
//      for(unsigned int row = 0; row < col; row++){
//        Z(row, col) = tanh_entries(iterator);
//        iterator ++;
//      }
//    }
//
//    // construct upper Cholesky from Z
//    Eigen::MatrixXd U(Q,Q);
//
//    for(unsigned int col = 0; col < Q; col++){
//      for(unsigned int row = 0; row < Q; row++){
//        if(row > col) U(row, col) = 0;
//        else if(row == col & row == 0) U(row, col) = 1;
//        else if(row > 0 & row == col) {
//          double prod = 1;
//          for(unsigned int row1 = 0; row1 < row; row1++){
//            prod *= pow(1-pow(Z(row1, col), 2),.5);
//          }
//          U(row, col) = prod;
//        }
//        else if(row == 0 & row < col) U(row, col) = Z(row, col);
//        else if(row > 0 & row < col) {
//          double prod = Z(row, col);
//          for(unsigned int row1 = 0; row1 < row; row1++){
//            prod *= pow(1-pow(Z(row1, col), 2),.5);
//          }
//          U(row, col) = prod;
//
//        }
//      }
//    }
//    // latent variable covariance matrix
//    Eigen::MatrixXd S = U.transpose()*U;
//
//    return S;
//  }



//' Get transformed parameters from latent correlation matrix
//'
//' Use Lewandowski-Kurowicka-Joe (2009) transformation on latent correlations
//' with Choleski decomposition.
//'
//' @param S Latent correlation matrix.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_par_from_S(const Eigen::MatrixXd &S){
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


// /* LATENT COVARIANCE */
// // Get latent correlation matrix from theta using
// // Lewandowski-Kurowicka-Joe (2009) transform
//  // [[Rcpp::export]]
// Eigen::MatrixXd grad_S(Eigen::Map<Eigen::MatrixXd> A,
//                               const Eigen::VectorXd &THETA,
//                               const unsigned int IDX
//  ){
//    unsigned int q = A.cols(); // number of latents
//    unsigned int d = THETA.size(); // number of parameters
//
//    // extract uncostrained parameters related to latent correlations
//    Eigen::VectorXd transformed_rhos = THETA.segment(d-q*(q-1)/2, q*(q-1)/2);
//
//    // inverse of Fisher's transformation
//    Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*transformed_rhos)).array().exp() - 1)/
//      ((Eigen::VectorXd(2*transformed_rhos)).array().exp() + 1);
//
//    // place tanh_entries in upper triangular matrix
//    // with unit diagonal
//    Eigen::MatrixXd Z(q,q); Z.setIdentity();
//    unsigned int iterator = 0; unsigned int i,j;
//    for(unsigned int col = 1; col < q; col ++){
//      for(unsigned int row = 0; row < col; row++){
//        if(iterator==IDX){i = row; j = col;}
//        Z(row, col) = tanh_entries(iterator);
//        iterator ++;
//      }
//    }
//
//    // construct upper Cholesky from Z
//    Eigen::MatrixXd U(q,q);
//    // for(unsigned int row = 0; row < q; row++){
//    //     for(unsigned int col = 0; col < q; col++){
//    //         if(row > col) U(row, col) = 0;
//    //         else if(row == col & row == 0) U(row, col) = 1;
//    //         else if(row == 0 & row < col) U(row, col) = Z(row, col);
//    //         else if(row > 0 & row <= col) U(row, col) = Z(row, col)/Z(row-1, col)*U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5);
//    //     }
//    // }
//    for( int col = 0; col < q; col++){
//      for( int row = 0; row < q; row++){
//        if(row > col) U(row, col) = 0;
//        else if(row == col & row == 0) U(row, col) = 1;
//        else if(row > 0 & row == col) {
//          double prod = 1;
//          for(unsigned int row1 = 0; row1 < row; row1++){
//            prod *= pow(1-pow(Z(row1, col), 2),.5);
//          }
//          U(row, col) = prod;
//        }
//        else if(row == 0 & row < col) U(row, col) = Z(row, col);
//        else if(row > 0 & row < col) {
//          double prod = Z(row, col);
//          for( int row1 = 0; row1 < row; row1++){
//            prod *= pow(1-pow(Z(row1, col), 2),.5);
//          }
//          U(row, col) = prod;
//
//        }
//      }
//    }
//    Eigen::MatrixXd dU(q,q); dU.setZero();
//    // for(unsigned int row = 0; row < q; row++){
//    //     if(row > j) dU(row, j) = 0;
//    //     else if(row == j & row == 0) dU(row, j) = 0;
//    //     else if(row == 0 & row < j & row == i ) dU(row, j) = 1;
//    //     else if(row > 0 & row <= j & row == i) dU(row, j) = 1/Z(row-1, j)*U(row-1, j)*pow(1-pow(Z(row-1, j),2), .5);
//    //     else if(row > 0 & row <= j & i == (row-1)) dU(row, j) = -Z(row, j)*pow(Z(row-1, j),-2)*U(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) + Z(row, j)/Z(row-1, j)*(dU(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) - U(row-1, j)* pow(1-pow(Z(row-1, j),2), -.5)*Z(row-1, j));
//    //     else if(row > 0 & row <= j & (row - i) > 1) dU(row, j) = dU(row-1, j)*Z(row, j)/Z(row-1, j)*pow(1-pow(Z(row-1, j),2), .5);
//    //     else if(row > 0 & row <= j) dU(row, j) = 2;
//    // }
//
//    for( int row = 0; row < q; row++){
//      if(row > j) dU(row, j) = 0;
//      else if(row == j & row == 0) dU(row, j) = 0;
//      else if(row > 0 & row == j & row == i)dU(row, j) = 0;
//      else if(row > 0 & row == j & i < row){
//        dU(row, j) = - U(row, j) * Z(i, j)/ (1-pow(Z(i,j), 2));
//      }
//      else if(row == 0 & row < j & row == i) dU(row, j) = 1;
//      else if(row > 0 & row < j & row == i ){
//        double prod = 1;
//        for( int row1 = 0; row1 < row; row1++){
//          prod *= pow(1-pow(Z(row1, j), 2),.5);
//        }
//        dU(row, j) = prod;
//      }
//      else if(row > 0 & row < j & i < row ){
//        double prod = Z(row, j);
//        for( int row1 = 0; row1 < row; row1++){
//          prod *= pow(1-pow(Z(row1, j), 2),.5);
//        }
//        dU(row, j) = - prod * Z(i, j) / (1-pow(Z(i,j), 2));
//      }
//    }
//
//    double dz = pow(cosh(transformed_rhos(IDX)),-2);
//    dU = dU*dz;
//
//    Eigen::MatrixXd dS = dU.transpose()*U + U.transpose()*dU;
//
//    return dS;
//  }

//' Get latent correlation matrix from theta
//'
//' get_S() extracts the latent correlation matrix from theta assuming
//' theta elements to be reparametrised following the
//' Lewandowski-Kurowicka-Joe (2009) transform.
//'
//' @param THETA Numerical vector of parameters.
//' @param Q Number of latent variables.
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_S(const Eigen::VectorXd &THETA,
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


// [[Rcpp::export]]
Eigen::MatrixXd grad_S(Eigen::Map<Eigen::MatrixXd> A,
                        const Eigen::VectorXd &THETA,
                        const unsigned int IDX
 ){
   unsigned int q = A.cols(); // number of latents
   unsigned int d = THETA.size(); // number of parameters

   // extract uncostrained parameters related to latent correlations
   Eigen::VectorXd transformed_rhos = THETA.segment(d-q*(q-1)/2, q*(q-1)/2);

   // inverse of Fisher's transformation
   Eigen::VectorXd tanh_entries = ((Eigen::VectorXd(2*transformed_rhos)).array().exp() - 1)/
     ((Eigen::VectorXd(2*transformed_rhos)).array().exp() + 1);

   // place tanh_entries in upper triangular matrix
   // with unit diagonal
   Eigen::MatrixXd Z(q,q); Z.setIdentity();
   unsigned int iterator = 0; unsigned int i,j;
   for(unsigned int col = 1; col < q; col ++){
     for(unsigned int row = 0; row < col; row++){
       if(iterator==IDX){i = row; j = col;}
       Z(row, col) = tanh_entries(iterator);
       iterator ++;
     }
   }

   // construct upper Cholesky from Z
   Eigen::MatrixXd U(q,q);
   for( int col = 0; col < q; col++){
     for( int row = 0; row < q; row++){
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
   // for(unsigned int col = 0; col < q; col++){
   //   for(unsigned int row = 0; row < q; row++){
   //
   //     if(row > col) U(row, col) = 0;
   //     else if(row == col & row == 0) U(row, col) = 1;
   //     else if(row == 0 & row < col) U(row, col) = Z(row, col);
   //     else if(row > 0 & row < col) {
   //       U(row, col) = Z(row, col)/Z(row-1, col)*U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5);
   //     }
   //     else if(row > 0 & row == col) {
   //       U(row, col) = U(row-1, col)*pow(1-pow(Z(row-1, col),2), .5)/Z(row-1, col);
   //     }
   //
   //   }
   // }
   Eigen::MatrixXd dU(q,q); dU.setZero();
   // for(unsigned int row = 0; row < q; row++){
   //     if(row > j) dU(row, j) = 0;
   //     else if(row == j & row == 0) dU(row, j) = 0;
   //     else if(row == 0 & row < j & row == i ) dU(row, j) = 1;
   //     else if(row > 0 & row <= j & row == i) dU(row, j) = 1/Z(row-1, j)*U(row-1, j)*pow(1-pow(Z(row-1, j),2), .5);
   //     else if(row > 0 & row <= j & i == (row-1)) dU(row, j) = -Z(row, j)*pow(Z(row-1, j),-2)*U(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) + Z(row, j)/Z(row-1, j)*(dU(row-1, j)*pow(1-pow(Z(row-1, j),2), .5) - U(row-1, j)* pow(1-pow(Z(row-1, j),2), -.5)*Z(row-1, j));
   //     else if(row > 0 & row <= j & (row - i) > 1) dU(row, j) = dU(row-1, j)*Z(row, j)/Z(row-1, j)*pow(1-pow(Z(row-1, j),2), .5);
   //     else if(row > 0 & row <= j) dU(row, j) = 2;
   // }

   for( int row = 0; row < q; row++){
     if(row > j) dU(row, j) = 0;
     else if(row == j & row == 0) dU(row, j) = 0;
     else if(row == 0 & row < j & row == i) dU(row, j) = 1;
     else if(row > 0 & row < j & row == i ){
       if(Z(row-1, j)==0){
         double prod = 1;
         for( int row1 = 0; row1 < row; row1++){
           prod *= pow(1-pow(Z(row1, j), 2),.5);
         }
         dU(row, j) = prod;
       }else{
         dU(row, j) = U(row-1,j)*pow(1-pow(Z(row-1, j),2), .5)/Z(row-1, j);
       }}
     else if(row > 0 & row < j & i < row ){
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
     else if(row > 0 & row == j & row == i) dU(row, j) = 0;
     else if(row > 0 & row == j & i < row){
       dU(row, j) = - U(row, j) * Z(i, j)/ (1-pow(Z(i,j), 2));
     }
   }

   double dz = pow(cosh(transformed_rhos(IDX)),-2);
   dU = dU*dz;

   Eigen::MatrixXd dS = dU.transpose()*U + U.transpose()*dU;

   return dS;
 }

#endif
