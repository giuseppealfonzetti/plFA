#ifndef latCorr_H
#define latCorr_H

namespace latcorr{
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

    Eigen::MatrixXd dU(q,q); dU.setZero();

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
}
#endif
