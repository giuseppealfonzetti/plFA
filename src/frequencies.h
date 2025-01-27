#ifndef frequencies_H
#define frequencies_H

//' Compute pairwise frequencies
//'
//' @param Y Integer matrix of dimension \eqn{n*p}, where \eqn{n} is the sample size
//' and \eqn{p} is the number of items considered. Categories must be coded starting from zero.
//' For example, an item with three categories can only accept values contained in
//' \eqn{\{0, 1, 2\}}.
//' @param C_VEC Integer vector indicating how many possible categories are associated to
//' each item in 'Y'.
//'
//' @return
//' It returns a 5-rows matrix with each combination of items and categories as columns.
//' Row0: item k, Row1: item l, Row2; category item k, Row3: category item l, Row4: freq
//' It is computed just once, before the optimization of the complete pairwise
//'
//' @export
// [[Rcpp::export]]
 Eigen::MatrixXd pairs_freq(
     Eigen::Map<Eigen::MatrixXd> Y,
     Eigen::Map<Eigen::VectorXd> C_VEC
 ){

   const unsigned int n = Y.rows(); // number of units
   const unsigned int p = Y.cols(); // number of items
   unsigned int iter; // helper iterator

   Eigen::MatrixXd freq(5,1);

   // Find how many possible pairs and setup freq matrix
   iter = 0;
   for(unsigned int k = 1; k < p; k++){
     const unsigned int ck = C_VEC(k);
     for(unsigned int l = 0; l < k; l ++){
       const unsigned int cl = C_VEC(l);
       for(unsigned int sk = 0; sk < ck; sk++){
         for(unsigned int sl = 0; sl < cl; sl ++){
           freq.conservativeResize(5, iter + 1);
           Eigen::VectorXd freq_coord(5);
           freq_coord << k, l, sk, sl, 0;
           freq.col(iter) = freq_coord;
           iter++;
         }
       }
     }
   }

   auto task = [&freq, &n, &Y](unsigned int r){

     // read data
     const unsigned int k  = freq(0,r);
     const unsigned int l  = freq(1,r);
     const unsigned int sk = freq(2,r);
     const unsigned int sl = freq(3,r);
     Eigen::MatrixXd obs_resp(n,2);
     obs_resp.col(0) = Y.col(k); obs_resp.col(1) = Y.col(l);

     // compute frequency
     unsigned int n_sksl = 0;
     for(unsigned int i = 0; i < n; i++){
       if(sk == obs_resp(i,0) && sl == obs_resp(i,1)) {
         n_sksl++;
       }
     }

     // update
     freq(4, r) = n_sksl;
   };

   // Loop to compute frequencies
   for(unsigned int r = 0; r < freq.cols(); r++){
     task(r);
   }


   return freq;
 }

#endif
