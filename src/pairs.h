#ifndef pairs_H
#define pairs_H

#include "genericUtils.h"
#include "bivariateProbs.h"

// Single pair contribution to
// 1. Log-likelihood
// 2. Gradient
void pair_contribution(
    // Parameters
    Eigen::Map<Eigen::MatrixXd> A,
    Eigen::Map<Eigen::VectorXd> C_VEC,
    const Eigen::VectorXd &THETA,
    const int CORRFLAG,

    // Input:
    const unsigned int k,
    const unsigned int l,
    const Eigen::MatrixXd &PAIRS_TABLE,

    // Options:
    const unsigned int SILENTFLAG,
    const unsigned int GRADFLAG,

    // Outputs:
    double &ll,
    Eigen::VectorXd &gradient
){
  const unsigned int p = A.rows();
  const unsigned int q = A.cols();
  const unsigned int d = THETA.size();
  const unsigned int c = C_VEC.sum();
  const unsigned int nthr = c-p;
  const unsigned int ncorr = q*(q-1)/2;
  const unsigned int nload = d-nthr-ncorr;

  double tmp_ll = 0;
  Eigen::VectorXd tmp_gradient = Eigen::VectorXd::Zero(d);

  // rearrange parameters
  Eigen::MatrixXd Lam            = get_Lam(A, c, THETA);
  Eigen::MatrixXd Sigma_u        = get_S(THETA, q);
  Eigen::VectorXd tau            = THETA.segment(0,c-p);
  Eigen::VectorXd transformed_rhos=THETA.segment(nthr+nload, ncorr);

  // Identifies quantities related to pair (k,l)
  const unsigned int ck = C_VEC(k);
  const unsigned int cl = C_VEC(l);
  const Eigen::VectorXd lambdak = Lam.row(k);
  const Eigen::VectorXd lambdal = Lam.row(l);
  const double rho_kl = lambdak.transpose() * Sigma_u * lambdal;
  Eigen::MatrixXd pairs_tab = PAIRS_TABLE;
  pairs_tab.conservativeResize(PAIRS_TABLE.rows() + 1, Eigen::NoChange_t() );
  // identify column index in freq table
  // i1: starting index item k
  unsigned int i1 = 0;
  if(k > 1){
    for(unsigned int u = 1; u < k; u++){
      const unsigned int cu = C_VEC(u);
      //if(SILENTFLAG == 0)Rcpp::Rcout << "u: " << u << ", cu: "<< cu << "\n";
      i1 += cu * C_VEC.segment(0,u).sum();
    }
  }

  // i2 starting index from i1 for item l
  unsigned int i2 = 0;
  if(l > 0){
    i2 = C_VEC.segment(0,l).sum() * C_VEC(k);
  }

  ////////////////////////////
  /* LIKELIHOOD COMPUTATION */
  ////////////////////////////
  for(unsigned int sk = 0; sk < ck; sk ++){

    // i3: starting index from i2 for cat sk
    const unsigned int i3 = sk * cl;

    for(unsigned int sl = 0; sl < cl; sl ++){

      // final column index for pairs_tab. Print to check
      const unsigned int r = i1 + i2 + i3 + sl;

      // read frequency
      const unsigned int n_sksl = PAIRS_TABLE(4, r);

      // identify thresholds
      const Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);

      // compute pi
      const double pi_sksl = compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
      if(SILENTFLAG == 0)Rcpp::Rcout << "("<<k<<","<<l<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
      pairs_tab(5,r) = pi_sksl;

      // update ll
      tmp_ll += n_sksl * log(pi_sksl+1e-8);
    }
  }

  //////////////////////////
  /* GRADIENT COMPUTATION */
  /////////////////////////

  if(GRADFLAG == 1){
    unsigned int iter = 0;

    /////////////////////////////////////////////////////
    // (k,l)-pair likelihood derivative wrt thresholds //
    /////////////////////////////////////////////////////
    if(SILENTFLAG == 0)Rcpp::Rcout << "- Gradient wrt thresholds: \n";

    // loop: iterate over elements of thresholds vector
    for(unsigned int s = 0; s < tau.size(); s++){
      double grs = 0; // temporary location for gradient related to s-th element of tau
      if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ gradient("<< s<< ")\n";

      // List three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
      if(s >= (C_VEC.segment(0, k).sum()) - (k) & s < C_VEC.segment(0, k + 1).sum() - (k + 1)){
        // [CASE 1]: threshold related to item k

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item k:\n";
        const unsigned int sk = s - (C_VEC.segment(0, k).sum()) + (k);

        // i3: starting index from i2 for cat sk and sk+1
        const unsigned int i3 = sk * cl;
        const unsigned int i3suc = (sk+1) * cl;
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ sk: " << sk << ". Summing over categories item l: ";

        // iterate over categories of item l
        for(unsigned int sl = 0; sl < cl; sl ++){
          if(SILENTFLAG == 0)Rcpp::Rcout << " ... cat" << sl ;

          // identify pairs_tab column for (sk,sl) and (sk+1, sl)
          const unsigned int r = i1 + i2 + i3 + sl;
          const unsigned int rsuc = i1 + i2 + i3suc + sl;

          // read frequences
          const unsigned int n_sksl = pairs_tab(4, r);
          const unsigned int n_sksucsl = pairs_tab(4, rsuc);

          // read probabilities
          const double pi_sksl = pairs_tab(5, r);
          const double pi_sksucsl = pairs_tab(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);

          // compute gradient
          const double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_sksucsl/(pi_sksucsl+1e-8)));
          const double tmp2 = R::dnorm(t_sk, 0, 1, 0);
          const double tmp3 = R::pnorm((t_sl-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          const double tmp4 = R::pnorm((t_sl_prev-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          grs += tmp1 * tmp2 * (tmp3 - tmp4);
        }
        if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

      }else if(s >= (C_VEC.segment(0, l).sum())-(l) & s<C_VEC.segment(0, l + 1).sum()-(l + 1)){
        // [CASE 2]: threshold related to item l

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item l\n";
        const unsigned int sl = s - (C_VEC.segment(0, l).sum()) + (l);

        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  sl: " << sl << ". Summing over categories item k: ";

        // iterate over categories item k
        for(unsigned int sk = 0; sk < ck; sk ++){

          // i3: starting index from i2 for cat sk
          const unsigned int i3 = sk * cl;

          // identify pairs_tab column for (sk,sl) and (sk, sl + 1)
          const unsigned int r = i1 + i2 + i3 + sl;
          const unsigned int rsuc = i1 + i2 + i3 + sl + 1;

          // read frequences
          const unsigned int n_sksl = pairs_tab(4, r);
          const unsigned int n_skslsuc = pairs_tab(4, rsuc);

          // read probabilities
          const double pi_sksl = pairs_tab(5, r);
          const double pi_skslsuc = pairs_tab(5, rsuc);

          // identify tau_sk, tau_sl, tau_sl-1
          const Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
          const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);


          if(SILENTFLAG == 0)Rcpp::Rcout<<"\n  |    |   |_ sk:"<< sk << ", r: "<< r<<", n_sksl:"
                                        << n_sksl<< ", n_sksl+1:" << n_skslsuc << ", pi_sksl:"
                                        << pi_sksl << ", pi_sksl+1:"<< pi_skslsuc << ", t_sk:"
                                        << t_sk<< ", t_sl:" << t_sl << "t_sk-1:"<< t_sk_prev;
          // compute gradient
          const double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_skslsuc/(pi_skslsuc+1e-8)));
          const double tmp2 = R::dnorm(t_sl, 0, 1, 0);
          const double tmp3 = R::pnorm((t_sk-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          const double tmp4 = R::pnorm((t_sk_prev-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
          if(SILENTFLAG == 0)Rcpp::Rcout<<" => out" << sk << ":" << tmp1 * tmp2 * (tmp3 - tmp4);
          grs += tmp1 * tmp2 * (tmp3 - tmp4);
        }
        if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

      }else{
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  tau of other item\n";
      }
      if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  Thresholds:: " << iter<< "/"<< tau.size()
                                     <<". Tot:: "<< iter << "/"<< d <<". Val ="<< grs <<"\n";
      tmp_gradient(iter) += grs;
      iter ++;
    }
    if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ Done. \n";

    ///////////////////////////////////////////////////////////
    // (k,l)-pair likelihood derivative wrt URV correlation: //
    // intermediate step for derivatives wrt                 //
    // loadings and factor correlations                      //
    ///////////////////////////////////////////////////////////
    if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Intermediate derivative for loadings and correlation: \n";
    double tmp_kl = 0; // temporary location of the gradient

    // double loop: iterate over each combination of categories of items k and l
    for(unsigned int sk = 0; sk < ck; sk ++){
      for(unsigned int sl = 0; sl < cl; sl ++){
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ sk: "<< sk << ", sl: " << sl << ": \n";

        // identify pairs_tab column for (sk,sl)
        const unsigned int i3 = sk * cl;
        const unsigned int r = i1 + i2 + i3 + sl;

        // read freq
        const unsigned int n_sksl = pairs_tab(4, r);

        // read prob
        const double pi_sksl = pairs_tab(5, r);

        // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
        const Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
        const double t_sk = pi_thresholds(0); const double t_sl = pi_thresholds(1); const double t_sk_prev = pi_thresholds(2); const double t_sl_prev = pi_thresholds(3);

        // phi(t_sk, t_sl; rho_kl)
        const double d1 = pbv_rcpp_dbvnorm0( t_sk, t_sl, rho_kl, 0);

        // phi(t_sk, t_sl-1; rho_kl)
        const double d2 = pbv_rcpp_dbvnorm0( t_sk, t_sl_prev, rho_kl, 0);

        // phi(t_sk-1, t_sl; rho_kl)
        const double d3 = pbv_rcpp_dbvnorm0( t_sk_prev, t_sl, rho_kl, 0);

        // phi(t_sk-1, t_sl-1; rho_kl)
        const double d4 = pbv_rcpp_dbvnorm0( t_sk_prev, t_sl_prev, rho_kl, 0);

        tmp_kl += (n_sksl/(pi_sksl+1e-8)) * ( d1 - d2 - d3 + d4);
      }
    }
    if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ tmp_kl:" << tmp_kl << "\n";

    ///////////////////////////////////////////////////
    // (k,l)-pair likelihood derivative wrt loadings //
    ///////////////////////////////////////////////////
    if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Gradient wrt loadings: \n";

    // double loop: iterate over elements of loadings matrix
    for(unsigned int j = 0; j < p; j++){
      for(unsigned int v = 0; v < q; v++){
        if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ visiting lambda_"<< j << v <<":\n";

        // elicit three cases: 1. free loading item k, 2. free loading l, 3. other
        if(j == k){
          if(A(j,v)!=0 ){
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item k, free loading:\n";
            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
            const double d_rho_kl = ev.transpose() * Sigma_u * lambdal;
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            tmp_gradient(iter) += tmp_kl * d_rho_kl;
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << "\n";

            iter ++;
          }
        }else if (j == l){
          if(A(j,v)!=0 ){
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item l, free loading:\n";
            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
            const double d_rho_kl = lambdak.transpose() * Sigma_u * ev;
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << "\n";
            tmp_gradient(iter) += tmp_kl * d_rho_kl;

            iter ++;
          }
        }else if(A(j,v)!=0){

          if(SILENTFLAG == 0)Rcpp::Rcout << "  |  |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << " [not included]\n";
          iter ++;
        }
      }
    }
    if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ Done. \n";

    //////////////////////////////////////////////////////////////
    // (k,l)-pair likelihood derivative wrt latent correlations //
    //////////////////////////////////////////////////////////////
    if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Gradient wrt correlations: \n";

    if(CORRFLAG == 1){

      for(int thro_idx = 0; thro_idx < ncorr; thro_idx ++){
        const Eigen::MatrixXd dSigma = grad_S(A, transformed_rhos, thro_idx);
        const double d_rho_kl = lambdak.transpose() * dSigma * lambdal;
        tmp_gradient(iter) += tmp_kl * d_rho_kl;
        iter++;
      }
    }

    if(SILENTFLAG == 0)Rcpp::Rcout << "\n=====> gradient r-th pair:\n" << gradient << "\n";
  }

  ll += tmp_ll;
  gradient += tmp_gradient;

}

//' Single pair contribution
//'
//' Wrapper of pair_contribution() used for unit tests
//' @export
// [[Rcpp::export]]
Rcpp::List compute_pair(
     Eigen::Map<Eigen::MatrixXd> A,
     Eigen::Map<Eigen::VectorXd> C_VEC,
     Eigen::VectorXd THETA,
     const int CORRFLAG,
     const unsigned int k,
     const unsigned int l,
     Eigen::MatrixXd PAIRS_TABLE,
     const unsigned int SILENTFLAG,
     const unsigned int GRADFLAG
 ){
   const unsigned int p = A.rows();
   const unsigned int q = A.cols();
   const unsigned int d = THETA.size();
   const unsigned int c = C_VEC.sum();
   const unsigned int nthr = c-p;
   const unsigned int ncorr = q*(q-1)/2;
   const unsigned int nload = d-nthr-ncorr;

   double ll = 0;
   Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);

   pair_contribution(A, C_VEC, THETA, CORRFLAG, k, l, PAIRS_TABLE, SILENTFLAG, GRADFLAG, ll, gradient);

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("ll") = -ll,
       Rcpp::Named("ngradient") = -gradient
     );
   return(output);

 }

// RcppParallel Worker to compute pair_contribution
// on multiple pairs in parallel
struct SubsetWorker : public RcppParallel::Worker{
  // Declaration parameters:
  //// Global:
  const Eigen::Map<Eigen::MatrixXd> constrmat;
  const Eigen::Map<Eigen::VectorXd> c_vec;
  const Eigen::MatrixXd &pairs_table;
  const Eigen::MatrixXd &items_pairs;
  const unsigned int corrFLAG;
  const unsigned int silentFLAG;
  const unsigned int gradFLAG;

  //// Iteration:
  const Eigen::VectorXd &theta;
  const std::vector<int> &index_vector;

  // Output quantities:
  double subset_ll;
  Eigen::VectorXd subset_gradient, subset_gradient2;

  // Constructor 1:
  SubsetWorker(
    const Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
    const Eigen::Map<Eigen::VectorXd> C_VEC,
    const Eigen::MatrixXd &PAIRS_TABLE,
    const Eigen::MatrixXd &ITEMS_PAIRS,
    const unsigned int CORRFLAG,
    const unsigned int SILENTFLAG,
    const unsigned int GRADFLAG,
    const Eigen::VectorXd &THETA,
    const std::vector<int> &INDEX_VECTOR
  ):
    constrmat(CONSTRMAT), c_vec(C_VEC), pairs_table(PAIRS_TABLE), items_pairs(ITEMS_PAIRS),
    corrFLAG(CORRFLAG), silentFLAG(SILENTFLAG), gradFLAG(GRADFLAG),
    theta(THETA), index_vector(INDEX_VECTOR), subset_ll(0.0),
    subset_gradient(Eigen::VectorXd::Zero(THETA.size())){}

  // Constructor 2:
  SubsetWorker(const SubsetWorker &OBJ, RcppParallel::Split):
    constrmat(OBJ.constrmat), c_vec(OBJ.c_vec), pairs_table(OBJ.pairs_table), items_pairs(OBJ.items_pairs),
    corrFLAG(OBJ.corrFLAG), silentFLAG(OBJ.silentFLAG), gradFLAG(OBJ.gradFLAG),
    theta(OBJ.theta), index_vector(OBJ.index_vector), subset_ll(0.0),
    subset_gradient(Eigen::VectorXd::Zero(theta.size())){}

  // MEMBER FUNCTIONS:
  //// Overload operator () for main computation:
  void operator()(std::size_t BEGIN, std::size_t END);
  void join(const SubsetWorker &RHS);

};

void SubsetWorker::operator()(std::size_t BEGIN, std::size_t END){
  // Run along the pairs identified by indexes in index_vector
  // computing their contribution to nll, gradient and eventually Hessian.
  unsigned int d = theta.size();
  //subset_gradient.resize(d); subset_gradient.setZero();
  for (unsigned int h = BEGIN; h < END; h++){

    // identify corresponding column in pairs_table
    unsigned int col = index_vector[h];

    // identify the pair
    unsigned int k = items_pairs(0, col);
    unsigned int l = items_pairs(1, col);

    // initialize empty pair-output
    double pair_ll = 0;
    Eigen::VectorXd pair_gradient = Eigen::VectorXd::Zero(d);

    // computation of log-likelihood, gradient
    pair_contribution(constrmat, c_vec, theta, corrFLAG, k, l, pairs_table, silentFLAG, gradFLAG, pair_ll, pair_gradient);

    // update
    {
      subset_ll += pair_ll;
      subset_gradient += pair_gradient;
    }
  }
}

void SubsetWorker::join(const SubsetWorker &RHS){
  subset_ll += RHS.subset_ll;
  for(unsigned i = 0; i < subset_gradient.size(); i++){
    subset_gradient(i) += RHS.subset_gradient(i);

  }
}


#endif
