// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pairs_freq
Eigen::MatrixXd pairs_freq(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC);
RcppExport SEXP _plFA_pairs_freq(SEXP YSEXP, SEXP C_VECSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    rcpp_result_gen = Rcpp::wrap(pairs_freq(Y, C_VEC));
    return rcpp_result_gen;
END_RCPP
}
// get_S
Eigen::MatrixXd get_S(const Eigen::VectorXd& THETA, const unsigned int Q);
RcppExport SEXP _plFA_get_S(SEXP THETASEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(get_S(THETA, Q));
    return rcpp_result_gen;
END_RCPP
}
// get_par_from_S
Eigen::VectorXd get_par_from_S(const Eigen::MatrixXd& S);
RcppExport SEXP _plFA_get_par_from_S(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(get_par_from_S(S));
    return rcpp_result_gen;
END_RCPP
}
// multiThread_completePairwise
Rcpp::List multiThread_completePairwise(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> TAU, Eigen::Map<Eigen::VectorXd> LAMBDA, Eigen::Map<Eigen::VectorXd> TRANSFORMED_RHOS, Eigen::Map<Eigen::MatrixXd> FREQ, int CORRFLAG, int GRFLAG, int SILENTFLAG);
RcppExport SEXP _plFA_multiThread_completePairwise(SEXP YSEXP, SEXP C_VECSEXP, SEXP ASEXP, SEXP TAUSEXP, SEXP LAMBDASEXP, SEXP TRANSFORMED_RHOSSEXP, SEXP FREQSEXP, SEXP CORRFLAGSEXP, SEXP GRFLAGSEXP, SEXP SILENTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TAU(TAUSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type LAMBDA(LAMBDASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TRANSFORMED_RHOS(TRANSFORMED_RHOSSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< int >::type GRFLAG(GRFLAGSEXP);
    Rcpp::traits::input_parameter< int >::type SILENTFLAG(SILENTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(multiThread_completePairwise(Y, C_VEC, A, TAU, LAMBDA, TRANSFORMED_RHOS, FREQ, CORRFLAG, GRFLAG, SILENTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// sampling_step
std::vector<int> sampling_step(std::vector<int>& FULL_POOL, const unsigned int METHODFLAG, const double PROB, const int PAIRS_PER_ITERATION, const unsigned int N_ITEMS, const unsigned int SEED, const int SILENTFLAG, const unsigned int ITER);
RcppExport SEXP _plFA_sampling_step(SEXP FULL_POOLSEXP, SEXP METHODFLAGSEXP, SEXP PROBSEXP, SEXP PAIRS_PER_ITERATIONSEXP, SEXP N_ITEMSSEXP, SEXP SEEDSEXP, SEXP SILENTFLAGSEXP, SEXP ITERSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int>& >::type FULL_POOL(FULL_POOLSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type METHODFLAG(METHODFLAGSEXP);
    Rcpp::traits::input_parameter< const double >::type PROB(PROBSEXP);
    Rcpp::traits::input_parameter< const int >::type PAIRS_PER_ITERATION(PAIRS_PER_ITERATIONSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_ITEMS(N_ITEMSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SEED(SEEDSEXP);
    Rcpp::traits::input_parameter< const int >::type SILENTFLAG(SILENTFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ITER(ITERSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_step(FULL_POOL, METHODFLAG, PROB, PAIRS_PER_ITERATION, N_ITEMS, SEED, SILENTFLAG, ITER));
    return rcpp_result_gen;
END_RCPP
}
// compute_pair
Rcpp::List compute_pair(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::VectorXd THETA, const int CORRFLAG, const unsigned int k, const unsigned int l, Eigen::MatrixXd PAIRS_TABLE, const unsigned int SILENTFLAG, const unsigned int GRADFLAG);
RcppExport SEXP _plFA_compute_pair(SEXP ASEXP, SEXP C_VECSEXP, SEXP THETASEXP, SEXP CORRFLAGSEXP, SEXP kSEXP, SEXP lSEXP, SEXP PAIRS_TABLESEXP, SEXP SILENTFLAGSEXP, SEXP GRADFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type l(lSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type PAIRS_TABLE(PAIRS_TABLESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SILENTFLAG(SILENTFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type GRADFLAG(GRADFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pair(A, C_VEC, THETA, CORRFLAG, k, l, PAIRS_TABLE, SILENTFLAG, GRADFLAG));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plFA_pairs_freq", (DL_FUNC) &_plFA_pairs_freq, 2},
    {"_plFA_get_S", (DL_FUNC) &_plFA_get_S, 2},
    {"_plFA_get_par_from_S", (DL_FUNC) &_plFA_get_par_from_S, 1},
    {"_plFA_multiThread_completePairwise", (DL_FUNC) &_plFA_multiThread_completePairwise, 10},
    {"_plFA_sampling_step", (DL_FUNC) &_plFA_sampling_step, 8},
    {"_plFA_compute_pair", (DL_FUNC) &_plFA_compute_pair, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_plFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
