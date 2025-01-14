// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_get_Lam
Eigen::MatrixXd cpp_get_Lam(Eigen::Map<Eigen::MatrixXd> A, const unsigned int C, Eigen::Map<Eigen::VectorXd> THETA);
RcppExport SEXP _plFA_cpp_get_Lam(SEXP ASEXP, SEXP CSEXP, SEXP THETASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type C(CSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_Lam(A, C, THETA));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_S
Eigen::MatrixXd cpp_get_S(Eigen::Map<Eigen::VectorXd> THETA, const unsigned int Q);
RcppExport SEXP _plFA_cpp_get_S(SEXP THETASEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_S(THETA, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_par_from_S
Eigen::VectorXd cpp_get_par_from_S(Eigen::Map<Eigen::MatrixXd> S);
RcppExport SEXP _plFA_cpp_get_par_from_S(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_par_from_S(S));
    return rcpp_result_gen;
END_RCPP
}
// cpp_compute_pair
Rcpp::List cpp_compute_pair(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::VectorXd THETA, const int CORRFLAG, const unsigned int k, const unsigned int l, Eigen::MatrixXd PAIRS_TABLE, const unsigned int SILENTFLAG, const unsigned int GRADFLAG);
RcppExport SEXP _plFA_cpp_compute_pair(SEXP ASEXP, SEXP C_VECSEXP, SEXP THETASEXP, SEXP CORRFLAGSEXP, SEXP kSEXP, SEXP lSEXP, SEXP PAIRS_TABLESEXP, SEXP SILENTFLAGSEXP, SEXP GRADFLAGSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cpp_compute_pair(A, C_VEC, THETA, CORRFLAG, k, l, PAIRS_TABLE, SILENTFLAG, GRADFLAG));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_thresholds_theta2vec
Eigen::VectorXd cpp_get_thresholds_theta2vec(Eigen::Map<Eigen::VectorXd> THETA, const unsigned int P, const unsigned int C);
RcppExport SEXP _plFA_cpp_get_thresholds_theta2vec(SEXP THETASEXP, SEXP PSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type P(PSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_thresholds_theta2vec(THETA, P, C));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_loadings_mat2vec
Eigen::VectorXd cpp_get_loadings_mat2vec(Eigen::Map<Eigen::MatrixXd> LOADINGS, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, const int NLOAD);
RcppExport SEXP _plFA_cpp_get_loadings_mat2vec(SEXP LOADINGSSEXP, SEXP CONSTRMATSEXP, SEXP NLOADSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type LOADINGS(LOADINGSSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_loadings_mat2vec(LOADINGS, CONSTRMAT, NLOAD));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_loadings_vec2mat
Eigen::MatrixXd cpp_get_loadings_vec2mat(Eigen::Map<Eigen::VectorXd> LOADINGS, Eigen::Map<Eigen::MatrixXd> CONSTRMAT);
RcppExport SEXP _plFA_cpp_get_loadings_vec2mat(SEXP LOADINGSSEXP, SEXP CONSTRMATSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type LOADINGS(LOADINGSSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_loadings_vec2mat(LOADINGS, CONSTRMAT));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_latvar_mat2vec
Eigen::VectorXd cpp_get_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S);
RcppExport SEXP _plFA_cpp_get_latvar_mat2vec(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_latvar_mat2vec(S));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_latvar_vec2mat
Eigen::MatrixXd cpp_get_latvar_vec2mat(Eigen::Map<Eigen::VectorXd> SVEC, const unsigned int Q);
RcppExport SEXP _plFA_cpp_get_latvar_vec2mat(SEXP SVECSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type SVEC(SVECSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_latvar_vec2mat(SVEC, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_loadings_theta2mat
Eigen::MatrixXd cpp_get_loadings_theta2mat(Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, const int P, const int Q, const int D, const int C);
RcppExport SEXP _plFA_cpp_get_loadings_theta2mat(SEXP THETASEXP, SEXP CONSTRMATSEXP, SEXP PSEXP, SEXP QSEXP, SEXP DSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< const int >::type P(PSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_loadings_theta2mat(THETA, CONSTRMAT, P, Q, D, C));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_latvar_theta2mat
Eigen::MatrixXd cpp_get_latvar_theta2mat(Eigen::Map<Eigen::VectorXd> THETA, const int Q, const int D);
RcppExport SEXP _plFA_cpp_get_latvar_theta2mat(SEXP THETASEXP, SEXP QSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_latvar_theta2mat(THETA, Q, D));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_latvar_theta2vec
Eigen::VectorXd cpp_get_latvar_theta2vec(Eigen::Map<Eigen::VectorXd> THETA, const unsigned int NTHR, const unsigned int NLOAD, const unsigned int NCORR);
RcppExport SEXP _plFA_cpp_get_latvar_theta2vec(SEXP THETASEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NCORR(NCORRSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_latvar_theta2vec(THETA, NTHR, NLOAD, NCORR));
    return rcpp_result_gen;
END_RCPP
}
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
// cpp_multiThread_completePairwise
Rcpp::List cpp_multiThread_completePairwise(const unsigned int N, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> FREQ, const unsigned int CORRFLAG, const unsigned int GRFLAG, const unsigned int SILENTFLAG);
RcppExport SEXP _plFA_cpp_multiThread_completePairwise(SEXP NSEXP, SEXP C_VECSEXP, SEXP CONSTRMATSEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP CORRFLAGSEXP, SEXP GRFLAGSEXP, SEXP SILENTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type GRFLAG(GRFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SILENTFLAG(SILENTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_multiThread_completePairwise(N, C_VEC, CONSTRMAT, THETA, FREQ, CORRFLAG, GRFLAG, SILENTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// multiThread_completePairwise_nll
double multiThread_completePairwise_nll(const unsigned int N, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::VectorXd THETA, Eigen::Map<Eigen::MatrixXd> FREQ, const unsigned int CORRFLAG, const unsigned int SILENTFLAG);
RcppExport SEXP _plFA_multiThread_completePairwise_nll(SEXP NSEXP, SEXP C_VECSEXP, SEXP CONSTRMATSEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP CORRFLAGSEXP, SEXP SILENTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SILENTFLAG(SILENTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(multiThread_completePairwise_nll(N, C_VEC, CONSTRMAT, THETA, FREQ, CORRFLAG, SILENTFLAG));
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
// estimate_H
Rcpp::List estimate_H(Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> FREQ, int N, int CORRFLAG);
RcppExport SEXP _plFA_estimate_H(SEXP C_VECSEXP, SEXP ASEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP NSEXP, SEXP CORRFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_H(C_VEC, A, THETA, FREQ, N, CORRFLAG));
    return rcpp_result_gen;
END_RCPP
}
// estimate_J
Rcpp::List estimate_J(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::VectorXd& THETA, int CORRFLAG);
RcppExport SEXP _plFA_estimate_J(SEXP YSEXP, SEXP C_VECSEXP, SEXP ASEXP, SEXP THETASEXP, SEXP CORRFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_J(Y, C_VEC, A, THETA, CORRFLAG));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plFA_cpp_get_Lam", (DL_FUNC) &_plFA_cpp_get_Lam, 3},
    {"_plFA_cpp_get_S", (DL_FUNC) &_plFA_cpp_get_S, 2},
    {"_plFA_cpp_get_par_from_S", (DL_FUNC) &_plFA_cpp_get_par_from_S, 1},
    {"_plFA_cpp_compute_pair", (DL_FUNC) &_plFA_cpp_compute_pair, 9},
    {"_plFA_cpp_get_thresholds_theta2vec", (DL_FUNC) &_plFA_cpp_get_thresholds_theta2vec, 3},
    {"_plFA_cpp_get_loadings_mat2vec", (DL_FUNC) &_plFA_cpp_get_loadings_mat2vec, 3},
    {"_plFA_cpp_get_loadings_vec2mat", (DL_FUNC) &_plFA_cpp_get_loadings_vec2mat, 2},
    {"_plFA_cpp_get_latvar_mat2vec", (DL_FUNC) &_plFA_cpp_get_latvar_mat2vec, 1},
    {"_plFA_cpp_get_latvar_vec2mat", (DL_FUNC) &_plFA_cpp_get_latvar_vec2mat, 2},
    {"_plFA_cpp_get_loadings_theta2mat", (DL_FUNC) &_plFA_cpp_get_loadings_theta2mat, 6},
    {"_plFA_cpp_get_latvar_theta2mat", (DL_FUNC) &_plFA_cpp_get_latvar_theta2mat, 3},
    {"_plFA_cpp_get_latvar_theta2vec", (DL_FUNC) &_plFA_cpp_get_latvar_theta2vec, 4},
    {"_plFA_pairs_freq", (DL_FUNC) &_plFA_pairs_freq, 2},
    {"_plFA_cpp_multiThread_completePairwise", (DL_FUNC) &_plFA_cpp_multiThread_completePairwise, 8},
    {"_plFA_multiThread_completePairwise_nll", (DL_FUNC) &_plFA_multiThread_completePairwise_nll, 7},
    {"_plFA_sampling_step", (DL_FUNC) &_plFA_sampling_step, 8},
    {"_plFA_estimate_H", (DL_FUNC) &_plFA_estimate_H, 6},
    {"_plFA_estimate_J", (DL_FUNC) &_plFA_estimate_J, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_plFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
