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

static const R_CallMethodDef CallEntries[] = {
    {"_plFA_pairs_freq", (DL_FUNC) &_plFA_pairs_freq, 2},
    {"_plFA_get_S", (DL_FUNC) &_plFA_get_S, 2},
    {"_plFA_get_par_from_S", (DL_FUNC) &_plFA_get_par_from_S, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_plFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}