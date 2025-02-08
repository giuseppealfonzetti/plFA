// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_compute_pair_ext
Rcpp::List cpp_compute_pair_ext(Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::VectorXd THETA, const int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const unsigned int K, const unsigned int L, Eigen::MatrixXd PAIRS_TABLE, const unsigned int SILENTFLAG, const unsigned int GRADFLAG, const int OPTION);
RcppExport SEXP _plFA_cpp_compute_pair_ext(SEXP CONSTRMATSEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP C_VECSEXP, SEXP THETASEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP KSEXP, SEXP LSEXP, SEXP PAIRS_TABLESEXP, SEXP SILENTFLAGSEXP, SEXP GRADFLAGSEXP, SEXP OPTIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type PAIRS_TABLE(PAIRS_TABLESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SILENTFLAG(SILENTFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type GRADFLAG(GRADFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type OPTION(OPTIONSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_compute_pair_ext(CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, THETA, CORRFLAG, NTHR, NLOAD, NCORR, NVAR, K, L, PAIRS_TABLE, SILENTFLAG, GRADFLAG, OPTION));
    return rcpp_result_gen;
END_RCPP
}
// cpp_loadings_theta2vec
Eigen::VectorXd cpp_loadings_theta2vec(Eigen::Map<Eigen::VectorXd> THETA, const int NTHR, const int NLOAD);
RcppExport SEXP _plFA_cpp_loadings_theta2vec(SEXP THETASEXP, SEXP NTHRSEXP, SEXP NLOADSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_loadings_theta2vec(THETA, NTHR, NLOAD));
    return rcpp_result_gen;
END_RCPP
}
// cpp_loadings_theta2mat
Eigen::MatrixXd cpp_loadings_theta2mat(Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, const std::vector<std::vector<std::vector<double>>> LLC, const int NTHR, const int NLOAD);
RcppExport SEXP _plFA_cpp_loadings_theta2mat(SEXP THETASEXP, SEXP CONSTRMATSEXP, SEXP LLCSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_loadings_theta2mat(THETA, CONSTRMAT, LLC, NTHR, NLOAD));
    return rcpp_result_gen;
END_RCPP
}
// cpp_loadings_mat2vec
Eigen::MatrixXd cpp_loadings_mat2vec(Eigen::Map<Eigen::MatrixXd> LOADINGS, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, const int NLOAD);
RcppExport SEXP _plFA_cpp_loadings_mat2vec(SEXP LOADINGSSEXP, SEXP CONSTRMATSEXP, SEXP NLOADSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type LOADINGS(LOADINGSSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_loadings_mat2vec(LOADINGS, CONSTRMAT, NLOAD));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_vec2cmat
Eigen::MatrixXd cpp_latvar_vec2cmat(Eigen::Map<Eigen::VectorXd> VEC, const int NCORR, const int Q);
RcppExport SEXP _plFA_cpp_latvar_vec2cmat(SEXP VECSEXP, SEXP NCORRSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type VEC(VECSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_vec2cmat(VEC, NCORR, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_vec2dmat
Eigen::MatrixXd cpp_latvar_vec2dmat(Eigen::Map<Eigen::VectorXd> VEC, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const int NCORR, const int NVAR, const int Q);
RcppExport SEXP _plFA_cpp_latvar_vec2dmat(SEXP VECSEXP, SEXP CONSTRLOGSDSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type VEC(VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_vec2dmat(VEC, CONSTRLOGSD, NCORR, NVAR, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_theta2cmat
Eigen::MatrixXd cpp_latvar_theta2cmat(Eigen::Map<Eigen::VectorXd> THETA, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const int Q);
RcppExport SEXP _plFA_cpp_latvar_theta2cmat(SEXP THETASEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_theta2dmat
Eigen::MatrixXd cpp_latvar_theta2dmat(Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const int Q);
RcppExport SEXP _plFA_cpp_latvar_theta2dmat(SEXP THETASEXP, SEXP CONSTRLOGSDSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_theta2mat
Eigen::MatrixXd cpp_latvar_theta2mat(Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const int Q);
RcppExport SEXP _plFA_cpp_latvar_theta2mat(SEXP THETASEXP, SEXP CONSTRLOGSDSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_theta2mat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_mat2vec
Eigen::VectorXd cpp_latvar_mat2vec(Eigen::Map<Eigen::MatrixXd> S, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const int NCORR, const int NVAR);
RcppExport SEXP _plFA_cpp_latvar_mat2vec(SEXP SSEXP, SEXP CONSTRLOGSDSEXP, SEXP NCORRSEXP, SEXP NVARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type S(SSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_mat2vec(S, CONSTRLOGSD, NCORR, NVAR));
    return rcpp_result_gen;
END_RCPP
}
// cpp_latvar_mat2cmat
Eigen::MatrixXd cpp_latvar_mat2cmat(Eigen::Map<Eigen::MatrixXd> S);
RcppExport SEXP _plFA_cpp_latvar_mat2cmat(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_latvar_mat2cmat(S));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sa_proj
Eigen::VectorXd cpp_sa_proj(Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> C_VEC, const int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR);
RcppExport SEXP _plFA_cpp_sa_proj(SEXP THETASEXP, SEXP CONSTRMATSEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP C_VECSEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< const int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sa_proj(THETA, CONSTRMAT, CONSTRLOGSD, LLC, C_VEC, CORRFLAG, NTHR, NLOAD, NCORR, NVAR));
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
Rcpp::List cpp_multiThread_completePairwise(const unsigned int N, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> FREQ, const int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const int GRFLAG, const int SILENTFLAG);
RcppExport SEXP _plFA_cpp_multiThread_completePairwise(SEXP NSEXP, SEXP C_VECSEXP, SEXP CONSTRMATSEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP GRFLAGSEXP, SEXP SILENTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< const int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type GRFLAG(GRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type SILENTFLAG(SILENTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_multiThread_completePairwise(N, C_VEC, CONSTRMAT, CONSTRLOGSD, LLC, THETA, FREQ, CORRFLAG, NTHR, NLOAD, NCORR, NVAR, GRFLAG, SILENTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// cpp_plSA
Rcpp::List cpp_plSA(Eigen::Map<Eigen::MatrixXd> FREQ, Eigen::Map<Eigen::MatrixXd> VALFREQ, const int N, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> THETA_INIT, Eigen::Map<Eigen::VectorXd> DIH, const int NTHR, const int NLOAD, const int NCORR, const int NVAR, const int SAMPLER, const int PAIRS_PER_ITERATION, const int SCHEDULE, const double STEP0, const double STEP1, const double STEP2, const double STEP3, const int BURN, const int MAXT, const int TOL_WINDOW, const double TOL_NLL, const int CHECK_TOL, const int CHECK_WINDOW, const int PATH_WINDOW, const int CLOCK_WINDOW, const int SEED, const bool VERBOSE);
RcppExport SEXP _plFA_cpp_plSA(SEXP FREQSEXP, SEXP VALFREQSEXP, SEXP NSEXP, SEXP C_VECSEXP, SEXP CONSTRMATSEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP THETA_INITSEXP, SEXP DIHSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP, SEXP SAMPLERSEXP, SEXP PAIRS_PER_ITERATIONSEXP, SEXP SCHEDULESEXP, SEXP STEP0SEXP, SEXP STEP1SEXP, SEXP STEP2SEXP, SEXP STEP3SEXP, SEXP BURNSEXP, SEXP MAXTSEXP, SEXP TOL_WINDOWSEXP, SEXP TOL_NLLSEXP, SEXP CHECK_TOLSEXP, SEXP CHECK_WINDOWSEXP, SEXP PATH_WINDOWSEXP, SEXP CLOCK_WINDOWSEXP, SEXP SEEDSEXP, SEXP VERBOSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type VALFREQ(VALFREQSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA_INIT(THETA_INITSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type DIH(DIHSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    Rcpp::traits::input_parameter< const int >::type SAMPLER(SAMPLERSEXP);
    Rcpp::traits::input_parameter< const int >::type PAIRS_PER_ITERATION(PAIRS_PER_ITERATIONSEXP);
    Rcpp::traits::input_parameter< const int >::type SCHEDULE(SCHEDULESEXP);
    Rcpp::traits::input_parameter< const double >::type STEP0(STEP0SEXP);
    Rcpp::traits::input_parameter< const double >::type STEP1(STEP1SEXP);
    Rcpp::traits::input_parameter< const double >::type STEP2(STEP2SEXP);
    Rcpp::traits::input_parameter< const double >::type STEP3(STEP3SEXP);
    Rcpp::traits::input_parameter< const int >::type BURN(BURNSEXP);
    Rcpp::traits::input_parameter< const int >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const int >::type TOL_WINDOW(TOL_WINDOWSEXP);
    Rcpp::traits::input_parameter< const double >::type TOL_NLL(TOL_NLLSEXP);
    Rcpp::traits::input_parameter< const int >::type CHECK_TOL(CHECK_TOLSEXP);
    Rcpp::traits::input_parameter< const int >::type CHECK_WINDOW(CHECK_WINDOWSEXP);
    Rcpp::traits::input_parameter< const int >::type PATH_WINDOW(PATH_WINDOWSEXP);
    Rcpp::traits::input_parameter< const int >::type CLOCK_WINDOW(CLOCK_WINDOWSEXP);
    Rcpp::traits::input_parameter< const int >::type SEED(SEEDSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSE(VERBOSESEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_plSA(FREQ, VALFREQ, N, C_VEC, CONSTRMAT, CONSTRLOGSD, LLC, THETA_INIT, DIH, NTHR, NLOAD, NCORR, NVAR, SAMPLER, PAIRS_PER_ITERATION, SCHEDULE, STEP0, STEP1, STEP2, STEP3, BURN, MAXT, TOL_WINDOW, TOL_NLL, CHECK_TOL, CHECK_WINDOW, PATH_WINDOW, CLOCK_WINDOW, SEED, VERBOSE));
    return rcpp_result_gen;
END_RCPP
}
// estimate_H
Rcpp::List estimate_H(Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> FREQ, int N, int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR);
RcppExport SEXP _plFA_estimate_H(SEXP C_VECSEXP, SEXP ASEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP NSEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_H(C_VEC, A, CONSTRLOGSD, LLC, THETA, FREQ, N, CORRFLAG, NTHR, NLOAD, NCORR, NVAR));
    return rcpp_result_gen;
END_RCPP
}
// estimate_J
Rcpp::List estimate_J(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::VectorXd& THETA, int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR);
RcppExport SEXP _plFA_estimate_J(SEXP YSEXP, SEXP C_VECSEXP, SEXP ASEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP THETASEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_J(Y, C_VEC, A, CONSTRLOGSD, LLC, THETA, CORRFLAG, NTHR, NLOAD, NCORR, NVAR));
    return rcpp_result_gen;
END_RCPP
}
// DH
Eigen::VectorXd DH(Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> CONSTRLOGSD, const std::vector<std::vector<std::vector<double>>> LLC, Eigen::Map<Eigen::VectorXd> THETA, Eigen::Map<Eigen::MatrixXd> FREQ, int N, int CORRFLAG, const int NTHR, const int NLOAD, const int NCORR, const int NVAR);
RcppExport SEXP _plFA_DH(SEXP C_VECSEXP, SEXP ASEXP, SEXP CONSTRLOGSDSEXP, SEXP LLCSEXP, SEXP THETASEXP, SEXP FREQSEXP, SEXP NSEXP, SEXP CORRFLAGSEXP, SEXP NTHRSEXP, SEXP NLOADSEXP, SEXP NCORRSEXP, SEXP NVARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type CONSTRLOGSD(CONSTRLOGSDSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::vector<double>>> >::type LLC(LLCSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const int >::type NTHR(NTHRSEXP);
    Rcpp::traits::input_parameter< const int >::type NLOAD(NLOADSEXP);
    Rcpp::traits::input_parameter< const int >::type NCORR(NCORRSEXP);
    Rcpp::traits::input_parameter< const int >::type NVAR(NVARSEXP);
    rcpp_result_gen = Rcpp::wrap(DH(C_VEC, A, CONSTRLOGSD, LLC, THETA, FREQ, N, CORRFLAG, NTHR, NLOAD, NCORR, NVAR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plFA_cpp_compute_pair_ext", (DL_FUNC) &_plFA_cpp_compute_pair_ext, 16},
    {"_plFA_cpp_loadings_theta2vec", (DL_FUNC) &_plFA_cpp_loadings_theta2vec, 3},
    {"_plFA_cpp_loadings_theta2mat", (DL_FUNC) &_plFA_cpp_loadings_theta2mat, 5},
    {"_plFA_cpp_loadings_mat2vec", (DL_FUNC) &_plFA_cpp_loadings_mat2vec, 3},
    {"_plFA_cpp_latvar_vec2cmat", (DL_FUNC) &_plFA_cpp_latvar_vec2cmat, 3},
    {"_plFA_cpp_latvar_vec2dmat", (DL_FUNC) &_plFA_cpp_latvar_vec2dmat, 5},
    {"_plFA_cpp_latvar_theta2cmat", (DL_FUNC) &_plFA_cpp_latvar_theta2cmat, 6},
    {"_plFA_cpp_latvar_theta2dmat", (DL_FUNC) &_plFA_cpp_latvar_theta2dmat, 7},
    {"_plFA_cpp_latvar_theta2mat", (DL_FUNC) &_plFA_cpp_latvar_theta2mat, 7},
    {"_plFA_cpp_latvar_mat2vec", (DL_FUNC) &_plFA_cpp_latvar_mat2vec, 4},
    {"_plFA_cpp_latvar_mat2cmat", (DL_FUNC) &_plFA_cpp_latvar_mat2cmat, 1},
    {"_plFA_cpp_sa_proj", (DL_FUNC) &_plFA_cpp_sa_proj, 10},
    {"_plFA_pairs_freq", (DL_FUNC) &_plFA_pairs_freq, 2},
    {"_plFA_cpp_multiThread_completePairwise", (DL_FUNC) &_plFA_cpp_multiThread_completePairwise, 14},
    {"_plFA_cpp_plSA", (DL_FUNC) &_plFA_cpp_plSA, 30},
    {"_plFA_estimate_H", (DL_FUNC) &_plFA_estimate_H, 12},
    {"_plFA_estimate_J", (DL_FUNC) &_plFA_estimate_J, 11},
    {"_plFA_DH", (DL_FUNC) &_plFA_DH, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_plFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
