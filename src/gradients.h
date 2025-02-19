#ifndef gradients_H
#define gradients_H
#include <RcppEigen.h>
#include "genericUtils.h"
#include "bivariateProbs.h"

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
  );

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
  );

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
  );

  Eigen::MatrixXd S(const Eigen::Ref<const Eigen::MatrixXd> A,
                    const Eigen::Ref<const Eigen::VectorXd> TRANSFORMED_RHOS,
                    const int Q,
                    const int IDX
  );

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
  );

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
  );

}
#endif
