#ifndef bivariateProbs_H
#define bivariateProbs_H

#include "bivariateNormal.h"
#include "gradients.h"

namespace biprobs{
  /*  PI (the probability of a bivariate response pattern) */
  // Compute specific pi_sksl
  double compute_pi(
      const Eigen::VectorXd c_vec,
      const Eigen::VectorXd pi_thresholds,
      double rho_kl,

      const unsigned int k,
      const unsigned int l,
      const unsigned int sk,
      const unsigned int sl
  );

}




#endif
