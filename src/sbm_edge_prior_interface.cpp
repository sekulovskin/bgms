#include <RcppArmadillo.h>
#include "math/explog_switch.h"

// ----------------------------------------------------------------------------|
// The c++ code below is based on the R code accompanying the paper:
//  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
//  Detection With Unknown Number of Communities, Journal of the American
//  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
arma::vec compute_Vn_mfm_sbm(arma::uword no_variables,
                             double dirichlet_alpha,
                             arma::uword t_max,
                             double lambda) {
  arma::vec log_Vn(t_max);
  double r;

  // Negative binomial parameters (overdispersed alternative to truncated Poisson)
  double nb_size = 100;
  double nb_prob = nb_size / (nb_size + lambda);

  for(arma::uword t = 0; t < t_max; t++) {
    r = -INFINITY;
    for(arma::uword k = t; k <= 500; k++){
      arma::vec b_linspace_1 = arma::linspace(k-t+1,k+1,t+1);
      arma::vec b_linspace_2 = arma::linspace((k+1)*dirichlet_alpha,(k+1)*dirichlet_alpha+no_variables-1, no_variables);

      double b = arma::accu(ARMA_MY_LOG(b_linspace_1)) -
        arma::accu(ARMA_MY_LOG(b_linspace_2)) +
        R::dnbinom((k+1)-1, nb_size, nb_prob, true);  // Match Poisson indexing

      double m = std::max(b,r);
      r = MY_LOG(MY_EXP(r-m) + MY_EXP(b-m)) + m;
    }
    log_Vn(t) = r;
  }
  return log_Vn;
}




