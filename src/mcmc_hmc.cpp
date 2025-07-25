#include <RcppArmadillo.h>
#include <functional>
#include "mcmc_hmc.h"
#include "mcmc_leapfrog.h"
#include "mcmc_utils.h"
using namespace Rcpp;



SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs
) {
  arma::vec theta = init_theta;
  arma::vec init_r = arma::randn(theta.n_elem);
  arma::vec r = init_r;

  std::tie(theta, r) = leapfrog(theta, r, step_size, grad, num_leapfrogs);

  // Hamiltonians
  double current_H = -log_post(init_theta) + 0.5 * arma::dot(init_r, init_r);
  double proposed_H = -log_post(theta) + 0.5 * arma::dot(r, r);
  double log_accept_prob = current_H - proposed_H;

  arma::vec state = (std::log(R::unif_rand()) < log_accept_prob) ? theta : init_theta;

  double accept_prob = std::min(1.0, std::exp(log_accept_prob));

  return {state, accept_prob};
}
