#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc_leapfrog.h"
#include "mcmc_utils.h"
using namespace Rcpp;



/**
 * Function: find_reasonable_initial_step_size
 *
 * Generic implementation of Algorithm 4 from:
 * Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn Sampler: Adaptively Setting
 * Path Lengths in Hamiltonian Monte Carlo. JMLR, 15, 1593–1623.
 *
 * Finds a step size that yields an acceptance probability close to the target
 * by simulating a single leapfrog step and adjusting step size heuristically.
 *
 * Inputs:
 *  - theta: Current parameter vector.
 *  - log_post: Function to compute log posterior.
 *  - grad: Function to compute gradient of log posterior.
 *  - target_acceptance: Target acceptance rate (default 0.65).
 *  - init_step: Initial step size to try (default 1.0).
 *  - max_attempts: Max number of doubling/halving attempts (default 20).
 *
 * Returns:
 *  - A step size epsilon resulting in log acceptance near -log(2).
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    double target_acceptance,
    double init_step,
    int max_attempts
) {
  double eps = init_step;
  arma::vec r = arma::randn(theta.n_elem);

  double logp0 = log_post(theta);
  double kin0 = 0.5 * arma::dot(r, r);

  // One leapfrog step
  arma::vec theta_new, r_new;
  std::tie(theta_new, r_new) = leapfrog(theta, r, eps, grad);

  double logp1 = log_post(theta_new);
  double kin1 = 0.5 * arma::dot(r_new, r_new);

  double H0 = logp0 - kin0;
  double H1 = logp1 - kin1;

  int direction = 2 * (H1 - H0 > std::log(0.5)) - 1;  // +1 or -1

  int attempts = 0;
  while (direction * (H1 - H0) > -direction * std::log(2.0) && attempts < max_attempts) {
    eps = (direction == 1) ? 2.0 * eps : 0.5 * eps;

    // One leapfrog step
    std::tie(theta_new, r_new) = leapfrog(theta, r, eps, grad);

    // Evaluate Hamiltonian
    logp1 = log_post(theta_new);
    kin1 = 0.5 * arma::dot(r_new, r_new);
    H1 = logp1 - kin1;

    attempts++;
  }

  return eps;
}
