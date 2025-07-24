#include "mcmc_nuts.h"
#include "mcmc_memoization.h"
#include <Rcpp.h>
using namespace Rcpp;



/**
 * Function: kinetic_energy
 *
 * Computes the kinetic energy of a momentum vector r, assuming a standard multivariate normal distribution.
 * This is used as part of the Hamiltonian energy in Hamiltonian Monte Carlo.
 *
 * Inputs:
 *  - r: The momentum vector.
 *
 * Returns:
 *  - The scalar kinetic energy value (0.5 * r^T * r).
 */
double kinetic_energy(const arma::vec& r) {
  return 0.5 * arma::dot(r, r);
}



/**
 * Function: leapfrog
 *
 * Performs a single leapfrog integration step using the standard gradient function.
 * Used to simulate Hamiltonian dynamics in HMC/NUTS.
 *
 * Inputs:
 *  - theta: Current position (parameter vector).
 *  - r: Current momentum vector.
 *  - eps: Step size for integration.
 *  - grad: Gradient function of the log posterior.
 *
 * Returns:
 *  - A pair containing:
 *      - Updated position vector.
 *      - Updated momentum vector.
 */
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  auto grad1 = grad(theta_new);
  r_half += 0.5 * eps * grad1;
  theta_new += eps * r_half;
  auto grad2 = grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}



/**
 * Function: leapfrog_memo
 *
 * Performs a single leapfrog step using memoization for the gradient evaluations.
 * This improves computational efficiency when the same positions are revisited.
 *
 * Inputs:
 *  - theta: Current position (parameter vector).
 *  - r: Current momentum vector.
 *  - eps: Step size for integration.
 *  - memo: Memoizer object caching gradient evaluations.
 *
 * Returns:
 *  - A pair containing:
 *      - Updated position vector.
 *      - Updated momentum vector.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  auto grad1 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad1;
  theta_new += eps * r_half;
  auto grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}



/**
 * Function: is_uturn
 *
 * Determines whether a U-turn has occurred in the trajectory of the NUTS algorithm.
 * This check ensures that the trajectory does not start to double back on itself.
 *
 * Inputs:
 *  - theta_min: Leftmost position in the trajectory.
 *  - theta_plus: Rightmost position in the trajectory.
 *  - r_min: Momentum at theta_min.
 *  - r_plus: Momentum at theta_plus.
 *
 * Returns:
 *  - true if a U-turn is detected; false otherwise.
 */
bool is_uturn(const arma::vec& theta_min,
              const arma::vec& theta_plus,
              const arma::vec& r_min,
              const arma::vec& r_plus) {
  arma::vec delta = theta_plus - theta_min;
  return arma::dot(delta, r_min) < 0 || arma::dot(delta, r_plus) < 0;
}



/**
 * Function: build_tree
 *
 * Recursively builds a binary tree of leapfrog steps in the NUTS algorithm.
 * This method explores forward or backward in time, evaluating trajectory termination criteria.
 *
 * This recursive tree-building procedure is based on Algorithm 6 in:
 *   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo.
 *   Journal of Machine Learning Research, 15(1), 1593–1623.
 *
 * Inputs:
 *  - theta: Current position at the base of the tree.
 *  - r: Current momentum at the base of the tree.
 *  - log_u: Log slice variable for accept/reject decision.
 *  - v: Direction of expansion (-1 for backward, +1 for forward).
 *  - j: Current tree depth.
 *  - step_size: Step size used in leapfrog integration.
 *  - theta_0: Initial position at the start of sampling.
 *  - r0: Initial momentum at the start of sampling.
 *  - logp0: Log posterior at initial position.
 *  - kin0: Kinetic energy at initial momentum.
 *  - memo: Memoizer object used for caching evaluations.
 *
 * Returns:
 *  - A BuildTreeResult struct containing updated position/momentum endpoints, candidate sample,
 *    subtree size and U-turn status, and average acceptance probability.
 */
BuildTreeResult build_tree(
    const arma::vec& theta,
    const arma::vec& r,
    double log_u,
    int v,
    int j,
    double step_size,
    const arma::vec& theta_0,
    const arma::vec& r0,
    const double logp0,
    const double kin0,
    Memoizer& memo
) {
  constexpr double Delta_max = 1000.0;
  if (j == 0) {

    arma::vec theta_new, r_new;
    std::tie(theta_new, r_new) = leapfrog_memo(theta, r, v * step_size, memo);

    auto logp = memo.cached_log_post(theta_new);
    double kin = kinetic_energy(r_new);
    int n_new = 1 * (log_u <= logp - kin);
    int s_new = 1 * (log_u <= Delta_max + logp - kin);
    double alpha = std::min(1.0, std::exp(logp - kin - logp0 + kin0));

    return {
      theta_new, r_new, theta_new, r_new, theta_new, n_new, s_new, alpha, 1
    };
  } else {
    BuildTreeResult result = build_tree(
      theta, r, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0, memo
    );

    arma::vec theta_min = result.theta_min;
    arma::vec r_min = result.r_min;
    arma::vec theta_plus = result.theta_plus;
    arma::vec r_plus = result.r_plus;
    arma::vec theta_prime = result.theta_prime;
    int n_prime = result.n_prime;
    int s_prime = result.s_prime;
    double alpha_prime = result.alpha;
    int n_alpha_prime = result.n_alpha;

    if (s_prime == 1) {
      if (v == -1) {
        result = build_tree(
          theta_min, r_min, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0, memo
        );
        theta_min = result.theta_min;
        r_min = result.r_min;
      } else {
        result = build_tree(
          theta_plus, r_plus, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0, memo
        );
        theta_plus = result.theta_plus;
        r_plus = result.r_plus;
      }

      arma::vec theta_double_prime = result.theta_prime;
      int n_double_prime = result.n_prime;
      int s_double_prime = result.s_prime;
      double alpha_double_prime = result.alpha;
      int n_alpha_double_prime = result.n_alpha;

      double denom = static_cast<double>(n_prime + n_double_prime);
      double prob = static_cast<double>(n_double_prime) / denom;

      if (R::unif_rand() < prob) {
        theta_prime = theta_double_prime;
      }
      alpha_prime += alpha_double_prime;
      n_alpha_prime += n_alpha_double_prime;
      s_prime = s_double_prime * !is_uturn(theta_min, theta_plus, r_min, r_plus);
      n_prime += n_double_prime;
    }

    return {
      theta_min, r_min, theta_plus, r_plus, theta_prime, n_prime, s_prime, alpha_prime, n_alpha_prime
    };
  }
}



/**
 * Function: nuts_sampler
 *
 * Runs the No-U-Turn Sampler (NUTS) using Hamiltonian Monte Carlo with dynamic path length.
 * Utilizes the recursive build_tree function to explore posterior space efficiently.
 *
 * This implementation is based on Algorithm 6 in:
 *   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo.
 *   Journal of Machine Learning Research, 15(1), 1593–1623.
 *
 * Inputs:
 *  - init_theta: Initial position (parameter vector).
 *  - step_size: Step size for leapfrog integration.
 *  - log_post: Log posterior function.
 *  - grad: Gradient of log posterior function.
 *  - max_depth: Maximum tree depth allowed for NUTS expansion (default = 10).
 *
 * Returns:
 *  - A SamplerResult struct containing:
 *      - Final sampled position.
 *      - Mean Metropolis acceptance probability.
 *      - Total number of proposals considered.
 */
SamplerResult nuts_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    int max_depth
) {
  // Here memo is created locally; terminates at end of nuts_sampler() call
  Memoizer memo(log_post, grad);

  arma::vec r0 = arma::randn(init_theta.n_elem);
  auto logp0 = memo.cached_log_post(init_theta);
  double kin0 = kinetic_energy(r0);
  double joint0 = logp0 - kin0;
  double log_u = log(R::unif_rand()) + joint0;
  arma::vec theta_min = init_theta, r_min = r0;
  arma::vec theta_plus = init_theta, r_plus = r0;
  arma::vec theta = init_theta;
  int j = 0;
  int n = 1, s = 1;

  double alpha = 0.5;
  int n_alpha = 1;

  while (s == 1 && j < max_depth) {
    int v = R::unif_rand() < 0.5 ? -1 : 1;

    BuildTreeResult result;
    if (v == -1) {
      result = build_tree(theta_min, r_min, log_u, v, j, step_size, init_theta, r0, logp0, kin0, memo);
      theta_min = result.theta_min;
      r_min = result.r_min;
    } else {
      result = build_tree(theta_plus, r_plus, log_u, v, j, step_size, init_theta, r0, logp0, kin0, memo);
      theta_plus = result.theta_plus;
      r_plus = result.r_plus;
    }

    alpha = result.alpha;
    n_alpha = result.n_alpha;

    if (result.s_prime == 1) {
      double prob = static_cast<double>(result.n_prime) / static_cast<double>(n);
      if (R::unif_rand() < prob) {
        theta = result.theta_prime;
      }
    }

    n += result.n_prime;
    s = result.s_prime == 1 ? !is_uturn(theta_min, theta_plus, r_min, r_plus) : 0;
    j++;
  }

  return {theta, alpha, n_alpha};
}