#ifndef MCMC_NUTS_H
#define MCMC_NUTS_H

#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc_memoization.h"

/**
 * Struct: BuildTreeResult
 *
 * Holds the return values of the recursive build_tree function used in the NUTS algorithm.
 * Each call to build_tree expands the sampling path and may return new candidate samples
 * or indicate when the trajectory should terminate.
 *
 * Fields:
 *  - theta_min: Leftmost position in the trajectory.
 *  - r_min: Corresponding momentum at theta_min.
 *  - theta_plus: Rightmost position in the trajectory.
 *  - r_plus: Corresponding momentum at theta_plus.
 *  - theta_prime: The current proposed sample (to possibly accept).
 *  - n_prime: Number of valid proposals from this subtree.
 *  - s_prime: Stop flag (1 = continue, 0 = stop expansion).
 *  - alpha: Sum of acceptance probabilities in the subtree.
 *  - n_alpha: Number of proposals contributing to alpha.
 */
struct BuildTreeResult {
  arma::vec theta_min;
  arma::vec r_min;
  arma::vec theta_plus;
  arma::vec r_plus;
  arma::vec theta_prime;
  int n_prime;
  int s_prime;
  double alpha;
  int n_alpha;
};



/**
 * Struct: SamplerResult
 *
 * Represents the final outcome of one iteration of the NUTS sampler.
 *
 * Fields:
 *  - theta: Final accepted position (parameter vector).
 *  - alpha: Average Metropolis acceptance probability across proposals.
 *  - n_alpha: Number of proposal steps contributing to alpha.
 */
struct SamplerResult {
  arma::vec theta;
  double alpha;
  int n_alpha;
};



/**
 * Function: leapfrog
 *
 * Performs a standard leapfrog integration step using the given gradient function.
 * This is a core component of Hamiltonian Monte Carlo.
 */
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad
);



/**
 * Function: leapfrog_memo
 *
 * Performs a leapfrog step using a memoization wrapper to avoid redundant gradient evaluations.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo
);


/**
 * Function: nuts_sampler
 *
 * Executes the No-U-Turn Sampler algorithm (NUTS).
 */
SamplerResult nuts_sampler(const arma::vec& init_theta,
                           double step_size,
                           const std::function<double(const arma::vec&)>& log_post,
                           const std::function<arma::vec(const arma::vec&)>& grad,
                           int max_depth = 10);

#endif