#include "mcmc_nuts.h"
#include "mcmc_memoization.h"
#include <Rcpp.h>
using namespace Rcpp;


// /**
//  * Struct: VecHash
//  *
//  * A hash function is a way of converting complex input data (like a vector of numbers)
//  * into a single fixed-size integer value. This hash acts like a digital fingerprint —
//  * it helps uniquely identify the input and allows fast lookup in data structures
//  * like hash tables.
//  *
//  * In C++, hash tables are implemented using containers like std::unordered_map,
//  * which require a hash function and an equality comparator for the key type.
//  * Since arma::vec (a vector of doubles) isn't natively supported as a key,
//  * this struct defines how to compute a hash for it.
//  *
//  * Why use hashing here?
//  *   - We cache (memoize) values like log-posterior or gradient at specific points.
//  *   - To retrieve them quickly later, we need a fast way to index by vector.
//  *   - Hashing gives constant-time access (on average) by turning the vector into a hash code.
//  *
//  * This implementation:
//  *   - Initializes a seed based on the vector size.
//  *   - Iteratively incorporates each element's hash into the seed using bitwise operations
//  *     and a constant based on the golden ratio to reduce collisions.
//  *
//  * The result is a unique and repeatable hash value that lets us efficiently
//  * store and retrieve cached evaluations.
//  *
//  * Provides a custom hash function for arma::vec so that it can be used as a key
//  * in standard hash-based containers like std::unordered_map.
//  *
//  * Since arma::vec is not hashable by default in C++, this functor allows you
//  * to compute a combined hash from the individual elements of the vector.
//  *
//  * Hashing Strategy:
//  *   - Starts with an initial seed equal to the vector's length.
//  *   - Iteratively combines the hash of each element into the seed using
//  *     a standard hashing recipe (XOR, bit-shift mix, and golden ratio).
//  *
//  * This ensures that vectors with similar but not identical elements
//  * produce distinct hash codes, reducing the risk of collisions.
//  */
// struct VecHash {
//   std::size_t operator()(const arma::vec& v) const {
//     std::size_t seed = v.n_elem;
//     for (size_t i = 0; i < v.n_elem; i++) {
//       seed ^= std::hash<double>{}(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//     }
//     return seed;
//   }
// };
//
//
//
// /**
//  * Struct: VecEqual
//  *
//  * Provides a custom equality comparator for arma::vec, to be used in conjunction
//  * with VecHash in unordered_map and similar containers.
//  *
//  * By default, arma::vec does not implement equality suitable for use as a key
//  * in hashed containers. This comparator solves that by relying on Armadillo's
//  * approx_equal function.
//  *
//  * Comparison Strategy:
//  *   - Uses "absdiff" mode to check if the absolute difference between vectors
//  *     is below a specified tolerance.
//  *   - The tolerance is hardcoded as 1e-12, assuming numerical precision errors
//  *     should not cause false inequality.
//  *
//  * This is crucial for numerical stability, especially when the same vector
//  * might be computed multiple times with slight rounding differences.
//  */
// struct VecEqual {
//   bool operator()(const arma::vec& a, const arma::vec& b) const {
//     return arma::approx_equal(a, b, "absdiff", 1e-12);
//   }
// };
//
//
//
// /**
//  * Class: Memoizer
//  *
//  * A utility class that caches (memoizes) evaluations of the log-posterior and its gradient
//  * at specific parameter values (theta) to avoid redundant and costly recomputation.
//  *
//  * Usage:
//  *   - Constructed with a log-posterior function and its gradient.
//  *   - Uses hash-based caching (via unordered_map) keyed by arma::vec inputs.
//  *     That means:
//  *       - Each parameter vector (theta) is used as a key.
//  *       - Results are stored in a lookup table so that repeated evaluations can be skipped.
//  *       - A custom VecHash function generates a hash code for arma::vec.
//  *       - A custom VecEqual comparator defines when two vectors are considered equal (with a small tolerance).
//  *     This enables fast retrieval (average constant time) of previously computed values for log_post and grad.
//  *
//  * Methods:
//  *  - cached_log_post(const arma::vec& theta):
//  *        Returns the cached log-posterior if available, otherwise computes, stores, and returns it.
//  *
//  *  - cached_grad(const arma::vec& theta):
//  *        Returns the cached gradient if available, otherwise computes, stores, and returns it.
//  *
//  * Internal Details:
//  *  - Relies on VecHash and VecEqual to use arma::vec as map keys.
//  *  - Assumes high-precision comparison for theta keys (absdiff tolerance 1e-12).
//  */
// class Memoizer {
// public:
//   std::function<double(const arma::vec&)> log_post;
//   std::function<arma::vec(const arma::vec&)> grad;
//
//   std::unordered_map<arma::vec, double, VecHash, VecEqual> logp_cache;
//   std::unordered_map<arma::vec, arma::vec, VecHash, VecEqual> grad_cache;
//
//   Memoizer(
//     const std::function<double(const arma::vec&)>& lp,
//     const std::function<arma::vec(const arma::vec&)>& gr
//   ) : log_post(lp), grad(gr) {}
//
//   double cached_log_post(const arma::vec& theta) {
//     auto it = logp_cache.find(theta);
//     if (it != logp_cache.end()) return it->second;
//     double val = log_post(theta);
//     logp_cache[theta] = val;
//     return val;
//   }
//
//   arma::vec cached_grad(const arma::vec& theta) {
//     auto it = grad_cache.find(theta);
//     if (it != grad_cache.end()) return it->second;
//     arma::vec val = grad(theta);
//     grad_cache[theta] = val;
//     return val;
//   }
// };



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