#include "nuts_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

static int leapfrog_counter = 0;

// Kinetic energy
double kinetic_energy(const arma::vec& r) {
  return 0.5 * arma::dot(r, r);
}

// Leapfrog integrator
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& logp_and_grad
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  auto [_, grad1] = logp_and_grad(theta_new);
  r_half += 0.5 * eps * grad1;
  theta_new += eps * r_half;
  auto [__, grad2] = logp_and_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}

// U-turn check
bool is_uturn(const arma::vec& theta_min,
              const arma::vec& theta_plus,
              const arma::vec& r_min,
              const arma::vec& r_plus) {
  arma::vec delta = theta_plus - theta_min;
  return arma::dot(delta, r_min) < 0 || arma::dot(delta, r_plus) < 0;
}


// Algorithm 6b in Hoffman & Gelman (2014)
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
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& logp_and_grad
) {
  constexpr double Delta_max = 1000.0;                                          // ~p. 1601
  if (j == 0) {
    leapfrog_counter++;

    // Base case: perform a single leapfrog step in direction v
    arma::vec theta_new, r_new;
    std::tie(theta_new, r_new) = leapfrog(theta, r, v * step_size, logp_and_grad);

    auto [logp, _] = logp_and_grad(theta_new);
    double kin = kinetic_energy(r_new);
    int n_new = 1 * (log_u <= logp - kin);
    int s_new = 1 * (log_u <= Delta_max + logp - kin);
    double alpha = std::min(1.0, std::exp(logp - kin - logp0 + kin0));

    return {
      theta_new, r_new, theta_new, r_new, theta_new, n_new, s_new, alpha, 1
    };
  } else {


    // Recursive case: build left and right subtrees
    BuildTreeResult result = build_tree (
      theta, r, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0,
      logp_and_grad
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

    if(s_prime == 1) {
      if(v == -1) {
        result = build_tree (
          theta_min, r_min, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0,
          logp_and_grad
        );

        theta_min = result.theta_min;
        r_min = result.r_min;
      } else {
        result = build_tree (
          theta_plus, r_plus, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0,
          logp_and_grad
        );
        theta_plus = result.theta_plus;
        r_plus = result.r_plus;
      }

      arma::vec theta_double_prime = result.theta_prime;
      int n_double_prime = result.n_prime;
      int s_double_prime = result.s_prime;
      double alpha_double_prime = result.alpha;
      int n_alpha_double_prime = result.n_alpha;

      double denom = static_cast<double>(n_prime) + static_cast<double>(n_double_prime);
      double prob = static_cast<double>(n_double_prime) / denom;

      if(R::unif_rand() < prob) {
        theta_prime = theta_double_prime;
      }
      alpha_prime += alpha_double_prime;
      n_alpha_prime += n_alpha_double_prime;
      s_prime = s_double_prime * !is_uturn(theta_min, theta_plus, r_min, r_plus);
      n_prime += n_double_prime;

      if (s_prime!=1) {
        Rcpp::Rcout << "[build_tree] U-turn detected at depth j = " << j << std::endl;
      }
    }

    return {
      theta_min, r_min, theta_plus, r_plus, theta_prime, n_prime, s_prime, alpha_prime, n_alpha_prime
    };
  }
}

// Algorithm 6a in Hoffman & Gelman (2014)
SamplerResult nuts_sampler(
    const arma::vec& init_theta,
    double step_size,
    std::function<std::pair<double, arma::vec>(const arma::vec&)> logp_and_grad,
    int max_depth
) {

  arma::vec r0 = arma::randn(init_theta.n_elem);
  auto [logp0, _] = logp_and_grad(init_theta);
  double kin0 = kinetic_energy(r0);
  double joint0 = logp0 - kin0;
  double log_u = log(R::unif_rand()) + joint0;
  arma::vec theta_min = init_theta, r_min = r0;
  arma::vec theta_plus = init_theta, r_plus = r0;
  arma::vec theta = init_theta;
  int j = 0;
  int n = 1, s = 1; // n: number of samples, s: success flag

  double alpha = 0.5;
  int n_alpha = 1;

  while (s == 1 && j < max_depth) {
    int v = R::unif_rand() < 0.5 ? -1 : 1;

    BuildTreeResult result;

    if (v == -1) {
      // Sample from the left subtree
      result = build_tree (
        theta_min, r_min, log_u, v, j, step_size, init_theta, r0, logp0, kin0,
        logp_and_grad
      );
      theta_min = result.theta_min;
      r_min = result.r_min;
    } else {
      // Sample from the right subtree
      result = build_tree (
        theta_plus, r_plus, log_u, v, j, step_size, init_theta, r0, logp0, kin0,
        logp_and_grad
      );
      theta_plus = result.theta_plus;
      r_plus = result.r_plus;
    }

    alpha = result.alpha;
    n_alpha = result.n_alpha;

    if(result.s_prime == 1) {
      double prob = static_cast<double>(result.n_prime) / static_cast<double>(n);
      if(R::unif_rand() < prob) {
        theta = result.theta_prime;
      }
    }

    n += result.n_prime;
    s = result.s_prime * !is_uturn(theta_min, theta_plus, r_min, r_plus);
    j++;
  }

  Rcpp::Rcout << "[nuts_sampler] alpha = " << alpha / n_alpha << std::endl;
  Rcpp::Rcout << "[nuts_sampler] total leapfrog steps = " << leapfrog_counter
              << ", max tree depth = " << j << std::endl;
  leapfrog_counter = 0;  // reset for next iteration

  return {
    theta,
    alpha,
    n_alpha
  };
}


