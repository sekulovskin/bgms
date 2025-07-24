#ifndef NUTS_FUNCTIONS_H
#define NUTS_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <functional>

// Struct to hold build_tree return value
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

// Struct to return sampler output
struct SamplerResult {
  arma::vec theta;
  double alpha;
  int n_alpha;
};

std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad
);

// NUTS sampler interface
SamplerResult nuts_sampler(const arma::vec& init_theta,
                           double step_size,
                           const std::function<double(const arma::vec&)>& log_post,
                           const std::function<arma::vec(const arma::vec&)>& grad,
                           int max_depth = 10);

#endif