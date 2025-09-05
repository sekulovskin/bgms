#pragma once

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <dqrng_generator.h>
#include <xoshiro.h>
#include <random>

// ============================================================
// Scalar RNG helpers
// ============================================================

// Uniform(0,1)
inline double runif(dqrng::xoshiro256plus& rng) {
  static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rng);
}

// Normal(mu, sigma)
inline double rnorm(dqrng::xoshiro256plus& rng,
                    double mu = 0.0, double sigma = 1.0) {
  static thread_local std::normal_distribution<double> dist(0.0, 1.0);
  return mu + sigma * dist(rng);
}

// Bernoulli(p)
inline int rbern(dqrng::xoshiro256plus& rng, double p) {
  return (runif(rng) < p) ? 1 : 0;
}

// Beta(a, b)
inline double rbeta(dqrng::xoshiro256plus& rng, double a, double b) {
  std::gamma_distribution<double> g1(a, 1.0);
  std::gamma_distribution<double> g2(b, 1.0);
  double x = g1(rng);
  double y = g2(rng);
  return x / (x + y);
}

// Exponential(lambda)
inline double rexp(dqrng::xoshiro256plus& rng, double lambda) {
  std::exponential_distribution<double> dist(lambda);
  return dist(rng);
}

// ============================================================
// Armadillo RNG helpers
// ============================================================

// arma::vec of N(mu, sigma)
inline arma::vec arma_rnorm_vec(dqrng::xoshiro256plus& rng,
                                arma::uword n,
                                double mu = 0.0, double sigma = 1.0) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = rnorm(rng, mu, sigma);
  return out;
}

// arma::mat of N(mu, sigma)
inline arma::mat arma_rnorm_mat(dqrng::xoshiro256plus& rng,
                                arma::uword nrow, arma::uword ncol,
                                double mu = 0.0, double sigma = 1.0) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = rnorm(rng, mu, sigma);
    }
  }
  return out;
}

// arma::vec of U(0,1)
inline arma::vec arma_runif_vec(dqrng::xoshiro256plus& rng, arma::uword n) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = runif(rng);
  return out;
}

// arma::mat of U(0,1)
inline arma::mat arma_runif_mat(dqrng::xoshiro256plus& rng,
                                arma::uword nrow, arma::uword ncol) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = runif(rng);
    }
  }
  return out;
}

// Random permutation (like arma::randperm)
inline arma::uvec arma_randperm(dqrng::xoshiro256plus& rng, arma::uword n) {
  arma::uvec out(n);
  std::iota(out.begin(), out.end(), 0);
  std::shuffle(out.begin(), out.end(), rng);
  return out;
}
