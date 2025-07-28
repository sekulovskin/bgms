// mcmc_adaptation.h
#pragma once

#include <armadillo>
#include <cmath>
#include <algorithm>
#include <vector>
#include "mcmc_utils.h"

class DualAveraging {
public:
  double log_step_size;
  double log_step_size_avg;
  double hbar;
  double mu;
  double gamma;
  double t0;
  double kappa;
  int t;

  DualAveraging(double initial_step_size)
    : log_step_size(std::log(initial_step_size)),
      log_step_size_avg(std::log(initial_step_size)),
      hbar(0.0),
      mu(std::log(10.0 * initial_step_size)),
      gamma(0.05),
      t0(10.0),
      kappa(0.75),
      t(1) {}

  void update(double accept_prob, double target_accept) {
    double eta = 1.0 / (t + t0);
    double error = target_accept - accept_prob;
    hbar = (1 - eta) * hbar + eta * error;
    log_step_size = mu - std::sqrt(t) / gamma * hbar;

    double weight = std::pow(t, -kappa);
    log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
    t++;
  }

  void restart(double new_step_size) {
    log_step_size = std::log(new_step_size);
    log_step_size_avg = std::log(new_step_size);
    mu = std::log(10.0 * new_step_size);
    hbar = 0.0;
    t = 1;
  }

  double current() const { return std::exp(log_step_size); }
  double averaged() const { return std::exp(log_step_size_avg); }
};



// === Diagonal Mass Matrix Estimator ===
class DiagMassMatrixAccumulator {
public:
  int count;
  arma::vec mean;
  arma::vec m2;

  DiagMassMatrixAccumulator(int dim)
    : count(0), mean(arma::zeros(dim)), m2(arma::zeros(dim)) {}

  void update(const arma::vec& sample) {
    count++;
    arma::vec delta = sample - mean;
    mean += delta / count;
    arma::vec delta2 = sample - mean;
    m2 += delta % delta2;
  }

  arma::vec variance() const {
    static constexpr double prior_weight = 5.0;
    static constexpr double prior_variance = 1.0;//1e-3;
    double n = static_cast<double>(count);

    arma::vec empirical = m2 / std::max(1.0, n - 1.0);
    arma::vec prior = arma::ones(empirical.n_elem) * prior_variance;
    arma::vec var = (n / (n + prior_weight)) * empirical
    + (prior_weight / (n + prior_weight)) * prior;
    return var;
  }

  arma::vec inverse_mass() const {
    return 1.0 / variance();
  }

  void reset() {
    count = 0;
    mean.zeros();
    m2.zeros();
  }
};


// === Stan-style Dynamic Warmup Schedule with Adaptive Windows ===
struct WarmupSchedule {
  int stage1_end;
  std::vector<int> window_ends;
  int stage3_start;
  int total_burnin;
  bool selection_enabled_globally;


  WarmupSchedule(int burnin, bool enable_selection)
    : total_burnin(burnin),
      selection_enabled_globally(enable_selection)
  {
    // 15% init, 10% term
    stage1_end  = int(0.15 * burnin);
    stage3_start= burnin - int(0.10 * burnin);

    // carve phase II into doubling windows
    int current  = stage1_end;
    int win_size = 25;  // Stan’s default base window
    while (current < stage3_start) {
      // clamp window so we don’t overshoot phase II
      int this_win = std::min(win_size, stage3_start - current);
      window_ends.push_back(current + this_win);
      current     += this_win;
      win_size    = std::min(win_size * 2, stage3_start - current);
    }
  }

  int current_window(int iteration) const {
    for (size_t i = 0; i < window_ends.size(); ++i) {
      if (iteration < window_ends[i]) return static_cast<int>(i);
    }
    return -1;
  }

  bool in_stage1(int iter) const { return iter < stage1_end; }
  bool in_stage2(int iter) const { return iter >= stage1_end && iter < stage3_start; }
  bool in_stage3(int iter) const { return iter >= stage3_start && iter < total_burnin; }

  bool selection_enabled(int iter) const {
    if (!selection_enabled_globally) return false;
    int win = current_window(iter);
    return window_ends.size() > 3 && win >= 3;
  }
};



// === HMC/NUTS Adaptation Controller ===
class HMCAdaptationController {
public:
  HMCAdaptationController(int dim,
                          double initial_step_size,
                          double target_accept,
                          WarmupSchedule& schedule_ref,
                          bool learn_mass_matrix = true)
    : schedule(schedule_ref),
      learn_mass_matrix_(learn_mass_matrix),
      mass_accumulator(dim),
      step_adapter(initial_step_size),
      inv_mass_(arma::ones<arma::vec>(dim)),
      step_size_(initial_step_size),
      target_accept_(target_accept),
      finalized_mass_(false) {}

  void update(const arma::vec& theta,
              double accept_prob,
              int iteration /* zero‐based: iteration=0 on first warmup draw */) {
    // 1) Step‐size adaptation on EVERY warmup iteration (phases I–III)
    if (iteration + 1 <= schedule.total_burnin) {
      step_adapter.update(accept_prob, target_accept_);
      step_size_ = step_adapter.current();
    }

    // 2) Mass‐matrix adaptation ONLY in phase II (slow windows)
    if (schedule.in_stage2(iteration) && learn_mass_matrix_) {
      // accumulate θ
      mass_accumulator.update(theta);

      // if we’ve just reached the end of the current window…
      int w = schedule.current_window(iteration);
      if (iteration + 1 == schedule.window_ends[w]) {
        // compute shrunk sample variance → diag mass matrix
        arma::vec var    = mass_accumulator.variance();
        inv_mass_        = 1.0 / var;

        // reset for next window
        mass_accumulator.reset();

        // restart dual‐averaging under the new metric
        step_adapter.restart(step_size_);
      }
    }

    // 3) At the last warmup iteration, freeze ε to its running average
    if (iteration + 1 == schedule.total_burnin) {
      step_size_ = step_adapter.averaged();
    }
  }

  double current_step_size() const { return step_size_; }
  double final_step_size() const { return step_adapter.averaged(); }
  const arma::vec& inv_mass_diag() const { return inv_mass_; }
  bool has_fixed_mass_matrix() const { return finalized_mass_; }

private:
  WarmupSchedule& schedule;
  bool learn_mass_matrix_;
  DiagMassMatrixAccumulator mass_accumulator;
  DualAveraging step_adapter;
  arma::vec inv_mass_;
  double step_size_;
  double target_accept_;
  bool finalized_mass_;
};



// RWMAdaptationController
// Handles Robbins-Monro adaptation for random walk Metropolis proposals
// Automatically adjusts step sizes during warmup using acceptance probabilities
class RWMAdaptationController {
public:
  arma::mat& proposal_sd;
  const int total_burnin;
  const double target_accept;

  RWMAdaptationController(arma::mat& proposal_sd_matrix,
                          const WarmupSchedule& warmup,
                          double target_accept_rate = 0.44)
    : proposal_sd(proposal_sd_matrix),
      total_burnin(warmup.total_burnin),
      target_accept(target_accept_rate) {}

  void update(const arma::umat& index_mask,
              const arma::mat& accept_prob_matrix,
              int iteration) {

    if (iteration >= total_burnin || iteration < 1)
      return;

    const double rm_decay_rate = 0.75;
    const double rm_weight = std::pow(iteration, -rm_decay_rate);

    for (arma::uword i = 0; i < proposal_sd.n_rows; ++i) {
      for (arma::uword j = 0; j < proposal_sd.n_cols; ++j) {
        if (index_mask(i, j) == 1) {
          const double accept_log = std::log(accept_prob_matrix(i, j));
          proposal_sd(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd(i, j),
            accept_log,
            rm_weight,
            target_accept
          );
        }
      }
    }
  }
};