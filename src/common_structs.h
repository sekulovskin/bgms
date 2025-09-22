#pragma once
#include <RcppArmadillo.h>
#include <string>

// Parameters updated during sampling
struct ModelParameters {
  arma::mat main_effects;
  arma::mat pairwise_effects;
  arma::imat inclusion_indicator;
};

// Fixed model structure
struct ModelStructure {
  arma::ivec num_categories;
  arma::uvec is_ordinal_variable;
  arma::ivec baseline_category;
  arma::imat pairwise_effect_indices;
  arma::imat main_effect_indices; // unused in bgm, used in bgmCompare
};

// Priors
struct PriorSettings {
  double main_alpha;
  double main_beta;
  double pairwise_scale;
  double difference_scale;        // 0.0 if unused
  std::string edge_prior;         // empty if unused
  std::string difference_prior;   // empty if unused
  arma::mat inclusion_probability;
  double beta_bernoulli_alpha;
  double beta_bernoulli_beta;
  double dirichlet_alpha;
  double lambda;
};

// Sampler configuration
struct SamplerSettings {
  int iter;
  int burnin;
  double target_accept;
  int hmc_num_leapfrogs;
  int nuts_max_depth;
  bool learn_mass_matrix;
  std::string update_method;
  bool na_impute;
  bool edge_selection;
  bool difference_selection;
};
