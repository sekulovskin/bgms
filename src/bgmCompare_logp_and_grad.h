#pragma once

#include <RcppArmadillo.h>



double log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const Rcpp::List& num_obs_categories_group,
    const Rcpp::List& sufficient_blume_capel_group,
    const Rcpp::List& sufficient_pairwise_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const double difference_scale
);

arma::vec gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const Rcpp::List& num_obs_categories_group,
    const Rcpp::List& sufficient_blume_capel_group,
    const Rcpp::List& sufficient_pairwise_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const double difference_scale,
    const arma::imat main_index,
    const arma::imat pair_index
);