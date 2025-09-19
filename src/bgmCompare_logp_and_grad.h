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
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>& sufficient_pairwise,
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
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>& sufficient_pairwise,
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

double log_pseudoposterior_main_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& num_obs_categories_group,
    const std::vector<arma::imat>& sufficient_blume_capel_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double difference_scale,
    int variable,
    int category, // for ordinal variables only
    int par, // for Blume-Capel variables only
    int h // Overall = 0, differences are 1, ....
);

double log_pseudoposterior_pair_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>&  sufficient_pairwise_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double interaction_scale,
    const double difference_scale,
    int variable1,
    int variable2,
    int h // Overall = 0, differences are 1, ....
);