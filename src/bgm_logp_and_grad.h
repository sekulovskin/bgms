#ifndef BGM_LOGP_AND_GRAD_H
#define BGM_LOGP_AND_GRAD_H

#include <RcppArmadillo.h>

// Log posterior for main effect parameters only
double log_pseudoposterior_thresholds(
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta
);

// Log posterior for a single component of main effect
double log_pseudoposterior_thresholds_component(
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const int variable,
    const int category,
    const int parameter
);

// Gradient for thresholds
arma::vec gradient_log_pseudoposterior_thresholds(
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta
);

// Gradient for all interactions
arma::vec gradient_log_pseudoposterior_interactions(
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
);

// Gradient for a single interaction parameter
double gradient_log_pseudoposterior_interaction_single(
    int var1,
    int var2,
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise
);

// Log posterior over all interactions
double log_pseudoposterior_interactions(
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise
);

// Pseudolikelihood ratio for a single variable
double compute_log_likelihood_ratio_for_variable(
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state,
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::mat& rest_matrix,
    const arma::imat& observations,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
);

// Pseudolikelihood ratio for an interaction
double log_pseudolikelihood_ratio_interaction(
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int num_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& rest_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const arma::imat& sufficient_pairwise
);

// Full gradient
arma::vec gradient_log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
);

// Full log posterior
double log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
);

#endif // BGM_LOGP_AND_GRAD_H
