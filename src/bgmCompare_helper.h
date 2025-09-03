#pragma once

#include <RcppArmadillo.h>



arma::vec compute_group_main_effects(
    const int variable,
    const int num_groups,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::vec& proj_group
);

double compute_group_pairwise_effects(
    const int var1,
    const int var2,
    const int num_groups,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::vec& proj_group
);

arma::vec vectorize_model_parameters(
    const arma::mat& main_effects,                 // [n_main_rows × G]
    const arma::mat& pairwise_effects,             // [n_pair_rows × G]
    const arma::imat& inclusion_indicator,         // [V × V]
    const arma::imat& main_effect_indices,         // [V × 2], inclusive [start,end]
    const arma::imat& pairwise_effect_indices,      // [V × V]
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

void unvectorize_model_parameters(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,                 // [n_main_rows × G]
    arma::mat& pairwise_effects_out,             // [n_pair_rows × G]
    const arma::imat& inclusion_indicator,       // [V × V]
    const arma::imat& main_effect_indices,       // [V × 2], inclusive [start,end]
    const arma::imat& pairwise_effect_indices,   // [V × V]
    const int num_groups,                          // G
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

std::pair<arma::imat, arma::imat> build_index_maps(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

arma::vec inv_mass_active(
    const arma::vec& inv_diag,
    const arma::imat& inclusion_indicator,
    const int num_groups,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const bool& selection
);