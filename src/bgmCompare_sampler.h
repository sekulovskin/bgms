#pragma once
#include <dqrng.h>
#include <xoshiro.h>
#include <RcppArmadillo.h>
#include "sampler_output.h"

SamplerOutput run_gibbs_sampler_for_bgmCompare(
    int chain_id,
    arma::imat observations,
    int num_groups,
    std::vector<arma::imat> num_obs_categories,
    std::vector<arma::imat> sufficient_blume_capel,
    std::vector<arma::mat> sufficient_pairwise,
    const arma::ivec& num_categories,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    double difference_scale,
    double difference_selection_alpha,
    double difference_selection_beta,
    const std::string& difference_prior,
    int iter,
    int burnin,
    bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    arma::mat inclusion_probability,
    dqrng::xoshiro256plus& rng
);