#ifndef BGM_HELPER_H
#define BGM_HELPER_H

#include <RcppArmadillo.h>

// Count the number of main effect parameters
int count_num_main_effects(
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Vectorize threshold matrix
arma::vec vectorize_thresholds(
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Unvectorize thresholds into matrix
arma::mat unvectorize_thresholds(
    const arma::vec& threshold_vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Vectorize model parameters (main + interaction effects)
arma::vec vectorize_model_parameters(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Unvectorize model parameters back into matrices
void unvectorize_model_parameters(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,
    arma::mat& pairwise_effects_out,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

#endif // BGM_HELPER_H
