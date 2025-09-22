#pragma once
#include <RcppArmadillo.h>



/**
 * Count the total number of main-effect parameters across all variables.
 *
 * Each variable contributes:
 *  - Ordinal variable: `num_categories[v]` parameters (one per category).
 *  - Blume–Capel variable: always 2 parameters (linear and quadratic terms).
 *
 * Inputs:
 *  - num_categories: Number of categories per variable [V].
 *  - is_ordinal_variable: Indicator of whether each variable is ordinal [V].
 *
 * Returns:
 *  - Total number of main-effect parameters (int).
 */
inline int count_num_main_effects(const arma::ivec& num_categories,
                                  const arma::uvec& is_ordinal_variable) {
  int n_params = 0;
  for (arma::uword i = 0; i < num_categories.n_elem; ++i) {
    n_params += is_ordinal_variable[i] ? num_categories[i] : 2;
  }
  return n_params;
}
