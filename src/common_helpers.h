#pragma once
#include <RcppArmadillo.h>

/**
 * Function: count_num_main_effects
 *
 * Computes the total number of main effect (threshold) parameters across variables.
 *
 * Ordinal variables contribute one parameter per category.
 * Blume-Capel variables contribute exactly two parameters.
 *
 * Inputs:
 *  - num_categories: Vector of category counts per variable.
 *  - is_ordinal_variable: Logical vector (0 or 1) indicating which variables are ordinal.
 *
 * Returns:
 *  - Total number of main effect parameters to estimate.
 */
inline int count_num_main_effects(const arma::ivec& num_categories,
                                  const arma::uvec& is_ordinal_variable) {
  int n_params = 0;
  for (arma::uword i = 0; i < num_categories.n_elem; ++i) {
    n_params += is_ordinal_variable[i] ? num_categories[i] : 2;
  }
  return n_params;
}
