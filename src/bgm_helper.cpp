#include <RcppArmadillo.h>
#include "bgm_helper.h"
using namespace Rcpp;

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
int count_num_main_effects (
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  int n_params = 0;
  for (int i = 0; i < num_categories.n_elem; i++) {
    n_params += is_ordinal_variable[i] ? num_categories[i] : 2;
  }
  return n_params;
}



/**
 * Function: vectorize_thresholds
 *
 * Converts a matrix of main effect parameters into a flat vector for optimization.
 * Respects the structure of ordinal vs. Blume-Capel variables.
 *
 * Inputs:
 *  - main_effects: Matrix of main effect parameters.
 *  - num_categories: Category count per variable.
 *  - is_ordinal_variable: Logical vector indicating ordinal (1) or Blume-Capel (0).
 *
 * Returns:
 *  - A flat vector of all parameters, in row-major order.
 *    Ordinal: one parameter per category;
 *    Blume-Capel: two parameters (linear and quadratic).
 */
arma::vec vectorize_thresholds (
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_parameters = count_num_main_effects (num_categories, is_ordinal_variable);
  arma::vec threshold_vector (num_parameters);
  int offset = 0;

  for (int variable = 0; variable < main_effects.n_rows; variable++) {
    const int num_pars = is_ordinal_variable (variable) ? num_categories(variable) : 2;
    threshold_vector.subvec (offset, offset + num_pars - 1) =
      main_effects.row (variable).cols (0, num_pars - 1).t ();
    offset += num_pars;
  }

  return threshold_vector;
}



/**
 * Function: unvectorize_thresholds
 *
 * Reconstructs a threshold matrix from a flat vector of parameters,
 * reversing the operation of `vectorize_thresholds()`.
 *
 * Inputs:
 *  - vector: Flattened vector of threshold parameters.
 *  - num_categories: Vector of category counts per variable.
 *  - is_ordinal_variable: Logical vector (1 for ordinal, 0 for Blume-Capel).
 *
 * Returns:
 *  - A matrix of main effect parameters, with each row corresponding to a variable.
 *    Entries beyond the actual number of parameters per variable are zero-filled.
 */
arma::mat unvectorize_thresholds (
    const arma::vec& threshold_vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int max_categories = num_categories.max ();

  arma::mat matrix (num_variables, max_categories, arma::fill::zeros);

  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_pars = is_ordinal_variable[variable] ? num_categories[variable] : 2;
    matrix.row (variable).cols (0, num_pars - 1) =
      threshold_vector.subvec (offset, offset + num_pars - 1).t ();
    offset += num_pars;
  }

  return matrix;
}



/**
 * Flattens threshold and interaction parameters into a single vector.
 *
 * This utility function creates a joint parameter vector used for
 * Fisher-preconditioned MALA updates. Threshold (main) effects are added
 * first in variable order, followed by interaction terms in upper-triangle order.
 * Interaction values are only included if marked as active in the inclusion indicator.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters [V × max_categories].
 *  - pairwise_effects: Symmetric matrix of interaction strengths [V × V].
 *  - inclusion_indicator: Symmetric binary matrix of active interactions [V × V].
 *  - num_categories: Vector of category counts per variable [V].
 *  - is_ordinal_variable: Binary indicator for ordinal variables [V].
 *
 * Returns:
 *  - Parameter vector (main effects + active interactions).
 */
arma::vec vectorize_model_parameters(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = main_effects.n_rows;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  arma::vec param_vec(num_main + num_interactions, arma::fill::zeros);
  int offset = 0;

  // Threshold parameters
  for (int v = 0; v < num_variables; ++v) {
    if (is_ordinal_variable(v)) {
      int num_cats = num_categories(v);
      for (int c = 0; c < num_cats; ++c) {
        param_vec(offset++) = main_effects(v, c);
      }
    } else {
      // Blume-Capel: always two params
      param_vec(offset++) = main_effects(v, 0);  // linear
      param_vec(offset++) = main_effects(v, 1);  // quadratic
    }
  }

  // Interaction parameters (upper triangle)
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      int idx = num_main + (v1 * num_variables + v2 - ((v1 + 2) * (v1 + 1)) / 2);
      if (inclusion_indicator(v1, v2) == 1) {
        param_vec(idx) = pairwise_effects(v1, v2);
      }
    }
  }

  return param_vec;
}



/**
 * Reconstructs threshold and interaction matrices from joint parameter vector.
 *
 * This function unpacks a flattened parameter vector into two matrices:
 * - main_effects: Threshold parameters [V × max_categories]
 * - pairwise_effects: Symmetric matrix of interaction strengths [V × V]
 *
 * Thresholds are unpacked first based on variable type (ordinal vs. Blume-Capel),
 * followed by all upper-triangle interactions (assumed symmetric).
 *
 * Inputs:
 *  - param_vec: Flattened parameter vector (main + interactions).
 *  - num_categories: Vector of number of categories per variable [V].
 *  - is_ordinal_variable: Indicator for ordinal (vs. Blume-Capel) variables [V].
 *
 * Outputs:
 *  - main_effects_out: Matrix of thresholds (overwritten in-place).
 *  - pairwise_effects_out: Symmetric matrix of pairwise effects (overwritten in-place).
 */
void unvectorize_model_parameters(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,
    arma::mat& pairwise_effects_out,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  const int max_num_categories = num_categories.max();

  main_effects_out.set_size(num_variables, max_num_categories);
  pairwise_effects_out.set_size(num_variables, num_variables);
  main_effects_out.fill(0.0);
  pairwise_effects_out.fill(0.0);

  int offset = 0;

  // Thresholds
  for (int v = 0; v < num_variables; ++v) {
    if (is_ordinal_variable(v)) {
      int num_cats = num_categories(v);
      for (int c = 0; c < num_cats; ++c) {
        main_effects_out(v, c) = param_vec(offset++);
      }
    } else {
      main_effects_out(v, 0) = param_vec(offset++);
      main_effects_out(v, 1) = param_vec(offset++);
    }
  }

  // Interactions (upper triangle)
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      int idx = num_main + (v1 * num_variables + v2 - ((v1 + 2) * (v1 + 1)) / 2);
      double val = param_vec(idx);
      pairwise_effects_out(v1, v2) = val;
      pairwise_effects_out(v2, v1) = val;
    }
  }
}