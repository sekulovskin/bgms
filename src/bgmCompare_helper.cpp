#include <RcppArmadillo.h>
#include <cmath>
#include "bgmCompare_helper.h"
#include "common_helpers.h"

using namespace Rcpp;



/**
 * Function: compute_group_main_effects
 * Purpose: Computes main_effects for a specific variable in a given group.
 * Inputs:
 *  - variable: Index of the variable for which main_effects are computed.
 *  - group: Index of the group for which main_effects are computed.
 *  - num_groups: Total number of groups in the analysis.
 *  - main_effects: Matrix of main effects across variables and groups.
 *  - main_effect_indices: Indices for main effect parameters.
 *  - projection: Projection matrix for group-specific scaling.
 *  - independent_main_effects: Whether main_effects are modeled independently.
 * Outputs:
 *  - A `arma::vec` containing the computed main_effects for the variable in the group.
 */
arma::vec compute_group_main_effects(
    const int variable,
    const int num_groups,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::vec& proj_group
) {
  // Base index for accessing main effects for this variable
  int base_category_index = main_effect_indices(variable, 0);
  int last_category_index = main_effect_indices(variable, 1);

  arma::vec Groupmain_effects =
    arma::conv_to<arma::vec>::from(
      main_effects.rows(base_category_index, last_category_index).col(0));
  Groupmain_effects += main_effects.rows(
    base_category_index, last_category_index
  ).cols(
      1, num_groups - 1
  ) *
    proj_group;
  return Groupmain_effects;
}



// Computes the pairwise effect between var1 and var2 for a given group.
//
// pairwise_effects: rows index pairs; columns are [shared | group-contrasts...]
// pairwise_effect_indices: maps (var1, var2) -> row index in pairwise_effects
// inclusion_indicator(var1,var2): 0 means use shared-only (no group contrast), nonzero includes contrast
//
// Returns the scalar interaction for that pair in that group.
double compute_group_pairwise_effects(
    const int var1,
    const int var2,
    const int num_groups,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::vec& proj_group
) {
  // Look up the row index for this pair; assume indices are provided for var1<var2
  // If your index matrix is symmetric, either order works; otherwise enforce (min,max).
  const int i = pairwise_effect_indices(var1, var2);

  // Shared/base effect
  double w = pairwise_effects(i, 0);

  // Optional group contrast
  if (inclusion_indicator(var1, var2) != 0) {
    // cols(1..G-1) * projection[group]^T  → scalar
    w += arma::as_scalar(pairwise_effects.row(i).cols(1, num_groups - 1) *
      proj_group);
  }
  return w;
}




/**
 * Flattens overall and (selected) difference parameters for bgmCompare.
 *
 * Order:
 *   1) MAIN overall main_effects (column 0), in variable order and within-variable row order.
 *   2) PAIRWISE overall interactions (column 0), in upper-triangle order (v1 < v2).
 *   3) MAIN group-difference main_effects (columns 1..G-1) for variables with inclusion_indicator(v,v) == 1,
 *      emitted in variable order, within-variable row order, and column-wise g = 1..G-1.
 *   4) PAIRWISE group-difference interactions (columns 1..G-1) for pairs with inclusion_indicator(v1,v2) == 1,
 *      emitted in upper-triangle order and column-wise g = 1..G-1.
 *
 * Returns:
 *   arma::vec of size:
 *     n_main_rows                              // all MAIN overall
 *   + n_pair_rows                              // all PAIRWISE overall
 *   + sum_v[ 1(v included) * rows_v * (num_groups-1) ]  // selected MAIN diffs
 *   + sum_{v1<v2}[ 1(pair included) * (num_groups-1) ]  // selected PAIR diffs
 */
arma::vec vectorize_model_parameters(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_groups = main_effects.n_cols;
  const int n_main_rows = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int n_pair_rows = num_variables * (num_variables - 1) / 2;

  // total length
  int total_len = 0;
  total_len += n_main_rows;      // MAIN overall (col 0)
  total_len += n_pair_rows;      // PAIR overall (col 0)

  // MAIN differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 1) {
      const int r0 = main_effect_indices(v, 0);
      const int r1 = main_effect_indices(v, 1);
      total_len += static_cast<int>(r1 - r0 + 1) * (num_groups - 1);
    }
  }
  // PAIRWISE differences
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) total_len += (num_groups - 1);
    }
  }

  arma::vec param_vec(total_len, arma::fill::zeros);
  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  param_vec.subvec(off, off + n_main_rows - 1) = main_effects.col(0);
  off += n_main_rows;

  // 2) PAIRWISE overall (col 0) — vectorized
  // (Relies on rows being in the same upper-triangle order as constructed in R.)
  param_vec.subvec(off, off + n_pair_rows - 1) = pairwise_effects.col(0);
  off += n_pair_rows;

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) != 1) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        param_vec(off++) = main_effects(r, g);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) != 1) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        param_vec(off++) = pairwise_effects(row, g);
      }
    }
  }

  return param_vec;
}



/**
 * Reconstructs bgmCompare matrices from a flattened parameter vector.
 *
 * Shapes are derived from the indices and num_groups (G). Overall effects
 * (column 0) are always filled. Difference columns (1..G-1) are filled only
 * where inclusion_indicator == 1 (diag for MAIN, off-diagonal for PAIRWISE).
 * Not-selected difference entries remain zero.
 *
 * Inputs:
 *   - param_vec: flattened vector produced by vectorize_model_parameters_bgmCompare
 *   - inclusion_indicator: [V×V] symmetric selection matrix
 *   - main_effect_indices: [V×2] inclusive [start,end] row indices per variable
 *   - pairwise_effect_indices: [V×V] row index per pair (v1,v2), symmetric
 *   - num_groups: G (total columns = 1 overall + (G-1) differences)
 *
 * Outputs (overwritten):
 *   - main_effects_out: [n_main_rows × G]
 *   - pairwise_effects_out: [n_pair_rows × G]
 */
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
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_main = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  main_effects_out.set_size(num_main, num_groups);
  pairwise_effects_out.set_size(num_pair, num_groups);
  main_effects_out.zeros();
  pairwise_effects_out.zeros();

  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  main_effects_out.col(0) = param_vec.subvec(off, off + num_main - 1);
  off += num_main;

  // 2) PAIRWISE overall (col 0) — vectorized
  pairwise_effects_out.col(0) = param_vec.subvec(off, off + num_pair - 1);
  off += num_pair;

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; v++) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        main_effects_out(r, g) = param_vec(off++);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 0) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        pairwise_effects_out(row, g) = param_vec(off++);
      }
    }
  }
}



// Build index maps for param_vec layout.
// Returns two imats with the same dims as main_effects and pairwise_effects.
// Entry (i,j) = position in param_vec, or -1 if that entry is never stored.
std::pair<arma::imat, arma::imat> build_index_maps(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_groups    = main_effects.n_cols;
  const int num_main = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  arma::imat main_index(num_main, num_groups, arma::fill::value(-1));
  arma::imat pair_index(num_main, num_groups, arma::fill::value(-1));

  int off = 0;

  // 1) main overall (col 0)
  for (int r = 0; r < num_main; ++r) {
    main_index(r,0) = off++;
  }

  // 2) pair overall (col 0)
  for (int r = 0; r < num_pair; ++r) {
    pair_index(r,0) = off++;
  }

  // 3) main differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v,v) == 0) continue;
    const int r0 = main_effect_indices(v,0);
    const int r1 = main_effect_indices(v,1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        main_index(r,g) = off++;
      }
    }
  }

  // 4) pairwise differences
  for (int v1 = 0; v1 < num_variables-1; ++v1) {
    for (int v2 = v1+1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1,v2) == 0) continue;
      const int row = pairwise_effect_indices(v1,v2);
      for (int g = 1; g < num_groups; ++g) {
        pair_index(row,g) = off++;
      }
    }
  }

  return {main_index, pair_index};
}



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
) {
  if(selection == false)
    return inv_diag;

  const int num_variables = inclusion_indicator.n_rows;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  const int num_pair = num_variables * (num_variables - 1) / 2;

  // total length
  int total_len = 0;
  total_len += num_main;      // MAIN overall (col 0)
  total_len += num_pair;      // PAIR overall (col 0)

  // MAIN differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    total_len += static_cast<int>(r1 - r0 + 1) * (num_groups - 1);
  }

  // PAIRWISE differences
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) total_len += (num_groups - 1);
    }
  }

  arma::vec active_inv_diag(total_len, arma::fill::zeros);

  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  active_inv_diag.subvec(off, off + num_main - 1) = inv_diag.subvec(off, off + num_main - 1);
  off += num_main;

  // 2) PAIRWISE overall (col 0) — vectorized
  // (Relies on rows being in the same upper-triangle order as constructed in R.)
  active_inv_diag.subvec(off, off + num_pair - 1) = inv_diag.subvec(off, off + num_pair - 1);
  off += num_pair;

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        int idx = main_index(r, g);
        active_inv_diag(off++) = inv_diag(idx);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) != 1) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        int idx = pair_index(row, g);
        active_inv_diag(off++) = inv_diag(idx);
      }
    }
  }

  return active_inv_diag;
}