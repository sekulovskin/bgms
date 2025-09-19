#include <RcppArmadillo.h>
#include "bgmCompare_helper.h"
#include "bgmCompare_logp_and_grad.h"
#include <cmath>

using namespace Rcpp;



double log_pseudoposterior(
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
    const std::vector<arma::mat>&  sufficient_pairwise_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const double difference_scale
) {
  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::imat num_obs_categories = num_obs_categories_group[group];
    const arma::imat sufficient_blume_capel = sufficient_blume_capel_group[group];

    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    for (int v = 0; v < num_variables; ++v) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );

      // store into row v (padded with zeros if variable has < max_num_categories params)
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      // upper triangle incl. base value; mirror to keep symmetry
      for (int u = v; u < num_variables; ++u) { // Combines with loop over v
        if(u == v) continue;
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }

    // ---- data contribution pseudolikelihood (linear terms) ----
      const int num_cats = num_categories(v);
      if (is_ordinal_variable(v)) {
        // use group-specific main_effects
        for (int c = 0; c < num_cats; ++c) {
          const double val = main_group(v, c);
          log_pp += static_cast<double>(num_obs_categories(c, v)) * val;
        }
      } else {
        // two sufficient stats for binary-ish coding? keep original shape
        log_pp += static_cast<double>(sufficient_blume_capel(0, v)) * main_group(v, 0);
        log_pp += static_cast<double>(sufficient_blume_capel(1, v)) * main_group(v, 1);
      }
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat sufficient_pairwise = sufficient_pairwise_group[group];

    log_pp += arma::accu(pairwise_group % sufficient_pairwise); // trace(X' * W * X) = sum(W %*% (X'X))

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const arma::mat rest_matrix = obs * pairwise_group;
    for (int v = 0; v < num_variables; ++v) {
      const int num_cats = num_categories(v);
      const arma::vec rest_score = rest_matrix.col(v);

      // bound to stabilize exp; use group-specific params consistently
      arma::vec bound = num_cats * rest_score;
      bound = arma::clamp(bound, 0.0, arma::datum::inf);

      arma::vec denom(rest_score.n_elem, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        // base term exp(-bound)
        denom = arma::exp(-bound);
        // main_effects from main_group
        for (int c = 0; c < num_cats; ++c) {
          const double th = main_group(v, c);
          const arma::vec exponent = th + (c + 1) * rest_score - bound;
          denom += arma::exp(exponent);
        }
      } else {
        // linear/quadratic main effects from main_group
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        const int ref = baseline_category(v);
        for (int c = 0; c <= num_cats; ++c) {
          const int centered = c - ref;
          const double quad = quad_effect * centered * centered;
          const double lin  = lin_effect * c;
          const arma::vec exponent = lin + quad + c * rest_score - bound;
          denom += arma::exp(exponent);
        }
      }
      // - sum_i [ bound_i + log denom_i ]
      log_pp -= arma::accu(bound + arma::log(denom));
    }
  }

  // ---- priors ----
  auto log_beta_prior = [&](double x) {
    return x * main_alpha - std::log1p(std::exp(x)) * (main_alpha + main_beta);
  };

  // Main effects prior
  for (int v = 0; v < num_variables; ++v) {
    const int row0 = main_effect_indices(v, 0);
    const int row1 = main_effect_indices(v, 1);
    for (int r = row0; r <= row1; ++r) {
      log_pp += log_beta_prior(main_effects(r, 0));

      if (inclusion_indicator(v, v) == 0) continue;
      for (int eff = 1; eff < num_groups; ++eff) {
        log_pp += R::dcauchy(main_effects(r, eff), 0.0, difference_scale, true);
      }
    }
  }

  // Pairwise effects prior
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      const int idx = pairwise_effect_indices(v1, v2);
      log_pp += R::dcauchy(pairwise_effects(idx, 0), 0.0, interaction_scale, true);

      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int eff = 1; eff < num_groups; ++eff) {
        log_pp += R::dcauchy(pairwise_effects(idx, eff), 0.0, difference_scale, true);
      }
    }
  }

  return log_pp;
}



arma::vec gradient(
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
    const std::vector<arma::mat>&  sufficient_pairwise_group,
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
) {
  const int num_variables      = observations.n_cols;
  const int max_num_categories = num_categories.max();
  const int n_main_rows        = main_effects.n_rows;
  const int n_pair_rows        = pairwise_effects.n_rows;

  // total length (must match vectorize_model_parameters)
  long long total_len = 0;
  total_len += n_main_rows;      // main col 0
  total_len += n_pair_rows;      // pair col 0
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 1) {
      const int r0 = main_effect_indices(v, 0);
      const int r1 = main_effect_indices(v, 1);
      total_len += static_cast<long long>(r1 - r0 + 1) * (num_groups - 1);
    }
  }
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) total_len += (num_groups - 1);
    }
  }

  arma::vec grad(total_len, arma::fill::zeros);
  int off;

  // -------------------------------
  // Observed sufficient statistics
  // -------------------------------
  for (int g = 0; g < num_groups; ++g) {
    // list access
    arma::imat num_obs_categories      = num_obs_categories_group[g];
    arma::imat sufficient_blume_capel  = sufficient_blume_capel_group[g];

    // Main effects
    for (int v = 0; v < num_variables; ++v) {
      const int base     = main_effect_indices(v, 0);
      const int num_cats = num_categories(v);

      if (is_ordinal_variable(v)) {
        for (int c = 0; c < num_cats; ++c) {
          // overall
          off = main_index(base + c, 0);
          grad(off) += num_obs_categories(c, v);

          // diffs
          if (inclusion_indicator(v, v) != 0) {
            for (int k = 1; k < num_groups; ++k) {
              off = main_index(base + c, k);
              grad(off) += num_obs_categories(c, v) * projection(g, k - 1);
            }
          }
        }
      } else {
        // overall (2 stats)
        off = main_index(base, 0);
        grad(off) += sufficient_blume_capel(0, v);

        off = main_index(base + 1, 0);
        grad(off) += sufficient_blume_capel(1, v);

        // diffs
        if (inclusion_indicator(v, v) != 0) {
          for (int k = 1; k < num_groups; ++k) {
            off = main_index(base, k);
            grad(off) += sufficient_blume_capel(0, v) * projection(g, k - 1);

            off = main_index(base + 1, k);
            grad(off) += sufficient_blume_capel(1, v) * projection(g, k - 1);
          }
        }
      }
    }

    // Pairwise (observed)
    arma::mat sufficient_pairwise = sufficient_pairwise_group[g];
    for (int v1 = 0; v1 < num_variables - 1; ++v1) {
      for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
        const int row = pairwise_effect_indices(v1, v2);

        off = pair_index(row, 0);
        grad(off) += 2.0 * sufficient_pairwise(v1, v2); // upper tri counted once

        if (inclusion_indicator(v1, v2) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = pair_index(row, k);
          grad(off) += 2.0 * sufficient_pairwise(v1, v2) * projection(g, k - 1);
        }
      }
    }
  }

  // --------------------------------
  // Expected sufficient statistics
  // --------------------------------
  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);

    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);
    const arma::vec proj_g = projection.row(g).t(); // length = num_groups-1

    // build group-specific params
    for (int v = 0; v < num_variables; ++v) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      for (int u = v; u < num_variables; ++u) { // Combines with loop over v
        if(u == v) continue;
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }

    // group slice and rest matrix
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat rest_matrix = obs * pairwise_group;
    const int num_group_obs = obs.n_rows;

    for (int v = 0; v < num_variables; ++v) {
      const int K   = num_categories(v);
      const int ref = baseline_category(v);

      arma::vec rest_score = rest_matrix.col(v);
      arma::vec bound      = K * rest_score;
      bound = arma::clamp(bound, 0.0, arma::datum::inf);

      arma::mat exponents(num_group_obs, K + 1, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        exponents.col(0) -= bound;
        arma::vec main_param = main_group.row(v).cols(0, K - 1).t();
        for (int j = 0; j < K; j++) {
          exponents.col(j+1) = main_param(j) + (j + 1) * rest_score - bound;
        }
      } else {
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        for (int s = 0; s <= K; ++s) {
          const int centered = s - ref;
          const double lin  = lin_effect * s;
          const double quad = quad_effect * centered * centered;
          exponents.col(s) = lin + quad + s * rest_score - bound;
        }
      }

      arma::mat probs = arma::exp(exponents);
      arma::vec denom = arma::sum(probs, 1); // base term
      probs.each_col() /= denom;

      // ---- MAIN expected ----
      const int base = main_effect_indices(v, 0);

      if (is_ordinal_variable(v)) {
        for (int s = 1; s <= K; ++s) {
          const int j = s - 1;
          off = main_index(base + j, 0);
          grad(off) -= arma::accu(probs.col(s));

          if (inclusion_indicator(v, v) == 0) continue;
          for (int k = 1; k < num_groups; ++k) {
            off = main_index(base + j, k);
            grad(off) -= projection(g, k - 1) * arma::accu(probs.col(s));
          }
        }
      } else {
        arma::vec lin_score  = arma::regspace<arma::vec>(0, K);          // length K+1
        arma::vec quad_score = arma::square(lin_score - ref);

        off = main_index(base, 0);
        grad(off) -= arma::accu(probs * lin_score);

        off = main_index(base + 1, 0);
        grad(off) -= arma::accu(probs * quad_score);

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = main_index(base, k);
          grad(off) -= projection(g, k - 1) * arma::accu(probs * lin_score);
          off = main_index(base + 1, k);
          grad(off) -= projection(g, k - 1) * arma::accu(probs * quad_score);
        }
      }

      // ---- PAIRWISE expected ----
      for (int v2 = 0; v2 < num_variables; ++v2) {
        if (v == v2) continue;

        const int row = (v < v2) ? pairwise_effect_indices(v, v2)
          : pairwise_effect_indices(v2, v);

        arma::vec expected_value(num_group_obs, arma::fill::zeros);

        for (int s = 1; s <= K; ++s) {
          expected_value += s * probs.col(s) % obs.col(v2);
        }

        off = pair_index(row, 0);
        grad(off) -= arma::accu(expected_value);

        if (inclusion_indicator(v, v2) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = pair_index(row, k);
          grad(off) -= projection(g, k - 1) * arma::accu(expected_value);
        }
      }
    }
  }

  // -------------------------------
  // Priors
  // -------------------------------
  // Main
  for (int v = 0; v < num_variables; ++v) {
    const int base     = main_effect_indices(v, 0);
    const int num_cats = num_categories(v);

    if (is_ordinal_variable(v)) {
      for (int c = 0; c < num_cats; ++c) {
        off = main_index(base + c, 0);
        double value = main_effects(base + c, 0);
        const double p = 1.0 / (1.0 + std::exp(-value));
        grad(off) += main_alpha - (main_alpha + main_beta) * p;

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = main_index(base + c, k);
          double value = main_effects(base + c, k);
          grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
        }

      }
    } else {
      off = main_index(base, 0);
      double value = main_effects(base, 0);
      double p = 1.0 / (1.0 + std::exp(-value));
      grad(off) += main_alpha - (main_alpha + main_beta) * p;

      off = main_index(base + 1, 0);
      value = main_effects(base + 1, 0);
      p = 1.0 / (1.0 + std::exp(-value));
      grad(off) += main_alpha - (main_alpha + main_beta) * p;


      if (inclusion_indicator(v, v) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = main_index(base, k);
        double value = main_effects(base, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);

        off = main_index(base + 1, k);
        value = main_effects(base + 1, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
      }
    }
  }

  // Pairwise (Cauchy)
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      const int row = pairwise_effect_indices(v1, v2);

      // overall uses interaction_scale
      off = pair_index(row, 0);
      double value = pairwise_effects(row, 0);
      grad(off) -= 2.0 * value / (value * value + interaction_scale * interaction_scale);


      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = pair_index(row, k);
        double value = pairwise_effects(row, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
      }
    }
  }

  return grad;
}



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
) {
  if(h > 0 && inclusion_indicator(variable, variable) == 0) {
    return 0.0; // No contribution if differences not included
  }

  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::imat num_obs_categories = num_obs_categories_group[group];
    const arma::imat sufficient_blume_capel = sufficient_blume_capel_group[group];

    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    arma::vec me = compute_group_main_effects(
      variable, num_groups, main_effects, main_effect_indices, proj_g
    );

    // store into row v (padded with zeros if variable has < max_num_categories params)
    main_group(variable, arma::span(0, me.n_elem - 1)) = me.t();

    // upper triangle incl. base value; mirror to keep symmetry
    for (int u = 0; u < num_variables; u++) {
      if(u == variable) continue;
      double w = compute_group_pairwise_effects(
        variable, u, num_groups, pairwise_effects, pairwise_effect_indices,
        inclusion_indicator, proj_g
      );
      pairwise_group(variable, u) = w;
      pairwise_group(u, variable) = w;
    }

    // ---- data contribution pseudolikelihood (linear terms) ----
    if (is_ordinal_variable(variable)) {
      const double val = main_group(variable, category);
      log_pp += static_cast<double>(num_obs_categories(category, variable)) *
        val;
    } else {
      log_pp += static_cast<double>(sufficient_blume_capel(par, variable)) *
        main_group(variable, par);
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const arma::vec rest_score = obs * pairwise_group.col(variable);
    const int num_cats = num_categories(variable);

    // bound to stabilize exp; use group-specific params consistently
    arma::vec bound = num_cats * rest_score;
    bound = arma::clamp(bound, 0.0, arma::datum::inf);

    arma::vec denom(rest_score.n_elem, arma::fill::zeros);
    if (is_ordinal_variable(variable)) {
      // base term exp(-bound)
      denom = arma::exp(-bound);
      // main_effects from main_group
      for (int cat = 0; cat < num_cats; cat++) {
        const double th = main_group(variable, cat);
        const arma::vec exponent = th + (cat + 1) * rest_score - bound;
        denom += arma::exp(exponent);
      }
    } else {
      // linear/quadratic main effects from main_group
      const double lin_effect  = main_group(variable, 0);
      const double quad_effect = main_group(variable, 1);
      const int ref = baseline_category(variable);
      for (int cat = 0; cat <= num_cats; cat++) {
        const int centered = cat - ref;
        const double quad = quad_effect * centered * centered;
        const double lin  = lin_effect * cat;
        const arma::vec exponent = lin + quad + cat * rest_score - bound;
        denom += arma::exp(exponent);
      }
    }
    // - sum_i [ bound_i + log denom_i ]
    log_pp -= arma::accu(bound + arma::log(denom));
  }

  // ---- priors ----
  if (h == 0) {
    // overall
    auto log_beta_prior = [&](double x) {
      return x * main_alpha - std::log1p(std::exp(x)) * (main_alpha + main_beta);
    };

    // Main effects prior
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += log_beta_prior(main_effects(r, 0));
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += log_beta_prior(main_effects(r, 0));
    }
  } else {
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += R::dcauchy(main_effects(r, h), 0.0, difference_scale, true);
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += R::dcauchy(main_effects(r, h), 0.0, difference_scale, true);
    }
  }

  return log_pp;
}



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
) {
  if(h > 0 && inclusion_indicator(variable1, variable2) == 0) {
    return 0.0; // No contribution if differences not included
  }

  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;
  int idx = pairwise_effect_indices(variable1, variable2);

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    for (int v : {variable1, variable2}) { // Only populate two rows
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();
    }
    for (int v = 0; v < num_variables; v++) {
      for (int u = v; u < num_variables; ++u) { // Combines with loop over v
        if(u == v) continue;
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }

    // ---- data contribution pseudolikelihood ----
    const arma::mat sufficient_pairwise = sufficient_pairwise_group[group];
    const double suff_pair = sufficient_pairwise(variable1, variable2);

    if(h == 0) {
      log_pp += 2.0 * suff_pair * pairwise_effects(idx, h);
    } else {
      log_pp += 2.0 * suff_pair * proj_g(h-1) * pairwise_effects(idx, h);
    }

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat rest_matrix = obs * pairwise_group;

    for (int v : {variable1, variable2}) {
      const int num_cats = num_categories(v);
      const arma::vec rest_score = rest_matrix.col(v);

      // bound to stabilize exp; use group-specific params consistently
      arma::vec bound = num_cats * rest_score;
      bound = arma::clamp(bound, 0.0, arma::datum::inf);

      arma::vec denom(rest_score.n_elem, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        // base term exp(-bound)
        denom = arma::exp(-bound);
        // main_effects from main_group
        for (int c = 0; c < num_cats; ++c) {
          const double th = main_group(v, c);
          const arma::vec exponent = th + (c + 1) * rest_score - bound;
          denom += arma::exp(exponent);
        }
      } else {
        // linear/quadratic main effects from main_group
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        const int ref = baseline_category(v);
        for (int c = 0; c <= num_cats; ++c) {
          const int centered = c - ref;
          const double quad = quad_effect * centered * centered;
          const double lin  = lin_effect * c;
          const arma::vec exponent = lin + quad + c * rest_score - bound;
          denom += arma::exp(exponent);
        }
      }
      // - sum_i [ bound_i + log denom_i ]
      log_pp -= arma::accu(bound + arma::log(denom));
    }
  }

  // ---- priors ----
  if (h == 0) {
    log_pp += R::dcauchy(pairwise_effects(idx, 0), 0.0, interaction_scale, true);
  } else {
    log_pp += R::dcauchy(pairwise_effects(idx, h), 0.0, difference_scale, true);
  }
  return log_pp;
}