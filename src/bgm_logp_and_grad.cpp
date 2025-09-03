#include <RcppArmadillo.h>
#include "bgm_helper.h"
#include "bgm_logp_and_grad.h"
#include "common_helpers.h"
using namespace Rcpp;



/**
 * Function: log_pseudoposterior_thresholds_component
 *
 * Computes the log pseudo-posterior for one of the main effect parameters.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - rest_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - main_alpha, main_beta: Prior hyperparameters.
 *  - variable: Which variable to compute the log pseudo-posterior for
 *  - category: If ordinal, which category to compute the log pseudo-posterior for
 *  - parameter: If Blume-Capel, which parameter to compute the log pseudo-posterior for (0 = linear, 1 = quadratic)
 *
 * Returns:
 *  - Scalar log pseudo-posterior.
 */
double log_pseudoposterior_thresholds_component (
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const int variable,
    const int category,
    const int parameter
) {
  const int num_persons = rest_matrix.n_rows;
  double log_posterior = 0.0;

  auto log_beta_prior = [&](double threshold_param) {
    return threshold_param * main_alpha - std::log1p (std::exp (threshold_param)) * (main_alpha + main_beta);
  };

  const int num_cats = num_categories(variable);

  if (is_ordinal_variable(variable)) {
    // Prior contribution + sufficient statistic
    const double value = main_effects(variable, category);
    log_posterior += value * num_obs_categories(category + 1, variable);
    log_posterior += log_beta_prior (value);

    // Vectorized likelihood contribution
    // For each person, we compute the unnormalized log-likelihood denominator:
    //   denom = exp (-bound) + sum_c exp (threshold_param_c + (c+1) * rest_score - bound)
    // Where:
    //   - rest_score is the summed interaction score excluding the variable itself
    //   - bound = num_cats * rest_score (for numerical stability)
    //   - threshold_param_c is the threshold parameter for category c (0-based)
    arma::vec rest_score = rest_matrix.col (variable);                     // rest scores for all persons
    arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
    arma::vec denom = arma::exp (-bound);                                      // initialize with base term
    arma::vec threshold_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // threshold parameters

    for (int cat = 0; cat < num_cats; cat++) {
      arma::vec exponent = threshold_param(cat) + (cat + 1) * rest_score - bound;       // exponent per person
      denom += arma::exp (exponent);                                           // accumulate exp terms
    }

    // We then compute the total log-likelihood contribution as:
    //   log_posterior -= bound + log (denom), summed over all persons
    log_posterior -= arma::accu (bound + arma::log (denom));                    // total contribution
  } else {
    const double value = main_effects(variable, parameter);
    const double linear_threshold = main_effects(variable, 0);
    const double quadratic_threshold = main_effects(variable, 1);
    const int ref = reference_category(variable);

    // Prior contribution + sufficient statistic
    log_posterior += value * sufficient_blume_capel(parameter, variable);
    log_posterior += log_beta_prior(value);

    // Vectorized likelihood contribution
    // For each person, we compute the unnormalized log-likelihood denominator:
    //   denom = sum_c exp (θ_lin * c + θ_quad * (c - ref)^2 + c * rest_score - bound)
    // Where:
    //   - θ_lin, θ_quad are linear and quadratic thresholds
    //   - ref is the reference category (used for centering)
    //   - bound = num_cats * rest_score (stabilizes exponentials)
    arma::vec rest_score = rest_matrix.col(variable);                     // rest scores for all persons
    arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
    arma::vec denom(num_persons, arma::fill::zeros);                          // initialize denominator
    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;                                               // centered category
      double quad_term = quadratic_threshold * centered * centered;                    // precompute quadratic term
      double lin_term = linear_threshold * cat;                                      // precompute linear term

      arma::vec exponent = lin_term + quad_term + cat * rest_score - bound;
      denom += arma::exp (exponent);                                           // accumulate over categories
    }

    // The final log-likelihood contribution is then:
    //   log_posterior -= bound + log (denom), summed over all persons
    log_posterior -= arma::accu (bound + arma::log (denom));                    // total contribution
  }

  return log_posterior;
}



/**
 * Function: log_pseudoposterior_interactions_component
 *
 * Computes the log pseudo-posterior contribution from a single pairwise interaction.
 *
 * This includes both the data likelihood contribution and the prior (Cauchy) term.
 * The likelihood is evaluated over both variables jointly and accounts for
 * their variable types (ordinal or Blume-Capel).
 *
 * Inputs:
 *  - pairwise_effects: Matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of main effect parameters.
 *  - observations: Integer matrix of observed category scores.
 *  - num_categories: Vector of category counts per variable.
 *  - inclusion_indicator: Symmetric matrix indicating active interactions.
 *  - is_ordinal_variable: Logical vector indicating which variables are ordinal.
 *  - reference_category: Reference categories for Blume-Capel variables.
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - sufficient_pairwise: Sufficient statistics for each interaction pair.
 *  - var1, var2: Indices of the two variables involved in the interaction.
 *
 * Returns:
 *  - Scalar log pseudo-posterior for the given interaction term.
 */
double log_pseudoposterior_interactions_component (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const int var1,
    const int var2
) {
  const int num_observations = observations.n_rows;

  double log_pseudo_posterior = 2.0 * pairwise_effects(var1, var2) * sufficient_pairwise(var1, var2);

  for (int var : {var1, var2}) {
    int num_categories_var = num_categories (var);

    // Compute rest score: contribution from other variables
    arma::vec rest_scores = observations * pairwise_effects.col (var);
    arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var;
    arma::vec denominator = arma::zeros (num_observations);

    if (is_ordinal_variable (var)) {
      // Ordinal variable: denominator includes exp (-bounds)

      denominator += arma::exp (-bounds);
      for (int category = 0; category < num_categories_var; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * rest_scores - bounds;
        denominator += arma::exp(exponent);
      }

    } else {
      // Binary/categorical variable: quadratic + linear term
      const int ref_cat = reference_category (var);
      for (int category = 0; category <= num_categories_var; category++) {
        int centered_cat = category - ref_cat;
        double lin_term = main_effects (var, 0) * category;
        double quad_term = main_effects (var, 1) * centered_cat * centered_cat;
        arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
        denominator += arma::exp (exponent);
      }
    }

    // Subtract log partition function and bounds adjustment
    log_pseudo_posterior -= arma::accu (arma::log (denominator));
    log_pseudo_posterior -= arma::accu (bounds);
  }

  // Add Cauchy prior terms for included pairwise effects
  if (inclusion_indicator (var1, var2) == 1) {
    log_pseudo_posterior += R::dcauchy (pairwise_effects (var1, var2), 0.0, interaction_scale, true);
  }

  return log_pseudo_posterior;
}



/**
 * Function: log_pseudoposterior
 *
 * Computes the full log pseudo-posterior of the model, including contributions
 * from main effects, pairwise interactions, and all relevant priors.
 *
 * This combines observed sufficient statistics with expected log-likelihood terms,
 * using vectorized likelihood calculations for efficiency.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - pairwise_effects: Matrix of interaction weights.
 *  - inclusion_indicator: Symmetric matrix indicating active pairwise terms.
 *  - observations: Matrix of observed variable scores.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Matrix of observed counts per category and variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel main effects.
 *  - reference_category: Reference categories per variable.
 *  - is_ordinal_variable: Logical vector indicating ordinal variables.
 *  - main_alpha, main_beta: Prior hyperparameters for main effects.
 *  - interaction_scale: Scale parameter for Cauchy prior over interactions.
 *  - sufficient_pairwise: Sufficient statistics for pairwise interactions.
 *  - rest_matrix: Matrix of residual predictors.
 *
 * Returns:
 *  - Scalar total log pseudo-posterior value.
 */
double log_pseudoposterior (
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
) {

  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;

  double log_pseudoposterior = 0.0;

  // Calculate the contribution from the data and the prior
  auto log_beta_prior = [&](double threshold_param) {
    return threshold_param * main_alpha - std::log1p (std::exp (threshold_param)) * (main_alpha + main_beta);
  };

  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        double value = main_effects(variable, cat);
        log_pseudoposterior += num_obs_categories(cat + 1, variable) * value;
        log_pseudoposterior += log_beta_prior(value);
      }
    } else {
      double value = main_effects(variable, 0);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += sufficient_blume_capel(0, variable) * value;

      value = main_effects(variable, 1);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += sufficient_blume_capel(1, variable) * value;
    }
  }
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0) continue;

      double value = pairwise_effects(var1, var2);
      log_pseudoposterior += 2.0 * sufficient_pairwise(var1, var2) * value;
      log_pseudoposterior += R::dcauchy(value, 0.0, interaction_scale, true); // Cauchy prior
    }
  }

  // Calculate the log denominators
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec rest_score = rest_matrix.col (variable);                    // rest scores for all persons
    arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
    bound = arma::clamp(bound, 0.0, arma::datum::inf);                        // only positive bounds to prevent overflow

    arma::vec denom;
    if (is_ordinal_variable(variable)) {
      denom = arma::exp (-bound);                                     // initialize with base term
      arma::vec threshold_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // threshold parameters for variable
      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec exponent = threshold_param(cat) + (cat + 1) * rest_score - bound; // exponent per person
        denom += arma::exp (exponent);                                          // accumulate exp terms
      }
    } else {
      const double lin_effect = main_effects(variable, 0);
      const double quad_effect = main_effects(variable, 1);
      const int ref = reference_category(variable);

      denom.zeros(num_persons);
      for (int cat = 0; cat <= num_cats; cat++) {
        int centered = cat - ref;                                               // centered category
        double quad = quad_effect * centered * centered;                    // precompute quadratic term
        double lin = lin_effect * cat;                                      // precompute linear term
        arma::vec exponent = lin + quad + cat * rest_score - bound;
        denom += arma::exp (exponent);                                           // accumulate over categories
      }
    }

    log_pseudoposterior -= arma::accu (bound + arma::log (denom));                    // total contribution
  }

  return log_pseudoposterior;
}



/**
 * Function: gradient_log_pseudoposterior
 *
 * Computes the full gradient of the log pseudo-posterior with respect to
 * all model parameters (main and pairwise effects).
 *
 * The gradient is computed as the difference between observed and expected
 * sufficient statistics, plus the gradient from prior terms.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - pairwise_effects: Matrix of interaction weights.
 *  - inclusion_indicator: Binary matrix indicating active pairwise terms.
 *  - observations: Matrix of observed scores per variable.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Observed counts of each category per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel main effects.
 *  - reference_category: Reference categories per variable (BC variables).
 *  - is_ordinal_variable: Logical vector for ordinal variables.
 *  - main_alpha, main_beta: Beta prior parameters for main effects.
 *  - interaction_scale: Scale of the Cauchy prior on interactions.
 *  - sufficient_pairwise: Sufficient statistics for pairwise effects.
 *  - rest_matrix: Residual predictor matrix for likelihood computation.
 *
 * Returns:
 *  - Gradient vector corresponding to all parameters (main + interactions).
 */
arma::vec gradient_log_pseudoposterior (
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
) {
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;
  const int num_main = count_num_main_effects (num_categories, is_ordinal_variable);
  const int num_pairwise = num_variables * (num_variables - 1) / 2;

  arma::umat index_matrix(num_variables, num_variables);
  int counter = num_main; // Start after the main effects
  for(int var1 = 0; var1 < num_variables-1; var1++) {
    for(int var2 = var1 + 1; var2 < num_variables; var2++) {
      index_matrix(var1, var2) = counter++;
    }
  }

  arma::vec gradient (num_main + num_pairwise, arma::fill::zeros);

  // Gradients are built up as O - E + gradient_prior
  // - O is the observed value of the sufficient statistic
  // - E is the expected value of the sufficient statistic
  // - gradient_prior is the gradient of the prior distribution

  // Calculate the observed sufficient statistics
  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }
      offset += num_cats;
    } else {
      gradient(offset) = sufficient_blume_capel(0, variable);
      gradient(offset + 1) = sufficient_blume_capel(1, variable);
      offset += 2;
    }
  }
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0)
        continue;
      int location = index_matrix(var1, var2);
      gradient(location) = 2.0 * sufficient_pairwise(var1, var2);
    }
  }

  // Calculate the expected sufficient statistics
  offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec rest_score = rest_matrix.col(variable);
    arma::vec bound = num_cats * rest_score;
    bound = arma::clamp(bound, 0.0, arma::datum::inf);

    if (is_ordinal_variable(variable)) {
      arma::vec threshold_param = main_effects.row(variable).cols(0, num_cats - 1).t();
      bound += threshold_param.max();

      arma::mat exponents(num_persons, num_cats);
      for (int cat = 0; cat < num_cats; cat++) {
        exponents.col(cat) = threshold_param(cat) + (cat + 1) * rest_score - bound;
      }

      arma::mat probs = arma::exp (exponents);
      arma::vec denom = arma::sum(probs, 1) + arma::exp (-bound);
      probs.each_col() /= denom;

      // Expected sufficient statistics main effects
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) -= arma::accu (probs.col(cat));
      }

      // Expected sufficient statistics pairwise effects
      for (int var2 = 0; var2 < num_variables; var2++) {
        if (inclusion_indicator(variable, var2) == 0 || variable == var2)
          continue;

        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat < num_cats; cat++) {
          expected_value += (cat + 1) * probs.col(cat) % observations.col(var2);
        }
        int location = (variable < var2) ? index_matrix(variable, var2) : index_matrix(var2, variable);
        gradient(location) -= arma::accu(expected_value);
      }

      offset += num_cats;
    } else {
      const int ref = reference_category(variable);
      const double lin_effect = main_effects(variable, 0);
      const double quad_effect = main_effects(variable, 1);

      arma::mat exponents(num_persons, num_cats + 1);
      for (int cat = 0; cat <= num_cats; cat++) {
        int score = cat;
        int centered = score - ref;
        double lin = lin_effect * score;
        double quad = quad_effect * centered * centered;
        exponents.col(cat) = lin + quad + score * rest_score - bound;
      }
      arma::mat probs = arma::exp (exponents);
      arma::vec denom = arma::sum(probs, 1);
      probs.each_col() /= denom;

      arma::ivec lin_score =  arma::regspace<arma::ivec>(0, num_cats);
      arma::ivec quad_score = arma::square(lin_score - ref);

      // Expected sufficient statistics main effects
      gradient(offset) -= arma::accu(probs * lin_score);
      gradient(offset + 1) -= arma::accu(probs * quad_score);;

      // Expected sufficient statistics pairwise effects
      for (int var2 = 0; var2 < num_variables; var2++) {
        if (inclusion_indicator(variable, var2) == 0 || variable == var2)
          continue;

        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat < num_cats; cat++) {
          expected_value += (cat + 1) * probs.col(cat + 1) % observations.col(var2); // Here the zero category score is in probs, so we skip it
        }
        int location = (variable < var2) ? index_matrix(variable, var2) : index_matrix(var2, variable);
        gradient(location) -= arma::accu(expected_value);
      }
      offset += 2;
    }
  }

  // Calculate the gradient contribution from the prior distribution
  offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        const double p = 1.0 / (1.0 + std::exp (-main_effects(variable, cat)));
        gradient(offset + cat) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += num_cats;
    } else {
      for (int i = 0; i < 2; i++) {
        const double threshold_param = main_effects(variable, i);
        const double p = 1.0 / (1.0 + std::exp (-threshold_param));
        gradient(offset + i) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += 2;
    }
  }
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0)
        continue;

      // ---- Gradient contribution from Cauchy prior
      int location = index_matrix(var1, var2);
      const double effect = pairwise_effects (var1, var2);
      gradient (location) -= 2.0 * effect / (effect * effect + interaction_scale * interaction_scale);
    }
  }

  return gradient;
}



/**
 * Function: gradient_log_pseudoposterior_active
 *
 * Computes a compressed gradient of the log pseudo-posterior, containing only
 * the active (included) pairwise effects and all main effects.
 *
 * This function extracts the relevant subset of the full gradient for use
 * in optimization or sampling where inactive parameters are not updated.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - pairwise_effects: Matrix of interaction weights.
 *  - inclusion_indicator: Matrix indicating which interactions are active.
 *  - observations: Matrix of categorical data.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category index per variable.
 *  - is_ordinal_variable: Logical vector of ordinal variable flags.
 *  - main_alpha, main_beta: Prior parameters for main effects.
 *  - interaction_scale: Cauchy scale for interaction prior.
 *  - sufficient_pairwise: Sufficient statistics for interactions.
 *  - rest_matrix: Matrix of residual predictor terms.
 *
 * Returns:
 *  - Vector of active parameter gradients (main effects + included interactions).
 */
arma::vec gradient_log_pseudoposterior_active (
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
) {

  // Step 1: Reconstruct full parameter matrices
  arma::vec full_gradient = gradient_log_pseudoposterior(
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, main_alpha, main_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  // Step 2: Extract compressed gradient (same layout as new vectorizer)
  const int num_variables = num_categories.n_elem;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);

  int num_active = 0;
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 1) num_active++;
    }
  }

  arma::vec reduced_gradient(num_main + num_active);
  int offset = 0;

  // Copy main effect gradient directly
  for (int i = 0; i < num_main; ++i) {
    reduced_gradient[offset++] = full_gradient[i];
  }

  // Copy active interaction gradient values
  int full_idx = num_main;
  for (int i = 0; i < num_variables - 1; ++i) {
    for (int j = i + 1; j < num_variables; ++j) {
      if (inclusion_indicator(i, j) == 1) {
        reduced_gradient[offset++] = full_gradient[full_idx];
      }
      full_idx++;
    }
  }

  return reduced_gradient;
}



/**
 * Function: compute_log_likelihood_ratio_for_variable
 *
 * Computes the log pseudo-likelihood ratio contribution for a single variable,
 * comparing a proposed vs. current interaction value. This is used to evaluate
 * Metropolis-Hastings updates to a pairwise interaction parameter.
 *
 * The function is vectorized over persons and supports both ordinal and
 * Blume-Capel variables.
 *
 * Inputs:
 *  - variable: Index of the variable whose likelihood contribution is evaluated.
 *  - interacting_score: Vector of category scores for the interacting variable (one per person).
 *  - proposed_state: Proposed interaction value.
 *  - current_state: Current interaction value.
 *  - main_effects: Matrix of threshold parameters [variables × categories].
 *  - num_categories: Vector with number of categories per variable.
 *  - rest_matrix: Current matrix of residual predictors (one column per variable).
 *  - observations: Data matrix of categorical scores (only used for row count).
 *  - is_ordinal_variable: Logical vector (1 if ordinal, 0 if Blume-Capel).
 *  - reference_category: Reference category per variable (for BC variables).
 *
 * Returns:
 *  - The total log pseudo-likelihood ratio for the given variable, summed over all persons.
 */
double compute_log_likelihood_ratio_for_variable (
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
) {
  // Convert interaction score vector to double precision
  arma::vec interaction = arma::conv_to<arma::vec>::from (interacting_score);

  const int num_persons = rest_matrix.n_rows;
  const int num_categories_var = num_categories (variable);

  // Compute adjusted linear predictors without the current interaction
  arma::vec rest_scores = rest_matrix.col (variable) - interaction * current_state;

  // Stability bound for softmax (scaled by number of categories)
  arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_persons)) * num_categories_var;

  arma::vec denom_current = arma::zeros (num_persons);
  arma::vec denom_proposed = arma::zeros (num_persons);

  if (is_ordinal_variable (variable)) {
    denom_current += arma::exp(-bounds);
    denom_proposed += arma::exp(-bounds);

    for (int category = 0; category < num_categories_var; category++) {
      const double main = main_effects(variable, category);
      const int score = category + 1;

      for (int person = 0; person < num_persons; person++) {
        const double base = main + score * rest_scores[person] - bounds[person];

        const double exp_current = std::exp(base + score * interaction[person] * current_state);
        const double exp_proposed = std::exp(base + score * interaction[person] * proposed_state);

        denom_current[person] += exp_current;
        denom_proposed[person] += exp_proposed;
      }
    }

    // Vectorized version (= slower?):
    // denom_current += arma::exp (-bounds);
    // denom_proposed += arma::exp (-bounds);
    //
    // for (int category = 0; category < num_categories_var; category++) {
    //   arma::vec exponent = main_effects (variable, category) + (category + 1) * rest_scores;
    //   denom_current += arma::exp (exponent + (category + 1) * interaction * current_state - bounds);
    //   denom_proposed += arma::exp (exponent + (category + 1) * interaction * proposed_state - bounds);
    // }

    // ----- Version with minimal exp calls -----
    // arma::vec base_current  = arma::exp(interaction * current_state);
    // arma::vec base_proposed = arma::exp(interaction * proposed_state);
    //
    // arma::vec e_accum_current  = base_current;
    // arma::vec e_accum_proposed = base_proposed;
    //
    // arma::vec exp_bounds = arma::exp(-bounds);
    // arma::vec exp_rest   = arma::exp(rest_scores);
    // arma::vec er_accum   = exp_rest;
    //
    // denom_current  += exp_bounds;
    // denom_proposed += exp_bounds;
    //
    // for (int category = 0; category < num_categories_var; category++) {
    //   const double em    = std::exp(main_effects(variable, category));
    //   for (int person = 0; person < num_persons; person++) {
    //     double w_common = em * er_accum[person] * exp_bounds[person];
    //     denom_current[person]  += w_common * e_accum_current[person];
    //     denom_proposed[person] += w_common * e_accum_proposed[person];
    //     er_accum[person]        *= exp_rest[person];
    //     e_accum_current[person] *= base_current[person];
    //     e_accum_proposed[person] *= base_proposed[person];
    //   }
    // }

  } else {
    // Binary or categorical variable: linear + quadratic score
    const int ref_cat = reference_category (variable);

    for (int category = 0; category <= num_categories_var; category++) {
      int centered = category - ref_cat;
      double lin_term = main_effects (variable, 0) * category;
      double quad_term = main_effects (variable, 1) * centered * centered;
      arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;

      denom_current += arma::exp (exponent + category * interaction * current_state);
      denom_proposed += arma::exp (exponent + category * interaction * proposed_state);
    }
  }

  // Accumulated log-likelihood difference across persons
  return arma::accu (arma::log (denom_current) - arma::log (denom_proposed));
}



/**
 * Function: log_pseudolikelihood_ratio_interaction
 *
 * Computes the change in log pseudo-likelihood when a pairwise interaction parameter
 * between two variables is updated from a current value to a proposed value.
 *
 * This function evaluates:
 *  1. The direct contribution of the interaction to the joint linear predictor.
 *  2. The change in pseudo-likelihood for each affected variable (via softmax terms),
 *     accounting for the influence of the interaction on their respective rest scores.
 *
 * Inputs:
 *  - pairwise_effects: Current matrix of pairwise interaction parameters [V × V].
 *  - main_effects: Matrix of main effect (threshold) parameters [V × max_categories].
 *  - observations: Integer matrix of observed category scores [N × V].
 *  - num_categories: Vector of category counts per variable [V].
 *  - num_persons: Number of individuals (rows in the data).
 *  - variable1, variable2: Indices of the two interacting variables (0-based).
 *  - proposed_state: Proposed new interaction weight.
 *  - current_state: Current interaction weight.
 *  - rest_matrix: Matrix of residual linear predictors [N × V].
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal (1) or BC (0).
 *  - reference_category: Vector of reference categories per variable (used for BC variables).
 *
 * Returns:
 *  - The log pseudo-likelihood ratio:
 *      log p(y | β_proposed) - log p(y | β_current)
 */
double log_pseudolikelihood_ratio_interaction (
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
) {
  double log_ratio = 0.0;
  const double delta = proposed_state - current_state;

  // Extract score vectors for both variables across all persons
  arma::ivec score1 = observations.col(variable1);
  arma::ivec score2 = observations.col(variable2);

  // (1) Direct interaction contribution to the linear predictor:
  //     Δβ × ∑(score1_i × score2_i) for all persons i
  log_ratio += 2.0 * sufficient_pairwise(variable1, variable2) * delta;

  // (2) Change in pseudo-likelihood for variable1 due to the update in its interaction with variable2
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable1, score2, proposed_state, current_state, main_effects,
    num_categories, rest_matrix, observations, is_ordinal_variable,
    reference_category
  );

  // (3) Symmetric change for variable2 due to its interaction with variable1
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable2, score1, proposed_state, current_state, main_effects,
    num_categories, rest_matrix, observations, is_ordinal_variable,
    reference_category
  );

  return log_ratio;
}