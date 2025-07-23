//#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gibbs_functions_edge_prior.h"
#include "nuts_functions.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;


/**
 * Adapts the log step size using dual averaging during MCMC burn-in.
 *
 * Implements the dual averaging algorithm to adaptively update the MALA step size
 * toward a target acceptance rate. Only used during burn-in.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - state: Vector of length 3, modified in-place:
 *      * state[0] = current log step size
 *      * state[1] = running average of log step size
 *      * state[2] = running average of acceptance error
 *
 * Modifies:
 *  - `state` vector in-place to update the log step size parameters.
 */
inline void update_step_size_with_dual_averaging (
    const double initial_step_size,
    const double acceptance_probability,
    const int iteration,
    arma::vec& state,
    const double target_acceptance
) {
  const double target_log_step_size = std::log (10.0 * initial_step_size);
  constexpr int stabilization_offset = 10;

  double& log_step_size = state[0];
  double& log_step_size_avg = state[1];
  double& acceptance_error_avg = state[2];

  const double adjusted_iter = iteration + stabilization_offset;
  const double error = target_acceptance - acceptance_probability;

  acceptance_error_avg = (1.0 - 1.0 / adjusted_iter) * acceptance_error_avg +
    (1.0 / adjusted_iter) * error;

  log_step_size = target_log_step_size - std::sqrt (static_cast<double> (iteration)) / 0.05 * acceptance_error_avg;

  const double weight = std::pow (static_cast<double> (iteration), -0.75);
  log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
}



/**
 * Function: update_step_size_with_robbins_monro
 *
 * Performs Robbins-Monro adaptation of the step size on the log scale.
 * Applies exponential decay to gradually reduce adaptation influence over time.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - step_size_mala: Step size to update (in-place).
 *
 * Notes:
 *  - Decay rate is 0.75 for stable adaptation.
 */
inline void update_step_size_with_robbins_monro (
    const double acceptance_probability,
    const int iteration,
    double& step_size_mala,
    const double target_acceptance
) {
  constexpr double decay_rate = 0.75;

  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow (static_cast<double> (iteration), -decay_rate);

  double log_step_size = std::log (step_size_mala);
  log_step_size += error * decay;
  step_size_mala = std::exp (log_step_size);
}



/**
 * Function: update_proposal_sd_with_robbins_monro
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 *
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - rm_weight: Robbins-Monro adaptation weight (e.g. iteration^{-0.75}).
 *
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double rm_weight,
    const double target_acceptance
) {
  constexpr double rm_lower_bound = 0.001;
  constexpr double rm_upper_bound = 2.0;

  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = std::exp (observed_log_acceptance_probability);
  }

  // Robbins-Monro update step
  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  // Handle NaNs robustly
  if (std::isnan (updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp (updated_sd, rm_lower_bound, rm_upper_bound);
}



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
inline int count_num_main_effects (
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



/**
 * Function: impute_missing_values_for_graphical_model
 *
 * Imputes missing values in the observation matrix using the current model parameters.
 * For each missing entry, a category is sampled from the posterior predictive distribution,
 * and if the value changes, all dependent structures are updated accordingly.
 *
 * Inputs:
 *  - pairwise_effects: matrix of interaction weights.
 *  - main_effects: matrix of threshold parameters (main effects).
 *  - missing_index: matrix of (person, variable) indices for missing entries.
 *  - num_categories: vector of category counts per variable.
 *  - is_ordinal_variable: logical vector (0/1) indicating ordinal vs. Blume-Capel.
 *  - reference_category: vector of reference categories per variable.
 *
 * Updates (in-place):
 *  - observations: matrix of imputed categorical scores.
 *  - num_obs_categories: counts of observed categories per variable.
 *  - sufficient_blume_capel: sufficient stats for Blume-Capel thresholds.
 *  - rest_matrix: updated linear predictors after imputation.
 */
void impute_missing_values_for_graphical_model (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& observations,
    arma::imat& num_obs_categories,
    arma::imat& sufficient_blume_capel,
    const arma::ivec& num_categories,
    arma::mat& rest_matrix,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    arma::imat& sufficient_pairwise
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_index.n_rows;
  const int max_num_categories = num_categories.max ();

  arma::vec category_probabilities (max_num_categories + 1);

  for (int miss = 0; miss < num_missings; miss++) {
    const int person = missing_index (miss, 0);
    const int variable = missing_index (miss, 1);

    const double rest_score = rest_matrix (person, variable);
    const int num_cats = num_categories (variable);
    const bool is_ordinal = is_ordinal_variable (variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const double exponent = main_effects (variable, cat) + score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const int ref = reference_category (variable);

      cumsum = std::exp (main_effects (variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const int centered = score - ref;
        const double exponent =
          main_effects (variable, 0) * score +
          main_effects (variable, 1) * centered * centered +
          score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    }

    // Sample from categorical distribution via inverse transform
    const double u = R::unif_rand () * cumsum;
    int sampled_score = 0;
    while (u > category_probabilities[sampled_score]) {
      sampled_score++;
    }

    const int new_value = sampled_score;
    const int old_value = observations(person, variable);

    if (new_value != old_value) {
      // Update observation matrix
      observations(person, variable) = new_value;

      if (is_ordinal) {
        num_obs_categories(old_value, variable)--;
        num_obs_categories(new_value, variable)++;
      } else {
        const int ref = reference_category(variable);
        const int delta = new_value - old_value;
        const int delta_sq =
          (new_value - ref) * (new_value - ref) -
          (old_value - ref) * (old_value - ref);

        sufficient_blume_capel(0, variable) += delta;
        sufficient_blume_capel(1, variable) += delta_sq;
      }
      // Update residuals across all variables
      for (int var = 0; var < num_variables; var++) {
        const double delta_score = (new_value - old_value) * pairwise_effects(var, variable);
        rest_matrix(person, var) += delta_score;
      }
    }
  }

  // Update sufficient statistics for pairwise effects
  sufficient_pairwise = observations.t() * observations;
  // This could be done more elegantly....
}



/**
 * Function: log_pseudoposterior_thresholds
 *
 * Computes the log pseudo-posterior for all main effect parameters.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - rest_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Returns:
 *  - Scalar log pseudo-posterior.
 */
double log_pseudoposterior_thresholds (
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta
) {
  const int num_variables = rest_matrix.n_cols;
  const int num_persons = rest_matrix.n_rows;
  double log_posterior = 0.0;

  auto log_beta_prior = [&](double threshold_param) {
    return threshold_param * threshold_alpha - std::log1p (std::exp (threshold_param)) * (threshold_alpha + threshold_beta);
  };

  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Prior contribution + sufficient statistic
      for (int cat = 0; cat < num_cats; cat++) {
        const double threshold_param = main_effects(variable, cat);
        log_posterior += threshold_param * num_obs_categories(cat + 1, variable);
        log_posterior += log_beta_prior (threshold_param);
      }

      // Vectorized likelihood contribution
      // For each person, we compute the unnormalized log-likelihood denominator:
      //   denom = exp (-bound) + sum_c exp (threshold_param_c + (c+1) * rest_score - bound)
      // Where:
      //   - rest_score is the summed interaction score excluding the variable itself
      //   - bound = num_cats * rest_score (for numerical stability)
      //   - threshold_param_c is the threshold parameter for category c (0-based)
      arma::vec rest_score = rest_matrix.col (variable);                    // rest scores for all persons
      arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
      bound = arma::clamp(bound, 0.0, arma::datum::inf);                        // only positive bounds to prevent overflow
      arma::vec denom = arma::exp (-bound);                                     // initialize with base term
      arma::vec threshold_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // threshold parameters for variable

      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec exponent = threshold_param(cat) + (cat + 1) * rest_score - bound; // exponent per person
        denom += arma::exp (exponent);                                          // accumulate exp terms
      }

      // We then compute the total log-likelihood contribution as:
      //   log_posterior -= bound + log (denom), summed over all persons
      log_posterior -= arma::accu (bound + arma::log (denom));                    // total contribution
    } else {
      const double linear_threshold = main_effects(variable, 0);
      const double quadratic_threshold = main_effects(variable, 1);
      const int ref = reference_category(variable);

      // Prior contribution + sufficient statistic
      log_posterior += linear_threshold * sufficient_blume_capel(0, variable);
      log_posterior += log_beta_prior(linear_threshold);
      log_posterior += quadratic_threshold * sufficient_blume_capel(1, variable);
      log_posterior += log_beta_prior(quadratic_threshold);;

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
  }

  return log_posterior;
}



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
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
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
    const double threshold_alpha,
    const double threshold_beta,
    const int variable,
    const int category,
    const int parameter
) {
  const int num_persons = rest_matrix.n_rows;
  double log_posterior = 0.0;

  auto log_beta_prior = [&](double threshold_param) {
    return threshold_param * threshold_alpha - std::log1p (std::exp (threshold_param)) * (threshold_alpha + threshold_beta);
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
 * Function: gradient_log_pseudoposterior_thresholds
 *
 * Computes the gradient of the log pseudo-posterior for threshold parameters.
 * Includes both data likelihood and logistic-Beta prior components.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - rest_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Returns:
 *  - A flat gradient vector in the same order as vectorized thresholds.
 */
arma::vec gradient_log_pseudoposterior_thresholds (
    const arma::mat& main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha = 1.0,
    const double threshold_beta = 1.0
) {
  const int num_variables = rest_matrix.n_cols;
  const int num_persons = rest_matrix.n_rows;
  const int num_parameters = count_num_main_effects(num_categories, is_ordinal_variable);

  arma::vec gradient(num_parameters, arma::fill::zeros);
  int offset = 0;

  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }

      // Vectorized computation of expected category counts
      //
      // For each person, we compute probabilities over categories:
      //   probs[p, c] = exp (θ_c + (c+1) * rest_score_p - bound_p)
      // where:
      //   - θ_c is the threshold for category c
      //   - bound_p = max(θ) + num_cats * rest_score_p (numerical stabilization)
      // We normalize across categories and subtract expected counts from the gradient.

      arma::vec rest_score = rest_matrix.col(variable);                     // rest scores per person
      arma::vec threshold_param = main_effects.row(variable).cols(0, num_cats - 1).t();   // thresholds for current variable
      arma::vec bound = threshold_param.max() + num_cats * rest_score;          // vector of bounds per person
      bound = arma::clamp(bound, 0.0, arma::datum::inf); //only positive bounds

      arma::mat exponents(num_persons, num_cats);                               // log unnormalized probabilities
      for (int cat = 0; cat < num_cats; cat++) {
        exponents.col(cat) = threshold_param(cat) + (cat + 1) * rest_score - bound;
      }

      arma::mat probs = arma::exp (exponents);                                   // unnormalized probabilities
      arma::vec denom = arma::sum(probs, 1) + arma::exp (-bound);                // normalization constants per person

      // Accumulate gradient contributions by subtracting expected counts
      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec normalized = probs.col(cat) / denom;                          // normalized prob for category
        gradient(offset + cat) -= arma::accu (normalized);                      // accumulate gradient
      }

      // Compute prior contribution to gradient
      for (int cat = 0; cat < num_cats; cat++) {
        const double p = 1.0 / (1.0 + std::exp (-threshold_param(cat)));
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += num_cats;
    } else {
      // Blume-Capel variable
      const int ref = reference_category(variable);
      const double linear_threshold = main_effects(variable, 0);
      const double quadratic_threshold = main_effects(variable, 1);

      gradient(offset) = sufficient_blume_capel(0, variable);
      gradient(offset + 1) = sufficient_blume_capel(1, variable);

      // Vectorized computation of expected statistics
      //
      // For each person, compute the expected sufficient statistics:
      //   - sum_lin  = E[score]
      //   - sum_quad = E[(score - ref)^2]
      // where probabilities are softmax-like terms derived from θ_lin and θ_quad.
      //
      // This replaces the nested loop with vectorized accumulation over categories.
      arma::vec rest_score = rest_matrix.col(variable);                     // Residuals per person
      arma::vec bound = num_cats * rest_score;                                  // Stabilization bound
      arma::vec denom = arma::exp (quadratic_threshold * ref * ref - bound);    // Initial term at score = 0

      arma::vec sum_lin(num_persons, arma::fill::zeros);                        // E[score]
      arma::vec sum_quad = ref * ref * denom;                                   // E[(score - ref)^2], starts at score = 0

      for (int cat = 0; cat < num_cats; cat++) {
        int score = cat + 1;
        int centered = score - ref;

        double lin_term = linear_threshold * score;
        double quad_term = quadratic_threshold * centered * centered;

        arma::vec exponent = lin_term + quad_term + score * rest_score - bound;
        arma::vec weight = arma::exp (exponent);                                 // Unnormalized probabilities

        sum_lin += weight * score;                                              // Accumulate score-weighted terms
        sum_quad += weight * centered * centered;                               // Accumulate centered^2-weighted terms
        denom += weight;                                                        // Update normalization constant
      }

      // Finalize the gradient updates
      gradient(offset) -= arma::accu (sum_lin / denom);                          // Gradient for θ_lin
      gradient(offset + 1) -= arma::accu (sum_quad / denom);                     // Gradient for θ_quad


      // Compute prior contribution to gradient
      for (int i = 0; i < 2; i++) {
        const double threshold_param = main_effects(variable, i);
        const double p = 1.0 / (1.0 + std::exp (-threshold_param));
        gradient(offset + i) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += 2;
    }
  }

  return gradient;
}



/**
 * Function: find_reasonable_initial_step_size_thresholds
 *
 * Heuristically selects an initial step size for MALA updates of the threshold parameters
 * by adapting the acceptance rate toward a target value (~0.574). The procedure
 * increases or decreases the log step size until the Metropolis-Hastings acceptance
 * crosses the threshold, using a method adapted from Algorithm 4 in:
 *
 *   Hoffman, M.D. & Gelman, A. (2014).
 *   The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.
 *   Journal of Machine Learning Research, 15, 1593–1623.
 *
 * Inputs:
 *  - main_effects: Initial threshold matrix.
 *  - rest_matrix: Current residual predictor matrix.
 *  - num_categories, num_obs_categories: Category structure for each variable.
 *  - sufficient_blume_capel, reference_category: Data model components.
 *  - is_ordinal_variable: Logical vector indicating ordinal predictors.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters for thresholds.
 *
 * Returns:
 *  - A positive scalar step size tuned for MALA with target acceptance ~0.574.
 *
 * Notes:
 *  - Uses log-space doubling/halving; returns 0.01 if the search fails.
 *  - Assumes Euclidean (identity) preconditioning.
 */
double find_reasonable_initial_step_size_thresholds (
    arma::mat main_effects,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double target_acceptance
) {
  constexpr double initial_step_size = 0.1;
  constexpr int max_attempts = 20;
  constexpr double max_log_step = 10.0;

  double log_step_size = std::log (initial_step_size);

  // Current state and log posterior
  const arma::vec current_state = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
  const arma::vec current_grad = gradient_log_pseudoposterior_thresholds(
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  const double current_log_post = log_pseudoposterior_thresholds(
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  int direction = 0;
  double accept_prob = 0.0;

  // Propose initial MALA step: θ' = θ + ½ε ∇log p(θ) + √ε * N(0, I)
  {
    const double step_size = std::exp (log_step_size);
    const double sqrt_step = std::sqrt(step_size);
    const arma::vec proposed_state = current_state + 0.5 * step_size * current_grad +
      sqrt_step * arma::randn(current_state.n_elem);
    const arma::mat proposed_main_effects = unvectorize_thresholds(
      proposed_state, num_categories, is_ordinal_variable
    );

    const double proposed_log_post = log_pseudoposterior_thresholds(
      proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );
    const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds(
      proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );

    const arma::vec forward = current_state + 0.5 * step_size * current_grad;
    const arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

    const double log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
    const double log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
    const double log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

    accept_prob = std::min(1.0, std::exp (log_accept));
    if (std::abs(accept_prob - target_acceptance) < 0.1) {
      return std::exp (log_step_size);
    }
    direction = (accept_prob > target_acceptance) ? 1 : -1;
  }

  // Log-scale search for step size that brackets the target acceptance rate
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    log_step_size += direction;
    const double step_size = std::exp (log_step_size);
    const double sqrt_step = std::sqrt(step_size);

    const arma::vec proposed_state = current_state + 0.5 * step_size * current_grad +
      sqrt_step * arma::randn(current_state.n_elem);
    const arma::mat proposed_main_effects = unvectorize_thresholds(
      proposed_state, num_categories, is_ordinal_variable
    );

    const double proposed_log_post = log_pseudoposterior_thresholds(
      proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );
    const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds(
      proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );

    const arma::vec forward = current_state + 0.5 * step_size * current_grad;
    const arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

    const double log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
    const double log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
    const double log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

    const double new_accept_prob = std::min(1.0, std::exp (log_accept));

    // Exit if acceptance flips across the target
    if ((direction == 1 && new_accept_prob < target_acceptance) ||
        (direction == -1 && new_accept_prob > target_acceptance)) {
      break;
    }

    accept_prob = new_accept_prob;

    if (std::abs(log_step_size) > max_log_step) {
      Rcpp::Rcout << "Warning: Step size search failed. Falling back to default (0.01)." << std::endl;
      return 0.01;
    }
  }

  return std::exp (log_step_size);
}



/**
 * Function: initialize_fisher_preconditioner
 *
 * Initializes the square root of the inverse Fisher information matrix,
 * based on the outer product of the gradient. This serves as the starting
 * preconditioner for Fisher-MALA after burn-in.
 *
 * Based on:
 *   Titsias, M.K. (2023).
 *   Gradient-Based MCMC Using Preconditioned Langevin Dynamics.
 *   JMLR, 24(216):1–40.
 *
 * Inputs:
 *  - grad: Gradient vector at the initial state.
 *
 * Returns:
 *  - Square root of the inverse Fisher approximation (d × d matrix).
 *
 * Notes:
 *  - A damping parameter is used to regularize the preconditioner.
 */
inline arma::mat initialize_fisher_preconditioner(
    const arma::vec& grad
) {
  constexpr double damping_par = 10.0;

  const int dim = grad.n_elem;
  const double inner = arma::dot (grad, grad);
  const arma::mat outer = grad * grad.t();

  // Shrinkage ratio: balances curvature and damping
  const double shrinkage = 1.0 / (1.0 + std::sqrt(damping_par / (damping_par + inner)));

  // Fisher inverse root approximation (Titsias 2023, Eq. 12 form)
  arma::mat sqrt_inv_fisher =
    (1.0 / std::sqrt(damping_par)) *
    (arma::eye(dim, dim) - shrinkage * outer / (damping_par + inner));

  return sqrt_inv_fisher;
}



/**
 * Function: update_fisher_preconditioner
 *
 * Performs a rank-1 update of the inverse Fisher matrix approximation,
 * as described in:
 *
 *   Titsias, M.K. (2023).
 *   Gradient-Based MCMC Using Preconditioned Langevin Dynamics.
 *   JMLR, 24(216):1–40.
 *
 * Inputs:
 *  - sqrt_inv_fisher: Current square root inverse Fisher matrix (updated in-place).
 *  - score_diff: Rao-Blackwellized score difference vector (∝ proposed_grad - current_grad).
 */
inline void update_fisher_preconditioner (
    arma::mat& sqrt_inv_fisher,
    const arma::vec& score_diff
) {
  // Transform score difference into scaled preconditioned space
  const arma::vec phi = sqrt_inv_fisher.t() * score_diff;
  const double inner = arma::dot (phi, phi);
  const arma::mat outer = phi * phi.t();

  // Compute the Titsias update scaling factor
  const double r = 1.0 / (1.0 + std::sqrt(1.0 / (1.0 + inner)));

  // Apply rank-1 update to the inverse Fisher root
  sqrt_inv_fisher -= r * sqrt_inv_fisher * outer / (1.0 + inner);
}



/**
 * Function: update_thresholds_with_fisher_mala
 *
 * Performs a Fisher-preconditioned MALA update of the threshold parameters.
 * The proposal uses the inverse Fisher matrix (either identity during burn-in or
 * adapted afterward) to scale the Langevin drift and noise.
 *
 * During burn-in, the step size is adapted using dual averaging. After burn-in,
 * the step size is adapted using Robbins-Monro, and the Fisher matrix is updated
 * with a rank-1 outer product of the score difference between the proposed and
 * current states (Titsias, 2023).
 *
 * The step size is trace-scaled to match the average proposal magnitude across
 * dimensions, following Proposition 1 from Titsias (2023).
 *
 * Modifies:
 *  - main_effects
 *  - step_size
 *  - dual_averaging_state
 *  - sqrt_inv_fisher
 */
void update_thresholds_with_fisher_mala (
    arma::mat& main_effects,
    double& step_size,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int iteration,
    arma::vec& dual_averaging_state,
    arma::mat& sqrt_inv_fisher,
    const double threshold_alpha,
    const double threshold_beta,
    const double initial_step_size,
    const double target_accept_thresholds,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end
) {
  // --- Compute current parameter vector and its gradient ---
  const arma::vec current_state = vectorize_thresholds (
    main_effects, num_categories, is_ordinal_variable
  );

  const arma::vec current_grad = gradient_log_pseudoposterior_thresholds (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post = log_pseudoposterior_thresholds (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const int dim = current_state.n_elem;

  // --- If burn-in just ended, initialize Fisher matrix using gradient
  if (iteration == warmup_stageI_end) {
    sqrt_inv_fisher = initialize_fisher_preconditioner (current_grad);
  }

  // --- Set inverse Fisher matrix: identity during warm-up, adapted afterward
  arma::mat inv_fisher;
  if (iteration >= warmup_stageI_end) {
    inv_fisher = sqrt_inv_fisher * sqrt_inv_fisher.t();
  } else {
    inv_fisher.eye(dim, dim);
  }

  // --- Construct Fisher-preconditioned MALA proposal ---
  const double trace_inv_fisher = arma::trace(inv_fisher);
  const double scaled_step_size = step_size / (trace_inv_fisher / dim);
  const double sqrt_step = std::sqrt(scaled_step_size);

  // Drift (mean shift) and stochastic noise
  const arma::vec proposal_drift = 0.5 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec noise = sqrt_inv_fisher * arma::randn(dim);

  const arma::vec proposed_state = current_state + proposal_drift + sqrt_step * noise;

  // --- Map vector back to main_effects matrix
  const arma::mat proposed_main_effects = unvectorize_thresholds (
    proposed_state, num_categories, is_ordinal_variable
  );
  // --- Evaluate proposed log posterior and gradient
  const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds (
    proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post_prop = log_pseudoposterior_thresholds (
    proposed_main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Compute MALA acceptance correction (Titsias 2023, Prop. 1) ---
  const arma::vec forward_proposal_residual =
    proposed_state - current_state -
    0.25 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec reverse_proposal_residual =
    current_state - proposed_state -
    0.25 * scaled_step_size * inv_fisher * proposed_grad;

  const double log_forward = 0.5 * arma::dot (forward_proposal_residual, current_grad);
  const double log_backward = 0.5 * arma::dot (reverse_proposal_residual, proposed_grad);

  const double log_accept = log_post_prop - log_post + (log_backward - log_forward);
  const double accept_prob = std::min(1.0, std::exp (log_accept));

  // --- Accept or reject proposed move ---
  if (std::log (R::unif_rand()) < log_accept) {
    main_effects = proposed_main_effects;
  }

  // --- Update step size and Fisher matrix ---
  if (iteration < warmup_stageI_end) {
    // During warm-up stage I: dual averaging adaptation
    update_step_size_with_dual_averaging (
        initial_step_size, accept_prob, iteration + 1, dual_averaging_state,
        target_accept_thresholds
    );
    step_size = std::exp (dual_averaging_state[1]);
  } else if (iteration < warmup_stageII_end) {
    // After warm-up stage I: Robbins-Monro + Fisher preconditioner update
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageI_end + 1, step_size,
        target_accept_thresholds
    );

    const arma::vec score_diff =
      std::sqrt(accept_prob) * (proposed_grad - current_grad);

    update_fisher_preconditioner(sqrt_inv_fisher, score_diff);
  } else if (iteration < warmup_stageIII_end) {
    // After warm-up stages I & II: Robbins-Monro
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageII_end + 1, step_size,
        target_accept_thresholds
    );
  }
}



/**
 * Function: update_thresholds_with_adaptive_metropolis
 *
 * Performs an adaptive Metropolis update of all threshold parameters, with
 * Robbins-Monro tuning.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - num_persons: Number of observations.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - rest_matrix: Residual scores.
 *  - proposal_sd_main: Matrix of proposal SDs for each variable (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *
 * Modifies:
 *  - main_effects (for the given variable)
 *  - proposal_sd_main
 */
void update_thresholds_with_adaptive_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int num_persons,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& rest_matrix,
    arma::mat& proposal_sd_main,
    const double exp_neg_log_t_rm_adaptation_rate,
    const double target_accept_thresholds,
    const int iteration,
    const int total_burnin
) {
  const int num_vars = observations.n_cols;

  for(int variable = 0; variable < num_vars; variable++) {
    const int num_cats = num_categories(variable);
    if(is_ordinal_variable[variable] == true) {
      const int parameter = -1;                                                 /* unused */
for (int category = 0; category < num_cats; category++) {

  double proposal_sd = proposal_sd_main(variable, category);
  double current = main_effects(variable, category);
  double proposed = R::rnorm(current, proposal_sd);

  // Compute log posterior at proposed state
  main_effects(variable, category) = proposed;
  double log_accept = log_pseudoposterior_thresholds_component (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta, variable, category, parameter
  );

  // Compute log posterior at current state
  main_effects(variable, category) = current;
  log_accept -= log_pseudoposterior_thresholds_component (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta, variable, category, parameter
  );

  if (std::log (R::unif_rand()) < log_accept) {
    main_effects(variable, category) = proposed;
  }

  double updated_proposal_sd = update_proposal_sd_with_robbins_monro (
    proposal_sd, log_accept, exp_neg_log_t_rm_adaptation_rate,
    target_accept_thresholds
  );

  proposal_sd_main(variable, category) = updated_proposal_sd;
}
    } else {
      const int category = -1;                                                  /* unused */
for (int parameter = 0; parameter < 2; parameter++) {

  double proposal_sd = proposal_sd_main(variable, parameter);
  double current = main_effects(variable, parameter);
  double proposed = R::rnorm(current, proposal_sd);

  // Compute log posterior at proposed state
  main_effects(variable, parameter) = proposed;
  double log_accept = log_pseudoposterior_thresholds_component (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta, variable, category, parameter
  );

  // Compute log posterior at current state
  main_effects(variable, category) = current;
  log_accept -= log_pseudoposterior_thresholds_component (
    main_effects, rest_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta, variable, category, parameter
  );

  if (std::log (R::unif_rand()) < log_accept) {
    main_effects(variable, parameter) = proposed;
  }

  if(iteration < total_burnin) {
    double updated_proposal_sd = update_proposal_sd_with_robbins_monro (
      proposal_sd, log_accept, exp_neg_log_t_rm_adaptation_rate,
      target_accept_thresholds
    );

    proposal_sd_main(variable, parameter) = updated_proposal_sd;
  }
}
    }
  }
}



/**
 * Function: gradient_log_pseudoposterior_interactions
 *
 * Computes the gradient of the log pseudoposterior with respect to the
 * active interaction parameters. This is used in MALA updates of the
 * pairwise interaction matrix.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters [V × V].
 *  - main_effects: Matrix of main effect (threshold) parameters [V × max_categories].
 *  - observations: Matrix of ordinal and Blume-Capel scores (row = individual, col = variable).
 *  - num_categories: Vector of number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix indicating which interactions are active.
 *  - is_ordinal_variable: Logical vector indicating ordinal (1) or Blume-Capel (0) variables.
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Cauchy prior scale on interaction weights.
 *
 * Returns:
 *  - Gradient vector for all pairwise interactions (in upper-triangle order).
 */
arma::vec gradient_log_pseudoposterior_interactions (
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
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  arma::umat index_matrix(num_variables, num_variables);
  int counter = 0;
  for(int var1 = 0; var1 < num_variables-1; var1++) {
    for(int var2 = var1 + 1; var2 < num_variables; var2++) {
      index_matrix(var1, var2) = counter;
      counter++;
    }
  }

  arma::vec gradient (num_interactions, arma::fill::zeros);

  for (int var = 0; var < num_variables; var++) {
    int num_cats = num_categories (var);
    arma::mat score_weights (num_observations, num_cats, arma::fill::zeros);    //First column would be zero

    arma::vec rest_scores = rest_matrix.col (var);
    arma::vec denominator = arma::zeros (num_observations);
    arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_cats;

    if (is_ordinal_variable (var)) {
      denominator += arma::exp (-bounds);
      for (int category = 0; category < num_cats; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * rest_scores - bounds;
        arma::vec weight = arma::exp (exponent);
        denominator += weight;
        score_weights.col(category) = (category + 1) * weight;
      }

      // ----- Version with minimal exp calls -----
      // This fails
      // arma::vec exp_bounds = arma::exp(-bounds);
      // arma::vec exp_rest = arma::exp(rest_scores);
      // arma::vec er_accum = arma::ones<arma::vec>(rest_scores.n_elem);  // start at exp(0)
      // denominator += exp_bounds;  // initialize directly
      //
      // for (int category = 0; category < num_cats; category++) {
      //   double main = main_effects(var, category);
      //   double score = category + 1;
      //   er_accum %= exp_rest;  // exp((category + 1) * rest_scores)
      //   arma::vec weight = std::exp(main) * er_accum % exp_bounds;
      //   denominator += weight;
      //   score_weights.col(category) = score * weight;
      // }

    } else {
      const int ref_cat = reference_category (var);
      // Zero category
      double quad_term = main_effects (var, 1) * ref_cat * ref_cat;
      arma::vec exponent = quad_term - bounds;
      denominator = arma::exp (exponent);
      for (int category = 1; category <= num_cats; category++) {
        int centered_cat = category - ref_cat;
        double lin_term = main_effects (var, 0) * category;
        double quad_term = main_effects (var, 1) * centered_cat * centered_cat;
        arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
        arma::vec weight = arma::exp (exponent);
        denominator += weight;
        score_weights.col(category - 1) = category * weight;
      }
    }
    score_weights.each_col() /= denominator;

    for(int var2 = 0; var2 < num_variables; var2++) {
      if (inclusion_indicator (var, var2) == 0)
        continue;

      arma::vec expected_value(num_observations, arma::fill::zeros);
      if(var == var2)
        continue;

      arma::ivec xv = observations.col(var2);

      for (int category = 0; category < num_cats; category++) {
        expected_value += score_weights.col(category) % xv;
      }

      int location = (var < var2) ? index_matrix(var, var2) : index_matrix(var2, var);

      gradient(location) -= arma::accu(expected_value);
    }
  }

  for(int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator (var1, var2) == 0)
        continue;

      int location = index_matrix(var1, var2);

      gradient (location) += 2.0 * sufficient_pairwise(var1, var2);
      // ---- Gradient contribution from Cauchy prior
      const double effect = pairwise_effects (var1, var2);
      gradient (location) -= 2.0 * effect / (effect * effect + interaction_scale * interaction_scale);
    }
  }

  return gradient;
}



/**
 * Function: gradient_log_pseudoposterior_interaction_single
 *
 * Computes the gradient of the log pseudoposterior with respect to a single
 * interaction parameter (i, j), assuming symmetric pairwise_effects matrix.
 *
 * Inputs:
 *  - i, j: Variable indices for the interaction pair (i < j).
 *  - pairwise_effects: Matrix of interaction parameters.
 *  - main_effects: Matrix of main effect parameters.
 *  - observations: Matrix of categorical data [N × V].
 *  - num_categories: Vector with number of categories per variable.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Scale of Cauchy prior.
 *
 * Returns:
 *  - A single double value: the gradient with respect to β_ij.
 */
double gradient_log_pseudoposterior_interaction_single (
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
) {
  const int num_persons = observations.n_rows;
  double gradient = 2.0 * sufficient_pairwise(var1, var2);

  // ---------- Variable 1 ----------
  arma::ivec x = observations.col(var2);
  arma::vec rest = rest_matrix.col(var1);
  int num_cats = num_categories(var1);
  arma::vec denom(num_persons, arma::fill::zeros);
  arma::vec numer(num_persons, arma::fill::zeros);
  arma::vec bounds = arma::clamp(rest, 0.0, arma::datum::inf) * num_cats;

  if (is_ordinal_variable(var1)) {
    denom += arma::exp(-bounds);
    for (int cat = 0; cat < num_cats; cat++) {
      const double meff = main_effects(var1, cat);
      const int score = cat + 1;
      for (int person = 0; person < num_persons; person++) {
        double exp_term = meff + score * rest[person] - bounds[person];
        double w = std::exp(exp_term);
        denom[person] += w;
        numer[person] += score * x[person] * w;
      }
    }
    // Vectorized version (= slower?)
    // denom += arma::exp (-bounds);
    // for (int cat = 0; cat < num_cats; cat++) {
    //   arma::vec exponent = main_effects (var1, cat) + (cat + 1) * rest - bounds;
    //   arma::vec weight = arma::exp (exponent);
    //   denom += weight;
    //   numer += (cat + 1) * x % weight;
    // }

    // ----- Version with minimal exp calls -----
    // denom += arma::exp(-bounds);
    //
    // arma::vec exp_rest = arma::exp(rest);
    // arma::vec exp_bounds = arma::exp(-bounds);
    // arma::vec er_accum = exp_rest;
    //
    // for (int cat = 0; cat < num_cats; cat++) {
    //   const double em = std::exp(main_effects(var1, cat));
    //   const int score = cat + 1;
    //
    //   for (int person = 0; person < num_persons; person++) {
    //     double w = em * er_accum[person] * exp_bounds[person];
    //     denom[person] += w;
    //     numer[person] += score * x[person] * w;
    //
    //     // Accumulate er^(score+1) for next iteration
    //     er_accum[person] *= exp_rest[person];
    //   }
    // }
  } else {
    const int ref = reference_category(var1);
    const double lin_coef = main_effects(var1, 0);
    const double quad_coef = main_effects(var1, 1);

    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;
      double lin_term = lin_coef * cat;
      double quad_term = quad_coef * centered * centered;

      for (int person = 0; person < num_persons; person++) {
        double exp_term = lin_term + quad_term + cat * rest[person] - bounds[person];
        double w = std::exp(exp_term);
        denom[person] += w;
        numer[person] += cat * x[person] * w;
      }
    }
  }

  gradient -= arma::accu(numer / denom);

  // ---------- Variable 2 ----------
  num_cats = num_categories(var2);

  denom.zeros();
  numer.zeros();
  rest = rest_matrix.col(var2);
  x = observations.col(var1);

  bounds = arma::clamp(rest, 0.0, arma::datum::inf) * num_cats;

  if (is_ordinal_variable(var2)) {
    denom += arma::exp(-bounds);
    for (int cat = 0; cat < num_cats; cat++) {
      const double meff = main_effects(var2, cat);
      const int score = cat + 1;
      for (int person = 0; person < num_persons; person++) {
        double exp_term = meff + score * rest[person] - bounds[person];
        double w = std::exp(exp_term);
        denom[person] += w;
        numer[person] += score * x[person] * w;
      }
    }
    // Vectorized version (= slower?)
    // denom += arma::exp (-bounds);
    // for (int cat = 0; cat < num_cats; cat++) {
    //   arma::vec exponent = main_effects (var2, cat) + (cat + 1) * rest - bounds;
    //   arma::vec weight = arma::exp (exponent);
    //   denom += weight;
    //   numer += (cat + 1) * x % weight;
    // }

    // ----- Version with minimal exp calls -----
    // denom += arma::exp(-bounds);
    //
    // arma::vec exp_rest = arma::exp(rest);
    // arma::vec exp_bounds = arma::exp(-bounds);
    // arma::vec er_accum = exp_rest;
    //
    // for (int cat = 0; cat < num_cats; cat++) {
    //   const double em = std::exp(main_effects(var2, cat));
    //   const int score = cat + 1;
    //
    //   for (int person = 0; person < num_persons; person++) {
    //     double w = em * er_accum[person] * exp_bounds[person];
    //     denom[person] += w;
    //     numer[person] += score * x[person] * w;
    //
    //     // Accumulate er^(score+1) for next iteration
    //     er_accum[person] *= exp_rest[person];
    //   }
    // }
  } else {
    const int ref = reference_category(var2);
    const double lin_coef = main_effects(var2, 0);
    const double quad_coef = main_effects(var2, 1);

    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;
      double lin_term = lin_coef * cat;
      double quad_term = quad_coef * centered * centered;

      for (int person = 0; person < num_persons; person++) {
        double exp_term = lin_term + quad_term + cat * rest[person] - bounds[person];
        double w = std::exp(exp_term);
        denom[person] += w;
        numer[person] += cat * x[person] * w;
      }
    }
  }

  gradient -= arma::accu(numer / denom);

  // ---------- Prior ----------
  const double val = pairwise_effects(var1, var2);
  const double val2 = val * val;
  const double scale2 = interaction_scale * interaction_scale;

  gradient -= 2.0 * val / (val2 + scale2);

  return gradient;
}


/**
 * Function: log_pseudoposterior_interactions
 *
 * Computes the log of the pseudoposterior distribution over interaction parameters.
 * Includes:
 *  - The pairwise quadratic sufficient statistic: trace(X * B * Xᵗ)
 *  - The (unnormalized) nodewise log-likelihood contributions
 *  - The log Cauchy prior for included interaction terms
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters.
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of ordinal and BC scores.
 *  - num_categories: Vector of number of categories for each variable.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - is_ordinal_variable: Logical vector: 1 = ordinal, 0 = Blume-Capel.
 *  - reference_category: Reference category index per variable (for BC).
 *  - interaction_scale: Cauchy prior scale for interaction terms.
 *
 * Returns:
 *  - Scalar log pseudoposterior (double)
 */
double log_pseudoposterior_interactions (
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
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;

  // Convert to double matrix for trace calculation
  //arma::mat real_observations = arma::conv_to<arma::mat>::from (observations);

  // Leading term: trace(X * B * X^T)
  //double log_pseudo_posterior = arma::trace (real_observations * pairwise_effects * real_observations.t ());

  double log_pseudo_posterior = arma::accu(pairwise_effects % arma::conv_to<arma::mat>::from(sufficient_pairwise));

  for (int var = 0; var < num_variables; var++) {
    int num_categories_var = num_categories (var);

    // Compute rest score: contribution from other variables
    arma::vec rest_scores = rest_matrix.col (var);
    arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var;
    arma::vec denominator = arma::zeros (num_observations);

    if (is_ordinal_variable (var)) {
      // Ordinal variable: denominator includes exp (-bounds) + exp over categories

      denominator += arma::exp (-bounds);
      for (int category = 0; category < num_categories_var; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * rest_scores - bounds;
        denominator += arma::exp(exponent);
      }

      // ----- Version with minimal exp calls -----
      // This ran into underflow issues due to large bounds when we started with exp_bounds as er_accum
      // That is resolved now, but it still fails
      // arma::vec exp_bounds = arma::exp(-bounds);
      // arma::vec exp_rest = arma::exp(rest_scores);
      // arma::vec er_accum = arma::ones<arma::vec>(num_observations);
      //
      // denominator += exp_bounds;
      //
      // for (int category = 0; category < num_categories_var; category++) {
      //   er_accum %= exp_rest;
      //   double main = main_effects(var, category);
      //   denominator += std::exp(main) * er_accum % exp_bounds;
      // }
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
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator (var1, var2) == 1) {
        log_pseudo_posterior += R::dcauchy (pairwise_effects (var1, var2), 0.0, interaction_scale, true);
      }
    }
  }

  return log_pseudo_posterior;
}



/**
 * Function: find_reasonable_initial_step_size_interactions
 *
 * Heuristically selects an initial step size for MALA updates of the interaction parameters
 * by adapting the acceptance rate toward a target value (~0.574).
 *
 * Inputs:
 *  - pairwise_effects: Matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix of active interactions.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Vector of reference categories (for BC variables).
 *  - interaction_scale: Prior scale for Cauchy prior on interactions.
 *
 * Returns:
 *  - A positive step size value that gives roughly 0.574 acceptance rate.
 */
double find_reasonable_initial_step_size_interactions (
    arma::mat pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const double target_acceptance,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
) {
  const double initial_step_size = 0.1;
  const int max_attempts = 20;

  const int num_variables = pairwise_effects.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  double log_step_size = std::log (initial_step_size);

  // --- Step 1: Extract current interaction state as vector
  arma::vec current_state (num_interactions, arma::fill::zeros);
  int interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator (var1, var2) == 1) {
        current_state (interaction_index) = pairwise_effects (var1, var2);
      }
    }
  }

  // --- Step 2: Evaluate gradient and posterior at current state
  arma::vec gradient = gradient_log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  double log_post_current = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, rest_matrix, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise
  );

  int direction = 0;
  double accept_prob = 0.0;

  // --- Step 3: Exponential step-size search loop
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    double step_size = std::exp (log_step_size);

    // Generate Langevin proposal
    arma::vec noise = arma::randn<arma::vec> (num_interactions);
    arma::vec proposal = current_state +
      0.5 * step_size * gradient + std::sqrt (step_size) * noise;

    // Rebuild symmetric proposal matrix
    arma::mat proposal_matrix = pairwise_effects;
    interaction_index = -1;
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
      for (int var2 = var1 + 1; var2 < num_variables; var2++) {
        interaction_index++;
        if (inclusion_indicator (var1, var2) == 1) {
          proposal_matrix (var1, var2) = proposal (interaction_index);
          proposal_matrix (var2, var1) = proposal (interaction_index);
        }
      }
    }

    arma::mat proposal_rest = observations * proposal_matrix;

    // Evaluate posterior and gradient at proposed state
    double log_post_proposal = log_pseudoposterior_interactions (
      proposal_matrix, main_effects, proposal_rest, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale, sufficient_pairwise
    );

    arma::vec gradient_prop = gradient_log_pseudoposterior_interactions (
      proposal_matrix, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale, sufficient_pairwise, proposal_rest
    );

    // Compute log proposal densities
    arma::vec forward = current_state + 0.5 * step_size * gradient;
    arma::vec backward = proposal + 0.5 * step_size * gradient_prop;

    double log_q_forward = -0.5 / step_size * arma::accu (arma::square (proposal - forward));
    double log_q_reverse = -0.5 / step_size * arma::accu (arma::square (current_state - backward));

    double log_accept = log_post_proposal + log_q_reverse - log_post_current - log_q_forward;
    accept_prob = std::min (1.0, std::exp (log_accept));

    // --- Step 4: Decide direction based on first attempt
    if (attempt == 0) {
      direction = (accept_prob > target_acceptance) ? 1 : -1;
    } else {
      if ((direction == 1 && accept_prob < target_acceptance) ||
          (direction == -1 && accept_prob > target_acceptance)) {
        break;
      }
    }

    // Adjust log step size in chosen direction
    log_step_size += direction;
  }

  return std::exp (log_step_size);
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



/**
 * Function: update_interactions_with_adaptive_metropolis
 *
 * Performs adaptive Metropolis-Hastings updates for all active pairwise interactions.
 *
 * For each pair (i, j) with inclusion_indicator(i, j) == 1, proposes a new interaction strength.
 * If accepted, updates the interaction matrix and residual matrix.
 *
 * Inputs:
 *  - pairwise_effects: Current matrix of interaction parameters (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd_pairwise_effects: Matrix of proposal standard deviations (updated in-place).
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - rest_matrix: Matrix of residual scores (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Modifies:
 *  - pairwise_effects
 *  - rest_matrix
 *  - proposal_sd_pairwise_effects
 */
void update_interactions_with_adaptive_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    arma::mat& proposal_sd_pairwise_effects,
    const double interaction_scale,
    const int num_persons,
    const int num_variables,
    arma::mat& rest_matrix,
    const double exp_neg_log_t_rm_adaptation_rate,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double target_accept_interactions,
    const int iteration,
    const int total_burnin,
    const arma::imat& sufficient_pairwise
) {
  for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
    for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        // Sample proposal from Gaussian centered at current state
        const double current_state = pairwise_effects(variable1, variable2);
        const double proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(variable1, variable2));

        // Compute log acceptance ratio: data + prior
        double log_acceptance = log_pseudolikelihood_ratio_interaction(
          pairwise_effects, main_effects, observations, num_categories, num_persons,
          variable1, variable2, proposed_state, current_state,
          rest_matrix, is_ordinal_variable, reference_category, sufficient_pairwise
        );

        // Add symmetric Cauchy prior ratio
        log_acceptance += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_acceptance -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        // Accept proposal with MH step
        if (std::log (R::unif_rand()) < log_acceptance) {
          const double delta = proposed_state - current_state;

          // Update effect matrix symmetrically
          pairwise_effects(variable1, variable2) = proposed_state;
          pairwise_effects(variable2, variable1) = proposed_state;

          // Vectorized update of residual matrix for both variables
          rest_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
          rest_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
        }

        if(iteration < total_burnin) {
          // Robbins-Monro adaptation of proposal SD
          proposal_sd_pairwise_effects(variable1, variable2) =
            update_proposal_sd_with_robbins_monro (
                proposal_sd_pairwise_effects(variable1, variable2), log_acceptance,
                exp_neg_log_t_rm_adaptation_rate,
                target_accept_interactions
            );
        }
      }
    }
  }
}



/**
 * Function: update_indicator_interaction_pair_with_metropolis
 *
 * Metropolis-Hastings update of pairwise inclusion indicators for a predefined set of edges.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of interaction weights (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - indicator: Matrix of edge inclusion flags (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd: Matrix of proposal standard deviations for pairwise effects.
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - index: List of interaction pairs to update.
 *  - num_interactions: Number of interaction pairs.
 *  - num_persons: Number of observations.
 *  - rest_matrix: Residual scores matrix (updated in-place).
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - is_ordinal_variable: Logical vector indicating variable type.
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Modifies:
 *  - indicator
 *  - pairwise_effects
 *  - rest_matrix
 */
void update_indicator_interaction_pair_with_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const arma::imat& sufficient_pairwise
) {
  for (int cntr = 0; cntr < num_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    const double current_state = pairwise_effects(variable1, variable2);

    // Propose a new state: either add a new edge or remove an existing one
    const bool proposing_addition = (indicator(variable1, variable2) == 0);
    const double proposed_state = proposing_addition ? R::rnorm(current_state, proposal_sd(variable1, variable2)) : 0.0;

    // Compute log pseudo-likelihood ratio
    double log_accept = log_pseudolikelihood_ratio_interaction (
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      variable1, variable2, proposed_state, current_state, rest_matrix,
      is_ordinal_variable, reference_category, sufficient_pairwise
    );

    // Add prior ratio and proposal correction
    const double inclusion_probability_ij = inclusion_probability(variable1, variable2);
    const double sd = proposal_sd(variable1, variable2);

    if (proposing_addition) {
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept -= R::dnorm(proposed_state, current_state, sd, true);
      log_accept += std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    } else {
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_accept += R::dnorm(current_state, proposed_state, sd, true);
      log_accept -= std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    }

    // Metropolis-Hastings accept step
    if (std::log (R::unif_rand()) < log_accept) {
      const int updated_indicator = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = updated_indicator;
      indicator(variable2, variable1) = updated_indicator;

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double delta = proposed_state - current_state;

      // Vectorized residual update
      rest_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      rest_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}



/**
 * Function: update_interactions_with_fisher_mala
 *
 * Updates all included pairwise interaction effects using Fisher-preconditioned MALA.
 *
 * This function performs a joint update of the interaction weights using a Fisher-preconditioned
 * Metropolis-adjusted Langevin algorithm (MALA). The inverse Fisher matrix is either identity
 * (during burn-in) or estimated and adapted post-burn-in using rank-one score difference updates.
 *
 * The step size is normalized using the trace of the inverse Fisher matrix (Titsias, 2023), ensuring
 * consistent average proposal scale.
 *
 * During burn-in, the step size is adapted using dual averaging. After burn-in, Robbins-Monro is used
 * along with Fisher matrix updates.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of interaction parameters (symmetric, updated in-place).
 *  - rest_matrix: Person-by-variable linear predictor matrix (updated if proposal is accepted).
 *  - main_effects: Matrix of current threshold (main effect) parameters.
 *  - observations: Person-by-variable observation matrix.
 *  - num_categories: Number of categories per variable.
 *  - inclusion_indicator: Binary indicator matrix showing active interaction pairs.
 *  - is_ordinal_variable: Flags indicating ordinal variables.
 *  - reference_category: Reference category index per variable.
 *  - interaction_scale: Scale parameter for the Cauchy prior on interaction weights.
 *  - step_size_pairwise: Current MALA step size (updated in-place).
 *  - initial_step_size_pairwise: Initial step size for dual averaging.
 *  - iteration: Current MCMC iteration.
 *  - total_burnin: Total burn-in length.
 *  - dual_averaging_state: Dual averaging parameters (updated during burn-in).
 *  - sqrt_inv_fisher_pairwise: Square root of inverse Fisher information matrix (updated post-burn-in).
 *
 * Modifies:
 *  - pairwise_effects (on accept)
 *  - rest_matrix (on accept)
 *  - step_size_pairwise (if adaptive)
 *  - dual_averaging_state (during burn-in)
 *  - sqrt_inv_fisher_pairwise (post-burn-in)
 */
void update_interactions_with_fisher_mala (
    arma::mat& pairwise_effects,
    arma::mat& rest_matrix,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    double& step_size_pairwise,
    const double initial_step_size_pairwise,
    const int iteration,
    arma::vec& dual_averaging_state,
    arma::mat& sqrt_inv_fisher_pairwise,
    const double target_accept_interactions,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end,
    const arma::imat& sufficient_pairwise
) {
  const int num_variables = pairwise_effects.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  // --- Flatten current pairwise interactions into vector format
  arma::vec current_state(num_interactions, arma::fill::zeros);
  int interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator(var1, var2) == 1) {
        current_state(interaction_index) = pairwise_effects(var1, var2);
      }
    }
  }

  // --- Compute gradient and log-posterior at current state
  arma::vec current_grad = gradient_log_pseudoposterior_interactions(
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  double current_log_post = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, rest_matrix, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise
  );

  // --- Initialize Fisher matrix after burn-in
  if (iteration == warmup_stageI_end) {
    sqrt_inv_fisher_pairwise = initialize_fisher_preconditioner(current_grad);
  }

  const int dim = current_state.n_elem;

  // --- Construct inverse Fisher matrix (identity during burn-in)
  arma::mat inv_fisher;
  if (iteration >= warmup_stageI_end) {
    inv_fisher = sqrt_inv_fisher_pairwise * sqrt_inv_fisher_pairwise.t();
  } else {
    inv_fisher.eye(dim, dim);
  }

  // --- Rescale step size based on trace-normalization
  double trace_inv_fisher = arma::trace(inv_fisher);
  double scaled_step_size = step_size_pairwise / (trace_inv_fisher / dim);
  const double sqrt_step = std::sqrt(scaled_step_size);

  // --- Compute proposal: drift + noise
  const arma::vec proposal_drift = 0.5 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec noise = sqrt_inv_fisher_pairwise * arma::randn(dim);
  const arma::vec proposed_state = current_state + proposal_drift + sqrt_step * noise;

  // --- Reconstruct proposed interaction matrix
  arma::mat proposed_matrix = pairwise_effects;
  interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator(var1, var2) == 1) {
        proposed_matrix(var1, var2) = proposed_state(interaction_index);
        proposed_matrix(var2, var1) = proposed_state(interaction_index);
      }
    }
  }

  // --- Compute gradient and posterior at proposed state
  arma::mat proposed_rest = observations * proposed_matrix;

  double proposed_log_post = log_pseudoposterior_interactions(
    proposed_matrix, main_effects, proposed_rest, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise
  );

  arma::vec proposed_grad = gradient_log_pseudoposterior_interactions(
    proposed_matrix, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale, sufficient_pairwise, proposed_rest
  );

  // --- Compute forward/reverse correction terms (Titsias 2023, Prop. 1)
  const arma::vec forward_proposal_residual =
    proposed_state - current_state -
    0.25 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec reverse_proposal_residual =
    current_state - proposed_state -
    0.25 * scaled_step_size * inv_fisher * proposed_grad;

  const double log_forward = 0.5 * arma::dot (forward_proposal_residual, current_grad);
  const double log_backward = 0.5 * arma::dot (reverse_proposal_residual, proposed_grad);

  const double log_accept = proposed_log_post - current_log_post + (log_backward - log_forward);
  const double accept_prob = std::min(1.0, std::exp (log_accept));

  // --- Accept/reject step
  if (std::log(R::unif_rand()) < log_accept) {
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
      for (int var2 = var1 + 1; var2 < num_variables; var2++) {
        if (inclusion_indicator (var1, var2) == 1) {
          double delta = proposed_matrix (var1, var2) - pairwise_effects (var1, var2);
          pairwise_effects (var1, var2) = proposed_matrix (var1, var2);
          pairwise_effects (var2, var1) = proposed_matrix (var1, var2);
          rest_matrix.col (var1) += arma::conv_to<arma::vec>::from (observations.col (var2)) * delta;
          rest_matrix.col (var2) += arma::conv_to<arma::vec>::from (observations.col (var1)) * delta;
        }
      }
    }
  }

  // --- Step size and Fisher adaptation
  if (iteration < warmup_stageI_end) {
    // Stage I warmup: Use dual averaging
    update_step_size_with_dual_averaging (
        initial_step_size_pairwise, accept_prob, iteration,
        dual_averaging_state, target_accept_interactions
    );
    step_size_pairwise = std::exp (dual_averaging_state[0]);
  } else if (iteration < warmup_stageII_end) {
    // Stage II warmup: Robbins-Monro and Fisher update
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageI_end + 1, step_size_pairwise,
        target_accept_interactions
    );

    // Scaled score difference for rank-1 Fisher update
    arma::vec score_diff = std::sqrt (accept_prob) * (proposed_grad - current_grad);
    update_fisher_preconditioner(
      sqrt_inv_fisher_pairwise, score_diff
    );
  } else if (iteration < warmup_stageIII_end) {
    // Stage III warmup: Robbins-Monro
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageII_end + 1, step_size_pairwise,
        target_accept_interactions
    );
  }
}



/**
 * Function: update_indicator_interaction_pair_with_mala
 *
 * Updates the inclusion indicators and associated interaction weights using MALA.
 *
 * This function iterates over all candidate interaction pairs and proposes either:
 *   - The addition of an edge with a Fisher-preconditioned Langevin step
 *   - The removal of an existing edge by proposing a zero interaction
 *
 * Proposals are evaluated using a Metropolis-Hastings step that includes:
 *   - Pseudolikelihood difference (pairwise-only)
 *   - Cauchy prior for interaction weights
 *   - Bernoulli prior for inclusion
 *   - Langevin forward/reverse proposal densities
 *
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters (modified in-place).
 *  - main_effects: Matrix of threshold (main effect) parameters.
 *  - indicator: Symmetric matrix of inclusion indicators (modified in-place).
 *  - observations: Person-by-variable matrix of observed scores.
 *  - num_categories: Number of categories per variable.
 *  - step_size_pairwise: Global MALA step size (scaled internally).
 *  - interaction_scale: Scale parameter for Cauchy prior on interaction weights.
 *  - index: Matrix listing candidate interactions: [index, var1, var2].
 *  - num_persons: Number of observations (rows in observations matrix).
 *  - rest_matrix: Linear predictor matrix (updated in-place on accept).
 *  - inclusion_probability: Prior inclusion probabilities for interaction pairs.
 *  - is_ordinal_variable: Indicator of ordinal variables.
 *  - reference_category: Reference category for each variable.
 *  - num_pairwise: Total number of candidate interactions.
 *  - iteration: Current MCMC iteration.
 *  - total_burnin: Number of warm-up iterations.
 *
 * Modifies:
 *  - indicator
 *  - pairwise_effects (on accept)
 *  - rest_matrix (on accept)
 */
void update_indicator_interaction_pair_with_mala (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double step_size,
    const double interaction_scale,
    const arma::imat& index,
    const int num_persons,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const int num_pairwise,
    const arma::imat& sufficient_pairwise
) {
  const double sd = std::sqrt(step_size);

  for (int pair_index = 0; pair_index < num_pairwise; ++pair_index) {
    const int var1 = index(pair_index, 1);
    const int var2 = index(pair_index, 2);

    const double current_state = pairwise_effects(var1, var2);
    const double inclusion_prob = inclusion_probability(var1, var2);

    double proposed_state = 0.0;
    double log_accept = 0.0;

    if (indicator(var1, var2) == 0) {
      const double grad = gradient_log_pseudoposterior_interaction_single(
        var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      const double drift = 0.5 * step_size * grad;
      const double forward_mean = current_state + drift;

      proposed_state = forward_mean + R::rnorm(0.0, sd);
      log_accept -= R::dnorm(proposed_state, forward_mean, sd, true);
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept += std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);
    } else {
      // Change in place
      const double proposed_state = 0.0;

      const double tmp = pairwise_effects(var1, var2);  // cache for restore
      pairwise_effects(var1, var2) = proposed_state;
      pairwise_effects(var2, var1) = proposed_state;

      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      double delta = proposed_state - current_state;

      rest_matrix.col(var1) += obs_var2 * delta;
      rest_matrix.col(var2) += obs_var1 * delta;

      const double grad = gradient_log_pseudoposterior_interaction_single(
        var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      const double drift = 0.5 * step_size * grad;
      const double backward_mean = proposed_state + drift;
      log_accept += R::dnorm(current_state, backward_mean, sd, true);
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_accept -= std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);

      // Restore
      pairwise_effects(var1, var2) = tmp;
      pairwise_effects(var2, var1) = tmp;
      rest_matrix.col(var1) -= obs_var2 * delta;
      rest_matrix.col(var2) -= obs_var1 * delta;
    }

    log_accept += log_pseudolikelihood_ratio_interaction(
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      var1, var2, proposed_state, current_state, rest_matrix,
      is_ordinal_variable, reference_category, sufficient_pairwise
    );

    if (std::log(R::unif_rand()) < log_accept) {
      const int new_value = 1 - indicator(var1, var2);
      indicator(var1, var2) = new_value;
      indicator(var2, var1) = new_value;

      const double delta = proposed_state - current_state;
      pairwise_effects(var1, var2) = proposed_state;
      pairwise_effects(var2, var1) = proposed_state;

      // Pre-fetch as double vectors
      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      rest_matrix.col(var1) += obs_var2 * delta;
      rest_matrix.col(var2) += obs_var1 * delta;
    }
  }
}




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
    const double threshold_alpha,
    const double threshold_beta,
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
      index_matrix(var1, var2) = counter;
      counter++;
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
      arma::vec rest_score = rest_matrix.col(variable);
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

      arma::ivec lin_score = arma::conv_to<arma::ivec>::from(arma::regspace<arma::vec>(0, num_cats));
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
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }
      offset += num_cats;
    } else {
      for (int i = 0; i < 2; i++) {
        const double threshold_param = main_effects(variable, i);
        const double p = 1.0 / (1.0 + std::exp (-threshold_param));
        gradient(offset + i) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
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
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    const arma::mat& rest_matrix
) {

  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;

  double log_pseudoposterior = 0.0;

  // Calculate the contribution from the data and the prior
  auto log_beta_prior = [&](double threshold_param) {
    return threshold_param * threshold_alpha - std::log1p (std::exp (threshold_param)) * (threshold_alpha + threshold_beta);
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

double find_reasonable_initial_step_size (
    const arma::mat main_effects,
    const arma::mat pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const double target_acceptance,
    const arma::imat& sufficient_pairwise
) {
  constexpr double initial_step_size = 0.1;
  constexpr int max_attempts = 20;
  constexpr double max_log_step = 10.0;

  double log_step_size = std::log (initial_step_size);

  // Current state and log posterior
  const arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator, num_categories,
    is_ordinal_variable);

  const arma::vec current_grad = gradient_log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator,
    observations, num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  const double current_log_post = log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  int direction = 0;
  double accept_prob = 0.0;

  // Propose initial MALA step: θ' = θ + ½ε ∇log p(θ) + √ε * N(0, I)

  double step_size = std::exp (log_step_size);
  double sqrt_step = std::sqrt(step_size);
  arma::vec proposed_state = current_state + 0.5 * step_size * current_grad +
    sqrt_step * arma::randn(current_state.n_elem);

  arma::mat main_proposed = arma::mat(main_effects.n_rows, main_effects.n_cols,
                                      arma::fill::zeros);
  arma::mat pairwise_proposed = arma::mat(pairwise_effects.n_rows, pairwise_effects.n_cols,
                                          arma::fill::zeros);
  unvectorize_model_parameters(
    proposed_state, main_proposed, pairwise_proposed, num_categories,
    is_ordinal_variable
  );

  arma::mat proposed_rest = observations * pairwise_proposed;

  double proposed_log_post = log_pseudoposterior (
    main_proposed, pairwise_proposed, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, proposed_rest
  );

  arma::vec proposed_grad = gradient_log_pseudoposterior (
    main_proposed, pairwise_proposed, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, proposed_rest
  );

  arma::vec forward = current_state + 0.5 * step_size * current_grad;
  arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

  double log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
  double log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
  double log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

  accept_prob = std::min(1.0, std::exp (log_accept));
  if (std::abs(accept_prob - target_acceptance) < 0.1) {
    return std::exp (log_step_size);
  }
  direction = (accept_prob > target_acceptance) ? 1 : -1;


  // Log-scale search for step size that brackets the target acceptance rate
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    log_step_size += direction;
    step_size = std::exp (log_step_size);
    sqrt_step = std::sqrt(step_size);

    proposed_state = current_state + 0.5 * step_size * current_grad +
      sqrt_step * arma::randn(current_state.n_elem);

    unvectorize_model_parameters(
      proposed_state, main_proposed, pairwise_proposed, num_categories,
      is_ordinal_variable
    );

    proposed_rest = observations * pairwise_proposed;

    proposed_log_post = log_pseudoposterior (
      main_proposed, pairwise_proposed, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, proposed_rest
    );

    proposed_grad = gradient_log_pseudoposterior (
      main_proposed, pairwise_proposed, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, proposed_rest
    );

    arma::vec forward = current_state + 0.5 * step_size * current_grad;
    arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

    log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
    log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
    log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

    double new_accept_prob = std::min(1.0, std::exp (log_accept));

    // Exit if acceptance flips across the target
    if ((direction == 1 && new_accept_prob < target_acceptance) ||
        (direction == -1 && new_accept_prob > target_acceptance)) {
      break;
    }

    accept_prob = new_accept_prob;

    if (std::abs(log_step_size) > max_log_step) {
      Rcpp::Rcout << "Warning: Step size search failed. Falling back to default (0.01)." << std::endl;
      return 0.01;
    }
  }

  return std::exp (log_step_size);
}

void update_parameters_with_fisher_mala (
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
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
    double& step_size,
    const int iteration,
    arma::vec& dual_averaging_state,
    arma::mat& sqrt_inv_fisher,
    const double initial_step_size,
    const double target_accept,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end,
    const arma::imat& sufficient_pairwise,
    arma::mat& rest_matrix
) {
  // --- Compute current parameter vector and its gradient ---
  const arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator, num_categories,
    is_ordinal_variable);

  const arma::vec current_grad = gradient_log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator,
    observations, num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  const double current_log_post = log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  const int dim = current_state.n_elem;

  // --- If burn-in just ended, initialize Fisher matrix using gradient
  if (iteration == warmup_stageI_end) {
    sqrt_inv_fisher = initialize_fisher_preconditioner (current_grad);
  }

  // --- Set inverse Fisher matrix: identity during warm-up, adapted afterward
  arma::mat inv_fisher;
  if (iteration >= warmup_stageI_end) {
    inv_fisher = sqrt_inv_fisher * sqrt_inv_fisher.t();
  } else {
    inv_fisher.eye(dim, dim);
  }

  // --- Construct Fisher-preconditioned MALA proposal ---
  const double trace_inv_fisher = arma::trace(inv_fisher);
  const double scaled_step_size = step_size / (trace_inv_fisher / dim);
  const double sqrt_step = std::sqrt(scaled_step_size);

  // Drift (mean shift) and stochastic noise
  const arma::vec proposal_drift = 0.5 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec noise = sqrt_inv_fisher * arma::randn(dim);

  const arma::vec proposed_state = current_state + proposal_drift + sqrt_step * noise;

  arma::mat main_proposed = arma::mat(
    main_effects.n_rows, main_effects.n_cols, arma::fill::zeros
  );
  arma::mat pairwise_proposed = arma::mat(
    pairwise_effects.n_rows, pairwise_effects.n_cols, arma::fill::zeros
  );
  unvectorize_model_parameters(
    proposed_state, main_proposed, pairwise_proposed, num_categories,
    is_ordinal_variable
  );

  arma::mat proposed_rest = observations * pairwise_proposed;

  const double proposed_log_post = log_pseudoposterior (
    main_proposed, pairwise_proposed, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, proposed_rest
  );

  const arma::vec proposed_grad = gradient_log_pseudoposterior (
    main_proposed, pairwise_proposed, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, proposed_rest
  );

  // --- Compute MALA acceptance correction (Titsias 2023, Prop. 1) ---
  const arma::vec forward_proposal_residual =
    proposed_state - current_state -
    0.25 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec reverse_proposal_residual =
    current_state - proposed_state -
    0.25 * scaled_step_size * inv_fisher * proposed_grad;

  const double log_forward = 0.5 * arma::dot (forward_proposal_residual, current_grad);
  const double log_backward = 0.5 * arma::dot (reverse_proposal_residual, proposed_grad);

  double log_accept = proposed_log_post - current_log_post;
  log_accept += (log_backward - log_forward);

  const double accept_prob = std::min(1.0, std::exp (log_accept));

  // --- Accept or reject proposed move ---
  if (std::log (R::unif_rand()) < log_accept) {
    main_effects = main_proposed;
    pairwise_effects = pairwise_proposed;
    rest_matrix = proposed_rest;
  }

  // --- Update step size and Fisher matrix ---
  if (iteration < warmup_stageI_end) {
    // During warm-up stage I: dual averaging adaptation
    update_step_size_with_dual_averaging (
        initial_step_size, accept_prob, iteration + 1, dual_averaging_state,
        target_accept
    );
    step_size = std::exp (dual_averaging_state[1]);
  } else if (iteration < warmup_stageII_end) {

    // After warm-up stage I: Robbins-Monro + Fisher preconditioner update
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageI_end + 1, step_size,
        target_accept
    );

    const arma::vec score_diff =
      std::sqrt(accept_prob) * (proposed_grad - current_grad);

    update_fisher_preconditioner(sqrt_inv_fisher, score_diff);
  } else if (iteration < warmup_stageIII_end) {
    // After warm-up stages I & II: Robbins-Monro
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageII_end + 1, step_size,
        target_accept
    );
  }
}





inline void update_diag_mass_matrix(
    int count,
    arma::vec& mean,
    arma::vec& M2,
    const arma::vec& x
) {
  arma::vec delta = x - mean;
  mean += delta / count;
  M2 += delta % (x - mean);  // element-wise square
}

double find_reasonable_initial_step_size_hmc(
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
    const double target_acceptance,
    const arma::imat& sufficient_pairwise,
    arma::mat& rest_matrix,
    arma::vec& mass_inv) {

  constexpr double initial_step_size = 0.1;
  constexpr int max_attempts = 20;
  constexpr double max_log_step = 10.0;

  arma::vec theta = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator, num_categories,
    is_ordinal_variable
  );

  double log_step_size = std::log(initial_step_size);

  // Compute potential and gradient at current state
  double current_U = -log_pseudoposterior(
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  arma::vec grad = -gradient_log_pseudoposterior(
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  // Sample initial momentum
  arma::vec p = arma::randn(theta.n_elem) / arma::sqrt(mass_inv);

  // Start with 1 leapfrog step
  double accept_prob = 0.0;
  int direction = 0;

  for (int attempt = 0; attempt < max_attempts; attempt++) {
    double step_size = std::exp(log_step_size);

    arma::vec theta_prop = theta;
    arma::vec p_prop = p;
    arma::vec grad_prop = grad;

    // Half step for momentum
    p_prop -= 0.5 * step_size * grad_prop;

    // Full step for position
    theta_prop += step_size * (mass_inv % p_prop);

    // Unvectorize proposed state
    arma::mat main_prop = arma::zeros<arma::mat>(main_effects.n_rows, main_effects.n_cols);
    arma::mat pairwise_prop = arma::zeros<arma::mat>(pairwise_effects.n_rows, pairwise_effects.n_cols);
    unvectorize_model_parameters(theta_prop, main_prop, pairwise_prop, num_categories, is_ordinal_variable);
    arma::mat rest_prop = observations * pairwise_prop;

    // Compute new gradient
    grad_prop = -gradient_log_pseudoposterior(
      main_prop, pairwise_prop, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rest_prop
    );

    // Final half-step for momentum
    p_prop -= 0.5 * step_size * grad_prop;

    // Compute proposed Hamiltonian
    double proposed_U = -log_pseudoposterior(
      main_prop, pairwise_prop, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rest_prop
    );
    double current_K = 0.5 * arma::dot(p % mass_inv, p);
    double proposed_K = 0.5 * arma::dot(p_prop % mass_inv, p_prop);

    double log_accept = (current_U + current_K) - (proposed_U + proposed_K);
    accept_prob = std::min(1.0, std::exp(log_accept));

    if (std::abs(accept_prob - target_acceptance) < 0.1) {
      return step_size;
    }

    direction = (accept_prob > target_acceptance) ? 1 : -1;
    log_step_size += direction;

    if (std::abs(log_step_size) > max_log_step) {
      Rcpp::Rcout << "Warning: Failed to find reasonable initial step size. Falling back to 0.01.\n";
      return 0.01;
    }
  }

  return std::exp(log_step_size);
}

void update_parameters_with_hmc(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
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
    double& step_size,
    const int iteration,
    arma::vec& dual_averaging_state,
    const double initial_step_size,
    const double target_accept,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end,
    const arma::imat& sufficient_pairwise,
    arma::mat& rest_matrix,
    const int L,
    arma::vec& mass_inv,
    arma::vec mean,
    arma::vec M2
) {
  arma::mat main_proposed = arma::mat(
    main_effects.n_rows, main_effects.n_cols, arma::fill::zeros
  );
  arma::mat pairwise_proposed = arma::mat(
    pairwise_effects.n_rows, pairwise_effects.n_cols, arma::fill::zeros
  );
  arma::mat rest_proposed;

  // Deep copy of current theta
  arma::vec theta = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator, num_categories,
    is_ordinal_variable);

  // Sample momentum p ~ N(0, M)
  arma::vec p = arma::randn(theta.n_elem) / arma::sqrt(mass_inv);
  arma::vec p_new = p; // Copy for later

  // Compute gradient of potential energy (negative log posterior)
  arma::vec grad = -gradient_log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );

  // --- Leapfrog integration ---
  p_new -= 0.5 * step_size * grad;

  for (int i = 0; i < L; i++) {
    theta += step_size * (mass_inv % p_new); // Full step for position

    // Recompute gradient
    unvectorize_model_parameters(
      theta, main_proposed, pairwise_proposed, num_categories,
      is_ordinal_variable
    );
    rest_proposed = observations * pairwise_proposed;

    grad = -gradient_log_pseudoposterior (
      main_proposed, pairwise_proposed, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rest_proposed
    );

    if (i != L - 1) {
      p_new -= step_size * grad; // Full step for momentum (except last)
    }
  }
  p_new -= 0.5 * step_size * grad; // Final half-step for momentum

  // Compute Hamiltonians
  double current_U = -log_pseudoposterior (
    main_effects, pairwise_effects, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_matrix
  );
  double current_K = 0.5 * arma::dot(p % mass_inv, p);
  double proposed_U = -log_pseudoposterior (
    main_proposed, pairwise_proposed, inclusion_indicator, observations,
    num_categories, num_obs_categories, sufficient_blume_capel,
    reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
    interaction_scale, sufficient_pairwise, rest_proposed
  );
  double proposed_K = 0.5 * arma::dot(p_new % mass_inv, p_new);

  double log_accept_prob = (current_U + current_K) - (proposed_U + proposed_K);
  if (std::log(R::runif(0, 1)) < log_accept_prob) {
    main_effects = main_proposed;
    pairwise_effects = pairwise_proposed;
    rest_matrix = rest_proposed;
  }

  // --- Update step size and Mass matrix ---
  const double accept_prob = std::min(1.0, std::exp (log_accept_prob));
  if (iteration < warmup_stageI_end) {
    // During warm-up stage I: dual averaging adaptation
    update_step_size_with_dual_averaging (
        initial_step_size, accept_prob, iteration + 1, dual_averaging_state,
        target_accept
    );
    step_size = std::exp (dual_averaging_state[1]);
  } else if (iteration < warmup_stageII_end) {

    // After warm-up stage I: Robbins-Monro + Fisher preconditioner update
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageI_end + 1, step_size,
        target_accept
    );

    // At each Stage II iteration:
    int stageII_iter = iteration - warmup_stageI_end + 2;
    update_diag_mass_matrix(stageII_iter, mean, M2, theta);
    //    mass_inv = 1.0 / (M2 / (stageII_iter - 1));

  } else if (iteration < warmup_stageIII_end) {
    // After warm-up stages I & II: Robbins-Monro
    update_step_size_with_robbins_monro (
        accept_prob, iteration - warmup_stageII_end + 1, step_size,
        target_accept
    );
  }
}

void update_indicator_interaction_pair_with_hmc(
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double step_size,
    const double interaction_scale,
    const arma::imat& index,
    const int num_persons,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const int num_pairwise,
    const arma::imat& sufficient_pairwise,
    const arma::vec& mass_inv,
    const int L) {

  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);

  for (int pair_index = 0; pair_index < num_pairwise; ++pair_index) {
    const int var1 = index(pair_index, 1);
    const int var2 = index(pair_index, 2);
    const int param_index = num_main + index(pair_index, 0);

    const double current_state = pairwise_effects(var1, var2);
    const double inclusion_prob = inclusion_probability(var1, var2);

    double proposed_state = 0.0;
    double log_accept = 0.0;

    // --- HMC step on proposed_state ---
    double theta = current_state;
    double grad = gradient_log_pseudoposterior_interaction_single(
      var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
      num_categories, is_ordinal_variable, reference_category,
      interaction_scale, sufficient_pairwise);

    double p = R::rnorm(0.0, std::sqrt(1.0 / mass_inv[param_index]));
    double p_prop = p;

    // Leapfrog integration
    p_prop -= 0.5 * step_size * grad;
    for (int i = 0; i < L; ++i) {
      theta += step_size * (mass_inv[param_index] * p_prop);

      pairwise_effects(var1, var2) = theta;
      pairwise_effects(var2, var1) = theta;

      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      double delta = theta - current_state;
      rest_matrix.col(var1) += obs_var2 * delta;
      rest_matrix.col(var2) += obs_var1 * delta;

      grad = gradient_log_pseudoposterior_interaction_single(
        var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise);

      if (i != L - 1) {
        p_prop -= step_size * grad;
      }
    }
    p_prop -= 0.5 * step_size * grad;

    proposed_state = theta;

    double current_K = 0.5 * mass_inv[param_index] * p * p;
    double proposed_K = 0.5 * mass_inv[param_index] * p_prop * p_prop;

    double log_prior_ratio = R::dcauchy(proposed_state, 0.0, interaction_scale, true)
      - R::dcauchy(current_state, 0.0, interaction_scale, true);

    double log_indicator_ratio = std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);

    log_accept = -proposed_K + current_K + log_prior_ratio + log_indicator_ratio;

    log_accept += log_pseudolikelihood_ratio_interaction(
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      var1, var2, proposed_state, current_state, rest_matrix,
      is_ordinal_variable, reference_category, sufficient_pairwise);

    if (std::log(R::unif_rand()) < log_accept) {
      const int new_value = 1 - indicator(var1, var2);
      indicator(var1, var2) = new_value;
      indicator(var2, var1) = new_value;
      pairwise_effects(var1, var2) = proposed_state;
      pairwise_effects(var2, var1) = proposed_state;
    } else {
      indicator(var1, var2) = indicator(var1, var2);
      indicator(var2, var1) = indicator(var2, var1);
      pairwise_effects(var1, var2) = current_state;
      pairwise_effects(var2, var1) = current_state;

      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      double delta = proposed_state - current_state;
      rest_matrix.col(var1) -= obs_var2 * delta;
      rest_matrix.col(var2) -= obs_var1 * delta;
    }
  }
}


// Algorithm 4 in // "NUTS: The No-U-Turn Sampler" by Hoffman and Gelman (2014)
double find_reasonable_initial_step_size_nuts(
    arma::mat main_effects,
    arma::mat pairwise_effects,
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
    const double target_acceptance,
    const arma::imat& sufficient_pairwise
) {
  constexpr double initial_step_size = 1;
  constexpr int max_attempts = 20;
  double eps = initial_step_size;

  auto logp_and_grad = [&](const arma::vec& theta_vec) {
    arma::mat m, p;
    unvectorize_model_parameters(theta_vec, m, p, num_categories, is_ordinal_variable);
    arma::mat rest = observations * p;
    double logp = log_pseudoposterior(m, p, inclusion_indicator, observations,
                                      num_categories, num_obs_categories, sufficient_blume_capel,
                                      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
                                      interaction_scale, sufficient_pairwise, rest);
    arma::vec grad = gradient_log_pseudoposterior(m, p, inclusion_indicator, observations,
                                                  num_categories, num_obs_categories, sufficient_blume_capel,
                                                  reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
                                                  interaction_scale, sufficient_pairwise, rest);
    return std::make_pair(logp, (-grad).eval());
  };


  arma::vec theta = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );
  arma::vec r = 0.01 * arma::randn(theta.n_elem);
  auto [logp0, dummy0] = logp_and_grad(theta);
  double kin0 = 0.5 * arma::dot(r, r);

  arma::vec theta_new, r_new;
  std::tie(theta_new, r_new) = leapfrog(
    theta, r, eps, logp_and_grad
  );
  auto [logp, dummy1] = logp_and_grad(theta_new);
  double kin = 0.5 * arma::dot(r_new, r_new);

  double H0 = logp0 - kin0;
  double H = logp - kin;

  int a = 2 * (H - H0 > std::log(0.5)) - 1;

  int attempts = 0;
  while(a * (H - H0) > - a * log(2.0) && attempts < max_attempts) {
    if(a == 1) {
      eps *= 2.0;
    } else {
      eps /= 2.0;
    }

    std::tie(theta_new, r_new) = leapfrog(
      theta, r, eps, logp_and_grad
    );
    std::tie(logp, std::ignore) = logp_and_grad(theta_new);
    kin = 0.5 * arma::dot(r_new, r_new);
    H = logp - kin;
    attempts++;
  }
  return eps;
}


void update_parameters_with_nuts(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
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
    double& step_size,
    const double initial_step_size,
    arma::vec& dual_averaging_state,
    const double target_acceptance,
    const int iteration,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end,
    const arma::imat& sufficient_pairwise,
    arma::mat& rest_matrix
) {


  // --- Vectorized starting point
  arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  // Working memory for model effects (to avoid mutating input inside logp)
  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  // --- Log-probability and gradient wrapper
  auto logp_and_grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair,
                                 num_categories, is_ordinal_variable);
    arma::mat rmat = observations * current_pair;

    double logp = log_pseudoposterior(current_main, current_pair, inclusion_indicator,
                                      observations, num_categories, num_obs_categories,
                                      sufficient_blume_capel, reference_category, is_ordinal_variable,
                                      threshold_alpha, threshold_beta, interaction_scale,
                                      sufficient_pairwise, rmat);

    arma::vec grad = gradient_log_pseudoposterior(current_main, current_pair, inclusion_indicator,
                                                  observations, num_categories, num_obs_categories,
                                                  sufficient_blume_capel, reference_category, is_ordinal_variable,
                                                  threshold_alpha, threshold_beta, interaction_scale,
                                                  sufficient_pairwise, rmat);

    return std::make_pair(logp, (-grad).eval());
  };

  // --- Run NUTS
  SamplerResult out = nuts_sampler(current_state, step_size, logp_and_grad, 5);

  // --- Update output effects
  current_state = out.theta;
  unvectorize_model_parameters(
    current_state, main_effects, pairwise_effects,
    num_categories, is_ordinal_variable);
  rest_matrix = observations * pairwise_effects;

  double accept_prob = out.alpha / static_cast<double>(out.n_alpha);

  // --- Step size adaptation
  if (iteration < warmup_stageII_end) {
    update_step_size_with_dual_averaging(
      initial_step_size, accept_prob, iteration + 1,
      dual_averaging_state, target_acceptance
    );
    step_size = std::exp(dual_averaging_state[1]);
  } else if (iteration < warmup_stageIII_end) {
    update_step_size_with_robbins_monro(
      accept_prob, iteration - warmup_stageII_end + 1, step_size,
      target_acceptance
    );
  }
}



void update_indicator_interaction_pair_with_nuts (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    double step_size,
    const double interaction_scale,
    const arma::imat& index,
    const int num_persons,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const int num_pairwise,
    const arma::imat& sufficient_pairwise
) {
  for (int pair_index = 0; pair_index < num_pairwise; ++pair_index) {
    const int var1 = index(pair_index, 1);
    const int var2 = index(pair_index, 2);

    const double current_val = pairwise_effects(var1, var2);
    const bool included = indicator(var1, var2);
    const double inclusion_prob = inclusion_probability(var1, var2);

    // Define 1D logp + grad for this edge
    auto logp_and_grad_1d = [&](const arma::vec& x) {
      arma::mat tmp = pairwise_effects;
      tmp(var1, var2) = tmp(var2, var1) = x[0];

      arma::mat rest = observations * tmp;

      double lp = log_pseudoposterior_interactions(
        tmp, main_effects, rest, observations, num_categories,
        indicator, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      double grad = -gradient_log_pseudoposterior_interaction_single(
        var1, var2, tmp, main_effects, observations, rest,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      return std::make_pair(lp, arma::vec({grad}));
    };

    SamplerResult out = nuts_sampler(arma::vec({current_val}), step_size, logp_and_grad_1d);
    double proposed_val = out.theta[0];

    // Compute MH acceptance ratio (as if toggling indicator)
    double log_accept = 0.0;

    if (!included) {
      log_accept += R::dcauchy(proposed_val, 0.0, interaction_scale, true);
      log_accept += std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);
    } else {
      log_accept -= R::dcauchy(current_val, 0.0, interaction_scale, true);
      log_accept -= std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);
    }

    log_accept += log_pseudolikelihood_ratio_interaction(
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      var1, var2, proposed_val, current_val, rest_matrix,
      is_ordinal_variable, reference_category, sufficient_pairwise
    );

    if (std::log(R::unif_rand()) < log_accept) {
      // Accept: toggle indicator and update state
      const int new_value = 1 - indicator(var1, var2);
      indicator(var1, var2) = new_value;
      indicator(var2, var1) = new_value;

      const double delta = proposed_val - current_val;
      pairwise_effects(var1, var2) = proposed_val;
      pairwise_effects(var2, var1) = proposed_val;

      const arma::vec obs1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      rest_matrix.col(var1) += obs2 * delta;
      rest_matrix.col(var2) += obs1 * delta;
    }
  }
}



/**
 * Performs a single iteration of the Gibbs sampler for graphical model parameters.
 *
 * This function performs a full Gibbs update sweep over:
 *   1. Inclusion indicators for pairwise interactions (if selection is enabled)
 *   2. Pairwise interaction coefficients (with MALA or adaptive Metropolis)
 *   3. Main effect (threshold) parameters (with MALA, Fisher-MALA, or Metropolis)
 *
 * The interaction updates support optional Fisher preconditioning, and the step sizes
 * are adapted using either dual averaging (during burn-in) or Robbins-Monro.
 *
 * Inputs:
 *  - observations: Matrix of observed categorical scores (persons × variables).
 *  - num_categories: Number of categories for each variable.
 *  - interaction_scale: Cauchy prior scale for pairwise interaction coefficients.
 *  - proposal_sd_pairwise: Proposal SDs for interaction updates (adaptive Metropolis).
 *  - proposal_sd_main: Proposal SDs for threshold updates (Blume-Capel variables).
 *  - index: List of candidate interaction pairs.
 *  - num_obs_categories: Number of observations per category per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Hyperparameters for main effect priors.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - num_pairwise: Number of candidate interaction pairs.
 *  - num_main: Number of main effect parameters.
 *  - inclusion_indicator: Symmetric binary matrix of active interactions (updated).
 *  - pairwise_effects: Symmetric matrix of interaction strengths (updated).
 *  - main_effects: Matrix of threshold parameters (updated).
 *  - rest_matrix: Linear predictor matrix (updated if interaction/main effects change).
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - rm_decay_rate: Robbins-Monro decay rate (e.g. 0.75).
 *  - is_ordinal_variable: Indicator vector for ordinal variables (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference categories for Blume-Capel variables.
 *  - edge_selection: Whether to update inclusion indicators this iteration.
 *  - step_size_main: Step size for MALA threshold updates (updated).
 *  - iteration: Current iteration number (starts at 0).
 *  - dual_averaging_main: Dual averaging state vector for main effect step size (updated).
 *  - total_burnin: Total number of burn-in iterations.
 *  - use_mala: Whether to use MALA (vs. Metropolis) for updates.
 *  - initial_step_size_main: Initial step size for MALA threshold updates.
 *  - sqrt_inv_fisher_main: Square root inverse Fisher matrix for threshold parameters (updated).
 *  - step_size_pairwise: Step size for interaction MALA updates (updated).
 *  - dual_averaging_pairwise: Dual averaging state for interaction step size (updated).
 *  - initial_step_size_pairwise: Initial step size for interaction MALA updates.
 *  - use_fisher_for_interactions: Whether to use Fisher-preconditioned MALA for interactions.
 *  - sqrt_inv_fisher_pairwise: Square root inverse Fisher matrix for interactions (updated).
 *
 * Updates (in-place):
 *  - inclusion_indicator
 *  - pairwise_effects
 *  - main_effects
 *  - rest_matrix
 *  - step_size_main
 *  - step_size_pairwise
 *  - dual_averaging_main
 *  - dual_averaging_pairwise
 *  - proposal_sd_main
 *  - proposal_sd_pairwise
 *  - sqrt_inv_fisher_main
 *  - sqrt_inv_fisher_pairwise
 */
void gibbs_update_step_for_graphical_model_parameters (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd_pairwise,
    arma::mat& proposal_sd_main,
    const arma::imat& index,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const int num_persons,
    const int num_variables,
    const int num_pairwise,
    const int num_main,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const double rm_decay_rate,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection,
    double& step_size_main,
    const int iteration,
    arma::vec& dual_averaging_main,
    const int total_burnin,
    const double initial_step_size_main,
    arma::mat& sqrt_inv_fisher_main,
    double& step_size_pairwise,
    arma::vec& dual_averaging_pairwise,
    const double initial_step_size_pairwise,
    arma::mat& sqrt_inv_fisher_pairwise,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end,
    arma::imat& sufficient_pairwise,
    double& step_size_joint,
    arma::vec& dual_averaging_joint,
    const double initial_step_size_joint,
    arma::mat& sqrt_inv_fisher_joint,
    const int L,
    arma::vec& mass_inv,
    arma::vec mean,
    arma::vec M2
) {
  // --- Robbins-Monro weight for adaptive Metropolis updates
  const double exp_neg_log_t_rm_adaptation_rate =
    std::exp (-std::log (static_cast<double>(iteration)) * rm_decay_rate);

  // Step 1: Edge selection via MH indicator updates (if enabled)
  if (edge_selection) {
    if (update_method == "fisher-mala-block") {
      // Use Fisher-preconditioned MALA for inclusion indicators
      update_indicator_interaction_pair_with_mala (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, step_size_pairwise, interaction_scale, index,
          num_persons, rest_matrix, inclusion_probability, is_ordinal_variable,
          reference_category, num_pairwise, sufficient_pairwise
      );
    } else if (update_method == "adaptive-metropolis") {
      // Use standard Metropolis-Hastings for indicator updates
      update_indicator_interaction_pair_with_metropolis (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, proposal_sd_pairwise, interaction_scale, index,
          num_pairwise, num_persons, rest_matrix, inclusion_probability,
          is_ordinal_variable, reference_category, sufficient_pairwise
      );
    } // else if (update_method == "fisher-mala-joint") {
    //   // Use Fisher-preconditioned MALA for inclusion indicators
    //   update_indicator_interaction_pair_with_mala (
    //       pairwise_effects, main_effects, inclusion_indicator, observations,
    //       num_categories, step_size_pairwise, interaction_scale, index,
    //       num_persons, rest_matrix, inclusion_probability, is_ordinal_variable,
    //       reference_category, num_pairwise, sufficient_pairwise
    //   );
    // } else if (update_method == "hamiltonian-mc") {
    //   // Use Fisher-preconditioned MALA for inclusion indicators
    //   // update_indicator_interaction_pair_with_mala (
    //   //     pairwise_effects, main_effects, inclusion_indicator, observations,
    //   //     num_categories, step_size_pairwise, interaction_scale, index,
    //   //     num_persons, rest_matrix, inclusion_probability, is_ordinal_variable,
    //   //     reference_category, num_pairwise, sufficient_pairwise
    //   // );
    //   //
    //   update_indicator_interaction_pair_with_hmc(
    //     pairwise_effects, main_effects, indicator, observations, num_categories,
    //     step_size_joint, interaction_scale, index, num_persons, rest_matrix,
    //     inclusion_probability, is_ordinal_variable, reference_category,
    //     num_pairwise, sufficient_pairwise, mass_inv, L
    //   );
    // } else if (update_method == "nuts") {
    // update_indicator_interaction_pair_with_nuts (
    //     pairwise_effects, main_effects, indicator, observations, num_categories,
    //     step_size, interaction_scale, index, num_persons, rest_matrix,
    //     inclusion_probability, is_ordinal_variable, reference_category,
    //     num_pairwise, sufficient_pairwise
    // );
    //}
  }

  // Step 2: Update interaction weights for active edges
  if (update_method == "fisher-mala-block") {
    update_interactions_with_fisher_mala (
        pairwise_effects, rest_matrix, main_effects, observations,
        num_categories, inclusion_indicator, is_ordinal_variable,
        reference_category, interaction_scale, step_size_pairwise,
        initial_step_size_pairwise, iteration, dual_averaging_pairwise,
        sqrt_inv_fisher_pairwise, target_accept, warmup_stageI_end,
        warmup_stageII_end, warmup_stageIII_end, sufficient_pairwise
    );
  } else if (update_method == "adaptive-metropolis") {
    update_interactions_with_adaptive_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, interaction_scale,
        num_persons, num_variables, rest_matrix,
        exp_neg_log_t_rm_adaptation_rate, is_ordinal_variable,
        reference_category, target_accept, iteration,
        warmup_stageIII_end, sufficient_pairwise
    );
  }

  // Step 3: Update main effect (threshold) parameters
  if (update_method == "fisher-mala-block"){
    // Update using Fisher-preconditioned MALA
    update_thresholds_with_fisher_mala (
        main_effects, step_size_main, rest_matrix, num_categories,
        num_obs_categories, sufficient_blume_capel, reference_category,
        is_ordinal_variable, iteration,
        dual_averaging_main,
        sqrt_inv_fisher_main, threshold_alpha, threshold_beta,
        initial_step_size_main, target_accept,
        warmup_stageI_end,
        warmup_stageII_end,
        warmup_stageIII_end
    );
  } else if (update_method == "adaptive-metropolis") {
    // Update using adaptive Metropolis
    update_thresholds_with_adaptive_metropolis (
        main_effects, observations, num_categories, num_obs_categories,
        sufficient_blume_capel, reference_category, is_ordinal_variable,
        num_persons, threshold_alpha, threshold_beta, rest_matrix,
        proposal_sd_main, exp_neg_log_t_rm_adaptation_rate,
        target_accept, iteration, warmup_stageIII_end
    );
  }

  // Step 4: Update joint parameters if applicable
  if (update_method == "fisher-mala-joint") {
    update_parameters_with_fisher_mala (
        main_effects, pairwise_effects, inclusion_indicator, observations,
        num_categories, num_obs_categories, sufficient_blume_capel,
        reference_category, is_ordinal_variable, threshold_alpha,
        threshold_beta, interaction_scale, step_size_joint, iteration,
        dual_averaging_joint, sqrt_inv_fisher_joint,
        initial_step_size_joint, target_accept, warmup_stageI_end,
        warmup_stageII_end, warmup_stageIII_end, sufficient_pairwise,
        rest_matrix
    );
  } else if (update_method == "hamiltonian-mc") {
    update_parameters_with_hmc(
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, step_size_joint, iteration, dual_averaging_joint,
      initial_step_size_joint, target_accept, warmup_stageI_end, warmup_stageII_end,
      warmup_stageIII_end, sufficient_pairwise, rest_matrix, L, mass_inv, mean,
      M2
    );
  } else if (update_method == "nuts") {
    update_parameters_with_nuts(
      main_effects, pairwise_effects, inclusion_indicator,
      observations, num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, step_size_joint, initial_step_size_joint,
      dual_averaging_joint, target_accept, iteration,
      warmup_stageI_end, warmup_stageII_end, warmup_stageIII_end,
      sufficient_pairwise, rest_matrix
    );
  }
}



/**
 * Function: run_gibbs_sampler_for_bgm
 *
 * Runs the full Gibbs sampler for a graphical model with ordinal and/or Blume-Capel variables.
 * Optionally performs edge selection using a Beta-Bernoulli or Stochastic Block prior.
 *
 * During each iteration, the algorithm:
 *  - Updates missing values (if imputation enabled)
 *  - Updates main effects and interactions using adaptive Metropolis, MALA, or Fisher-MALA
 *  - Updates inclusion indicators if edge selection is active
 *  - Adapts proposal variances using Robbins-Monro
 *  - Saves MCMC samples or updates running averages
 *
 * Inputs:
 *  - observations: Imputed categorical data.
 *  - num_categories: Number of categories per variable.
 *  - interaction_scale: Scale for Cauchy prior on pairwise interactions.
 *  - edge_prior: Type of edge prior ("Beta-Bernoulli" or "Stochastic-Block").
 *  - inclusion_probability: Matrix of edge inclusion probabilities (updated if SBM/Beta-Bernoulli).
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Beta prior parameters for edge inclusion.
 *  - dirichlet_alpha, lambda: SBM prior parameters.
 *  - interaction_index_matrix: Matrix of pairwise edge indices (E × 3).
 *  - iter: Number of post-burn-in iterations.
 *  - burnin: Number of burn-in iterations.
 *  - num_obs_categories: Category counts for each variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Prior parameters for logistic-Beta.
 *  - na_impute: Whether to impute missing values.
 *  - missing_index: List of missing value locations.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference category for centering (for Blume-Capel).
 *  - save_main, save_pairwise, save_indicator: Whether to save MCMC samples.
 *  - display_progress: Show progress bar during sampling.
 *  - edge_selection: Whether to enable edge selection during burn-in.
 *  - update_method: One of "adaptive-metropolis", "adaptive-mala", or "fisher-mala".
 *    Controls which sampler is used for main effects and interactions.
 *
 * Returns:
 *  - List containing:
 *    - "main": Posterior mean of main effects
 *    - "pairwise": Posterior mean of interaction effects
 *    - "inclusion_indicator": Posterior mean of inclusion matrix
 *    - (optional) "main_samples", "pairwise_samples", "inclusion_indicator_samples"
 *    - (optional) "allocations": Cluster allocations (if SBM used)
 */
// [[Rcpp::export]]
List run_gibbs_sampler_for_bgm (
    arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    const String& edge_prior,
    arma::mat& inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    arma::imat& num_obs_categories,
    arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool save_main,
    const bool save_pairwise,
    const bool save_indicator,
    const bool display_progress,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    arma::imat& sufficient_pairwise,
    const int L
) {
  // --- Setup: dimensions and storage structures
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;
  const int max_num_categories = num_categories.max();
  const int num_pairwise = interaction_index_matrix.n_rows;

  // Initialize model parameter matrices
  arma::mat main_effects(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_effects(num_variables, num_variables, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Posterior mean accumulators
  arma::mat posterior_mean_main(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat posterior_mean_pairwise(num_variables, num_variables, arma::fill::zeros);
  arma::mat posterior_mean_indicator(num_variables, num_variables, arma::fill::zeros);

  // Residuals used in pseudo-likelihood computation
  arma::mat rest_matrix(num_persons, num_variables, arma::fill::zeros);

  // Allocate optional storage for MCMC samples
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat* main_effect_samples = nullptr;
  arma::mat* pairwise_effect_samples = nullptr;
  arma::imat* indicator_samples = nullptr;

  if (save_main) main_effect_samples = new arma::mat(iter, num_main);
  if (save_pairwise) pairwise_effect_samples = new arma::mat(iter, num_pairwise);
  if (save_indicator) indicator_samples = new arma::imat(iter, num_pairwise);

  // Edge update shuffling setup
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pairwise - 1);
  arma::uvec order(num_pairwise);
  arma::imat index(num_pairwise, 3);

  // SBM-specific structures
  arma::uvec K_values;
  arma::uvec cluster_allocations(num_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);
  arma::imat out_allocations(iter, num_variables);

  // --- Initialize SBM prior if applicable
  if (edge_prior == "Stochastic-Block") {
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for (int i = 2; i < num_variables; i++) {
      cluster_allocations[i] = (R::unif_rand() > 0.5) ? 1 : 0;
    }

    cluster_prob = block_probs_mfm_sbm(
      cluster_allocations, arma::conv_to<arma::umat>::from(inclusion_indicator),
      num_variables, beta_bernoulli_alpha, beta_bernoulli_beta
    );

    for (int i = 0; i < num_variables - 1; i++) {
      for (int j = i + 1; j < num_variables; j++) {
        inclusion_probability(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
        inclusion_probability(j, i) = inclusion_probability(i, j);
      }
    }

    log_Vn = compute_Vn_mfm_sbm(num_variables, dirichlet_alpha, num_variables + 10, lambda);
  }

  // Initialize proposal SDs, MALA and HMC tracking
  arma::mat proposal_sd_main(num_main, max_num_categories, arma::fill::ones);
  arma::mat proposal_sd_pairwise(num_variables, num_variables, arma::fill::ones);
  arma::vec mass_inv(num_main + num_pairwise, arma::fill::ones);
  arma::vec mean = arma::zeros(num_main + num_pairwise);
  arma::vec M2 = arma::zeros(num_main + num_pairwise);

  double step_size_main = 0.01;
  double step_size_pairwise = 0.01;
  double step_size_joint = 0.01;
  arma::vec dual_averaging_main(3, arma::fill::zeros);
  arma::vec dual_averaging_pairwise(3, arma::fill::zeros);
  arma::vec dual_averaging_joint(3, arma::fill::zeros);
  double initial_step_size_main = 0.01;
  double initial_step_size_pairwise = 0.01;
  double initial_step_size_joint = 0.01;
  const double rm_decay_rate = 0.75;

  // --- Optional MALA warmup stage (step size tuning only)
  arma::mat sqrt_inv_fisher_main(num_main, num_main, arma::fill::eye);  // Identity at start
  arma::mat sqrt_inv_fisher_pairwise(num_pairwise, num_pairwise, arma::fill::eye);
  arma::mat sqrt_inv_fisher_joint(num_main + num_pairwise, num_main + num_pairwise, arma::fill::eye);
  if(update_method == "fisher-mala-block") {
    // Warmup phase for MALA: find initial step size (per Hoffman & Gelman, 2014)
    // Uses one MALA step to tune the log step size before dual averaging.
    initial_step_size_pairwise = find_reasonable_initial_step_size_interactions (
      pairwise_effects, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale, target_accept, sufficient_pairwise,
      rest_matrix
    );
    step_size_pairwise = initial_step_size_pairwise;
    dual_averaging_pairwise[0] = std::log (step_size_pairwise);

    initial_step_size_main = find_reasonable_initial_step_size_thresholds (
      main_effects, rest_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta, target_accept
    );
    step_size_main = initial_step_size_main;
    dual_averaging_main[0] = std::log (step_size_main);
  } else if (update_method == "fisher-mala-joint") {
    initial_step_size_joint = find_reasonable_initial_step_size (
      main_effects, pairwise_effects, inclusion_indicator, observations,
      rest_matrix, num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha,
      threshold_beta, interaction_scale, target_accept, sufficient_pairwise
    );
    step_size_joint = initial_step_size_joint;
    dual_averaging_joint[0] = std::log (step_size_joint);
  } else if (update_method == "hamiltonian-mc") {
    initial_step_size_joint = find_reasonable_initial_step_size_hmc (
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, target_accept, sufficient_pairwise, rest_matrix,
      mass_inv
    );
    step_size_joint = initial_step_size_joint;
    dual_averaging_joint[0] = std::log (step_size_joint);
  } else if (update_method == "nuts") {
    initial_step_size_joint = find_reasonable_initial_step_size_nuts (
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha,
      threshold_beta, interaction_scale, target_accept, sufficient_pairwise
    );
    step_size_joint = initial_step_size_joint;
    dual_averaging_joint[0] = std::log (step_size_joint);
  }

  // --- Set up total number of iterations (burn-in + sampling)
  bool enable_edge_selection = edge_selection;
  int warmup_stageI_end = static_cast<int>(std::round(0.1 * burnin));
  int warmup_stageII_len = burnin - static_cast<int>(std::round(0.2 * burnin));
  int warmup_stageII_end = warmup_stageI_end + warmup_stageII_len * (enable_edge_selection ? 2 : 1);
  int warmup_stageIII_end = burnin + warmup_stageII_len * (enable_edge_selection ? 1 : 0);// Add stage II if edge selection is enabled
  int total_burnin = warmup_stageIII_end;

  edge_selection = false;
  const int total_iter = total_burnin + iter;
  Progress p(total_iter, display_progress);

  // --- Main Gibbs sampling loop
  for (int iteration = 0; iteration < total_iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(
        Named("main") = posterior_mean_main,
        Named("pairwise") = posterior_mean_pairwise,
        Named("inclusion_indicator") = posterior_mean_indicator
      );
    }
    p.increment();

    // Re-enable edge selection halfway through burn-in
    if (enable_edge_selection && iteration == warmup_stageI_end + warmup_stageII_len)
      edge_selection = true;

    // Shuffle update order of edge indices
    order = arma::randperm(num_pairwise);
    for (int i = 0; i < num_pairwise; i++) {
      index.row(i) = interaction_index_matrix.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_values_for_graphical_model (
          pairwise_effects, main_effects, observations, num_obs_categories,
          sufficient_blume_capel, num_categories, rest_matrix,
          missing_index, is_ordinal_variable, reference_category,
          sufficient_pairwise
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_for_graphical_model_parameters (
        observations, num_categories, interaction_scale, proposal_sd_pairwise,
        proposal_sd_main, index, num_obs_categories, sufficient_blume_capel,
        threshold_alpha, threshold_beta, num_persons, num_variables, num_pairwise,
        num_main, inclusion_indicator, pairwise_effects, main_effects,
        rest_matrix, inclusion_probability, rm_decay_rate, is_ordinal_variable,
        reference_category, edge_selection, step_size_main, iteration,
        dual_averaging_main, total_burnin, initial_step_size_main,
        sqrt_inv_fisher_main, step_size_pairwise, dual_averaging_pairwise,
        initial_step_size_pairwise, sqrt_inv_fisher_pairwise,
        update_method, pairwise_effect_indices, target_accept,
        warmup_stageI_end, warmup_stageII_end, warmup_stageIII_end,
        sufficient_pairwise, step_size_joint, dual_averaging_joint,
        initial_step_size_joint, sqrt_inv_fisher_joint, L,
        mass_inv, mean, M2
    );

    // --- Update edge probabilities under the prior (if edge selection is active)
    if (edge_selection) {
      if (edge_prior == "Beta-Bernoulli") {
        int num_edges_included = 0;
        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            num_edges_included += inclusion_indicator(i, j);

        double prob = R::rbeta(
          beta_bernoulli_alpha + num_edges_included,
          beta_bernoulli_beta + num_pairwise - num_edges_included
        );

        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            inclusion_probability(i, j) = inclusion_probability(j, i) = prob;

      } else if (edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm(
          cluster_allocations, num_variables, log_Vn, cluster_prob,
          arma::conv_to<arma::umat>::from(inclusion_indicator), dirichlet_alpha,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        cluster_prob = block_probs_mfm_sbm(
          cluster_allocations,
          arma::conv_to<arma::umat>::from(inclusion_indicator), num_variables,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        for (int i = 0; i < num_variables - 1; i++) {
          for (int j = i + 1; j < num_variables; j++) {
            inclusion_probability(i, j) = inclusion_probability(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }

    // --- Save samples and update posterior means
    if (iteration >= total_burnin) {
      int iter_adj = iteration - total_burnin + 1;

      // Running posterior means
      posterior_mean_main = (posterior_mean_main * (iter_adj - 1) + main_effects) / iter_adj;
      posterior_mean_pairwise = (posterior_mean_pairwise * (iter_adj - 1) + pairwise_effects) / iter_adj;

      if (edge_selection) {
        posterior_mean_indicator = (posterior_mean_indicator * (iter_adj - 1) +
          arma::conv_to<arma::mat>::from(inclusion_indicator)) / iter_adj;
      }

      int sample_index = iteration - total_burnin;

      if (save_main) {
        arma::vec vectorized_main = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
        main_effect_samples->row(sample_index) = vectorized_main.t();
      }

      if (save_pairwise) {
        arma::vec vectorized_pairwise(num_pairwise);
        for (int i = 0; i < num_pairwise; i++) {
          vectorized_pairwise(i) = pairwise_effects(interaction_index_matrix(i, 1), interaction_index_matrix(i, 2));
        }
        pairwise_effect_samples->row(sample_index) = vectorized_pairwise.t();
      }

      if (save_indicator) {
        arma::ivec vectorized_indicator(num_pairwise);
        for (int i = 0; i < num_pairwise; i++) {
          vectorized_indicator(i) = inclusion_indicator(interaction_index_matrix(i, 1), interaction_index_matrix(i, 2));
        }
        indicator_samples->row(sample_index) = vectorized_indicator.t();
      }

      if (edge_prior == "Stochastic-Block") {
        for (int j = 0; j < num_variables; j++)
          out_allocations(sample_index, j) = cluster_allocations[j] + 1;
      }
    }
  }

  // --- Final output
  List out = List::create(
    Named("main") = posterior_mean_main,
    Named("pairwise") = posterior_mean_pairwise,
    Named("inclusion_indicator") = posterior_mean_indicator
  );

  if (save_main) {
    out["main_samples"] = *main_effect_samples;
    delete main_effect_samples;
  }
  if (save_pairwise) {
    out["pairwise_samples"] = *pairwise_effect_samples;
    delete pairwise_effect_samples;
  }
  if (save_indicator) {
    out["inclusion_indicator_samples"] = *indicator_samples;
    delete indicator_samples;
  }
  if (edge_prior == "Stochastic-Block") {
    out["allocations"] = out_allocations;
  }

  return out;
}