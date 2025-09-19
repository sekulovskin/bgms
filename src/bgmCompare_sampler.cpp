#include <RcppArmadillo.h>
#include "bgmCompare_helper.h"
#include "bgmCompare_logp_and_grad.h"
#include "bgmCompare_sampler.h"
#include "common_helpers.h"
#include "mcmc_adaptation.h"
#include "mcmc_hmc.h"
#include "mcmc_leapfrog.h"
#include "mcmc_nuts.h"
#include "mcmc_rwm.h"
#include "mcmc_utils.h"
#include "print_mutex.h"
#include "rng_utils.h"
#include "sampler_output.h"
#include <string>

using namespace Rcpp;


/**
 * Function: impute_missing_data_for_anova_model
 * Purpose: Imputes missing data for independent samples designs by generating new observations
 *          based on the model parameters and pseudo-likelihood.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise interaction effects between variables across groups.
 *  - main_effect_indices: Integer matrix mapping variable indices to main effect parameters.
 *  - pairwise_effect_indices: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling for effects.
 *  - observations: Integer matrix of observed data (individuals x variables), with missing data encoded.
 *  - num_groups: Number of groups in the analysis.
 *  - group_membership: Integer vector mapping individuals to their respective groups.
 *  - num_obs_categories: List of matrices, one per group, recording category frequencies per variable.
 *  - sufficient_blume_capel: List of matrices storing sufficient statistics for Blume-Capel variables.
 *  - num_categories: Integer matrix of category counts for each variable and group.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - missing_data_indices: Integer matrix of indices indicating missing observations (row x column pairs).
 *  - is_ordinal_variable: Logical vector indicating whether variables are ordinal.
 *  - baseline_category: Integer vector of reference categories for Blume-Capel variables.
 *
 * Outputs:
 *  - A List containing:
 *    - `observations`: Updated observation matrix with imputed values.
 *    - `num_obs_categories`: Updated list of category counts per group.
 *    - `sufficient_blume_capel`: Updated sufficient statistics for Blume-Capel variables.
 *    - `residual_matrix`: Updated residual effects matrix.
 */
void impute_missing_data_for_graphical_model(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    arma::imat& observations,
    const int num_groups,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    std::vector<arma::imat>& num_obs_categories,
    std::vector<arma::imat>& sufficient_blume_capel,
    std::vector<arma::mat>& sufficient_pairwise,
    const arma::imat& num_categories,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    SafeRNG& rng
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_data_indices.n_rows;
  const int max_num_categories = arma::max(arma::vectorise(num_categories));

  arma::vec category_response_probabilities(max_num_categories + 1);
  double exponent, cumsum, u;
  int score, person, variable, new_observation, old_observation, gr;

  //Impute missing data
  for(int missing = 0; missing < num_missings; missing++) {
    // Identify the observation to impute
    person = missing_data_indices(missing, 0);
    variable = missing_data_indices(missing, 1);
    gr = group_membership[person];

    const arma::vec proj_g = projection.row(gr).t();
    // Compute thresholds for the variable in the given group
    arma::vec GroupThresholds = compute_group_main_effects(
      variable, num_groups, main_effects,  main_effect_indices, proj_g);

    // Generate a new observation based on the model
    arma::mat GroupInteractions(num_variables, num_variables, arma::fill::zeros);
    for(int v1 = 0; v1 < num_variables-1; v1++) {
      for(int v2 = v1 + 1; v2 < num_variables; v2++) {
        double w = compute_group_pairwise_effects(
            v1, v2, num_groups, pairwise_effects, pairwise_effect_indices,
            inclusion_indicator, proj_g
        );
        GroupInteractions(v1, v2) = w;
        GroupInteractions(v2, v1) = w;
      }
    }

    double rest_score =
      arma::as_scalar(observations.row(person) * GroupInteractions.col(variable));
    if(is_ordinal_variable[variable] == true) {
      // For regular binary or ordinal variables
      cumsum = 1.0;
      category_response_probabilities[0] = 1.0;
      for(int category = 1; category <= num_categories(variable, gr); category++) {
        exponent = GroupThresholds(category - 1);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        category_response_probabilities[category] = cumsum;
      }
    } else {
      // For Blume-Capel variables
      cumsum = 0.0;
      for(int category = 0; category <= num_categories(variable, gr); category++) {
        exponent = GroupThresholds[0] * category;
        exponent += GroupThresholds[1] *
          (category - baseline_category[variable]) *
          (category - baseline_category[variable]);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        category_response_probabilities[category] = cumsum;
      }
    }

    // Sample a new value based on computed probabilities
    u = cumsum * runif(rng);
    score = 0;
    while (u > category_response_probabilities[score]) {
      score++;
    }
    new_observation = score;
    old_observation = observations(person, variable);

    if(old_observation != new_observation) {
      // Update raw observations
      observations(person, variable) = new_observation;

      // Update sufficient statistics for main effects
      if(is_ordinal_variable[variable] == true) {
        arma::imat num_obs_categories_gr = num_obs_categories[gr];
        if(old_observation > 0)
          num_obs_categories_gr(old_observation, variable)--;
        if(new_observation > 0)
          num_obs_categories_gr(new_observation, variable)++;
        num_obs_categories[gr] = num_obs_categories_gr;
      } else {
        arma::imat sufficient_blume_capel_gr = sufficient_blume_capel[gr];
        sufficient_blume_capel_gr(0, variable) -= old_observation;
        sufficient_blume_capel_gr(0, variable) += new_observation;
        sufficient_blume_capel_gr(1, variable) -=
          (old_observation - baseline_category[variable]) *
          (old_observation - baseline_category[variable]);
        sufficient_blume_capel_gr(1, variable) +=
          (new_observation - baseline_category[variable]) *
          (new_observation - baseline_category[variable]);
        sufficient_blume_capel[gr] = sufficient_blume_capel_gr;
      }

      // Update sufficient statistics for pairwise effects
      const int r0 = group_indices(gr, 0);
      const int r1 = group_indices(gr, 1);
      arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
      arma::mat sufficient_pairwise_gr = obs.t() * obs; // crossprod
      sufficient_pairwise[gr] = sufficient_pairwise_gr;
    }
  }

  return;
}



/**
 * Function: update_main_effects_with_metropolis_cmp
 */
void update_main_effects_with_metropolis_cmp (
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double difference_scale,
    const double main_alpha,
    const double main_beta,
    const int iteration,
    RWMAdaptationController& rwm_adapt,
    SafeRNG& rng,
    arma::mat& proposal_sd_main
) {
  const int num_vars = observations.n_cols;
  arma::umat index_mask_main = arma::zeros<arma::umat>(proposal_sd_main.n_rows,
                                                       proposal_sd_main.n_cols);
  arma::mat accept_prob_main = arma::zeros<arma::mat>(proposal_sd_main.n_rows,
                                                      proposal_sd_main.n_cols);

  // --- helper for one update ---
  auto do_update = [&](int row, int variable, int category, int par, int h) {
    index_mask_main(row, h) = 1;
    double& current = main_effects(row, h);
    double proposal_sd = proposal_sd_main(row, h);

    auto log_post = [&](double theta) {
      main_effects(row, h) = theta;
      return log_pseudoposterior_main_component(
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, projection, observations,
        group_indices, num_categories, num_obs_categories,
        sufficient_blume_capel, num_groups, inclusion_indicator,
        is_ordinal_variable, baseline_category,
        main_alpha, main_beta, difference_scale,
        variable, category, par, h
      );
    };

    SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);
    current = result.state[0];
    accept_prob_main(row, h) = result.accept_prob;
  };

  // --- loop over variables ---
  for (int variable = 0; variable < num_vars; ++variable) {
    int base_category_index = main_effect_indices(variable, 0);
    int num_cats = num_categories(variable);
    bool group_differences = (inclusion_indicator(variable, variable) == 1);

    if (is_ordinal_variable[variable]) {
      // ordinal: loop categories
      for (int category = 0; category < num_cats; ++category) {
        int row = base_category_index + category;
        int hmax = group_differences ? num_groups : 1;
        for (int h = 0; h < hmax; ++h)
          do_update(row, variable, category, -1, h);
      }
    } else {
      // non-ordinal: two parameters
      for (int par = 0; par < 2; ++par) {
        int row = base_category_index + par;
        int hmax = group_differences ? num_groups : 1;
        for (int h = 0; h < hmax; ++h)
          do_update(row, variable, -1, par, h);
      }
    }
  }

  rwm_adapt.update(index_mask_main, accept_prob_main, iteration);
}



/**
 * Function: update_pairwise_effects_with_metropolis_cmp
 */
void update_pairwise_effects_with_metropolis_cmp (
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::mat>& sufficient_pairwise,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const double difference_scale,
    const int iteration,
    RWMAdaptationController& rwm_adapt,
    SafeRNG& rng,
    arma::mat& proposal_sd_pair
) {
  int num_variables = observations.n_cols;
  int num_pairs = num_variables * (num_variables - 1) / 2;
  arma::mat accept_prob_pair = arma::zeros<arma::mat>(num_pairs,
                                                      num_groups);
  arma::umat index_mask_pair = arma::zeros<arma::umat>(num_pairs,
                                                       num_groups);

  // --- helper for one update ---
  auto do_update = [&](int var1, int var2, int h) {
    int idx = pairwise_effect_indices(var1, var2);
    index_mask_pair(idx, h) = 1;
    double& current = pairwise_effects(idx, h);
    double proposal_sd = proposal_sd_pair(idx, h);

    auto log_post = [&](double theta) {
      pairwise_effects(idx, h) = theta;
      return log_pseudoposterior_pair_component(
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, projection, observations, group_indices,
        num_categories, sufficient_pairwise, num_groups,
        inclusion_indicator, is_ordinal_variable, baseline_category,
        pairwise_scale, difference_scale, var1, var2, h
      );
    };

    SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);
    current = result.state[0];
    accept_prob_pair(idx, h) = result.accept_prob;
  };

  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      bool group_differences = (inclusion_indicator(var1, var2) == 1);
      int hmax = group_differences ? num_groups : 1;
      for (int h = 0; h < hmax; h++) {
        do_update(var1, var2, h);
      }
    }
  }

  rwm_adapt.update(index_mask_pair, accept_prob_pair, iteration);
}



/**
 * Function: find_reasonable_initial_step_size_cmp
 *
 * Heuristically finds a reasonable initial step size for leapfrog-based MCMC algorithms
 * (such as HMC and NUTS), following the procedure described in:
 *
 *   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn Sampler: Adaptively Setting
 *   Path Lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15, 1593–1623.
 *   [Algorithm 4: Heuristic for Choosing an Initial Value of ε]
 *
 * The algorithm simulates a single leapfrog step and compares the resulting change
 * in the log joint density (Hamiltonian). If the proposal is too likely (acceptance too high),
 * the step size is increased. If it is too unlikely (acceptance too low), it is decreased.
 * The process repeats until the log acceptance probability is approximately −log(2),
 * corresponding to a target acceptance probability of ~0.5.
 *
 * Inputs:
 *  - main_effects: Matrix of main (threshold) parameters.
 *  - pairwise_effects: Matrix of pairwise interaction parameters.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Observed category counts per variable (matrix).
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - baseline_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - main_alpha, main_beta: Prior hyperparameters for main effects.
 *  - pairwise_scale: Scale of the Cauchy prior on interactions.
 *  - target_acceptance: Desired acceptance probability (typically ~0.65).
 *  - sufficient_pairwise: Sufficient statistics for pairwise interactions.
 *
 * Returns:
 *  - A scalar step size ε that yields roughly the target acceptance probability
 *    under a single leapfrog step.
 *
 * Note:
 *  - This function is suitable for both NUTS and standard HMC algorithms.
 *  - It is typically called once before warm-up/adaptation.
 */
double find_reasonable_initial_step_size_cmp(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>& sufficient_pairwise,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const double difference_scale,
    const double main_alpha,
    const double main_beta,
    const double target_acceptance,
    SafeRNG& rng
) {
  arma::vec theta = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator, main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable

  );
  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects,
    inclusion_indicator,
    main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  auto& main_index = index_maps.first;
  auto& pair_index = index_maps.second;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale, main_index, pair_index
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return log_pseudoposterior(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale
    );
  };

  return heuristic_initial_step_size(theta, log_post, grad, rng, target_acceptance);
}



void update_parameters_with_hmc(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>& sufficient_pairwise,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const double difference_scale,
    const double main_alpha,
    const double main_beta,
    const int num_leapfrogs,
    const int iteration,
    HMCAdaptationController& hmc_adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_categories,
    is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects,
    inclusion_indicator,
    main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  auto& main_index = index_maps.first;
  auto& pair_index = index_maps.second;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale, main_index, pair_index
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return log_pseudoposterior(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale
    );
  };

  //adapt
  arma::vec active_inv_mass = inv_mass_active(
    hmc_adapt.inv_mass_diag(), inclusion_indicator, num_groups, num_categories,
    is_ordinal_variable, main_index, pair_index, main_effect_indices,
    pairwise_effect_indices, selection
  );

  SamplerResult result = hmc_sampler(
    current_state, hmc_adapt.current_step_size(), log_post, grad, num_leapfrogs,
    active_inv_mass, rng
  );

  current_state = result.state;
  unvectorize_model_parameters(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
    is_ordinal_variable
  );

  hmc_adapt.update(current_state, result.accept_prob, iteration);
}



/**
 * Function: update_parameters_with_nuts
 *
 * Performs one update of the main and pairwise effect parameters using
 * the No-U-Turn Sampler (NUTS), with centralized step size and mass matrix adaptation.
 *
 * Step size and mass matrix adaptation are handled via the HMCAdaptationController,
 * which manages warmup phases and dual averaging internally.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - pairwise_effects: Matrix of interaction parameters (updated in-place).
 *  - inclusion_indicator: Binary matrix indicating which interactions are active.
 *  - observations: Matrix of categorical observations.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Matrix of observed category counts.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel thresholds.
 *  - sufficient_pairwise: Sufficient statistics for pairwise terms.
 *  - baseline_category: Reference category per variable.
 *  - is_ordinal_variable: Logical vector indicating ordinal vs. Blume-Capel (1 = ordinal).
 *  - main_alpha, main_beta: Hyperparameters for main effect priors.
 *  - pairwise_scale: Scale parameter for the Cauchy prior on interactions.
 *  - sufficient_pairwise: Sufficient statistics for pairwise interactions.
 *  - rest_matrix: Matrix of residual scores (observations × variables), updated in-place.
 *  - nuts_max_depth: Maximum tree depth for the NUTS trajectory expansion.
 *  - iteration: Current iteration number.
 *  - adapt: Adaptation controller (step size + mass matrix).
 *
 * Modifies (in-place):
 *  - main_effects, pairwise_effects: Updated if NUTS proposal is accepted.
 *  - rest_matrix: Recomputed from the updated pairwise_effects.
 *  - adapt: Updated with step size and mass matrix changes if within warmup.
 */
SamplerResult update_parameters_with_nuts(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>& sufficient_pairwise,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const double difference_scale,
    const double main_alpha,
    const double main_beta,
    const int nuts_max_depth,
    const int iteration,
    HMCAdaptationController& hmc_adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_categories,
    is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects,
    inclusion_indicator,
    main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  auto& main_index = index_maps.first;
  auto& pair_index = index_maps.second;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale, main_index, pair_index
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return log_pseudoposterior(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, observations, group_indices, num_categories,
      num_obs_categories, sufficient_blume_capel,
      sufficient_pairwise, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, main_alpha, main_beta,
      pairwise_scale, difference_scale
    );
  };

  //adapt
  arma::vec active_inv_mass = inv_mass_active(
    hmc_adapt.inv_mass_diag(), inclusion_indicator, num_groups, num_categories,
    is_ordinal_variable, main_index, pair_index, main_effect_indices,
    pairwise_effect_indices, selection
  );

  SamplerResult result = nuts_sampler(
    current_state, hmc_adapt.current_step_size(), log_post, grad,
    active_inv_mass, rng, nuts_max_depth
  );

  current_state = result.state;
  unvectorize_model_parameters(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
    is_ordinal_variable
  );

  hmc_adapt.update(current_state, result.accept_prob, iteration);

  return result;
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
 *  - pairwise_scale: Cauchy prior scale for pairwise interaction coefficients.
 *  - proposal_sd_pairwise: Proposal SDs for interaction updates (adaptive Metropolis).
 *  - proposal_sd_main: Proposal SDs for threshold updates (Blume-Capel variables).
 *  - index: List of candidate interaction pairs.
 *  - num_obs_categories: Number of observations per category per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for main effect priors.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - num_pair: Number of candidate interaction pairs.
 *  - num_main: Number of main effect parameters.
 *  - inclusion_indicator: Symmetric binary matrix of active interactions (updated).
 *  - pairwise_effects: Symmetric matrix of interaction strengths (updated).
 *  - main_effects: Matrix of threshold parameters (updated).
 *  - rest_matrix: Linear predictor matrix (updated if interaction/main effects change).
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - rm_decay_rate: Robbins-Monro decay rate (e.g. 0.75).
 *  - is_ordinal_variable: Indicator vector for ordinal variables (1 = ordinal, 0 = Blume-Capel).
 *  - baseline_category: Reference categories for Blume-Capel variables.
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
void gibbs_update_step_for_graphical_model_parameters_cmp (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double pairwise_scale,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const double main_alpha,
    const double main_beta,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    const arma::imat& pairwise_effect_indices,
    const std::vector<arma::mat>& sufficient_pairwise,
    const int nuts_max_depth,
    HMCAdaptationController& hmc_adapt,
    RWMAdaptationController& rwm_adapt_main,
    RWMAdaptationController& rwm_adapt_pair,
    const bool learn_mass_matrix,
    WarmupSchedule const& schedule,
    arma::ivec& treedepth_samples,
    arma::ivec& divergent_samples,
    arma::vec& energy_samples,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const int num_groups,
    const arma::imat group_indices,
    double difference_scale,
    SafeRNG& rng,
    arma::mat& inclusion_probability,
    int hmc_nuts_leapfrogs,
    const std::string& update_method,
    arma::mat& proposal_sd_main,
    arma::mat& proposal_sd_pair
) {

  // Step 0: Initialise random graph structure when edge_selection = TRUE
  if (schedule.selection_enabled(iteration) && iteration == schedule.stage3c_start) {
    initialise_graph(
      inclusion_indicator, main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_probability, rng
    );
  }

  // Step 1: Difference selection via MH indicator updates (if enabled)
  // ....

  // Step 2: Update parameters
  if(update_method == "adaptive-metropolis") {
    update_main_effects_with_metropolis_cmp (
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, inclusion_indicator, projection,
        num_categories, observations, num_groups, group_indices,
        num_obs_categories, sufficient_blume_capel, is_ordinal_variable,
        baseline_category, difference_scale, main_alpha, main_beta, iteration,
        rwm_adapt_main, rng, proposal_sd_main
    );

    update_pairwise_effects_with_metropolis_cmp (
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, inclusion_indicator, projection,
        num_categories, observations, num_groups, group_indices,
        sufficient_pairwise, is_ordinal_variable, baseline_category,
        pairwise_scale, difference_scale, iteration, rwm_adapt_pair, rng,
        proposal_sd_pair
    );
  } else if (update_method == "hamiltonian-mc") {
    update_parameters_with_hmc(
      main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_indicator, projection, num_categories,
      observations, num_groups, group_indices, num_obs_categories,
      sufficient_blume_capel, sufficient_pairwise, is_ordinal_variable,
      baseline_category, pairwise_scale, difference_scale, main_alpha,
      main_beta, hmc_nuts_leapfrogs, iteration, hmc_adapt, learn_mass_matrix,
      schedule.selection_enabled(iteration), rng
    );
  } else if (update_method == "nuts") {
    SamplerResult result = update_parameters_with_nuts(
      main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_indicator, projection, num_categories,
      observations, num_groups, group_indices, num_obs_categories,
      sufficient_blume_capel, sufficient_pairwise, is_ordinal_variable,
      baseline_category, pairwise_scale, difference_scale, main_alpha,
      main_beta, nuts_max_depth, iteration, hmc_adapt, learn_mass_matrix,
      schedule.selection_enabled(iteration), rng
    );

    if (iteration >= schedule.total_burnin) {
      int sample_index = iteration - schedule.total_burnin;
      if (auto diag = std::dynamic_pointer_cast<NUTSDiagnostics>(result.diagnostics)) {
        treedepth_samples(sample_index) = diag->tree_depth;
        divergent_samples(sample_index) = diag->divergent ? 1 : 0;
        energy_samples(sample_index) = diag->energy;
      }
    }
  }

  /* --- 2b.  proposal-sd tuning during Stage-3b ------------------------------ */
  // ....

}



SamplerOutput run_gibbs_sampler_for_bgmCompare(
    int chain_id,
    arma::imat observations,
    const int num_groups,
    std::vector<arma::imat>& num_obs_categories,
    std::vector<arma::imat>& sufficient_blume_capel,
    std::vector<arma::mat>& sufficient_pairwise,
    const arma::ivec& num_categories,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const double difference_scale,//new
    const double difference_selection_alpha,//new
    const double difference_selection_beta,//new
    const std::string& difference_prior,//new
    const int iter,
    const int burnin,
    const bool na_impute,
    const arma::imat& missing_data_indices,//updated
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,//new
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat& projection,//new
    const arma::ivec& group_membership,//new
    const arma::imat& group_indices,//new
    const arma::imat& interaction_index_matrix,//new
    arma::mat inclusion_probability,//new
    SafeRNG& rng,
    const std::string& update_method,
    const int hmc_num_leapfrogs
) {
  // --- Setup: dimensions and storage structures
  const int num_variables = observations.n_cols;
  const int num_main = count_num_main_effects (
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  // Initialize model parameter matrices
  arma::mat main_effects(num_main, num_groups, arma::fill::zeros);
  arma::mat pairwise_effects(num_pair, num_groups, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Allocate optional storage for MCMC samples
  arma::mat main_effect_samples(iter, num_main * num_groups);
  arma::mat pairwise_effect_samples(iter, num_pair * num_groups);
  arma::imat indicator_samples;

  if (difference_selection) {
    indicator_samples.set_size(iter, num_pair + num_variables);
  }

  // For logging nuts performance
  arma::ivec treedepth_samples(iter, arma::fill::zeros);
  arma::ivec divergent_samples(iter, arma::fill::zeros);
  arma::vec energy_samples(iter, arma::fill::zeros);

  // Edge update shuffling setup
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pair - 1);
  arma::uvec order(num_pair);
  arma::imat index(num_pair, 3);

  // --- Initialize proposal SDs
  arma::mat proposal_sd_main(num_main, num_groups, arma::fill::ones);
  arma::mat proposal_sd_pair(num_pair, num_groups, arma::fill::ones);

  // --- Optional HMC/NUTS warmup stage
  double initial_step_size = 1.0;
  if (update_method == "hamiltonian-mc" || update_method == "nuts") {
    initial_step_size = find_reasonable_initial_step_size_cmp(
      main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_indicator, projection, num_categories,
      observations, num_groups, group_indices, num_obs_categories,
      sufficient_blume_capel, sufficient_pairwise, is_ordinal_variable,
      baseline_category, pairwise_scale, difference_scale, main_alpha, main_beta,
      target_accept, rng
    );
  }

  // --- Warmup scheduling + adaptation controller
  WarmupSchedule warmup_schedule(burnin, difference_selection, true);

  HMCAdaptationController hmc_adapt(
      (num_main + num_pair) * num_groups, initial_step_size, target_accept,
      warmup_schedule, learn_mass_matrix
  );

  RWMAdaptationController rwm_adapt_main(
      proposal_sd_main, warmup_schedule, target_accept
  );
  RWMAdaptationController rwm_adapt_pair(
      proposal_sd_pair, warmup_schedule, target_accept
  );

  const int total_iter = warmup_schedule.total_burnin + iter;
  const int print_every = std::max(1, total_iter / 10);

  // --- Main Gibbs sampling loop
  for (int iteration = 0; iteration < total_iter; iteration++) {
    if (iteration % print_every == 0) {
      tbb::mutex::scoped_lock lock(get_print_mutex());
      std::cout
      << "[bgm] chain " << chain_id
      << " iteration " << iteration
      << " / " << total_iter
      << std::endl;
    }

    // Shuffle update order of edge indices
    order = arma_randperm(rng, num_pair);
    for (int i = 0; i < num_pair; i++) {
      index.row(i) = interaction_index_matrix.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_data_for_graphical_model (
          main_effects, pairwise_effects, main_effect_indices,
          pairwise_effect_indices, inclusion_indicator, projection,
          observations, num_groups, group_membership, group_indices,
          num_obs_categories, sufficient_blume_capel, sufficient_pairwise,
          num_categories, missing_data_indices, is_ordinal_variable,
          baseline_category, rng
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_for_graphical_model_parameters_cmp (
        observations, num_categories, pairwise_scale, num_obs_categories,
        sufficient_blume_capel, main_alpha, main_beta, inclusion_indicator,
        pairwise_effects, main_effects, is_ordinal_variable, baseline_category,
        iteration, pairwise_effect_indices, sufficient_pairwise, nuts_max_depth,
        hmc_adapt, rwm_adapt_main, rwm_adapt_pair, learn_mass_matrix,
        warmup_schedule, treedepth_samples, divergent_samples, energy_samples,
        main_effect_indices, projection, num_groups, group_indices, difference_scale,//new line of args
        rng, inclusion_probability, hmc_num_leapfrogs, update_method,
        proposal_sd_main, proposal_sd_pair
    );

    // --- Update difference probabilities under the prior (if difference selection is active)
    if (warmup_schedule.selection_enabled(iteration)) {
      int sumG = 0;

      if (difference_prior == "Beta-Bernoulli") {
        // Update pairwise inclusion probabilities
        for (int i = 0; i < num_variables - 1; ++i) {
          for (int j = i + 1; j < num_variables; ++j) {
            sumG += inclusion_indicator(i, j);
          }
        }
        for(int i = 0; i < num_variables; i++) {
          sumG += inclusion_indicator(i, i);
        }
        double prob = rbeta(rng,
                            difference_selection_alpha + sumG,
                            difference_selection_beta + num_pair + num_variables - sumG);
        std::fill(inclusion_probability.begin(), inclusion_probability.end(), prob);
      }
    }

    // --- Store states
    if (iteration >= warmup_schedule.total_burnin) {
      int sample_index = iteration - warmup_schedule.total_burnin;


      int cntr = 0;
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_main; ++row) {
          main_effect_samples(sample_index, cntr) = main_effects(row, col);
          cntr++;
        }
      }

      cntr = 0;
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_pair; ++row) {
          pairwise_effect_samples(sample_index, cntr) = pairwise_effects(row, col);
          cntr++;
        }
      }

      if (difference_selection) {
        int cntr = 0;
        for (int i = 0; i < num_variables; ++i) {
          for (int j = i; j < num_variables; ++j) {
            indicator_samples(sample_index, cntr) = inclusion_indicator(i, j);
            cntr++;
          }
        }
      }
    }
  }

  SamplerOutput out;
  out.chain_id = chain_id;
  out.main_samples = main_effect_samples;
  out.pairwise_samples = pairwise_effect_samples;
  out.treedepth_samples = treedepth_samples;
  out.divergent_samples = divergent_samples;
  out.energy_samples = energy_samples;
  out.has_indicator = difference_selection;
  if (difference_selection) {
    out.indicator_samples = indicator_samples;
  } else {
    out.indicator_samples = arma::imat();
  }
  return out;
}