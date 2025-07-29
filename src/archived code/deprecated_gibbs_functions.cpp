// Deprecated: Former function to update main effect parameters of ordinal variables.
// Retained for reference and archival purposes only. Not included in build.
/**
 * Function: update_regular_thresholds_with_metropolis
 *
 * Performs a Metropolis-Hastings update for each threshold of an ordinal variable.
 * Uses a generalized beta-prime proposal and logistic-Beta prior.
 *
 * Inputs:
 *  - main_effects: Matrix of thresholds (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Vector of number of categories per variable.
 *  - num_obs_categories: Count matrix of observed scores.
 *  - num_persons: Number of individuals.
 *  - variable: Index of the variable being updated.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores for each observation and variable.
 *
 * Modifies:
 *  - main_effects (only for the specified variable)
 */
void update_regular_thresholds_with_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const int num_persons,
    const int variable,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix
) {
  arma::vec g(num_persons);
  arma::vec q(num_persons);
  const int num_cats = num_categories(variable);

  for (int category = 0; category < num_cats; category++) {
    double current = main_effects(variable, category);
    double exp_current = std::exp (current);
    double c = (threshold_alpha + threshold_beta) / (1.0 + exp_current);

    for (int person = 0; person < num_persons; person++) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;
      double numer = std::exp ((category + 1) * rest_score);

      for (int cat = 0; cat < num_cats; cat++) {
        if (cat != category) {
          denom += std::exp (main_effects(variable, cat) + (cat + 1) * rest_score);
        }
      }

      g(person) = denom;
      q(person) = numer;
      c += numer / (denom + numer * exp_current);
    }

    c /= (num_persons + threshold_alpha + threshold_beta - exp_current * c);

    // Sample from generalized beta-prime proposal
    double a = num_obs_categories(category + 1, variable) + threshold_alpha;
    double b = num_persons + threshold_beta - num_obs_categories(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed = std::log (tmp / (1.0 - tmp) / c);
    double exp_proposed = std::exp (proposed);

    // Compute MH acceptance probability
    double log_acceptance_probability = 0.0;
    for (int person = 0; person < num_persons; person++) {
      log_acceptance_probability += std::log (g(person) + q(person) * exp_current);
      log_acceptance_probability -= std::log (g(person) + q(person) * exp_proposed);
    }

    log_acceptance_probability -= (threshold_alpha + threshold_beta) * std::log1p(exp_proposed);
    log_acceptance_probability += (threshold_alpha + threshold_beta) * std::log1p(exp_current);
    log_acceptance_probability -= (a + b) * std::log1p(c * exp_current);
    log_acceptance_probability += (a + b) * std::log1p(c * exp_proposed);

    if (std::log (R::unif_rand()) < log_acceptance_probability) {
      main_effects(variable, category) = proposed;
    }
  }
}



// Deprecated: Former function to update main effect parameters of Blume-Capel variables.
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: update_blumecapel_thresholds_with_adaptive_metropolis
 *
 * Performs an adaptive Metropolis update of the Blume-Capel threshold parameters
 * (linear and quadratic) for a single variable, with Robbins-Monro tuning.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - num_persons: Number of observations.
 *  - variable: Index of the variable being updated.
 *  - reference_category: Reference category per variable.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores.
 *  - proposal_sd_blumecapel: Matrix of proposal SDs for each variable (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *
 * Modifies:
 *  - main_effects (for the given variable)
 *  - proposal_sd_blumecapel
 */
void update_blumecapel_thresholds_with_adaptive_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& sufficient_blume_capel,
    const int num_persons,
    const int variable,
    const arma::ivec& reference_category,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_blumecapel,
    const double exp_neg_log_t_rm_adaptation_rate
) {
  const int num_cats = num_categories(variable);
  const int ref = reference_category(variable);

  // --- Define helper for prior contribution
  auto log_beta_prior_diff = [&](double curr, double prop) {
    return (threshold_alpha + threshold_beta) *
      (std::log1p(std::exp (curr)) - std::log1p(std::exp (prop)));
  };

  // --- Update each threshold parameter: 0 = linear, 1 = quadratic
  for (int param = 0; param < 2; param++) {
    double& proposal_sd = proposal_sd_blumecapel(variable, param);
    double current = main_effects(variable, param);
    double proposed = R::rnorm(current, proposal_sd);
    double diff = proposed - current;

    arma::vec numer_current(num_cats + 1);
    arma::vec numer_proposed(num_cats + 1);

    // --- Step 1: Construct numerators for softmax (for all categories)
    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;
      if (param == 0) {
        // Linear update
        double quad = main_effects(variable, 1) * centered * centered;
        numer_current(cat) = current * cat + quad;
        numer_proposed(cat) = proposed * cat + quad;
      } else {
        // Quadratic update
        double lin = main_effects(variable, 0) * cat;
        numer_current(cat) = current * centered * centered + lin;
        numer_proposed(cat) = proposed * centered * centered + lin;
      }
    }

    // --- Step 2: Compute lbound for numerical stability
    double max_curr = numer_current.max();
    double max_prop = numer_proposed.max();
    double lbound = (max_curr > 0.0 || max_prop > 0.0) ? std::max(max_curr, max_prop) : 0.0;

    // --- Step 3: Likelihood ratio
    // Accumulate log acceptance probability based on change in pseudo-likelihood

    // Contribution from sufficient statistics and prior
    double log_accept = diff * (threshold_alpha + sufficient_blume_capel(param, variable));

    // Vectorized likelihood ratio for all persons:
    //
    // For each person, compute:
    //   log p(y_i | proposed) - log p(y_i | current)
    //   using softmax-style normalization for categorical probabilities.
    //
    // The bound stabilizes exponentials across categories and persons.
    arma::vec rest_score = residual_matrix.col(variable);                       // Person-wise residuals
    arma::vec bound = arma::max(rest_score, arma::zeros<arma::vec>(num_persons)) * num_cats + lbound;

    arma::vec denom_curr = arma::exp (numer_current(0) - bound);                 // Score = 0 contribution
    arma::vec denom_prop = arma::exp (numer_proposed(0) - bound);

    for (int cat = 0; cat < num_cats; cat++) {
      arma::vec score_term = (cat + 1) * rest_score - bound;

      // Compute exponentials for each category and add to denominator
      denom_curr += arma::exp (numer_current(cat + 1) + score_term);
      denom_prop += arma::exp (numer_proposed(cat + 1) + score_term);
    }

    // Accumulate the person-wise log ratio contributions
    log_accept += arma::accu (arma::log (denom_curr) - arma::log (denom_prop));

    // --- Step 4: Add prior ratio
    log_accept += log_beta_prior_diff(current, proposed);

    // --- Step 5: Metropolis accept/reject
    if (std::log (R::unif_rand()) < log_accept) {
      main_effects(variable, param) = proposed;
    }

    // --- Step 6: Robbins-Monro proposal adaptation
    proposal_sd = update_proposal_sd_with_robbins_monro (
      proposal_sd, log_accept, exp_neg_log_t_rm_adaptation_rate
    );
  }
}



// Deprecated: Former function to compute gradient vector for the interactions.
// Retained for reference and archival purposes only. Not included in build.

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
    const double interaction_scale
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  arma::vec gradient (num_interactions, arma::fill::zeros);
  int interaction_index = -1;

  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;

      if (inclusion_indicator (var1, var2) == 0)
        continue;

      // Convert observed scores to integer vectors
      const arma::ivec responses_var1 = arma::conv_to<arma::ivec>::from (observations.col (var1));
      const arma::ivec responses_var2 = arma::conv_to<arma::ivec>::from (observations.col (var2));


      // First-order gradient term from complete data
      gradient (interaction_index) = 2.0 * arma::dot (responses_var1, responses_var2);

      // --- Contribution from variable var1
      int num_categories_var1 = num_categories (var1);
      arma::vec rest_scores = observations * pairwise_effects.col (var1);
      arma::vec numerator = arma::zeros (num_observations);
      arma::vec denominator = arma::zeros (num_observations);
      arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var1;

      if (is_ordinal_variable (var1)) {
        denominator += arma::exp ( -bounds );
        for (int category = 0; category < num_categories_var1; category++) {
          arma::vec exponent = main_effects (var1, category) + (category + 1) * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += (category + 1) * responses_var2 % weight;
        }
      } else {
        const int ref_cat = reference_category (var1);
        for (int category = 0; category <= num_categories_var1; category++) {
          int centered_cat = category - ref_cat;
          double lin_term = main_effects (var1, 0) * category;
          double quad_term = main_effects (var1, 1) * centered_cat * centered_cat;
          arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += category * responses_var2 % weight;
        }
      }

      gradient (interaction_index) -= arma::accu (numerator / denominator);

      // --- Contribution from variable var2
      int num_categories_var2 = num_categories (var2);
      rest_scores = observations * pairwise_effects.col (var2);
      numerator.zeros ();
      denominator.zeros ();
      bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var2;

      if (is_ordinal_variable (var2)) {
        denominator += arma::exp ( -bounds );
        for (int category = 0; category < num_categories_var2; category++) {
          arma::vec exponent = main_effects (var2, category) + (category + 1) * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += (category + 1) * responses_var1 % weight;
        }
      } else {
        const int ref_cat = reference_category (var2);
        for (int category = 0; category <= num_categories_var2; category++) {
          int centered_cat = category - ref_cat;
          double lin_term = main_effects (var2, 0) * category;
          double quad_term = main_effects (var2, 1) * centered_cat * centered_cat;
          arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += category * responses_var1 % weight;
        }
      }

      gradient (interaction_index) -= arma::accu (numerator / denominator);

      // ---- Gradient contribution from Cauchy prior
      const double effect = pairwise_effects (var1, var2);
      gradient (interaction_index) -= 2.0 * effect / (effect * effect + interaction_scale * interaction_scale);
    }
  }

  return gradient;
}



// Deprecated: Former threshold update function from Gibbs sampler
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: update_thresholds_with_adaptive_mala
 *
 * Performs a MALA update of threshold parameters with adaptive step size tuning.
 * Applies dual averaging during burn-in and Robbins-Monro afterward.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - step_size_mala: MALA step size (updated in-place).
 *  - residual_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - iteration: Current iteration number.
 *  - burnin: Total number of burn-in iterations.
 *  - dual_averaging_state: Dual averaging state vector [log_eps, log_eps_avg, H_bar] (updated in-place).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Modifies:
 *  - main_effects
 *  - step_size_mala
 *  - dual_averaging_state
 */
void update_thresholds_with_adaptive_mala (
    arma::mat& main_effects,
    double& step_size_mala,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int iteration,
    const int burnin,
    arma::vec& dual_averaging_state,
    const double threshold_alpha,
    const double threshold_beta,
    const double initial_step_size_mala
) {
  // --- Step 1: Flatten current parameters and compute gradient & posterior
  arma::vec flat_theta = vectorize_thresholds (
    main_effects, num_categories, is_ordinal_variable
  );
  arma::vec grad = gradient_log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  const double log_post_current = log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Step 2: Propose new parameters using MALA
  const double sqrt_step = std::sqrt(step_size_mala);
  arma::vec proposal = flat_theta + 0.5 * step_size_mala * grad + sqrt_step * arma::randn(flat_theta.n_elem);
  arma::mat proposed_thresholds = reconstruct_threshold_matrix (
    proposal, num_categories, is_ordinal_variable
  );

  // --- Step 3: Evaluate proposed state
  const double log_post_proposal = log_pseudoposterior_thresholds (
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  arma::vec grad_proposal = gradient_log_pseudoposterior_thresholds (
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Step 4: Compute forward and backward proposal densities
  const arma::vec forward_mean = flat_theta + 0.5 * step_size_mala * grad;
  const arma::vec backward_mean = proposal + 0.5 * step_size_mala * grad_proposal;

  const double log_forward = -0.5 / step_size_mala * arma::accu(arma::square(proposal - forward_mean));
  const double log_backward = -0.5 / step_size_mala * arma::accu(arma::square(flat_theta - backward_mean));

  // --- Step 5: Accept/reject
  const double log_acceptance = log_post_proposal + log_backward - log_post_current - log_forward;
  if (std::log(R::unif_rand()) < log_acceptance) {
    main_effects = proposed_thresholds;
  }

  const double accept_prob = std::min(1.0, std::exp(log_acceptance));

  // --- Step 6: Adapt step size
  if (iteration <= burnin) {
    update_step_size_with_dual_averaging (
        initial_step_size_mala, accept_prob, iteration + 1, dual_averaging_state);
    step_size_mala = std::exp(dual_averaging_state[1]);
  } else {
    update_step_size_with_robbins_monro(accept_prob, iteration - burnin, step_size_mala);
  }
}

update_thresholds_with_adaptive_mala(
  main_effects, step_size_mala, residual_matrix, num_categories,
  num_obs_categories, sufficient_blume_capel, reference_category,
  is_ordinal_variable, iteration, total_burnin, dual_averaging_state,
  threshold_alpha, threshold_beta, initial_step_size_mala
);


// Deprecated: Former function to update pairwise effect parameters.
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: update_interactions_with_mala
 *
 * Performs a blockwise MALA update of the interaction parameters using a fixed
 * or adaptive step size. Uses full-length vectors over all pairwise interactions,
 * with inactive elements set to zero.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters (updated in-place).
 *  - residual_matrix: Residual matrix (updated in-place if proposal accepted).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of observed scores.
 *  - num_categories: Vector of number of categories per variable.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - is_ordinal_variable: Logical vector: 1 = ordinal, 0 = Blume-Capel.
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Scale parameter for Cauchy prior on interactions.
 *  - step_size_interactions: Current step size (updated if adaptive).
 *  - initial_step_size_interactions: Initial step size used during dual averaging.
 *  - iteration: Current MCMC iteration (0-based).
 *  - total_burnin: Total number of burn-in iterations.
 *  - dual_averaging_state: Vector (length 3) tracking dual averaging state.
 *
 * Modifies:
 *  - pairwise_effects
 *  - residual_matrix
 *  - step_size_interactions (if in burn-in or adaptation phase)
 *  - dual_averaging_state (during burn-in)
 */
void update_interactions_with_mala (
    arma::mat& pairwise_effects,
    arma::mat& residual_matrix,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    double& step_size_interactions,
    const double initial_step_size_interactions,
    const int iteration,
    const int total_burnin,
    arma::vec& dual_averaging_state,
    const double target_accept_interactions
) {
  const int num_variables = pairwise_effects.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  // --- Flatten current interaction matrix to vector
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

  // --- Compute gradient and posterior at current state
  arma::vec grad_current = gradient_log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  double log_post_current = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  // --- Generate Langevin proposal
  arma::vec noise = arma::randn<arma::vec> (num_interactions);
  arma::vec proposal = current_state +
    0.5 * step_size_interactions * grad_current +
    std::sqrt (step_size_interactions) * noise;

  // --- Build symmetric proposal matrix
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

  // --- Evaluate posterior and gradient at proposed state
  double log_post_proposal = log_pseudoposterior_interactions (
    proposal_matrix, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  arma::vec grad_proposed = gradient_log_pseudoposterior_interactions (
    proposal_matrix, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  // --- Compute MH log acceptance ratio
  arma::vec forward_mean = current_state + 0.5 * step_size_interactions * grad_current;
  arma::vec reverse_mean = proposal + 0.5 * step_size_interactions * grad_proposed;

  double log_q_forward = -0.5 / step_size_interactions * arma::accu (arma::square (proposal - forward_mean));
  double log_q_reverse = -0.5 / step_size_interactions * arma::accu (arma::square (current_state - reverse_mean));

  double log_accept = log_post_proposal + log_q_reverse - log_post_current - log_q_forward;

  double logu = std::log (R::runif (0.0, 1.0));
  bool accepted = (logu < log_accept);
  double acceptance_prob = std::min (1.0, std::exp (log_accept));

  // --- Adapt step size
  if (iteration < 500) {
    update_step_size_with_dual_averaging (
        initial_step_size_interactions,
        acceptance_prob,
        iteration + 1,
        dual_averaging_state,
        target_accept_interactions
    );
    step_size_interactions = std::exp (dual_averaging_state (1));
  } else if(iteration < total_burnin) {
    update_step_size_with_robbins_monro (
        acceptance_prob,
        iteration - 500 + 1,
        step_size_interactions,
        target_accept_interactions
    );
  }

  // --- Accept and update interaction matrix and residuals
  if (accepted) {
    interaction_index = -1;
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
      for (int var2 = var1 + 1; var2 < num_variables; var2++) {
        interaction_index++;
        if (inclusion_indicator (var1, var2) == 1) {
          double delta = proposal (interaction_index) - pairwise_effects (var1, var2);
          pairwise_effects (var1, var2) = proposal (interaction_index);
          pairwise_effects (var2, var1) = proposal (interaction_index);
          residual_matrix.col (var1) += arma::conv_to<arma::vec>::from (observations.col (var2)) * delta;
          residual_matrix.col (var2) += arma::conv_to<arma::vec>::from (observations.col (var1)) * delta;
        }
      }
    }
  }
}



// Deprecated: Former function to find step size for componentwise mala update pairwise effect parameters.
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function
 *
 * Inputs:
 *
 * Returns:
 */
double find_reasonable_initial_step_size_single_interaction(
    int var1, int var2,
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    double interaction_scale,
    const double target_acceptance
) {

  const double initial_step_size = 0.1;
  double log_step_size = std::log (initial_step_size);
  const int max_attempts = 20;

  int direction = 0;
  double accept_prob = 0.0;

  // --- Step 1: Extract current interaction state
  double current_state = pairwise_effects(var1, var2);

  // --- Step 2: Evaluate gradient and posterior at current state
  double grad = gradient_log_pseudoposterior_interaction_single(
    var1, var2, pairwise_effects, main_effects, observations,
    num_categories, is_ordinal_variable, reference_category, interaction_scale
  );
  double log_post_current = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale
  );


  // --- Step 3: Exponential step-size search loop
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    double step_size = std::exp (log_step_size);

    // Generate Langevin proposal
    double noise = R::rnorm(0.0, 1.0);
    double proposal = current_state + 0.5 * step_size * grad +
      std::sqrt(step_size) * noise;


    // Compute reverse proposal terms
    pairwise_effects(var1, var2) = proposal;
    pairwise_effects(var2, var1) = proposal;

    // Evaluate posterior and gradient at proposed state
    double log_post_prop = log_pseudoposterior_interactions (
      pairwise_effects, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale
    );

    double grad_prop = gradient_log_pseudoposterior_interaction_single(
      var1, var2, pairwise_effects, main_effects, observations,
      num_categories, is_ordinal_variable, reference_category, interaction_scale
    );

    pairwise_effects(var1, var2) = current_state;
    pairwise_effects(var2, var1) = current_state;

    // Compute log proposal densities
    double forward_mean = current_state + 0.5 * step_size * grad;
    double backward_mean = proposal + 0.5 * step_size * grad_prop;

    double log_q_forward = -0.5 / step_size * std::pow(proposal - forward_mean, 2);
    double log_q_reverse = -0.5 / step_size * std::pow(current_state - backward_mean, 2);

    double log_accept = log_post_prop + log_q_reverse - log_post_current - log_q_forward;

    accept_prob = std::min(1.0, std::exp(log_accept));

    // --- Step 4: Decide direction based on first attempt
    if (attempt == 0) {
      direction = (accept_prob > target_acceptance) ? 1 : -1;
    } else {
      if ((direction == 1 && accept_prob < target_acceptance) ||
          (direction == -1 && accept_prob > target_acceptance)) {
        break;
      }
    }

    log_step_size += direction;
  }

  return std::exp (log_step_size);
}



// Deprecated: Former function for componentwise mala update pairwise effect parameters.
// Retained for reference and archival purposes only. Not included in build.
/**
 * Function
 *
 * Inputs:
 *
 * Returns:
 */
void update_interactions_with_componentwise_mala (
    const arma::imat& pairwise_effect_indices,
    arma::mat& pairwise_effects,
    arma::mat& residual_matrix,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    arma::vec& component_wise_interactions_step_sizes,
    arma::mat& componentwise_dual_averaging_state,
    const double initial_step_size,
    const int iteration,
    const int total_burnin,
    const double target_accept_interactions
) {
  const int num_variables = pairwise_effects.n_rows;

  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0) continue;

      const int interaction_id = pairwise_effect_indices(var1, var2);

      double& beta = pairwise_effects(var1, var2);
      double grad = gradient_log_pseudoposterior_interaction_single(
        var1, var2, pairwise_effects, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );

      double& step_size = component_wise_interactions_step_sizes(interaction_id);
      double sqrt_step = std::sqrt(step_size);
      double noise = R::rnorm(0.0, 1.0);

      double proposal = beta + 0.5 * step_size * grad + sqrt_step * noise;

      double log_post_diff = log_pseudolikelihood_ratio_interaction(
        pairwise_effects, main_effects, observations, num_categories,
        observations.n_rows, var1, var2, proposal, beta, residual_matrix,
        is_ordinal_variable, reference_category
      );

      double grad_prop = gradient_log_pseudoposterior_interaction_single(
        var1, var2, pairwise_effects, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );

      double forward = beta + 0.5 * step_size * grad;
      double backward = proposal + 0.5 * step_size * grad_prop;
      double log_q_forward = -0.5 / step_size * std::pow(proposal - forward, 2);
      double log_q_reverse = -0.5 / step_size * std::pow(beta - backward, 2);
      double log_accept = log_post_diff + log_q_reverse - log_q_forward;

      if (std::log(R::unif_rand()) < log_accept) {
        double delta = proposal - beta;
        beta = proposal;
        pairwise_effects(var2, var1) = proposal;

        // Update residuals
        residual_matrix.col(var1) += arma::conv_to<arma::vec>::from(observations.col(var2)) * delta;
        residual_matrix.col(var2) += arma::conv_to<arma::vec>::from(observations.col(var1)) * delta;
      }

      double accept_prob = std::min(1.0, std::exp(log_accept));

      if (iteration < 500) {
        arma::rowvec state = componentwise_dual_averaging_state.row(interaction_id);
        arma::vec state_vec = state.t();  // Convert to column vec for dual averaging

        update_step_size_with_dual_averaging(
          initial_step_size, accept_prob, iteration + 1, state_vec,
          target_accept_interactions);

        step_size = std::exp(state_vec[1]);
        componentwise_dual_averaging_state.row(interaction_id) = state_vec.t();  // Store back as row
      } else if (iteration < total_burnin){
        update_step_size_with_robbins_monro(
          accept_prob, iteration - 500 + 1, step_size,
          target_accept_interactions);
      }

    }
  }
}




// Deprecated: Former function for mala update pairwise effect - indicator
// Retained for reference and archival purposes only. Not included in build.
/**
 * Function: update_indicator_interaction_pair_with_mala
 *
 * Proposes and accepts/rejects inclusion of pairwise interaction terms in a sparse graphical model
 * using a Metropolis-adjusted Langevin algorithm (MALA). Operates one pair at a time across a list
 * of possible interactions.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of pairwise interaction effects (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - indicator: Binary matrix indicating inclusion of pairwise interaction terms (updated in-place).
 *  - observations: Matrix of observed scores.
 *  - num_categories: Vector with number of categories per variable.
 *  - step_size_pairwise: Step size used for MALA proposals.
 *  - interaction_scale: Scale parameter of the Cauchy prior on interaction effects.
 *  - index: Matrix listing all candidate pairs for interaction updates.
 *  - num_interactions: Number of candidate interaction pairs to consider.
 *  - num_persons: Number of observations (individuals).
 *  - residual_matrix: Matrix of residuals (updated if proposal accepted).
 *  - inclusion_probability: Matrix of prior inclusion probabilities for pairwise interactions.
 *  - is_ordinal_variable: Indicator vector specifying which variables are ordinal.
 *  - reference_category: Vector of reference categories for binary/categorical variables.
 *
 * Modifies:
 *  - pairwise_effects
 *  - indicator
 *  - residual_matrix (only if proposal accepted)
 */
void update_indicator_interaction_pair_with_mala (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double step_size_pairwise,
    const double interaction_scale,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const arma::vec& component_wise_interactions_step_sizes,
    const std::string& update_method_interactions,
    const arma::imat& pairwise_effect_indices
) {
  for (int cntr = 0; cntr < num_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    double step_size = (update_method_interactions == "adaptive-mala")
      ? step_size_pairwise
    : component_wise_interactions_step_sizes(pairwise_effect_indices(variable1, variable2));

    // Determine if we are proposing to add (if currently absent)
    const bool proposing_addition = (indicator(variable1, variable2) == 0);
    const double current_state = pairwise_effects(variable1, variable2);
    double proposed_state = 0.0;
    double log_accept = 0.0;
    const double inclusion_probability_ij = inclusion_probability(variable1, variable2);

    if (proposing_addition) {
      // Compute gradient of log-pseudo-posterior for interaction term
      double grad = gradient_log_pseudoposterior_interaction_single (
        variable1, variable2, pairwise_effects, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );

      // MALA proposal: Langevin step forward
      double proposal_sd = std::sqrt(step_size);
      double noise = R::rnorm(0.0, proposal_sd);
      double forward_mean = current_state + 0.5 * step_size * grad;
      proposed_state = forward_mean + noise;
      log_accept -= R::dnorm(proposed_state, forward_mean, proposal_sd, true);

      // Cauchy prior on interaction effect
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);

      // Prior inclusion probability
      log_accept += std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    } else {

      // MALA proposal: Langevin step backward
      arma::mat proposed_matrix = pairwise_effects;
      proposed_matrix(variable1, variable2) = proposed_state;
      proposed_matrix(variable2, variable1) = proposed_state;
      double grad = gradient_log_pseudoposterior_interaction_single (
        variable1, variable2, proposed_matrix, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );
      double proposal_sd = std::sqrt(step_size);
      double backward_mean = proposed_state + 0.5 * step_size * grad;
      log_accept += R::dnorm(current_state, backward_mean, proposal_sd, true);

      // Cauchy prior on interaction effect
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      // Prior inclusion probability
      log_accept -= std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    }

    log_accept += log_pseudolikelihood_ratio_interaction (
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      variable1, variable2, proposed_state, current_state, residual_matrix,
      is_ordinal_variable, reference_category
    );

    // Metropolis-Hastings accept step
    if (std::log (R::unif_rand()) < log_accept) {
      const int new_value = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = new_value;
      indicator(variable2, variable1) = new_value;

      const double delta = proposed_state - current_state;
      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}

// Deprecated version of the Gibbs update step for graphical model parameters.
// Retained for reference and archival purposes only. Not included in build.
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
    arma::mat& residual_matrix,
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
    const std::string& update_method_interactions,
    const std::string& update_method_thresholds,
    arma::vec& component_wise_interactions_step_sizes,
    arma::mat& componentwise_dual_averaging_state,
    const arma::imat& pairwise_effect_indices,
    const double target_accept_interactions,
    const double target_accept_thresholds,
    const int warmup_stageI_end,
    const int warmup_stageII_end,
    const int warmup_stageIII_end
) {
  // --- Robbins-Monro weight for adaptive Metropolis updates
  const double exp_neg_log_t_rm_adaptation_rate =
    std::exp (-std::log (static_cast<double>(iteration)) * rm_decay_rate);

  // Step 1: Edge selection via MH indicator updates (if enabled)
  if (edge_selection) {
    if (update_method_interactions == "fisher-mala") {
      // Use Fisher-preconditioned MALA for inclusion indicators
      update_indicator_interaction_pair_with_fisher_mala (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, step_size_pairwise, interaction_scale, index,
          num_persons, residual_matrix, inclusion_probability, is_ordinal_variable,
          reference_category, sqrt_inv_fisher_pairwise, num_pairwise,
          iteration, total_burnin
      );
    } else if (update_method_interactions == "adaptive-mala") {
      // Use standard MALA for indicator updates
      update_indicator_interaction_pair_with_mala (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, step_size_pairwise, interaction_scale,
          index, num_pairwise, num_persons, residual_matrix,
          inclusion_probability, is_ordinal_variable, reference_category,
          component_wise_interactions_step_sizes, update_method_interactions,
          pairwise_effect_indices
      );
    } else if (update_method_interactions == "adaptive-componentwise-mala") {
      // Use standard MALA for indicator updates
      update_indicator_interaction_pair_with_mala (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, step_size_pairwise, interaction_scale,
          index, num_pairwise, num_persons, residual_matrix,
          inclusion_probability, is_ordinal_variable, reference_category,
          component_wise_interactions_step_sizes, update_method_interactions,
          pairwise_effect_indices
      );
    } else {
      // Use standard Metropolis-Hastings for indicator updates
      update_indicator_interaction_pair_with_metropolis (
          pairwise_effects, main_effects, inclusion_indicator, observations,
          num_categories, proposal_sd_pairwise, interaction_scale,
          index, num_pairwise, num_persons, residual_matrix,
          inclusion_probability, is_ordinal_variable, reference_category
      );
    }
  }

  // Step 2: Update interaction weights for active edges
  if (update_method_interactions == "fisher-mala") {
    update_interactions_with_fisher_mala (
        pairwise_effects, residual_matrix, main_effects, observations,
        num_categories, inclusion_indicator, is_ordinal_variable,
        reference_category, interaction_scale, step_size_pairwise,
        initial_step_size_pairwise, iteration, //total_burnin,
        dual_averaging_pairwise, sqrt_inv_fisher_pairwise,
        target_accept_interactions,
        warmup_stageI_end,
        warmup_stageII_end,
        warmup_stageIII_end
    );
  } else if (update_method_interactions == "adaptive-mala") {
    update_interactions_with_mala (
        pairwise_effects, residual_matrix, main_effects, observations,
        num_categories, inclusion_indicator, is_ordinal_variable,
        reference_category, interaction_scale,
        step_size_pairwise, initial_step_size_pairwise,
        iteration, total_burnin, dual_averaging_pairwise,
        target_accept_interactions
    );
  } else if (update_method_interactions == "adaptive-componentwise-mala") {
    update_interactions_with_componentwise_mala (
        pairwise_effect_indices, pairwise_effects, residual_matrix,
        main_effects, observations, num_categories, inclusion_indicator,
        is_ordinal_variable, reference_category, interaction_scale,
        component_wise_interactions_step_sizes,
        componentwise_dual_averaging_state, initial_step_size_pairwise,
        iteration, total_burnin, target_accept_interactions
    );
  } else {
    update_interactions_with_adaptive_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, interaction_scale,
        num_persons, num_variables, residual_matrix,
        exp_neg_log_t_rm_adaptation_rate, is_ordinal_variable,
        reference_category, target_accept_interactions, iteration, total_burnin
    );
  }

  // Step 3: Update main effect (threshold) parameters
  if (update_method_thresholds == "adaptive-metropolis") {
    // Update using adaptive Metropolis
    update_thresholds_with_adaptive_metropolis (
        main_effects, observations, num_categories, num_obs_categories,
        sufficient_blume_capel, reference_category, is_ordinal_variable,
        num_persons, threshold_alpha, threshold_beta, residual_matrix,
        proposal_sd_main, exp_neg_log_t_rm_adaptation_rate,
        target_accept_thresholds, iteration, total_burnin
    );
  } else {
    // Update using Fisher-preconditioned MALA
    update_thresholds_with_fisher_mala (
        main_effects, step_size_main, residual_matrix, num_categories,
        num_obs_categories, sufficient_blume_capel, reference_category,
        is_ordinal_variable, iteration, //total_burnin,
        dual_averaging_main,
        sqrt_inv_fisher_main, threshold_alpha, threshold_beta,
        initial_step_size_main, target_accept_thresholds,
        warmup_stageI_end,
        warmup_stageII_end,
        warmup_stageIII_end
    );
  }
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



/**
 * Function: update_indicator_interaction_pair_with_mala
 *
 * Updates the inclusion indicators and associated interaction weights using MALA.
 *
 * For each candidate interaction, this function proposes:
 *   - Adding an edge using a Fisher-preconditioned Langevin proposal
 *   - Removing an edge by proposing an effect value of zero
 *
 * The acceptance ratio includes:
 *   - Pseudolikelihood difference (pairwise-only)
 *   - Cauchy prior for interaction weights
 *   - Inclusion probability prior (Bernoulli)
 *   - Forward/reverse Langevin proposal correction
 *
 * Step size is preconditioned using an inverse mass term estimated during warmup.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters (modified in-place).
 *  - main_effects: Matrix of threshold (main effect) parameters.
 *  - indicator: Symmetric matrix of inclusion indicators (modified in-place).
 *  - observations: Person-by-variable matrix of observed scores.
 *  - num_categories: Number of categories per variable.
 *  - step_size: Global step size for MALA.
 *  - inv_mass_diag: Inverse diagonal mass vector (per interaction pair).
 *  - interaction_scale: Scale parameter for Cauchy prior on interaction weights.
 *  - index: Matrix listing candidate interactions: [index, var1, var2].
 *  - num_persons: Number of observations (rows in observations matrix).
 *  - rest_matrix: Linear predictor matrix (updated in-place on accept).
 *  - inclusion_probability: Prior inclusion probabilities for interaction pairs.
 *  - is_ordinal_variable: Indicator of ordinal variables.
 *  - reference_category: Reference category for each variable.
 *  - num_pairwise: Total number of candidate interactions.
 *  - sufficient_pairwise: Sufficient statistics for pairwise terms.
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
    const arma::vec& inv_mass_diag,
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
  for (int pair_index = 0; pair_index < num_pairwise; pair_index++) {
    const int var1 = index(pair_index, 1);
    const int var2 = index(pair_index, 2);
    const double current_state = pairwise_effects(var1, var2);
    const double inclusion_prob = inclusion_probability(var1, var2);

    const int true_index = index(pair_index, 0);
    const double inv_mass = inv_mass_diag(true_index);
    const double local_step_size = step_size * inv_mass;
    const double sd = std::sqrt(local_step_size);

    double proposed_state = 0.0;
    double log_accept = 0.0;

    if (indicator(var1, var2) == 0) {
      const double grad = gradient_log_pseudoposterior_interactions_component(
        var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      const double drift = 0.5 * local_step_size * grad;
      const double forward_mean = current_state + drift;

      proposed_state = forward_mean + R::rnorm(0.0, sd);
      log_accept -= R::dnorm(proposed_state, forward_mean, sd, true);
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept += std::log(inclusion_prob) - std::log(1.0 - inclusion_prob);
    } else {
      const double proposed_state = 0.0;
      const double tmp = pairwise_effects(var1, var2);

      pairwise_effects(var1, var2) = proposed_state;
      pairwise_effects(var2, var1) = proposed_state;

      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));
      const double delta = proposed_state - current_state;

      rest_matrix.col(var1) += obs_var2 * delta;
      rest_matrix.col(var2) += obs_var1 * delta;

      const double grad = gradient_log_pseudoposterior_interactions_component(
        var1, var2, pairwise_effects, main_effects, observations, rest_matrix,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale, sufficient_pairwise
      );

      const double drift = 0.5 * local_step_size * grad;
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

      const arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations.col(var1));
      const arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations.col(var2));

      rest_matrix.col(var1) += obs_var2 * delta;
      rest_matrix.col(var2) += obs_var1 * delta;
    }
  }
}