#' @export
bgmCompare2 = function(
    x,
    y,
    group_indicator,
    difference_selection = FALSE,
    variable_type = "ordinal",
    reference_category,
    difference_scale = 1,
    difference_prior = c("Bernoulli", "Beta-Bernoulli"),
    difference_probability = 0.5,
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1,
    interaction_scale = 2.5,
    threshold_alpha = 0.5,
    threshold_beta = 0.5,
    iter = 1e3,
    burnin = 1e3,
    na_action = c("listwise", "impute"),
    display_progress = TRUE,
    update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
    target_accept,
    hmc_num_leapfrogs = 100,
    nuts_max_depth = 10,
    learn_mass_matrix = FALSE,
    chains = 4,
    cores = parallel::detectCores(),
    seed = NULL
) {
  # Check update method
  update_method_input = update_method
  update_method = match.arg(update_method)

  # Check target acceptance rate
  if(hasArg(target_accept)) {
    target_accept = min(target_accept, 1 - sqrt(.Machine$double.eps))
    target_accept = max(target_accept, 0 + sqrt(.Machine$double.eps))
  } else {
    if(update_method == "adaptive-metropolis") {
      target_accept = 0.44
    } else if(update_method == "hamiltonian-mc") {
      target_accept = 0.65
    } else if(update_method == "nuts") {
      target_accept = 0.60
    }
  }

  # Check and preprocess data
  x = data_check(x, "x")
  if (hasArg(y)) {
    y = data_check(y, "y")
    if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }

  if(!hasArg(y) & !hasArg(group_indicator))
    stop(paste0("For multi-group designs, the bgmCompare function requires input for\n",
                "either y (group 2 data) or group_indicator (group indicator)."))

  # Validate group indicators
  if (!hasArg(y) && hasArg(group_indicator)) {
    group_indicator = as.vector(group_indicator)
    if (anyNA(group_indicator)) stop("group_indicator cannot contain missing values.")
    if (length(group_indicator) != nrow(x)) stop("Length of group_indicator must match number of rows in x.")
  }

  # Model and preprocessing
  if(!hasArg(y))
    y = NULL
  if(!hasArg(group_indicator))
    group_indicator = NULL

  model = check_compare2_model(
    x = x, y = y, g = group_indicator, difference_selection = difference_selection,
    variable_type = variable_type, reference_category = reference_category,
    difference_scale = difference_scale, difference_prior = difference_prior,
    difference_probability = difference_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    interaction_scale = interaction_scale, threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta
  )

  x = model$x
  group = model$group
  ordinal_variable = model$variable_bool
  reference_category = model$reference_category

  # Check Gibbs input
  check_positive_integer(iter, "iter")
  check_non_negative_integer(burnin, "burnin")

  # Check na_action
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if (inherits(na_action, "try-error")) {
    stop(sprintf("Invalid value for `na_action`. Expected 'listwise' or 'impute', got: %s", na_action_input))
  }

  # Check display_progress
  display_progress = check_logical(display_progress, "display_progress")

  ## Format data
  data = compare2_reformat_data(
    x = x, group = group,
    na_action = na_action,
    variable_bool = ordinal_variable,
    reference_category = reference_category
  )

  x = data$x
  group = data$group
  num_obs_groups = tabulate(group)
  missing_index = data$missing_index
  num_categories = data$num_categories

  na_impute = data$na_impute
  reference_category = data$reference_category
  num_variables = ncol(x)
  num_interactions = num_variables * (num_variables - 1) / 2

  # Compute `num_obs_categories`
  num_obs_categories = compute_num_obs_categories(
    x, num_categories, group
  )

  # Compute sufficient statistics for Blume-Capel variables
  sufficient_blume_capel = compute_sufficient_blume_capel(
    x, reference_category, ordinal_variable, group
  )

  # Compute sufficient statistics for pairwise interactions
  sufficient_pairwise = compute_sufficient_pairwise(
    x, group
  )


  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0, nrow = num_interactions, ncol = 3)
  counter = 0
  for(variable1 in 1:(num_variables - 1)) {
    for(variable2 in (variable1 + 1):num_variables) {
      counter =  counter + 1
      Index[counter, ] = c(counter, variable1 - 1, variable2 - 1)
    }
  }

  # Gibbs sampling
  # Prepare indices for main and pairwise effects
  main_effect_indices = matrix(NA, nrow = num_variables, ncol = 2)
  for (variable in seq_len(num_variables)) {
    if (variable > 1) {
      main_effect_indices[variable, 1] = 1 + main_effect_indices[variable - 1, 2]
    } else {
      main_effect_indices[variable, 1] = 0  # C++ starts at zero
    }
    if (ordinal_variable[variable]) {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + max(num_categories[variable, ]) - 1
    } else {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + 1
    }
  }

  pairwise_effect_indices = matrix(NA, nrow = num_variables, ncol = num_variables)
  tel = 0
  for (v1 in seq_len(num_variables - 1)) {
    for (v2 in seq((v1 + 1), num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1  # C++ starts at zero
    }
  }

  # Compute group-level data
  num_groups = length(unique(group))
  group_indices = matrix(NA, nrow = num_groups, ncol = 2)

  # Align observations with sorted group
  observations = x
  sorted_group = sort(group)
  for (g in unique(group)) {
    observations[which(sorted_group == g), ] = x[which(group == g), ]
    group_indices[g, 1] = min(which(sorted_group == g)) - 1  # C++ starts at zero
    group_indices[g, 2] = max(which(sorted_group == g)) - 1  # C++ starts at zero
  }

  # Compute projection matrix for group differences
  one = matrix(1, nrow = num_groups, ncol = num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups]
  if (num_groups == 2) {
    projection = matrix(projection, ncol = 1) / sqrt(2)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || any(is.na(seed)) || any(seed < 0)) {
      stop("Argument 'seed' must be a non-negative integer or vector of non-negative integers.")
    }
    # Force to integer type
    seed <- as.integer(seed)
    dqrng::dqset.seed(seed)
  }

  # Call the Rcpp function
  out = run_bgmCompare_parallel(
    observations = observations,
    num_groups = num_groups,
    num_obs_categories = num_obs_categories,
    sufficient_blume_capel = sufficient_blume_capel,
    sufficient_pairwise = sufficient_pairwise,
    num_categories = num_categories[, 1],
    main_alpha = threshold_alpha,
    main_beta = threshold_beta,
    pairwise_scale = interaction_scale,
    difference_scale = difference_scale,
    difference_selection_alpha = beta_bernoulli_alpha,
    difference_selection_beta = beta_bernoulli_beta,
    difference_prior = model$difference_prior, iter = iter, burnin = burnin,
    na_impute = na_impute, missing_data_indices = missing_index,
    is_ordinal_variable = ordinal_variable,
    baseline_category = reference_category,
    difference_selection = difference_selection,
    main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices,
    target_accept = target_accept,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    projection = projection,
    group_membership = sorted_group - 1,  ######################################
    group_indices = group_indices,
    interaction_index_matrix = Index,
    inclusion_probability = model$inclusion_probability_difference,
    num_chains = chains, nThreads = cores
  )

  # out = run_gibbs_sampler_for_bgmCompare(
  #   chain_id = 1,
  #   observations = observations,
  #   num_groups = num_groups,
  #   num_obs_categories = num_obs_categories,
  #   sufficient_blume_capel = sufficient_blume_capel,
  #   sufficient_pairwise = sufficient_pairwise,
  #   num_categories = num_categories[, 1],
  #   main_alpha = threshold_alpha,
  #   main_beta = threshold_beta,
  #   pairwise_scale = interaction_scale,
  #   difference_scale = difference_scale,
  #   difference_selection_alpha = beta_bernoulli_alpha,
  #   difference_selection_beta = beta_bernoulli_beta,
  #   difference_prior = model$difference_prior, iter = iter, burnin = burnin,
  #   na_impute = na_impute, missing_data_indices = missing_index,
  #   is_ordinal_variable = ordinal_variable,
  #   baseline_category = reference_category,
  #   difference_selection = difference_selection,
  #   main_effect_indices = main_effect_indices,
  #   pairwise_effect_indices = pairwise_effect_indices,
  #   target_accept = target_accept,
  #   nuts_max_depth = nuts_max_depth,
  #   learn_mass_matrix = learn_mass_matrix,
  #   projection = projection,
  #   group_membership = sorted_group - 1,  ######################################
  #   group_indices = group_indices,
  #   interaction_index_matrix = Index,
  #   inclusion_probability = model$inclusion_probability_difference)


  # Main output handler in the wrapper function
  # output = prepare_output_bgmCompare2(
  #   out = out, ...
  # )

  return(out)
}