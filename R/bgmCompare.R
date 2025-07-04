#' Bayesian Variable Selection or Bayesian Estimation for Differences in Markov Random Fields
#'
#' @description
#' The \code{bgmCompare} function estimates the pseudoposterior distribution of
#' the parameters of a Markov Random Field (MRF) model for mixed binary and ordinal
#' variables, as well as differences in pairwise interactions and category thresholds
#' across groups. Groups are assumed to be \code{G} independent samples.
#'
#' @details
#' This function models group differences in Markov Random Fields using a Bayesian
#' framework. It supports binary and ordinal variables, and includes options for
#' Bayesian variable selection on the differences in both pairwise interactions
#' and threshold parameters. Key components are described in the sections below.
#'
#' @section Pairwise Interactions:
#' Pairwise interactions between variables \code{i} and \code{j} are modeled as:
#' \deqn{\boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}{
#' \boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}
#' where:
#'
#' \itemize{
#'  \item \eqn{\boldsymbol{\theta}_{ij}}{\boldsymbol{\theta}_{ij}} is the vector of pairwise interaction parameters of length \code{G}.
#'  \item \eqn{\phi_{ij}}{\phi_{ij}} is the overall pairwise interaction (nuisance parameter).
#'  \item \eqn{\boldsymbol{\delta}_{ij}}{\boldsymbol{\delta}_{ij}} represents group-specific differences constrained to sum to zero for identification.
#' }
#'
#' @section Ordinal Variables:
#' The function supports two types of ordinal variables:
#'
#' \strong{Regular ordinal variables}:
#' Introduce a threshold parameter for each category except the lowest, modeled as:
#' \deqn{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}
#' where:
#' \itemize{
#'  \item \eqn{\tau_{ic}}{\tau_{ic}} denotes an overall effect (nuisance parameter).
#'  \item \eqn{\boldsymbol{\epsilon}_{ic}}{\boldsymbol{\epsilon}_{ic}} represents group-specific differences constrained to sum to zero.
#' }
#'
#' \strong{Blume-Capel ordinal variables}:
#' Assume a specific reference category and score responses based on distance to it:
#' \deqn{\boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}{
#' \boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}
#' where:
#' \itemize{
#'  \item \code{r} is the reference category.
#'  \item \eqn{\tau_{i1}}{\tau_{i1}} and \eqn{\tau_{i2}}{\tau_{i2}} are nuisance parameters.
#'  \item \eqn{\boldsymbol{\epsilon}_{i1}}{\boldsymbol{\epsilon}_{i1}} and \eqn{\boldsymbol{\epsilon}_{i2}}{\boldsymbol{\epsilon}_{i2}} represent group-specific differences.
#' }
#'
#' @section Variable Selection:
#' Bayesian variable selection enables testing of parameter differences or equivalence across groups. Independent spike-and-slab priors are applied to difference parameters:
#' \itemize{
#'  \item \strong{Bernoulli Model}: Assigns a fixed probability to parameter inclusion.
#'  \item \strong{Beta-Bernoulli Model}: Incorporates a beta prior to model inclusion probabilities.
#' }
#'
#' @section Gibbs Sampling:
#'
#' Parameters are estimated using a Metropolis-within-Gibbs sampling scheme.
#' When \code{difference_selection = TRUE}, the algorithm runs \code{2 * burnin} warmup iterations:
#' \itemize{
#'   \item First half without difference selection.
#'   \item Second half with edge selection enabled.
#' }
#' This warmup strategy improves stability of adaptive Metropolis-Hastings proposals and starting values.
#'
#' @section Saving Options:
#' Users can store sampled states for parameters (\code{main_effects}, \code{pairwise_effects}, \code{indicator}) during Gibbs sampling. Enabling these options (\code{save_main}, \code{save_pairwise}, \code{save_indicator}) increases output size and memory usage, so use them judiciously.
#'
#' @param x Data frame or matrix with binary and ordinal responses. Regular ordinal variables should be coded as integers starting from 0. Missing categories are collapsed for regular ordinal variables but retained for Blume-Capel variables.
#' @param y A data frame or matrix similar to \code{x}, used for two-group designs. \code{x} contains Group 1 data, and \code{y} contains Group 2 data. Ignored for multi-group designs.
#' @param g Group membership vector for rows in \code{x}. Required for multi-group designs; ignored if \code{y} is provided.
#' @param difference_selection Logical. If \code{TRUE}, the function models the inclusion or exclusion of parameter differences. Default: \code{TRUE}.
#' @param save_main,save_pairwise,save_indicator Logical. Enable saving sampled states for \code{main_effects}, \code{pairwise_effects}, and \code{indicator}, respectively. Default: \code{FALSE}.
#' @param main_difference_model Character. Specifies how to handle threshold differences when categories are unmatched. Options: \code{"Collapse"}, \code{"Free"}. The option "Collapse" collapses categories unobserved in one or more groups. The option "Free" option estimates thresholds separately for each group and does not model their difference. Default: \code{"Free"}.
#' @param variable_type Character or vector. Specifies the type of variables in \code{x} (\code{"ordinal"} or \code{"blume-capel"}). Default: \code{"ordinal"}.
#' @param reference_category Integer or vector. Reference category for Blume-Capel variables. Required if there is at least one Blume-Capel variable.
#' @param pairwise_difference_scale Double. Scale parameter for the Cauchy prior on pairwise differences. Default: \code{1}.
#' @param main_difference_scale Double. Scale parameter for the Cauchy prior on threshold differences. Default: \code{1}.
#' @param pairwise_difference_prior,main_difference_prior Character. Specifies the inclusion probability model (\code{"Bernoulli"} or \code{"Beta-Bernoulli"}). Default: \code{"Bernoulli"}.
#' @param pairwise_difference_probability A numeric value or a \eqn{p \times p} matrix specifying the prior inclusion probability of a pairwise difference in the Bernoulli model. A single value applies the same probability to all pairs, while a matrix allows for edge-specific probabilities. Default: 0.5 for equal prior probability for inclusion and exclusion.
#' @param main_difference_probability A numeric value or a length-\eqn{p} vector specifying the prior inclusion probability of a threshold difference in the Bernoulli model. A single value applies the same probability to all variables, while a vector allows for variable-specific probabilities. Default: 0.5 to indicate no prior preference.
#' @param iter,burnin Integer. Number of Gibbs iterations (\code{iter}) and burn-in iterations (\code{burnin}). Defaults: \code{iter = 1e4}, \code{burnin = 1e3}.
#' @param na_action Character. Specifies handling of missing data. \code{"listwise"} deletes rows with missing values; \code{"impute"} imputes values during Gibbs sampling. Default: \code{"listwise"}.
#' @param display_progress Logical. Show progress bar during computation. Default: \code{TRUE}.
#' @param threshold_alpha,threshold_beta Double. Shape parameters for the beta-prime prior on nuisance threshold parameters.
#' @param interaction_scale Double. Scale of the Cauchy prior for nuisance pairwise interactions. Default: \code{2.5}.
#' @param main_beta_bernoulli_alpha,main_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on threshold differences.
#' @param pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on pairwise differences.
#' @param save Logical. If true, sampled states for all parameters are returned. Deprecated.
#'
#' @return A list containing the posterior means and, optionally, sampled states based on the \code{save_*} options. The returned components include:
#' \itemize{
#'  \item \code{posterior_mean_main}, \code{posterior_mean_pairwise}, and \code{posterior_mean_indicator} for posterior means.
#'  \item If saving options are enabled, the list also includes:
#'    \itemize{
#'      \item \code{raw_samples_main} – sampled states of main effects.
#'      \item \code{raw_samples_pairwise} – sampled states of pairwise effects.
#'      \item \code{raw_samples_indicator} – sampled states of inclusion indicators.
#'    }
#' }
#' In addition to the results of the analysis, the output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @examples
#' \donttest{
#' # Store user par() settings
#' op <- par(no.readonly = TRUE)
#'
#' # Run bgmCompare on the Boredom dataset
#' # For publication-quality results, consider using at least 1e5 iterations
#' fit <- bgmCompare(x = Boredom[, -1], g = Boredom[, 1], iter = 1e4)
#'
#' #--- INCLUSION VS EDGE DIFFERENCE PLOT -------------------------------------
#' incl_probs <- fit$posterior_mean_indicator[lower.tri(fit$posterior_mean_indicator)]
#' edge_diffs <- 2 * fit$posterior_mean_pairwise[, 3]
#'
#' par(mar = c(5, 5, 1, 1) + 0.1, cex = 1.7)
#' plot(edge_diffs, incl_probs,
#'      pch = 21, bg = "gray", cex = 1.3,
#'      ylim = c(0, 1), axes = FALSE,
#'      xlab = "", ylab = "")
#' abline(h = c(0, 0.5, 1), lty = 2, col = "gray")
#' axis(1); axis(2, las = 1)
#' mtext("Posterior Mean Edge Difference", side = 1, line = 3, cex = 1.7)
#' mtext("Posterior Inclusion Probability", side = 2, line = 3, cex = 1.7)
#'
#' #--- EVIDENCE PLOT ----------------------------------------------------------
#' prior_odds <- 1
#' post_odds <- incl_probs / (1 - incl_probs)
#' log_bf <- log(post_odds / prior_odds)
#' log_bf <- pmin(log_bf, 5)  # cap extreme values
#'
#' plot(edge_diffs, log_bf,
#'      pch = 21, bg = "#bfbfbf", cex = 1.3,
#'      axes = FALSE, xlab = "", ylab = "",
#'      ylim = c(-5, 5.5), xlim = c(-0.3, 0.6))
#' axis(1); axis(2, las = 1)
#' abline(h = log(c(1/10, 10)), lwd = 2, col = "#bfbfbf")
#' text(0.4, log(1/10), "Evidence for Exclusion", pos = 1, cex = 1.1)
#' text(0.4, log(10),   "Evidence for Inclusion", pos = 3, cex = 1.1)
#' text(0.4, 0,         "Absence of Evidence", cex = 1.1)
#' mtext("Log-Inclusion Bayes Factor", side = 2, line = 3, cex = 1.7)
#' mtext("Posterior Mean Interaction Difference", side = 1, line = 3, cex = 1.7)
#'
#' #--- MEDIAN PROBABILITY DIFFERENCE NETWORK ----------------------------------
#' median_edges <- ifelse(incl_probs >= 0.5, edge_diffs, 0)
#' p <- ncol(Boredom) - 1
#' median_net <- matrix(0, nrow = p, ncol = p)
#' median_net[lower.tri(median_net)] <- median_edges
#' median_net <- median_net + t(median_net)
#'
#' node_labels <- colnames(Boredom[, -1])
#' dimnames(median_net) <- list(node_labels, node_labels)
#'
#' par(cex = 1)
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   qgraph::qgraph(median_net,
#'     theme = "TeamFortress", maximum = 0.5, fade = FALSE,
#'     color = "#f0ae0e", vsize = 10, repulsion = 0.9,
#'     label.cex = 0.9, label.scale = FALSE,
#'     labels = node_labels
#'   )
#' }
#'
#' # Restore user par() settings
#' par(op)
#' }
#'
#' @importFrom methods hasArg
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#'
#' @export
bgmCompare = function(x,
                      y,
                      g,
                      difference_selection = TRUE,
                      main_difference_model = c("Free", "Collapse", "Constrain"),
                      variable_type = "ordinal",
                      reference_category,
                      pairwise_difference_scale = 1,
                      main_difference_scale = 1,
                      pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                      main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                      pairwise_difference_probability = 0.5,
                      main_difference_probability = 0.5,
                      pairwise_beta_bernoulli_alpha = 1,
                      pairwise_beta_bernoulli_beta = 1,
                      main_beta_bernoulli_alpha = 1,
                      main_beta_bernoulli_beta = 1,
                      interaction_scale = 2.5,
                      threshold_alpha = 0.5,
                      threshold_beta = 0.5,
                      iter = 1e4,
                      burnin = 1e3,
                      na_action = c("listwise", "impute"),
                      save = FALSE,
                      save_main = FALSE,
                      save_pairwise = FALSE,
                      save_indicator = FALSE,
                      display_progress = TRUE) {

  # Deprecation warning for save parameter
  if(hasArg(save)) {
    warning("`save` is deprecated. Use `save_main`, `save_pairwise`, or `save_indicator` instead.")
    save_main = save_main || save
    save_pairwise = save_pairwise || save
    save_indicator = save_indicator || save
  }

  # Check and preprocess data
  x = data_check(x, "x")
  if (hasArg(y)) {
    y = data_check(y, "y")
    if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }

  if(!hasArg(y) & !hasArg(g))
    stop(paste0("For multi-group designs, the bgmCompare function requires input for\n",
                "either y (group 2 data) or g (group indicator)."))

  # Validate group indicators
  if (!hasArg(y) && hasArg(g)) {
    g = as.vector(g)
    if (anyNA(g)) stop("g cannot contain missing values.")
    if (length(g) != nrow(x)) stop("Length of g must match number of rows in x.")
  }

  # Model and preprocessing
  if(!hasArg(y))
    y = NULL
  if(!hasArg(g))
    g = NULL

  model = check_compare_model(
    x = x, y = y, g = g, difference_selection = difference_selection,
    variable_type = variable_type, reference_category = reference_category,
    pairwise_difference_scale = pairwise_difference_scale,
    main_difference_scale = main_difference_scale,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    pairwise_difference_probability = pairwise_difference_probability,
    main_difference_probability = main_difference_probability,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    interaction_scale = interaction_scale, threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    main_difference_model = main_difference_model
    )

  x = model$x
  group = model$group
  ordinal_variable = model$variable_bool
  reference_category = model$reference_category
  independent_thresholds = (model$main_difference_model == "Free")

  # Check Gibbs input
  check_positive_integer(iter, "iter")
  check_non_negative_integer(burnin, "burnin")

  # Check na_action
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if (inherits(na_action, "try-error")) {
    stop(sprintf("Invalid value for `na_action`. Expected 'listwise' or 'impute', got: %s", na_action_input))
  }

  # Check save options
  save_main = check_logical(save_main, "save_main")
  save_pairwise = check_logical(save_pairwise, "save_pairwise")
  save_indicator = check_logical(save_indicator, "save_indicator")

  # Check display_progress
  display_progress = check_logical(display_progress, "display_progress")

  ## Format data
  data = compare_reformat_data(
    x = x, group = group,
    na_action = na_action,
    variable_bool = ordinal_variable,
    reference_category = reference_category,
    main_difference_model = model$main_difference_model
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
  num_obs_categories = compute_num_obs_categories(x, num_categories, group)


  # Compute sufficient statistics for Blume-Capel variables
  sufficient_blume_capel = compute_sufficient_blume_capel(x, reference_category,ordinal_variable, group)


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

  # Call the Rcpp function
  out = compare_anova_gibbs_sampler(
    observations = observations, main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices, projection = projection,
    num_categories = num_categories, num_groups = num_groups,
    group_indices = group_indices, interaction_scale = interaction_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    main_difference_scale = main_difference_scale,
    pairwise_difference_prior = model$pairwise_difference_prior,
    main_difference_prior = model$main_difference_prior,
    inclusion_probability_difference = model$inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta, Index = Index,
    iter = iter, burnin = burnin, num_obs_categories = num_obs_categories,
    sufficient_blume_capel = sufficient_blume_capel,
    prior_threshold_alpha = threshold_alpha,
    prior_threshold_beta = threshold_beta,
    na_impute = na_impute, missing_data_indices = missing_index,
    is_ordinal_variable = ordinal_variable,
    baseline_category = reference_category,
    independent_thresholds = independent_thresholds, save_main = save_main,
    save_pairwise = save_pairwise, save_indicator = save_indicator,
    display_progress = display_progress,
    difference_selection = difference_selection)

  # Main output handler in the wrapper function
  output = prepare_output_bgmCompare(
    out = out, x = x,
    independent_thresholds = independent_thresholds, num_variables = num_variables,
    num_categories = num_categories,
    group = group, iter = iter,
    data_columnnames = if (is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
    save_options = list(save_main = save_main, save_pairwise = save_pairwise,
                        save_indicator = save_indicator),
    difference_selection = difference_selection, na_action = na_action,
    na_impute = na_impute, variable_type = variable_type, burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    main_difference_model = model$main_difference_model,
    pairwise_difference_prior = model$pairwise_difference_prior,
    main_difference_prior = model$main_difference_prior,
    inclusion_probability_difference = model$inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    main_difference_scale = main_difference_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    projection,
    is_ordinal_variable = ordinal_variable
  )

  return(output)
}
