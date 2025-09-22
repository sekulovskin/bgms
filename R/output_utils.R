prepare_output_bgm = function(
    out, x, num_categories, iter, data_columnnames, is_ordinal_variable,
    burnin, interaction_scale, threshold_alpha, threshold_beta,
    na_action, na_impute, edge_selection, edge_prior, inclusion_probability,
    beta_bernoulli_alpha, beta_bernoulli_beta, dirichlet_alpha, lambda,
    variable_type, update_method, target_accept, hmc_num_leapfrogs,
    nuts_max_depth, learn_mass_matrix, num_chains
) {
  arguments = list(
    prepared_data = x,
    num_variables = ncol(x),
    num_cases = nrow(x),
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda,
    na_action = na_action,
    version = packageVersion("bgms"),
    update_method = update_method,
    target_accept = target_accept,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    num_chains = num_chains
  )

  num_variables = ncol(x)
  results = list()

  # ======= Parameter name generation =======
  names_variable_categories = character()
  for (v in seq_len(num_variables)) {
    if (is_ordinal_variable[v]) {
      cats = seq_len(num_categories[v])
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], "(", cats, ")")
      )
    } else {
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], "(linear)"),
        paste0(data_columnnames[v], "(quadratic)")
      )
    }
  }

  edge_names = character()
  for (i in 1:(num_variables - 1)) {
    for (j in (i + 1):num_variables) {
      edge_names = c(edge_names, paste0(data_columnnames[i], "-", data_columnnames[j]))
    }
  }

  # ======= Summarize MCMC chains =======
  summary_list = summarize_fit(out, edge_selection = edge_selection)
  main_summary = summary_list$main[, -1]
  pairwise_summary = summary_list$pairwise[, -1]

  rownames(main_summary) = names_variable_categories
  rownames(pairwise_summary) = edge_names

  results$posterior_summary_main = main_summary
  results$posterior_summary_pairwise = pairwise_summary

  if (edge_selection) {
    indicator_summary = summarize_indicator(out, param_names = edge_names)[, -1]
    rownames(indicator_summary) = edge_names
    results$posterior_summary_indicator = indicator_summary
    if (identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) {
      # convergence diagnostics of the co-apperance of the nodes
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(out, `[[`, "allocations"),
        node_names = data_columnnames
      )
      # posterior summary of each pairwise cluster coclustering
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
    }
  }

  # ======= Posterior mean matrices (for legacy compatibility) =======
  results$posterior_mean_main = matrix(main_summary$mean, nrow = num_variables, byrow = TRUE)
  rownames(results$posterior_mean_main) = data_columnnames
  colnames(results$posterior_mean_main) = NULL

  results$posterior_mean_pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
  results$posterior_mean_pairwise[upper.tri(results$posterior_mean_pairwise)] = pairwise_summary$mean
  results$posterior_mean_pairwise[lower.tri(results$posterior_mean_pairwise)] =
    t(results$posterior_mean_pairwise)[lower.tri(results$posterior_mean_pairwise)]
  rownames(results$posterior_mean_pairwise) = data_columnnames
  colnames(results$posterior_mean_pairwise) = data_columnnames

  if (edge_selection) {
    indicator_means = indicator_summary$mean
    results$posterior_mean_indicator = matrix(0, nrow = num_variables, ncol = num_variables)
    results$posterior_mean_indicator[upper.tri(results$posterior_mean_indicator)] = indicator_means
    results$posterior_mean_indicator[lower.tri(results$posterior_mean_indicator)] =
      t(results$posterior_mean_indicator)[lower.tri(results$posterior_mean_indicator)]
    rownames(results$posterior_mean_indicator) = data_columnnames
    colnames(results$posterior_mean_indicator) = data_columnnames

    if (identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) {
      # convergence diagnostics of the co-apperance of the nodes
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(out, `[[`, "allocations"),
        node_names = data_columnnames
      )
      results$posterior_coclustering_matrix = sbm_convergence$co_occur_matrix
      # calculate the estimated clustering and block probabilities
      sbm_summary = posterior_summary_SBM(allocations = lapply(out, `[[`, "allocations"),
                                          arguments = arguments)  # check if only arguments would work
      # extract the posterior mean and median
      results$posterior_mean_allocations = sbm_summary$allocations_mean
      results$posterior_mode_allocations = sbm_summary$allocations_mode

      # extract the number of blocks and their estimated posterior probabilties
      results$posterior_num_blocks = sbm_summary$blocks
    }
  }


  results$arguments = arguments
  class(results) = "bgms"

  results$raw_samples = list(
    main      = lapply(out, function(chain) chain$main_samples),
    pairwise  = lapply(out, function(chain) chain$pairwise_samples),
    indicator = if (edge_selection) lapply(out, function(chain) chain$indicator_samples) else NULL,
    allocations = if (edge_selection && identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) lapply(out, `[[`, "allocations") else NULL,
    nchains   = length(out),
    niter     = nrow(out[[1]]$main_samples),
    parameter_names = list(
      main     = names_variable_categories,
      pairwise = edge_names,
      indicator = if (edge_selection) edge_names else NULL,
      allocations = if (identical(edge_prior, "Stochastic-Block")) edge_names else NULL
    )
  )

  return(results)
}



prepare_output_bgmCompare = function(out, x, independent_thresholds,
                                     num_variables, num_categories, group, iter,
                                     data_columnnames, save_options,
                                     difference_selection, na_action,
                                     na_impute, variable_type, burnin,
                                     interaction_scale, threshold_alpha,
                                     threshold_beta, main_difference_model,
                                     pairwise_difference_prior,
                                     main_difference_prior,
                                     inclusion_probability_difference,
                                     pairwise_beta_bernoulli_alpha,
                                     pairwise_beta_bernoulli_beta,
                                     main_beta_bernoulli_alpha,
                                     main_beta_bernoulli_beta,
                                     main_difference_scale,
                                     pairwise_difference_scale,
                                     projection, is_ordinal_variable) {

  save = any(c(save_options$save_main, save_options$save_pairwise, save_options$save_indicator))

  arguments = list(
    num_variables = num_variables, group = group, num_cases = tabulate(group),
    na_action = na_action, na_impute = na_impute, iter = iter, burnin = burnin,
    difference_selection = difference_selection,
    independent_thresholds = independent_thresholds,
    save_main = save_options$save_main,
    save_pairwise = save_options$save_pairwise,
    save_indicator = save_options$save_indicator, save = save,
    variable_type = variable_type, interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha, threshold_beta = threshold_beta,
    main_difference_model = main_difference_model,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    inclusion_probability_difference = inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    main_difference_scale = main_difference_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    version = packageVersion("bgms")
  )

  # Names for pairwise effects
  names_bycol = matrix(rep(data_columnnames, each = num_variables), ncol = num_variables)
  names_byrow = matrix(rep(data_columnnames, each = num_variables), ncol = num_variables, byrow = TRUE)
  names_comb = matrix(paste0(names_byrow, "-", names_bycol), ncol = num_variables)
  names_vec = names_comb[lower.tri(names_comb)]
  names_vec_t = names_comb[lower.tri(names_comb, diag = TRUE)]

  # Prepare output elements
  results = list()

  # Handle output from the new rcpp function (anova)
  num_groups = length(unique(group))
  num_main = nrow(out$posterior_mean_main)

  # Main effects
  tmp = out$posterior_mean_main
  if(independent_thresholds) {
    posterior_mean_main = matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp))
    for(g in 1:num_groups) {
      posterior_mean_main[, g] = tmp[, g]

      #This can probably be done prettier
      if(any(tmp[,g] == 0))
        posterior_mean_main[tmp[,g] == 0, g] = NA
    }
  } else {
    posterior_mean_main = matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp) + 1)
    posterior_mean_main[, 1] = tmp[, 1]
    for (row in 1:nrow(tmp)) {
      posterior_mean_main[row, -1] = projection %*% tmp[row, -1]
    }
  }
  results$posterior_mean_main = posterior_mean_main

  names_variable_categories = vector(length = nrow(tmp))
  counter = 0
  for (v in 1:num_variables) {
    if(is_ordinal_variable[v]) {
      for (c in 1:max(num_categories[v, ])) {
        counter = counter + 1
        names_variable_categories[counter] = paste0(data_columnnames[v], "(", c, ")")
      }
    } else {
      counter = counter + 1
      names_variable_categories[counter] = paste0(data_columnnames[v], "(linear)")
      counter = counter + 1
      names_variable_categories[counter] = paste0(data_columnnames[v], "(quadratic)")
    }

  }
  if(independent_thresholds) {
    dimnames(results$posterior_mean_main) = list(names_variable_categories, paste0("group_", 1:num_groups))
  } else {
    dimnames(results$posterior_mean_main) = list(names_variable_categories, c("overall", paste0("group_", 1:num_groups)))
  }

  # Pairwise effects
  tmp = out$posterior_mean_pairwise
  posterior_mean_pairwise = matrix(0, nrow = nrow(tmp), ncol = ncol(tmp) + 1)
  posterior_mean_pairwise[, 1] = tmp[, 1]
  for (row in 1:nrow(tmp)) {
    posterior_mean_pairwise[row, -1] = projection %*% tmp[row, -1]
  }

  results$posterior_mean_pairwise = posterior_mean_pairwise
  dimnames(results$posterior_mean_pairwise) = list(names_vec, c("overall", paste0("group_", 1:num_groups)))

  if (difference_selection && "posterior_mean_indicator" %in% names(out) ) {
    results$posterior_mean_indicator = out$posterior_mean_indicator
    dimnames(results$posterior_mean_indicator) = list(data_columnnames,
                                                      data_columnnames)
    if(main_difference_model == "Free"){
      diag(results$posterior_mean_indicator) = NA
    }
  }


  # Handle main_effect_samples
  if (save_options$save_main && "main_effect_samples" %in% names(out)) {
    main_effect_samples = out$main_effect_samples
    col_names = character()  # Vector to store column names

    if(independent_thresholds) {
      for (gr in 1:num_groups) {
        for (var in 1:num_variables) {
          if(is_ordinal_variable[var]) {
            for (cat in 1:max(num_categories[var, ])) {
              col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(", cat, ")"))
            }
          } else {
            col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(linear)"))
            col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(quadratic)"))
          }
        }
      }

    } else {

      for (gr in 1:num_groups) {
        if(gr == 1) {
          for (var in 1:num_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(num_categories[var, ])) {
                col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(", cat, ")"))
              }
            } else {
              col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(linear)"))
              col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(quadratic)"))
            }
          }
        } else {
          for (var in 1:num_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(num_categories[var, ])) {
                col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(", cat, ")"))
              }
            } else {
              col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(linear)"))
              col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(quadratic)"))
            }
          }
        }
      }
    }

    dimnames(main_effect_samples) = list(Iter. = 1:nrow(main_effect_samples), col_names)
    results$main_effect_samples = main_effect_samples
  }

  # Handle pairwise_effect_samples
  if (save_options$save_pairwise && "pairwise_effect_samples" %in% names(out)) {
    col_names = character()  # Vector to store column names
    for (v1 in 1:(num_variables - 1)) {
      for (v2 in (v1 + 1):num_variables) {
        col_names = c(col_names, paste0(data_columnnames[v1], "-", data_columnnames[v2]))
        for (gr in 2:num_groups) {
          col_names = c(col_names,
                        paste0("contrast_#",
                               gr-1,
                               "(",
                               data_columnnames[v1],
                               "-",
                               data_columnnames[v2],
                               ")"))
        }
      }
    }

    dimnames(out$pairwise_effect_samples) = list(Iter. = 1:nrow(out$pairwise_effect_samples), col_names)
    results$pairwise_effect_samples = out$pairwise_effect_samples
  }

  # Handle inclusion_indicator_samples
  if (difference_selection && save_options$save_indicator && "inclusion_indicator_samples" %in% names(out)) {
    if(independent_thresholds) {
      #
    } else {
      col_names = character()  # Vector to store column names
      for (v1 in 1:num_variables) {
        for (v2 in v1:num_variables) {
          col_names = c(col_names, paste0(data_columnnames[v1], "-", data_columnnames[v2]))
        }
      }
      dimnames(out$inclusion_indicator_samples) = list(Iter. = 1:nrow(out$inclusion_indicator_samples), col_names)
    }

    results$inclusion_indicator_samples = out$inclusion_indicator_samples
  }

  # Add arguments to the output
  results$arguments = arguments
  class(results) = c("bgmCompare")
  return(results)
}