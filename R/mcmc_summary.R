# Summary utilities for spike-and-slab MCMC output

# Combine MCMC chains into a 3D array [niter x nchains x nparam]
combine_chains = function(fit, component) {
  nchains = length(fit)
  samples_list = lapply(fit, function(x) x[[component]])
  niter = nrow(samples_list[[1]])
  nparam = ncol(samples_list[[1]])
  array3d = array(NA_real_, dim = c(niter, nchains, nparam))
  for (i in seq_len(nchains)) {
    array3d[, i, ] = samples_list[[i]]
  }
  array3d
}

# Compute effective sample size and Rhat (Gelman-Rubin diagnostic)
compute_rhat_ess = function(draws) {
  tryCatch({
    if (is.matrix(draws) && ncol(draws) > 1) {
      mcmc_list = coda::mcmc.list(
        lapply(seq_len(ncol(draws)), function(i) coda::mcmc(draws[, i]))
      )
      ess = coda::effectiveSize(mcmc_list)
      rhat = coda::gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1]
    } else {
      ess = coda::effectiveSize(draws)
      rhat = NA_real_
    }
    list(ess = ess, rhat = rhat)
  }, error = function(e) list(ess = NA_real_, rhat = NA_real_))
}

# Basic summarizer for continuous parameters
summarize_manual = function(fit, component = c("main_samples", "pairwise_samples"), param_names = NULL) {
  component = match.arg(component) # Add options later
  array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]

  result = matrix(NA, nparam, 5)
  colnames(result) = c("mean", "mcse", "sd", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec = as.vector(draws)
    result[j, "mean"] = mean(vec)
    result[j, "sd"] = sd(vec)
    est = compute_rhat_ess(draws)
    result[j, "mcse"] = sd(vec) / sqrt(est$ess)
    result[j, "n_eff"] = est$ess
    result[j, "Rhat"] = est$rhat
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("parameter [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }

}

# Summarize binary indicator variables
summarize_indicator = function(fit, component = c("indicator_samples"), param_names = NULL) {
  component = match.arg(component) # Add options later
  array3d = combine_chains(fit, component)

  nparam = dim(array3d)[3]
  nchains = dim(array3d)[2]
  niter = dim(array3d)[1]

  result = matrix(NA, nparam, 9)
  colnames(result) = c("mean", "sd", "mcse", "n0->0", "n0->1", "n1->0", "n1->1", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec = as.vector(draws)
    T = length(vec)
    g_next = vec[-1]
    g_curr = vec[-T]

    p_hat = mean(vec)
    sd = sqrt(p_hat * (1 - p_hat))
    n00 = sum(g_curr == 0 & g_next == 0)
    n01 = sum(g_curr == 0 & g_next == 1)
    n10 = sum(g_curr == 1 & g_next == 0)
    n11 = sum(g_curr == 1 & g_next == 1)

    if (n01 + n10 == 0) {
      n_eff = mcse = R = NA_real_
    } else {
      a = n01 / (n00 + n01)
      b = n10 / (n10 + n11)
      tau_int = (2 - (a + b)) / (a + b)
      n_eff = T / tau_int
      mcse = sd / sqrt(n_eff)
      est = compute_rhat_ess(draws)
      R = est$rhat
    }

    result[j, ] = c(p_hat, sd, mcse, n00, n01, n10, n11, n_eff, R)
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("indicator [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = paste0(param_names, "- indicator"),
               result, check.names = FALSE)
  }
}

# Summarize slab values where indicators are 1
summarize_slab = function(fit, component = c("pairwise_samples"), param_names = NULL) {
  component = match.arg(component) # Add options later
  array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]
  result = matrix(NA, nparam, 5)
  colnames(result) = c("mean", "sd", "mcse", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec = as.vector(draws)
    nonzero = vec != 0
    vec = vec[nonzero]
    T = length(vec)

    if (T > 10) {
      eap = mean(vec)
      sdev = sd(vec)
      est = compute_rhat_ess(draws)
      mcse = sdev / sqrt(est$ess)
      result[j, ] = c(eap, sdev, mcse, est$ess, est$rhat)
    }
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("weight [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = paste0(param_names, "- weight"),
               result, check.names = FALSE)
  }
}

# Combined summary for pairwise parameters with selection
summarize_pair = function(fit,
                           indicator_component = c("indicator_samples"),
                           slab_component = c("pairwise_samples"),
                           param_names = NULL
) {
  indicator_component = match.arg(indicator_component) # Add options later
  slab_component = match.arg(slab_component) # Add options later

  summ_ind = summarize_indicator(fit, component = indicator_component)
  summ_slab = summarize_slab(fit, component = slab_component)
  nparam = nrow(summ_ind)

  eap = summ_ind$mean * summ_slab$mean
  v = (summ_slab$mean^2 * summ_ind$sd^2) + (summ_ind$mean^2 * summ_slab$sd^2)
  mcse2 = (summ_slab$mean^2 * summ_ind$mcse^2) + (summ_ind$mean^2 * summ_slab$mcse^2)
  mcse = sqrt(mcse2)
  sd = sqrt(v)
  n_eff = v / mcse2

  rhat = rep(NA_real_, nparam)
  array3d_pw = combine_chains(fit, slab_component)
  array3d_id = combine_chains(fit, indicator_component)
  nchains = dim(array3d_pw)[2]
  T = prod(dim(array3d_pw)[1:2])

  for (j in seq_len(nparam)) {
    draws_pw = array3d_pw[, , j]
    draws_id = array3d_id[, , j]
    if (nchains > 1) {
      chain_means = numeric(nchains)
      chain_vars = numeric(nchains)
      for (chain in 1:nchains) {
        pi = mean(draws_id[, chain])
        tmp = draws_pw[, chain]
        e = mean(tmp[tmp != 0])
        v = var(tmp[tmp != 0])
        chain_means[chain] = pi * e
        chain_vars[chain] = pi * (v + (1 - pi) * e^2)
      }
      B = T * sum((chain_means - eap[j])^2) / (nchains - 1)
      W = mean(chain_vars)
      V = (T - 1) * W / T + B / T
      rhat[j] = sqrt(V / W)
    }
  }

  data.frame(parameter = paste0("V[", seq_len(nparam), "]"),
             mean = eap, sd = sd, mcse = mcse, n_eff = n_eff, Rhat = rhat,
             check.names = FALSE)

  if(is.null(param_names)) {
    data.frame(parameter = paste0("weight [", seq_len(nparam), "]"),
               mean = eap, sd = sd, mcse = mcse, n_eff = n_eff, Rhat = rhat,
               check.names = FALSE)
  } else {
    data.frame(parameter = paste0(param_names, "- weight"),
               mean = eap, sd = sd, mcse = mcse, n_eff = n_eff, Rhat = rhat,
               check.names = FALSE)
  }
}

# Unified summary dispatcher for either model type
summarize_fit = function(fit, edge_selection = FALSE) {
  main_summary = summarize_manual(fit, component = "main_samples")

  if (!edge_selection) {
    pair_summary = summarize_manual(fit, component = "pairwise_samples")
  } else {
    # Get indicators and slab draws
    ind_summary = summarize_indicator(fit, component = "indicator_samples")
    slab_summary = summarize_slab(fit, component = "pairwise_samples")

    all_selected = ind_summary$mean == 1

    # Use summarize_pair only where not always selected
    full_summary = summarize_pair(fit,
                                   indicator_component = "indicator_samples",
                                   slab_component = "pairwise_samples")
    manual_summary = summarize_manual(fit, component = "pairwise_samples")

    # Replace rows in full_summary with manual results for fully selected entries
    if (any(all_selected)) {
      full_summary[all_selected, ] = manual_summary[all_selected, ]
    }

    pair_summary = full_summary
  }

  list(main = main_summary, pairwise = pair_summary)
}


# summarize the SBM output -----------------------------------------------------

# Calculate convergence diagnostics on the pairwise cluster co-appearance values
summarize_alloc_pairs = function(allocations, node_names = NULL) {
  #stopifnot(is.list(allocations), length(allocations) >= 2)
  n_ch   = length(allocations)
  n_iter = nrow(allocations[[1]])
  no_variables    = ncol(allocations[[1]])
  for (c in seq_len(n_ch)) {
    stopifnot(nrow(allocations[[c]]) == n_iter, ncol(allocations[[c]]) == no_variables)
  }
  if (!is.null(node_names)) stopifnot(length(node_names) == no_variables)

  # all node pairs
  Pairs  = t(combn(seq_len(no_variables), 2))
  nparam = nrow(Pairs)

  result = matrix(NA, nparam, 5)
  colnames(result) = c("mean", "sd", "mcse", "n_eff", "Rhat")

  # helper to construct a "time-series"
  get_draws_pair = function(i, j) {
    out = matrix(NA, n_iter, n_ch)
    for (c in seq_len(n_ch)) {
      Zc = allocations[[c]]
      out[, c] = as.integer(Zc[, i] == Zc[, j])
    }
    out
  }

  for (p in seq_len(nparam)) {
    i = Pairs[p, 1]; j = Pairs[p, 2]
    draws = get_draws_pair(i, j)
    vec   = as.vector(draws)
    phat  = mean(vec)
    sdev  = sd(vec)

    est   = compute_rhat_ess(draws)
    n_eff = as.numeric(est$ess)
    Rhat  = as.numeric(est$rhat)

    mcse  = if (is.finite(n_eff) && n_eff > 0) sdev / sqrt(n_eff) else NA
    result[p, ] = c(phat, sdev, mcse, n_eff, Rhat)
  }
  if (is.null(node_names)) {
    rn = paste0(Pairs[,1], "-", Pairs[,2])
    dimn = as.character(seq_len(no_variables))
  } else {
    rn = paste0(node_names[Pairs[,1]], "-", node_names[Pairs[,2]])
    dimn = node_names
  }

  sbm_summary = as.data.frame(result, check.names = FALSE)
  rownames(sbm_summary) = rn

  # construct the co-appearance matrix
  co_occur_matrix = matrix(0, nrow = no_variables, ncol = no_variables,
                            dimnames = list(dimn, dimn))
  diag(co_occur_matrix) = 1
  for (p in seq_len(nparam)) {
    i = Pairs[p, 1]; j = Pairs[p, 2]
    m = sbm_summary[p, "mean"]
    co_occur_matrix[i, j] = m
    co_occur_matrix[j, i] = m
  }
  list(sbm_summary = sbm_summary, co_occur_matrix = co_occur_matrix)
}

# calculate a representative allocation vector using
# (1) the median of the posterior distribution of the cluster allocations
# (2) the mean which is based on Dahl's method: This part of the code
# was adapted from the R
# code accompanying the paper:
#  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
#  Detection With Unknown Number of Communities, Journal of the American
#  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
find_representative_clustering = function(cluster_matrix) {
  stopifnot(is.matrix(cluster_matrix))
  n_iter = nrow(cluster_matrix)
  p      = ncol(cluster_matrix)

  # Build co-clustering (membership) matrices for all iterations

  Ms = lapply(seq_len(n_iter), function(t) {
    z = cluster_matrix[t, ]
    (outer(z, z, FUN = "==")) * 1L
  })

  # Average (posterior similarity / co-clustering) matrix
  psm = Reduce(`+`, Ms) / n_iter

  # MEAN representative (Dahl's method)
  sqerr = vapply(Ms, function(M) sum((M - psm)^2), numeric(1))
  idx_dahl = which.min(sqerr)
  alloc_dahl = cluster_matrix[idx_dahl, , drop = TRUE]

  #  MODE representative
  hash_mat = function(M) paste(as.integer(t(M)), collapse = ",")
  keys = vapply(Ms, hash_mat, character(1))
  tab  = table(keys)
  key_mode = names(tab)[which.max(tab)]
  idx_mode = match(key_mode, keys)
  alloc_mode = cluster_matrix[idx_mode, , drop = TRUE]
  indicator_mode = matrix(
    as.integer(strsplit(key_mode, ",", fixed = TRUE)[[1]]),
    nrow = p, byrow = TRUE
  )
  p_dist = as.numeric(tab) / sum(tab)
  posterior_variance = (1 - sum(p_dist^2)) / (1 - 1 / length(p_dist))

  list(
    mean = alloc_dahl,
    mode = alloc_mode
  )
}

# Calculate the conditional probability of the number of blocks given the
# cardinality of a sampled allocation vector based on Equation (3.7) from
# Miller & Harrison (2018). Mixture Models With a Prior on the Number of
# blocks, Journal of the American Statistical Association, 113:521, 340-356,
# DOI:10.1080/01621459.2016.1255636
#' @importFrom stats dpois
compute_p_k_given_t = function(
    t,
    log_Vn,
    dirichlet_alpha,
    num_variables,
    lambda) {
  # Define the K_values
  K_values = as.numeric(1:num_variables)

  # Initialize vector for probabilities
  p_k_given_t = numeric(length(K_values))

  # Normalization constant for t
  log_vn_t = log_Vn[t]

  # Normalizing factor for the truncated Poisson distribution
  norm_factor = 1 - dpois(0, lambda)
  truncated_poisson_pmf = dpois(K_values, lambda) / norm_factor

  # Loop through each value of K
  for (i in seq_along(K_values)) {
    K = K_values[i]
    if (K >= t) {
      # Falling factorial
      falling_factorial = prod(K:(K - t + 1))
      # Rising factorial
      rising_factorial = prod((dirichlet_alpha * K) + 0:(num_variables - 1))
      # Compute log probability
      log_p_k = log(falling_factorial) - log(rising_factorial) +
        log(truncated_poisson_pmf[i]) - log_vn_t
      # Convert log probability to probability
      p_k_given_t[i] = exp(log_p_k)
    } else {
      p_k_given_t[i] = 0
    }
  }
  # Normalize probabilities
  p_k_given_t = p_k_given_t / sum(p_k_given_t)

  return(p_k_given_t)
}

# Wrapper function to compute the posterior summary for the Stochastic Block Model
posterior_summary_SBM = function(
    allocations,
    arguments) {

  # combine the allocations from the chains
  cluster_allocations = do.call(rbind, allocations)

  dirichlet_alpha = arguments$dirichlet_alpha
  lambda = arguments$lambda
  num_variables = ncol(cluster_allocations)

  # Pre-compute log_Vn for computing the cluster probabilities
  log_Vn = compute_Vn_mfm_sbm(
    num_variables, dirichlet_alpha, num_variables + 10, lambda)

  # Compute the number of unique clusters (t) for each iteration, i.e., the
  # cardinality  of the partition z
  clusters = apply(cluster_allocations, 1, function(row) length(unique(row)))

  # Compute the conditional probabilities of the number of clusters for each
  # row in clusters
  p_k_given_t = matrix(NA, nrow = length(clusters), ncol = num_variables)

  for (i in 1:length(clusters)) {
    p_k_given_t[i, ] = compute_p_k_given_t(
      clusters[i], log_Vn, dirichlet_alpha, num_variables, lambda)
  }

  # Average across all iterations
  p_k_given_t = colMeans(p_k_given_t)

  # Format the output
  #num_blocks = 1:num_variables
  blocks = cbind(p_k_given_t)
  colnames(blocks) = c("probability")

  # make blocks a data frame
  blocks = as.data.frame(blocks)

  # Compute the mean and mode of the allocations
  allocations = find_representative_clustering(cluster_allocations)

  return(list(blocks = blocks,
              allocations_mean = allocations$mean,
              allocations_mode = allocations$mode))
}


# Combine MCMC chains for bgmCompare into a 3D array [niter x nchains x nparam]
combine_chains_compare = function(fit, component) {
  nchains = length(fit)
  samples_list = lapply(fit, function(x) x[[component]])
  niter = nrow(samples_list[[1]])
  nparam = ncol(samples_list[[1]])
  array3d = array(NA_real_, dim = c(niter, nchains, nparam))
  for (i in seq_len(nchains)) {
    array3d[, i, ] = samples_list[[i]]
  }
  array3d
}


summarize_manual_compare = function(fit_or_array,
                                    component = c("main_samples", "pairwise_samples"),
                                    param_names = NULL) {
  component = match.arg(component)

  # allow either a fit list or a pre-combined 3D array
  if (is.array(fit_or_array)) {
    array3d = fit_or_array
  } else {
    array3d = combine_chains_compare(fit_or_array, component)
  }

  nparam = dim(array3d)[3]
  result = matrix(NA, nparam, 5)
  colnames(result) = c("mean", "mcse", "sd", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec   = as.vector(draws)
    result[j, "mean"] = mean(vec)
    result[j, "sd"]   = sd(vec)
    est = compute_rhat_ess(draws)
    result[j, "mcse"] = sd(vec) / sqrt(est$ess)
    result[j, "n_eff"] = est$ess
    result[j, "Rhat"]  = est$rhat
  }

  if (is.null(param_names)) {
    data.frame(parameter = paste0("param [", seq_len(nparam)), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }
}




summarize_indicator_compare = function(fit, component = "indicator_samples", param_names = NULL) {
  array3d = combine_chains_compare(fit, component)
  nparam = dim(array3d)[3]

  result = matrix(NA, nparam, 9)
  colnames(result) = c("mean", "sd", "mcse", "n0->0", "n0->1", "n1->0", "n1->1", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec = as.vector(draws)
    T = length(vec)
    g_next = vec[-1]
    g_curr = vec[-T]

    p_hat = mean(vec)
    sd = sqrt(p_hat * (1 - p_hat))
    n00 = sum(g_curr == 0 & g_next == 0)
    n01 = sum(g_curr == 0 & g_next == 1)
    n10 = sum(g_curr == 1 & g_next == 0)
    n11 = sum(g_curr == 1 & g_next == 1)

    if (n01 + n10 == 0) {
      n_eff = mcse = R = NA_real_
    } else {
      a = n01 / (n00 + n01)
      b = n10 / (n10 + n11)
      tau_int = (2 - (a + b)) / (a + b)
      n_eff = T / tau_int
      mcse = sd / sqrt(n_eff)
      est = compute_rhat_ess(draws)
      R = est$rhat
    }

    result[j, ] = c(p_hat, sd, mcse, n00, n01, n10, n11, n_eff, R)
  }

  if (is.null(param_names)) {
    data.frame(parameter = paste0("indicator [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }
}

summarize_main_diff_compare = function(
    fit,
    main_effect_indices,
    num_groups,
    param_names = NULL
) {
  main_effect_samples = combine_chains_compare(fit, "main_samples")
  indicator_samples   = combine_chains_compare(fit, "indicator_samples")

  nvars   = nrow(main_effect_indices)
  results = list()
  param_counter = 0

  for (v in seq_len(nvars)) {
    indicator_draws = indicator_samples[, , v]
    vec_id   = as.vector(indicator_draws)
    pi_hat   = mean(vec_id)

    start = main_effect_indices[v, 1] + 1L
    stop  = main_effect_indices[v, 2] + 1L

    for (row in start:stop) {
      category = row - start + 1
      for (h in 1:(num_groups - 1)) {
        param_counter = param_counter + 1
        col_index = h * nrow(main_effect_indices) + row
        effect_draws = main_effect_samples[, , col_index]
        vec_slab   = as.vector(effect_draws)

        eap_slab   = if (any(vec_id == 1)) mean(vec_slab[vec_id == 1]) else 0
        posterior_mean = pi_hat * eap_slab
        posterior_var  = (eap_slab^2 * var(vec_id)) + (pi_hat^2 * var(vec_slab))
        posterior_sd   = sqrt(posterior_var)

        est  = compute_rhat_ess(effect_draws)
        mcse_est = if (is.finite(est$ess) && est$ess > 0) posterior_sd / sqrt(est$ess) else NA

        results[[param_counter]] = data.frame(
          parameter = if (!is.null(param_names)) {
            param_names[param_counter]
          } else {
            paste0("var", v, "(diff", h, ", c", category, ")")
          },
          mean  = posterior_mean,
          sd    = posterior_sd,
          mcse  = mcse_est,
          n_eff = est$ess,
          Rhat  = est$rhat,
          check.names = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}



summarize_pairwise_diff_compare = function(
    fit,
    pairwise_effect_indices,
    num_variables,
    num_groups,
    param_names = NULL
) {
  pairwise_effect_samples = combine_chains_compare(fit, "pairwise_samples")
  indicator_samples = combine_chains_compare(fit, "indicator_samples")

  results = list()
  param_counter = 0
  ind_counter   = nrow(pairwise_effect_indices) # offset after diagonals

  for (v1 in 1:(num_variables - 1)) {
    for (v2 in (v1 + 1):num_variables) {
      ind_counter = ind_counter + 1
      indicator_draws = indicator_samples[, , ind_counter]
      vec_id   = as.vector(indicator_draws)
      pi_hat   = mean(vec_id)

      row = pairwise_effect_indices[v1, v2] + 1L
      for (h in 1:(num_groups - 1)) {
        param_counter = param_counter + 1
        col_index = h * (num_variables * (num_variables - 1) / 2) + row
        effect_draws = pairwise_effect_samples[, , col_index]
        vec_slab   = as.vector(effect_draws)

        eap_slab   = if (any(vec_id == 1)) mean(vec_slab[vec_id == 1]) else 0
        posterior_mean = pi_hat * eap_slab
        posterior_var  = (eap_slab^2 * var(vec_id)) + (pi_hat^2 * var(vec_slab))
        posterior_sd   = sqrt(posterior_var)

        est  = compute_rhat_ess(effect_draws)
        mcse_est = if (is.finite(est$ess) && est$ess > 0) posterior_sd / sqrt(est$ess) else NA

        results[[param_counter]] = data.frame(
          parameter = if (!is.null(param_names)) {
            param_names[param_counter]
          } else {
            paste0("V", v1, "-", "V", v2, "(diff", h, ")")
          },
          mean  = posterior_mean,
          sd    = posterior_sd,
          mcse  = mcse_est,
          n_eff = est$ess,
          Rhat  = est$rhat,
          check.names = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}




summarize_fit_compare = function(
    fit,
    main_effect_indices,
    pairwise_effect_indices,
    num_variables,
    num_groups,
    param_names_main = NULL,
    param_names_pairwise = NULL,
    param_names_main_diff = NULL,
    param_names_pairwise_diff = NULL,
    param_names_indicators = NULL
) {
  # Helper: extract baseline columns (g = 1) for each row
  extract_baseline = function(array3d, num_rows, num_groups) {
    idx = as.vector(sapply(0:(num_rows - 1), function(r) r * num_groups + 1))
    array3d[, , idx, drop = FALSE]
  }

  # --- main baseline
  array3d_main = combine_chains_compare(fit, "main_samples")
  n_main = nrow(main_effect_indices)  # number of rows in main_effects
  array3d_main_baseline = extract_baseline(array3d_main, n_main, num_groups)
  main_baseline = summarize_manual_compare(
    array3d_main_baseline,
    "main_samples",
    param_names = param_names_main
  )

  # --- pairwise baseline
  array3d_pair = combine_chains_compare(fit, "pairwise_samples")
  n_pair = num_variables * (num_variables - 1) / 2
  array3d_pair_baseline = extract_baseline(array3d_pair, n_pair, num_groups)
  pairwise_baseline = summarize_manual_compare(
    array3d_pair_baseline,
    "pairwise_samples",
    param_names = param_names_pairwise
  )

  # --- main differences
  main_differences = summarize_main_diff_compare(
    fit,
    main_effect_indices,
    num_groups,
    param_names = param_names_main_diff
  )

  # --- pairwise differences
  pairwise_differences = summarize_pairwise_diff_compare(
    fit,
    pairwise_effect_indices,
    num_variables,
    num_groups,
    param_names = param_names_pairwise_diff
  )

  # --- indicators
  indicators = summarize_indicator_compare(
    fit,
    "indicator_samples",
    param_names = param_names_indicators
  )

  list(
    main_baseline        = main_baseline,
    pairwise_baseline    = pairwise_baseline,
    main_differences     = main_differences,
    pairwise_differences = pairwise_differences,
    indicators           = indicators
  )
}


