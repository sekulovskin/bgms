# Summary utilities for spike-and-slab MCMC output

# Combine MCMC chains into a 3D array [niter x nchains x nparam]
combine_chains <- function(fit, component) {
  nchains <- length(fit)
  samples_list <- lapply(fit, function(x) x[[component]])
  niter <- nrow(samples_list[[1]])
  nparam <- ncol(samples_list[[1]])
  array3d <- array(NA_real_, dim = c(niter, nchains, nparam))
  for (i in seq_len(nchains)) {
    array3d[, i, ] <- samples_list[[i]]
  }
  array3d
}

# Compute effective sample size and Rhat (Gelman-Rubin diagnostic)
compute_rhat_ess <- function(draws) {
  tryCatch({
    nchains <- ncol(draws)
    if (nchains > 1) {
      #mcmc_list <- coda::mcmc.list(apply(draws, 2, coda::mcmc))
      mcmc_list <- coda::mcmc.list(lapply(1:ncol(draws), function(i) coda::mcmc(draws[, i])))
      ess <- coda::effectiveSize(mcmc_list)
      rhat <- coda::gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1]
    } else {
      ess <- coda::effectiveSize(draws)
      rhat <- NA_real_
    }
    list(ess = ess, rhat = rhat)
  }, error = function(e) list(ess = NA_real_, rhat = NA_real_))
}

# Basic summarizer for continuous parameters
summarize_manual <- function(fit, component = c("main_samples", "pairwise_samples"), param_names = NULL) {
  component <- match.arg(component) # Add options later
  array3d <- combine_chains(fit, component)
  nparam <- dim(array3d)[3]

  result <- matrix(NA, nparam, 5)
  colnames(result) <- c("mean", "mcse", "sd", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws <- array3d[, , j]
    vec <- as.vector(draws)
    result[j, "mean"] <- mean(vec)
    result[j, "sd"] <- sd(vec)
    est <- compute_rhat_ess(draws)
    result[j, "mcse"] <- sd(vec) / sqrt(est$ess)
    result[j, "n_eff"] <- est$ess
    result[j, "Rhat"] <- est$rhat
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("parameter [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }

}

# Summarize binary indicator variables
summarize_indicator <- function(fit, component = c("indicator_samples"), param_names = NULL) {
  component <- match.arg(component) # Add options later
  array3d <- combine_chains(fit, component)

  nparam <- dim(array3d)[3]
  nchains <- dim(array3d)[2]
  niter <- dim(array3d)[1]

  result <- matrix(NA, nparam, 9)
  colnames(result) <- c("mean", "sd", "mcse", "n0->1", "n0->0", "n1->0", "n1->1", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws <- array3d[, , j]
    vec <- as.vector(draws)
    T <- length(vec)
    g_next <- vec[-1]
    g_curr <- vec[-T]

    p_hat <- mean(vec)
    sd <- sqrt(p_hat * (1 - p_hat))
    n00 <- sum(g_curr == 0 & g_next == 0)
    n01 <- sum(g_curr == 0 & g_next == 1)
    n10 <- sum(g_curr == 1 & g_next == 0)
    n11 <- sum(g_curr == 1 & g_next == 1)

    if (any(c(n01, n10) == 0)) {
      n_eff <- mcse <- R <- NA_real_
    } else {
      a <- n01 / (n00 + n01)
      b <- n10 / (n10 + n11)
      tau_int <- (2 - (a + b)) / (a + b)
      n_eff <- T / tau_int
      mcse <- sd / sqrt(n_eff)
      est <- compute_rhat_ess(draws)
      R <- est$rhat
    }

    result[j, ] <- c(p_hat, sd, mcse, n01, n00, n10, n11, n_eff, R)
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("indicator [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = paste0(param_names, "- indicator"),
               result, check.names = FALSE)
  }
}

# Summarize slab values where indicators are 1
summarize_slab <- function(fit, component = c("pairwise_samples"), param_names = NULL) {
  component <- match.arg(component) # Add options later
  array3d <- combine_chains(fit, component)
  nparam <- dim(array3d)[3]
  result <- matrix(NA, nparam, 6)
  colnames(result) <- c("mean", "sd", "mcse", "n", "n_eff", "Rhat")

  for (j in seq_len(nparam)) {
    draws <- array3d[, , j]
    vec <- as.vector(draws)
    nonzero <- vec != 0
    vec <- vec[nonzero]
    T <- length(vec)

    if (T > 10) {
      eap <- mean(vec)
      sdev <- sd(vec)
      mcse <- sdev / sqrt(T)
      est <- compute_rhat_ess(draws)
      result[j, ] <- c(eap, sdev, mcse, T, est$ess, est$rhat)
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
summarize_pair <- function(fit,
                           indicator_component = c("indicator_samples"),
                           slab_component = c("pairwise_samples"),
                           param_names = NULL
) {
  indicator_component <- match.arg(indicator_component) # Add options later
  slab_component <- match.arg(slab_component) # Add options later

  summ_ind <- summarize_indicator(fit, component = indicator_component)
  summ_slab <- summarize_slab(fit, component = slab_component)
  nparam <- nrow(summ_ind)

  eap <- summ_ind$mean * summ_slab$mean
  v <- (summ_slab$mean^2 * summ_ind$sd^2) + (summ_ind$mean^2 * summ_slab$sd^2)
  mcse2 <- (summ_slab$mean^2 * summ_ind$mcse^2) + (summ_ind$mean^2 * summ_slab$mcse^2)
  mcse <- sqrt(mcse2)
  sd <- sqrt(v)
  n_eff <- v / mcse2

  rhat <- rep(NA_real_, nparam)
  array3d_pw <- combine_chains(fit, slab_component)
  array3d_id <- combine_chains(fit, indicator_component)
  nchains <- dim(array3d_pw)[2]
  T <- prod(dim(array3d_pw)[1:2])

  for (j in seq_len(nparam)) {
    draws_pw <- array3d_pw[, , j]
    draws_id <- array3d_id[, , j]
    if (nchains > 1) {
      chain_means <- numeric(nchains)
      chain_vars <- numeric(nchains)
      for (chain in 1:nchains) {
        pi <- mean(draws_id[, chain])
        tmp <- draws_pw[, chain]
        e <- mean(tmp[tmp != 0])
        v <- var(tmp[tmp != 0])
        chain_means[chain] <- pi * e
        chain_vars[chain] <- pi * (v + (1 - pi) * e^2)
      }
      eap_mean <- mean(chain_means)
      B <- T * sum((chain_means - eap_mean)^2) / (nchains - 1)
      W <- mean(chain_vars)
      V <- (T - 1) * W / T + B / T
      rhat[j] <- sqrt(V / W)
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
summarize_fit <- function(fit, edge_selection = FALSE) {
  main_summary <- summarize_manual(fit, component = "main_samples")

  if (!edge_selection) {
    pair_summary <- summarize_manual(fit, component = "pairwise_samples")
  } else {
    # Get indicators and slab draws
    ind_summary <- summarize_indicator(fit, component = "indicator_samples")
    slab_summary <- summarize_slab(fit, component = "pairwise_samples")

    all_selected <- ind_summary$mean == 1

    # Use summarize_pair only where not always selected
    full_summary <- summarize_pair(fit,
                                   indicator_component = "indicator_samples",
                                   slab_component = "pairwise_samples")
    manual_summary <- summarize_manual(fit, component = "pairwise_samples")

    # Replace rows in full_summary with manual results for fully selected entries
    if (any(all_selected)) {
      full_summary[all_selected, ] <- manual_summary[all_selected, ]
    }

    pair_summary <- full_summary
  }

  list(main = main_summary, pairwise = pair_summary)
}