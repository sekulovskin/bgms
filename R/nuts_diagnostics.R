summarize_nuts_diagnostics <- function(out, nuts_max_depth = 10, verbose = TRUE) {
  # Filter chains that include NUTS diagnostics
  nuts_chains <- Filter(function(chain) {
    all(c("treedepth__", "divergent__", "energy__") %in% names(chain))
  }, out)

  if(length(nuts_chains) == 0) {
    stop("No NUTS diagnostics found in output.")
  }

  # Combine per-field into matrices (chains × iterations)
  combine_diag <- function(field) {
    do.call(rbind, lapply(nuts_chains, function(chain) as.numeric(chain[[field]])))
  }

  treedepth_mat <- combine_diag("treedepth__")
  divergent_mat <- combine_diag("divergent__")
  energy_mat <- combine_diag("energy__")

  # Per-chain E-BFMI calculation
  compute_ebfmi <- function(energy) {
    diffs <- diff(energy)
    mean(diffs^2) / stats::var(energy)
  }
  ebfmi_per_chain <- apply(energy_mat, 1, compute_ebfmi)

  # Summaries
  n_total <- nrow(divergent_mat) * ncol(divergent_mat)
  total_divergences <- sum(divergent_mat)
  max_tree_depth_hits <- sum(treedepth_mat == nuts_max_depth)
  min_ebfmi <- min(ebfmi_per_chain)
  low_ebfmi_chains <- which(ebfmi_per_chain < 0.2)

  divergence_rate <- total_divergences / n_total
  depth_hit_rate <- max_tree_depth_hits / n_total

  if(verbose) {
    issues <- character(0)

    # Divergences
    if (total_divergences > 0) {
      if (divergence_rate > 0.001) {
        issues <- c(issues, sprintf(
          "Divergences: %d (%.2f%%) - increase target acceptance or use adaptive-metropolis",
          total_divergences, 100 * divergence_rate))
      } else {
        issues <- c(issues, sprintf(
          "Divergences: %d (%.3f%%) - check R-hat and ESS",
          total_divergences, 100 * divergence_rate))
      }
    }

    # Tree depth
    if (max_tree_depth_hits > 0) {
      if (depth_hit_rate > 0.01) {
        issues <- c(issues, sprintf(
          "Tree depth: %d hits (%.1f%%) - consider max_depth > %d",
          max_tree_depth_hits, 100 * depth_hit_rate, nuts_max_depth))
      } else {
        issues <- c(issues, sprintf(
          "Tree depth: %d hits (%.2f%%) - check ESS",
          max_tree_depth_hits, 100 * depth_hit_rate))
      }
    }

    # E-BFMI
    if (length(low_ebfmi_chains) > 0) {
      issues <- c(issues, sprintf(
        "E-BFMI: %.3f in chain%s %s - try adaptive-metropolis or more warmup",
        min_ebfmi,
        if(length(low_ebfmi_chains) > 1) "s" else "",
        paste(low_ebfmi_chains, collapse = ", ")))
    }

    # Only print if there are issues and verbose is enabled
    if (length(issues) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
      cat("NUTS issues:\n")
      for (issue in issues) {
        cat("  -", issue, "\n")
      }
    }
  }

  # Return structured summary
  invisible(list(
    treedepth = treedepth_mat,
    divergent = divergent_mat,
    energy = energy_mat,
    ebfmi = ebfmi_per_chain,
    summary = list(
      total_divergences = total_divergences,
      max_tree_depth_hits = max_tree_depth_hits,
      min_ebfmi = min_ebfmi
    )
  ))
}
