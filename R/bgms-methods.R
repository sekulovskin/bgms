#' @name print.bgms
#' @title Print method for `bgms` objects
#'
#' @description Minimal console output for `bgms` fit objects.
#'
#' @param x An object of class `bgms`.
#' @param ... Ignored.
#'
#' @export
print.bgms <- function(x, ...) {
  arguments <- extract_arguments(x)

  # Model type
  if (isTRUE(arguments$edge_selection)) {
    prior_msg <- switch(arguments$edge_prior,
                        "Bernoulli" = "Bayesian Edge Selection using a Bernoulli prior on edge inclusion",
                        "Beta-Bernoulli" = "Bayesian Edge Selection using a Beta-Bernoulli prior on edge inclusion",
                        "Stochastic-Block" = "Bayesian Edge Selection using a Stochastic Block prior on edge inclusion",
                        "Bayesian Edge Selection"
    )
    cat(prior_msg, "\n")
  } else {
    cat("Bayesian Estimation\n")
  }

  # Dataset info
  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if (isTRUE(arguments$na_impute)) {
    cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
  }

  # Iterations and chains
  if (!is.null(arguments$num_chains)) {
    total_iter <- arguments$iter * arguments$num_chains
    cat(paste0(" Number of post-burnin MCMC iterations: ", total_iter, "\n"))
    cat(paste0(" Number of MCMC chains: ", arguments$num_chains, "\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, "\n"))
  }

  cat("Use the `summary()` function for posterior summaries and chain diagnostics.\n")
  cat("See the `easybgm` package for summary and plotting tools.\n")
}

#' @name print.bgmCompare
#' @title  Print method for \code{bgms} objects
#'
#' @description Used to prevent bgms output cluttering the console.
#'
#' @param x An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @export
print.bgmCompare <- function(x, ...) {
  arguments = extract_arguments(x)
  if(arguments$difference_selection) {
    if(arguments$pairwise_difference_prior == "Bernoulli") {
      cat(paste0("Bayesian Variable Selection using a Bernoulli prior on the inclusion of \n",
                 "differences in pairwise interactions\n"))
    } else {
      cat(paste0("Bayesian Variable Selection using a Beta-Bernoulli prior on the inclusion of \n",
                 "differences in pairwise interactions\n"))
    }
    if(arguments$main_difference_model == "Free") {
      cat("Group specific category threshold parameters were estimated")
    } else {
      if(arguments$main_difference_prior == "Bernoulli") {
        cat(paste0("Bayesian Variable Selection using a Bernoulli prior on the inclusion of\n",
                   "differences in the category thresholds\n"))
      } else {
        cat(paste0("Bayesian Variable Selection using a Beta-Bernoulli prior on the inclusion of\n",
                   "differences in the category thresholds\n"))      }
    }
  } else {
    cat("Bayesian Estimation\n")
  }

  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  num_groups = length(unique(arguments$group))
  if(arguments$na_impute) {
    for(group in 1:num_groups) {
      cat(paste0(" Number of cases Group ", group,": ", arguments$num_cases[group], " (missings imputed)\n"))
    }
  } else {
    for(group in 1:num_groups) {
      cat(paste0(" Number of cases Group ", group,": ", arguments$num_cases[group],"\n"))
    }
  }
  if(arguments$save) {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (MCMC output saved)\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (posterior means saved)\n"))
  }
  cat("See the easybgm package for extensive summary and plotting functions \n")
}



#' @name summary.bgms
#' @title Summary method for `bgms` objects
#'
#' @description Returns posterior summaries and diagnostics for a fitted `bgms` model.
#'
#' @param object An object of class `bgms`.
#' @param ... Currently ignored.
#'
#' @return An object of class `summary.bgms` with posterior summaries.
#' @export
summary.bgms <- function(object, ...) {
  arguments <- extract_arguments(object)

  if (!is.null(object$posterior_summary_main) && !is.null(object$posterior_summary_pairwise)) {
    out <- list(
      main = object$posterior_summary_main,
      pairwise = object$posterior_summary_pairwise
    )
    class(out) <- "summary.bgms"
    return(out)
  }

  message("No summary statistics available for this model object.\n",
          "Try fitting the model again using the latest bgms version,\n",
          "or use the `easybgm` package for diagnostic summaries and plotting.")
  invisible(NULL)
}



#' @export
print.summary.bgms <- function(x, digits = 3, ...) {
  cat("Posterior summaries from Bayesian estimation:\n\n")

  if (!is.null(x$main)) {
    cat("Category thresholds:\n")
    print(round(head(x$main, 6), digits = digits))
    if (nrow(x$main) > 6) cat("... (use `summary(fit)$main` to see full output)\n")
    cat("\n")
  }

  if (!is.null(x$pairwise)) {
    cat("Pairwise interactions:\n")
    print(round(head(x$pairwise, 6), digits = digits))
    if (nrow(x$pairwise) > 6) cat("... (use `summary(fit)$pairwise` to see full output)\n")
    cat("\n")
  }

  if (!is.null(x$indicator)) {
    cat("Inclusion probabilities:\n")
    print(round(head(x$indicator, 6), digits = digits))
    if (nrow(x$indicator) > 6) cat("... (use `summary(fit)$indicator` to see full output)\n")
    cat("\n")
  }

  cat("Use `summary(fit)$<component>` to access full results.\n")
  cat("See the `easybgm` package for other summary and plotting tools.\n")
}



#' @title Extract Coefficients from a bgms Object
#' @name coef.bgms
#' @description Returns the posterior mean thresholds, pairwise effects, and edge inclusion indicators from a \code{bgms} model fit.
#'
#' @param object An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{main}{Posterior mean of the category threshold parameters.}
#'   \item{pairwise}{Posterior mean of the pairwise interaction matrix.}
#'   \item{indicator}{Posterior mean of the edge inclusion indicators (if available).}
#' }
#'
#' @export
coef.bgms <- function(object, ...) {
  out <- list(
    main = object$posterior_mean_main,
    pairwise = object$posterior_mean_pairwise
  )
  if (!is.null(object$posterior_mean_indicator)) {
    out$indicator <- object$posterior_mean_indicator
  }
  return(out)
}



#' @export
as_draws.bgms <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Install the 'posterior' package to use this method.")
  }

  arguments <- extract_arguments(x)

  # Helper to collapse list of chains into matrix and assign names
  make_matrix_block <- function(samples_list, param_names = NULL) {
    if (is.null(samples_list) || all(vapply(samples_list, is.null, logical(1)))) {
      return(NULL)
    }
    mat <- do.call(rbind, samples_list)
    if (!is.null(param_names)) colnames(mat) <- param_names
    return(mat)
  }

  # Helper to prefix column names
  prefix_colnames <- function(mat, prefix) {
    if (!is.null(mat)) {
      colnames(mat) <- paste0(prefix, colnames(mat))
    }
    mat
  }

  # Assemble all blocks
  main_mat <- make_matrix_block(x$raw_samples$main_samples, x$raw_samples$parameter_names$main)
  pairwise_mat <- make_matrix_block(x$raw_samples$pairwise_samples, x$raw_samples$parameter_names$pairwise)
  indicator_mat <- make_matrix_block(x$raw_samples$indicator, x$raw_samples$parameter_names$indicator)

  # Apply prefixes to avoid name duplication
  main_mat <- prefix_colnames(main_mat, "main_")
  pairwise_mat <- prefix_colnames(pairwise_mat, "pairwise_")
  indicator_mat <- prefix_colnames(indicator_mat, "indicator_")

  # Issue spike-and-slab warning (only once)
  if (isTRUE(arguments$edge_selection) && !getOption("bgms.as_draws_warning_given", FALSE)) {
    warning("This model includes spike-and-slab posteriors (e.g., pairwise effects or indicators).\n",
            "Posterior summaries and diagnostics for these require special treatment.\n",
            "See `vignette(\"spike_slab_manual\", package = \"bgms\")` for details.",
            call. = FALSE)
    options(bgms.as_draws_warning_given = TRUE)
  }

  # Combine all available samples
  all_samples <- cbind(main_mat, pairwise_mat, indicator_mat)

  # Return draws_df
  posterior::as_draws_df(all_samples)
}




.warning_issued <- FALSE
warning_once <- function(msg) {
  if (!.warning_issued) {
    warning(msg, call. = FALSE)
    .warning_issued <<- TRUE
  }
}
