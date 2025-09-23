##' Extractor Functions for bgms Objects
#'
#' These functions extract various components from objects returned by the `bgm()` function,
#' such as edge indicators, posterior inclusion probabilities, and parameter summaries.
#'
#' @section Functions:
#' - `extract_arguments()` – Extract model arguments
#' - `extract_indicators()` – Get sampled edge indicators
#' - `extract_posterior_inclusion_probabilities()` – Posterior edge inclusion probabilities
#' - `extract_pairwise_interactions()` – Posterior mean of pairwise interactions
#' - `extract_category_thresholds()` – Posterior mean of category thresholds
#' - `extract_indicator_priors()` – Prior structure used for edge indicators
#'
#' @name extractor_functions
#' @title Extractor Functions for bgms Objects
#' @keywords internal
NULL

#' @export
extract_arguments <- function(bgms_object) {
  UseMethod("extract_arguments")
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgms <- function(bgms_object) {
  if (!inherits(bgms_object, "bgms")) stop("Object must be of class 'bgms'.")
  if (is.null(bgms_object$arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }
  return(bgms_object$arguments)
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgmCompare <- function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(bgms_object$arguments)
  }
}

#' @rdname extractor_functions
#' @export
extract_indicators <- function(bgms_object) {
  UseMethod("extract_indicators")
}

#' @rdname extractor_functions
#' @export
extract_indicators.bgms <- function(bgms_object) {
  arguments <- extract_arguments(bgms_object)

  if (!isTRUE(arguments$edge_selection)) {
    stop("To access edge indicators, the model must be run with edge_selection = TRUE.")
  }

  # Resolve indicator samples
  indicators_list <- bgms_object$raw_samples$indicator
  if (is.null(indicators_list)) {
    if (!is.null(bgms_object$indicator)) {
      indicators_list <- bgms_object$indicator
    } else if (!is.null(bgms_object$gamma)) {
      indicators_list <- bgms_object$gamma
    } else {
      stop("No indicator samples found in this object.")
    }
  }

  # Combine all chains
  indicator_samples <- do.call(rbind, indicators_list)

  # Assign column names if available
  param_names <- bgms_object$raw_samples$parameter_names$indicator
  if (!is.null(param_names)) {
    colnames(indicator_samples) <- param_names
  }

  return(indicator_samples)
}

#' @rdname extractor_functions
#' @export
extract_indicators.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(arguments$difference_selection & arguments$save) {
    pairwise_difference_indicator = bgms_object$pairwise_difference_indicator
    if(arguments$independent_thresholds == FALSE) {
      main_difference_indicator = bgms_object$main_difference_indicator
    } else {
      main_difference_indicator = NULL
    }
    return(list(main_difference_indicator = main_difference_indicator,
                pairwise_difference_indicator = pairwise_difference_indicator))
  } else {
    stop(paste0("To access the sampled difference indicators the bgmCompare function needs to be run using \n",
                "difference_selection = TRUE and save = TRUE."))
  }
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities <- function(bgms_object) {
  UseMethod("extract_posterior_inclusion_probabilities")
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities.bgms <- function(bgms_object) {
  arguments <- extract_arguments(bgms_object)

  if (!isTRUE(arguments$edge_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgm() with edge_selection = TRUE.")
  }

  num_vars <- arguments$num_variables
  data_columnnames <- arguments$data_columnnames

  edge_means <- NULL
  # New format: use extract_indicators()
  if (!is.null(bgms_object$raw_samples$indicator)) {
    indicator_samples <- extract_indicators(bgms_object)
    edge_means <- colMeans(indicator_samples)
  } else if (!is.null(bgms_object$indicator)) {
    if (!is.null(arguments$save) && isTRUE(arguments$save)) {
      edge_means <- colMeans(bgms_object$indicator)
    } else {
      edge_means <- bgms_object$indicator
    }
  } else if (!is.null(bgms_object$gamma)) {
    if (!is.null(arguments$save) && isTRUE(arguments$save)) {
      edge_means <- colMeans(bgms_object$gamma)
    } else {
      edge_means <- bgms_object$gamma
    }
  } else {
    stop("No indicator data found to compute posterior inclusion probabilities.")
  }

  pip_matrix <- matrix(0, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] <- edge_means
  pip_matrix = pip_matrix + t(pip_matrix)

  colnames(pip_matrix) <- data_columnnames
  rownames(pip_matrix) <- data_columnnames

  return(pip_matrix)
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$difference_selection) {
    stop(paste0("To estimate the posterior inclusion probabilities for the between-group \n",
                "parameter differences , please run the bgmCompare function with \n",
                "difference_selection = TRUE."))
  }

  if(arguments$save) {
    pairwise_difference_means = colMeans(bgms_object$pairwise_difference_indicator)
    num_variables = arguments$num_variables

    posterior_inclusion_probabilities = matrix(0, num_variables, num_variables)
    posterior_inclusion_probabilities[lower.tri(posterior_inclusion_probabilities)] = pairwise_difference_means
    posterior_inclusion_probabilities = posterior_inclusion_probabilities +
      t(posterior_inclusion_probabilities)
    if(!arguments$independent_thresholds) {
      main_difference_means = colMeans(bgms_object$main_difference_indicator)
      diag(posterior_inclusion_probabilities) = main_difference_means
    }

    data_columnnames = arguments$data_columnnames
    colnames(posterior_inclusion_probabilities) = data_columnnames
    rownames(posterior_inclusion_probabilities) = data_columnnames

  } else {
    posterior_inclusion_probabilities = bgms_object$indicator
  }
  return(posterior_inclusion_probabilities)
}

#' @rdname extractor_functions
#' @export
extract_indicator_priors <- function(bgms_object) {
  UseMethod("extract_indicator_priors")
}

#' @rdname extractor_functions
#' @export
extract_indicator_priors.bgms <- function(bgms_object) {
  arguments <- extract_arguments(bgms_object)
  if (!isTRUE(arguments$edge_selection)) stop("No edge selection performed.")

  switch(
    arguments$edge_prior,
    "Bernoulli" = list(type = "Bernoulli", prior_inclusion_probability = arguments$inclusion_probability),
    "Beta-Bernoulli" = list(type = "Beta-Bernoulli", alpha = arguments$beta_bernoulli_alpha, beta = arguments$beta_bernoulli_beta),
    "Stochastic-Block" = list(
      type = "Stochastic-Block",
      beta_bernoulli_alpha = arguments$beta_bernoulli_alpha,
      beta_bernoulli_beta = arguments$beta_bernoulli_beta,
      dirichlet_alpha = arguments$dirichlet_alpha
    )
  )
}


#' @rdname extractor_functions
#' @export
extract_indicator_priors.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$difference_selection) {
    stop(paste0("The bgmCompare function did not perform selection on the between-group\n",
                "differences, so there are no indicator priors specified."))
  } else {
    if(arguments$pairwise_difference_prior == "Bernoulli") {
      difference_prior = list(pairwise_type = "Bernoulli",
                              prior_inclusion_probability = arguments$inclusion_probability_difference)
    } else {
      difference_prior = list(pairwise_type = "Beta-Bernoulli",
                              pairwise_alpha = arguments$pairwise_beta_bernoulli_alpha,
                              pairwise_beta = arguments$pairwise_beta_bernoulli_beta)
    }
    if(!arguments$independent_thresholds) {
      if(arguments$main_difference_prior == "Bernoulli") {
        difference_prior$main_type = "Bernoulli"
      } else {
        difference_prior$main_type = "Beta-Bernoulli"
        difference_prior$main_alpha = arguments$beta_bernoulli_alpha
        difference_prior$main_beta = arguments$beta_bernoulli_beta
      }
    }
  }
  return(difference_prior)
}



#' @rdname extractor_functions
#' @export
extract_pairwise_interactions <- function(bgms_object) {
  UseMethod("extract_pairwise_interactions")
}

#' @rdname extractor_functions
#' @export
extract_pairwise_interactions.bgms <- function(bgms_object) {
  arguments <- extract_arguments(bgms_object)
  num_vars <- arguments$num_variables
  var_names <- arguments$data_columnnames

  if(!is.null(bgms_object$raw_samples)) {
    nchains = length(bgms_object$raw_samples$pairwise)
    mat = NULL
    mats <- bgms_object$raw_samples$pairwise
    mat  <- do.call(rbind, mats)
  } else if (!is.null(bgms_object$posterior_summary_pairwise)) {
    vec <- bgms_object$posterior_summary_pairwise[, "mean"]
    mat <- matrix(0, nrow = num_vars, ncol = num_vars)
    mat[upper.tri(mat)] <- vec
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  } else if (!is.null(bgms_object$posterior_mean_pairwise)) {
    mat <- bgms_object$posterior_mean_pairwise
  } else if (!is.null(bgms_object$pairwise_effects)) {
    mat <- bgms_object$pairwise_effects
  } else {
    stop("No pairwise interaction effects found in the object.")
  }

  dimnames(mat) <- list(var_names, var_names)
  return(mat)
}


#' @rdname extractor_functions
#' @export
extract_category_thresholds <- function(bgms_object) {
  UseMethod("extract_category_thresholds")
}

#' @rdname extractor_functions
#' @export
extract_category_thresholds.bgms <- function(bgms_object) {
  arguments <- extract_arguments(bgms_object)
  var_names <- arguments$data_columnnames

  if (!is.null(bgms_object$posterior_summary_main)) {
    vec <- bgms_object$posterior_summary_main[, "mean"]
    num_vars <- arguments$num_variables
    variable_type <- arguments$variable_type
    if(length(variable_type) == 1) {
      variable_type <- rep(variable_type, num_vars)
    }
    num_cats <- arguments$num_categories
    max_cats <- max(num_cats)
    mat <- matrix(NA_real_, nrow = num_vars, ncol = max_cats)
    rownames(mat) <- var_names
    pos <- 1
    for (v in seq_len(num_vars)) {
      if (variable_type[v] == "ordinal") {
        k <- num_cats[v]
        mat[v, 1:k] <- vec[pos:(pos + k - 1)]
        pos <- pos + k
      } else {
        mat[v, 1:2] <- vec[pos:(pos + 1)]
        pos <- pos + 2
      }
    }
    return(mat)
  } else if (!is.null(bgms_object$posterior_mean_main)) {
    mat <- bgms_object$posterior_mean_main
  } else if (!is.null(bgms_object$main_effects)) {
    mat <- bgms_object$main_effects
  } else {
    stop("No threshold or main effects found in the object.")
  }

  rownames(mat) <- var_names
  return(mat)
}

#' @rdname extractor_functions
#' @export
extract_pairwise_difference.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$pairwise_difference)
}

#' @rdname extractor_functions
#' @export
extract_main_difference.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$main_difference)
}

#' @rdname extractor_functions
#' @export
extract_edge_indicators <- function(bgms_object) {
  warning(paste0("The ``extract_edge_indicators'' function is deprecated and will be removed in a \n",
                 "future release of bgms. Please use the ``extract_indicators'' function instead."))
  return(extract_indicators(bgms_object))
}

#' @rdname extractor_functions
#' @export
extract_pairwise_thresholds <- function(bgms_object) {
  warning(paste0("The ``extract_pairwise_thresholds'' function is deprecated and will be removed in a \n",
                 "future release of bgms. Please use the ``extract_category_thresholds'' function instead."))
  return(extract_category_thresholds(bgms_object))
}