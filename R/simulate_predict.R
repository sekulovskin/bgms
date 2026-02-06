# ==============================================================================
#   MRF Simulation and Prediction Functions
#
#   This file contains:
#     - simulate_mrf(): Standalone MRF data simulation
#     - mrfSampler(): Deprecated wrapper for simulate_mrf()
#     - simulate.bgms(): S3 method for simulating from fitted models
#     - predict.bgms(): S3 method for conditional probability prediction
# ==============================================================================


#' Simulate Observations from an Ordinal MRF
#'
#' @description
#' `simulate_mrf()` generates observations from an ordinal Markov Random Field
#' using Gibbs sampling with user-specified parameters.
#'
#' @details
#' The Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each variable from a logistic
#' model using the other variable states as predictor variables.
#'
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of baseline_category category, such as the neutral
#' category in a Likert scale. The Blume-Capel model specifies the following
#' quadratic model for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times (\text{c} - \text{r}) + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times (\text{c} - \text{r}) + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' baseline_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a
#' preference for responding in the baseline_category category.
#'
#' @param no_states The number of states of the ordinal MRF to be generated.
#'
#' @param no_variables The number of variables in the ordinal MRF.
#'
#' @param no_categories Either a positive integer or a vector of positive
#' integers of length \code{no_variables}. The number of response categories on top
#' of the base category: \code{no_categories = 1} generates binary states.
#'
#' @param interactions A symmetric \code{no_variables} by \code{no_variables} matrix of
#' pairwise interactions. Only its off-diagonal elements are used.
#'
#' @param thresholds A \code{no_variables} by \code{max(no_categories)} matrix of
#' category thresholds. The elements in row \code{i} indicate the thresholds of
#' variable \code{i}. If \code{no_categories} is a vector, only the first
#' \code{no_categories[i]} elements are used in row \code{i}. If the Blume-Capel
#' model is used for the category thresholds for variable \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#'
#' @param variable_type What kind of variables are simulated? Can be a single
#' character string specifying the variable type of all \code{p} variables at
#' once or a vector of character strings of length \code{p} specifying the type
#' for each variable separately. Currently, bgm supports ``ordinal'' and
#' ``blume-capel''. Binary variables are automatically treated as ``ordinal''.
#' Defaults to \code{variable_type = "ordinal"}.
#'
#' @param baseline_category An integer vector of length \code{no_variables} specifying the
#' baseline_category category that is used for the Blume-Capel model (details below).
#' Can be any integer value between \code{0} and \code{no_categories} (or
#' \code{no_categories[i]}).
#'
#' @param iter The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @param seed Optional integer seed for reproducibility. If \code{NULL},
#' a seed is generated from R's random number generator (so \code{set.seed()}
#' can be used before calling this function).
#'
#' @return A \code{no_states} by \code{no_variables} matrix of simulated states of
#' the ordinal MRF.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' no_variables = 5
#' no_categories = sample(1:5, size = no_variables, replace = TRUE)
#'
#' Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#' Thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))
#'
#' x = simulate_mrf(
#'   no_states = 1e3,
#'   no_variables = no_variables,
#'   no_categories = no_categories,
#'   interactions = Interactions,
#'   thresholds = Thresholds
#' )
#'
#' # Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
#' no_variables = 5
#' no_categories = 4
#'
#' Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#'
#' Thresholds = matrix(NA, no_variables, no_categories)
#' Thresholds[, 1] = -1
#' Thresholds[, 2] = -1
#' Thresholds[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#' Thresholds[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#'
#' x = simulate_mrf(
#'   no_states = 1e3,
#'   no_variables = no_variables,
#'   no_categories = no_categories,
#'   interactions = Interactions,
#'   thresholds = Thresholds,
#'   variable_type = c("b", "b", "o", "b", "o"),
#'   baseline_category = 2
#' )
#'
#' @seealso \code{\link{simulate.bgms}} for simulating from a fitted model.
#'
#' @export
simulate_mrf = function(no_states,
                      no_variables,
                      no_categories,
                      interactions,
                      thresholds,
                      variable_type = "ordinal",
                      baseline_category,
                      iter = 1e3,
                      seed = NULL) {
  # Check no_states, no_variables, iter --------------------------------------------
  if(no_states <= 0 ||
    abs(no_states - round(no_states)) > .Machine$double.eps) {
    stop("``no_states'' needs be a positive integer.")
  }
  if(no_variables <= 0 ||
    abs(no_variables - round(no_variables)) > .Machine$double.eps) {
    stop("``no_variables'' needs be a positive integer.")
  }
  if(iter <= 0 ||
    abs(iter - round(iter)) > .Machine$double.eps) {
    stop("``iter'' needs be a positive integer.")
  }

  # Check no_categories --------------------------------------------------------
  if(length(no_categories) == 1) {
    if(no_categories <= 0 ||
      abs(no_categories - round(no_categories)) > .Machine$double.eps) {
      stop("``no_categories'' needs be a (vector of) positive integer(s).")
    }
    no_categories = rep(no_categories, no_variables)
  } else {
    for(variable in 1:no_variables) {
      if(no_categories[variable] <= 0 ||
        abs(no_categories[variable] - round(no_categories[variable])) >
          .Machine$double.eps) {
        stop(paste("For variable", variable, "``no_categories'' was not a positive integer."))
      }
    }
  }

  # Check variable specification -----------------------------------------------
  if(length(variable_type) == 1) {
    variable_type = match.arg(
      arg = variable_type,
      choices = c("ordinal", "blume-capel")
    )

    if(variable_type == "blume-capel" && any(no_categories < 2)) {
      stop(paste0(
        "The Blume-Capel model only works for ordinal variables with more than two \n",
        "response options. But variables ", which(no_categories < 2), " are binary variables."
      ))
    }
    variable_type = rep(variable_type, no_variables)
  } else {
    if(length(variable_type) != no_variables) {
      stop(paste0(
        "The argument ``variable_type'' should be either a single character string or a \n",
        "vector of character strings of length ``no_variables''."
      ))
    } else {
      for(variable in 1:no_variables) {
        variable_type[variable] = match.arg(
          arg = variable_type[variable],
          choices = c("ordinal", "blume-capel")
        )
      }
      if(any(variable_type == "blume-capel" & no_categories < 2)) {
        stop(paste0(
          "The Blume-Capel model only works for ordinal variables with more than two \n",
          "response options. But variables ",
          which(variable_type == "blume-capel" & no_categories < 2),
          " are binary variables."
        ))
      }
    }
  }

  # Check the baseline_category for Blume-Capel variables ---------------------
  if(any(variable_type == "blume-capel")) {
    if(length(baseline_category) == 1) {
      baseline_category = rep(baseline_category, no_variables)
    }
    if(any(baseline_category < 0) || any(abs(baseline_category - round(baseline_category)) > .Machine$double.eps)) {
      stop(paste0("For variables ",
                  which(baseline_category < 0),
                  " ``baseline_category'' was either negative or not integer."))
    }
    if(any(baseline_category - no_categories > 0)) {
      stop(paste0("For variables ",
                  which(baseline_category - no_categories > 0),
                  " the ``baseline_category'' category was larger than the maximum category value."))
    }
  }

  # Check interactions ---------------------------------------------------------
  if(!inherits(interactions, what = "matrix")) {
    interactions = as.matrix(interactions)
  }
  if(!isSymmetric(interactions)) {
    stop("The matrix ``interactions'' needs to be symmetric.")
  }
  if(nrow(interactions) != no_variables) {
    stop("The matrix ``interactions'' needs to have ``no_variables'' rows and columns.")
  }

  # Check the threshold values -------------------------------------------------
  if(!inherits(thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(thresholds) == no_variables) {
        thresholds = matrix(thresholds, ncol = 1)
      } else {
        stop(paste0(
          "The matrix ``thresholds'' has ",
          length(thresholds),
          " elements, but requires",
          no_variables,
          "."
        ))
      }
    } else {
      stop("``Thresholds'' needs to be a matrix.")
    }
  }

  if(nrow(thresholds) != no_variables) {
    stop("The matrix ``thresholds'' needs to be have ``no_variables'' rows.")
  }

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      if(anyNA(thresholds[variable, 1:no_categories[variable]])) {
        tmp = which(is.na(thresholds[variable, 1:no_categories[variable]]))

        string = paste(tmp, sep = ",")

        stop(paste0(
          "The matrix ``thresholds'' contains NA(s) for variable ",
          variable,
          " in category \n",
          "(categories) ",
          paste(which(is.na(thresholds[variable, 1:no_categories[variable]])), collapse = ", "),
          ", where a numeric value is needed."
        ))
      }
      if(ncol(thresholds) > no_categories[variable]) {
        if(!anyNA(thresholds[variable, (no_categories[variable] + 1):ncol(thresholds)])) {
          warning(paste0(
            "The matrix ``thresholds'' contains numeric values for variable ",
            variable,
            " for category \n",
            "(categories, i.e., columns) exceding the maximum of ",
            no_categories[variable],
            ". These values will \n",
            "be ignored."
          ))
        }
      }
    } else {
      if(anyNA(thresholds[variable, 1:2])) {
        stop(paste0(
          "The Blume-Capel model is chosen for the category thresholds of variable ",
          variable,
          ". \n",
          "This model has two parameters that need to be placed in columns 1 and 2, row \n",
          variable,
          ", of the ``thresholds'' input matrix. Currently, there are NA(s) in these \n",
          "entries, where a numeric value is needed."
        ))
      }
      if(ncol(thresholds) > 2) {
        if(!anyNA(thresholds[variable, 3:ncol(thresholds)])) {
          warning(paste0(
            "The Blume-Capel model is chosen for the category thresholds of variable ",
            variable,
            ". \n",
            "This model has two parameters that need to be placed in columns 1 and 2, row \n",
            variable,
            ", of the ``thresholds'' input matrix. However, there are numeric values \n",
            "in higher categories. These values will be ignored."
          ))
        }
      }
    }
  }

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      for(category in 1:no_categories[variable]) {
        if(!is.finite(thresholds[variable, category])) {
          stop(paste(
            "The threshold parameter for variable", variable, "and category",
            category, "is NA or not finite."
          ))
        }
      }
    } else {
      if(!is.finite(thresholds[variable, 1])) {
        stop(paste0(
          "The alpha parameter for the Blume-Capel model for variable ",
          variable,
          " is NA \n",
          " or not finite."
        ))
      }
      if(!is.finite(thresholds[variable, 2])) {
        stop(paste0(
          "The beta parameter for the Blume-Capel model for variable",
          variable,
          "is NA \n",
          " or not finite."
        ))
      }
    }
  }

  # Handle seed ----------------------------------------------------------------
  if(is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  } else {
    if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
      stop("Argument 'seed' must be a single non-negative integer.")
    }
    seed <- as.integer(seed)
  }

  # The Gibbs sampler ----------------------------------------------------------
  if(!any(variable_type == "blume-capel")) {
    x <- sample_omrf_gibbs(
      no_states = no_states,
      no_variables = no_variables,
      no_categories = no_categories,
      interactions = interactions,
      thresholds = thresholds,
      iter = iter,
      seed = seed
    )
  } else {
    x <- sample_bcomrf_gibbs(
      no_states = no_states,
      no_variables = no_variables,
      no_categories = no_categories,
      interactions = interactions,
      thresholds = thresholds,
      variable_type_r = variable_type,
      baseline_category = baseline_category,
      iter = iter,
      seed = seed
    )
  }

  return(x)
}


# ==============================================================================
#   mrfSampler() - Deprecated Wrapper for simulate_mrf()
# ==============================================================================

#' Sample observations from the ordinal MRF
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `mrfSampler()` was renamed to [simulate_mrf()] as of bgms 0.1.6.3 to
#' follow the package's naming conventions.
#'
#' @inheritParams simulate_mrf
#'
#' @return A matrix of simulated observations (see [simulate_mrf()]).
#'
#' @seealso [simulate_mrf()] for the current function.
#'
#' @keywords internal
#' @export
mrfSampler = function(no_states,
                      no_variables,
                      no_categories,
                      interactions,
                      thresholds,
                      variable_type = "ordinal",
                      baseline_category,
                      iter = 1e3,
                      seed = NULL) {
  lifecycle::deprecate_warn("0.1.6.3", "mrfSampler()", "simulate_mrf()")

  simulate_mrf(
    no_states = no_states,
    no_variables = no_variables,
    no_categories = no_categories,
    interactions = interactions,
    thresholds = thresholds,
    variable_type = variable_type,
    baseline_category = baseline_category,
    iter = iter,
    seed = seed
  )
}


# ==============================================================================
#   simulate.bgms() - S3 Method for Simulating from Fitted Models
# ==============================================================================

#' Simulate Data from a Fitted bgms Model
#'
#' @description
#' Generates new observations from the Markov Random Field model using the
#' estimated parameters from a fitted \code{bgms} object.
#'
#' @param object An object of class \code{bgms}.
#' @param nsim Number of observations to simulate. Default: \code{500}.
#' @param seed Optional random seed for reproducibility.
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters (faster,
#'       single simulation).}
#'     \item{\code{"posterior-sample"}}{Sample from posterior draws, producing
#'       one dataset per draw (accounts for parameter uncertainty). This method
#'       uses parallel processing when \code{cores > 1}.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL}, uses all available draws.
#' @param iter Number of Gibbs iterations for equilibration before collecting
#'   samples. Default: \code{1000}.
#' @param cores Number of CPU cores for parallel execution when
#'   \code{method = "posterior-sample"}.
#'   Default: \code{parallel::detectCores()}.
#' @param display_progress Character string specifying the type of progress bar.
#'   Options: \code{"per-chain"}, \code{"total"}, \code{"none"}.
#'   Default: \code{"per-chain"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' If \code{method = "posterior-mean"}: A matrix with \code{nsim} rows and
#' \code{p} columns containing simulated observations.
#'
#' If \code{method = "posterior-sample"}: A list of matrices, one per posterior
#' draw, each with \code{nsim} rows and \code{p} columns.
#'
#' @details
#' This function uses the estimated interaction and threshold parameters to
#' generate new data via Gibbs sampling. When \code{method = "posterior-sample"},
#' parameter uncertainty is propagated to the simulated data by using different
#' posterior draws. Parallel processing is available for this method via the
#' \code{cores} argument.
#'
#' @seealso \code{\link{predict.bgms}} for computing conditional probabilities,
#'   \code{\link{simulate_mrf}} for simulation with user-specified parameters.
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit <- bgm(x = Wenchuan[, 1:5])
#'
#' # Simulate 100 new observations using posterior means
#' new_data <- simulate(fit, nsim = 100)
#'
#' # Simulate with parameter uncertainty (10 datasets)
#' new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample", ndraws = 10)
#'
#' # Use parallel processing for faster simulation
#' new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample",
#'                           ndraws = 100, cores = 4)
#' }
#'
#' @importFrom stats simulate
#' @export
simulate.bgms <- function(object,
                          nsim = 500,
                          seed = NULL,
                          method = c("posterior-mean", "posterior-sample"),
                          ndraws = NULL,
                          iter = 1000,
                          cores = parallel::detectCores(),
                          display_progress = c("per-chain", "total", "none"),
                          ...) {
  method <- match.arg(method)
  progress_type <- progress_type_from_display_progress(display_progress)

  # Validate cores
  if(!is.numeric(cores) || length(cores) != 1 || is.na(cores) || cores < 1) {
    stop("Argument 'cores' must be a positive integer.")
  }
  cores <- as.integer(cores)

  # Setting the seed
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
      stop("Argument 'seed' must be a single non-negative integer.")
    }
    seed <- as.integer(seed)
  } else {
    # Generate a seed for the C++ RNG
    seed <- sample.int(.Machine$integer.max, 1L)
  }

  # Extract model information

  arguments <- extract_arguments(object)
  num_variables <- arguments$num_variables
  num_categories <- arguments$num_categories
  variable_type <- arguments$variable_type
  data_columnnames <- arguments$data_columnnames

  # Handle variable_type

  if(length(variable_type) == 1) {
    variable_type <- rep(variable_type, num_variables)
  }

  # Get baseline_category (for Blume-Capel variables)
  baseline_category <- arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category <- rep(0L, num_variables)
  }


  if(method == "posterior-mean") {
    # Use posterior mean parameters
    interactions <- object$posterior_mean_pairwise
    thresholds <- object$posterior_mean_main

    # Set R's RNG for simulate_mrf
    if(!is.null(seed)) set.seed(seed)

    # Call simulate_mrf
    result <- simulate_mrf(
      no_states = nsim,
      no_variables = num_variables,
      no_categories = num_categories,
      interactions = interactions,
      thresholds = thresholds,
      variable_type = variable_type,
      baseline_category = baseline_category,
      iter = iter
    )

    colnames(result) <- data_columnnames
    return(result)

  } else {
    # Use posterior samples with parallel processing
    pairwise_samples <- do.call(rbind, object$raw_samples$pairwise)
    main_samples <- do.call(rbind, object$raw_samples$main)

    total_draws <- nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws <- total_draws
    }
    ndraws <- min(ndraws, total_draws)

    # Sample which draws to use
    if(!is.null(seed)) set.seed(seed)
    draw_indices <- sample.int(total_draws, ndraws)

    # Call parallel C++ function
    results <- run_simulation_parallel(
      pairwise_samples = pairwise_samples,
      main_samples = main_samples,
      draw_indices = as.integer(draw_indices),
      no_states = as.integer(nsim),
      no_variables = as.integer(num_variables),
      no_categories = as.integer(num_categories),
      variable_type_r = variable_type,
      baseline_category = as.integer(baseline_category),
      iter = as.integer(iter),
      nThreads = cores,
      seed = seed,
      progress_type = progress_type
    )

    # Add column names
    for(i in seq_along(results)) {
      colnames(results[[i]]) <- data_columnnames
    }

    return(results)
  }
}


# ==============================================================================
#   predict.bgms() - S3 Method for Conditional Probability Prediction
# ==============================================================================

#' Predict Conditional Probabilities from a Fitted bgms Model
#'
#' @description
#' Computes conditional probability distributions for one or more variables
#' given the observed values of other variables in the data.
#'
#' @param object An object of class \code{bgms}.
#' @param newdata A matrix or data frame with \code{n} rows and \code{p} columns
#'   containing the observed data. Must have the same variables (columns) as
#'   the original data used to fit the model.
#' @param variables Which variables to predict. Can be:
#'   \itemize{
#'     \item A character vector of variable names
#'     \item An integer vector of column indices
#'     \item \code{NULL} (default) to predict all variables
#'   }
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"probabilities"}}{Return the full conditional probability
#'       distribution for each variable and observation.}
#'     \item{\code{"response"}}{Return the predicted category (mode of the
#'       conditional distribution).}
#'   }
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters.}
#'     \item{\code{"posterior-sample"}}{Average predictions over posterior draws.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL}, uses all available draws.
#' @param seed Optional random seed for reproducibility when
#'   \code{method = "posterior-sample"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' \code{num_categories + 1} columns containing \eqn{P(X_j = c | X_{-j})}{P(X_j = c | X_-j)} for each
#' observation and category.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories.
#'
#' When \code{method = "posterior-sample"}, probabilities are averaged over
#' posterior draws, and an attribute \code{"sd"} is included containing the
#' standard deviation across draws.
#'
#' @details
#' For each observation, the function computes the conditional distribution
#' of the target variable(s) given the observed values of all other variables.
#' This is the same conditional distribution used internally by the Gibbs
#' sampler.
#'
#' @seealso \code{\link{simulate.bgms}} for generating new data from the model.
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit <- bgm(x = Wenchuan[, 1:5])
#'
#' # Compute conditional probabilities for all variables
#' probs <- predict(fit, newdata = Wenchuan[1:10, 1:5])
#'
#' # Predict the first variable only
#' probs_v1 <- predict(fit, newdata = Wenchuan[1:10, 1:5], variables = 1)
#'
#' # Get predicted categories
#' pred_class <- predict(fit, newdata = Wenchuan[1:10, 1:5], type = "response")
#' }
#'
#' @importFrom stats predict
#' @export
predict.bgms <- function(object,
                         newdata,
                         variables = NULL,
                         type = c("probabilities", "response"),
                         method = c("posterior-mean", "posterior-sample"),
                         ndraws = NULL,
                         seed = NULL,
                         ...) {
  type <- match.arg(type)
  method <- match.arg(method)

  # Setting the seed (for R's RNG used by sample.int for draw selection)
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
      stop("Argument 'seed' must be a single non-negative integer.")
    }
    seed <- as.integer(seed)
    set.seed(seed)
  }

  # Validate newdata
  if(missing(newdata)) {
    stop("Argument 'newdata' is required. Provide the data for which to compute predictions.")
  }

  if(!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
    stop("'newdata' must be a matrix or data frame.")
  }

  if(inherits(newdata, "data.frame")) {
    newdata <- data.matrix(newdata)
  }

  # Extract model information
  arguments <- extract_arguments(object)
  num_variables <- arguments$num_variables
  num_categories <- arguments$num_categories
  variable_type <- arguments$variable_type
  data_columnnames <- arguments$data_columnnames

  # Validate dimensions

  if(ncol(newdata) != num_variables) {
    stop(paste0(
      "'newdata' must have ", num_variables, " columns (same as fitted model), ",
      "but has ", ncol(newdata), "."
    ))
  }

  # Handle variable_type
  if(length(variable_type) == 1) {
    variable_type <- rep(variable_type, num_variables)
  }

  # Get baseline_category
  baseline_category <- arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category <- rep(0L, num_variables)
  }

  # Convert variable_type to is_ordinal logical vector
  is_ordinal <- variable_type != "blume-capel"

  # Determine which variables to predict
  if(is.null(variables)) {
    predict_vars <- seq_len(num_variables)
  } else if(is.character(variables)) {
    predict_vars <- match(variables, data_columnnames)
    if(anyNA(predict_vars)) {
      stop("Variable names not found: ", paste(variables[is.na(predict_vars)], collapse = ", "))
    }
  } else {
    predict_vars <- as.integer(variables)
    if(any(predict_vars < 1 | predict_vars > num_variables)) {
      stop("Variable indices must be between 1 and ", num_variables)
    }
  }

  # Recode data to 0-based integers (matching what bgm() does)
  newdata_recoded <- recode_data_for_prediction(newdata, num_categories, is_ordinal)

  if(method == "posterior-mean") {
    # Use posterior mean parameters
    interactions <- object$posterior_mean_pairwise
    thresholds <- object$posterior_mean_main

    probs <- compute_conditional_probs(
      observations = newdata_recoded,
      predict_vars = predict_vars - 1L,  # C++ uses 0-based indexing
      interactions = interactions,
      thresholds = thresholds,
      no_categories = num_categories,
      variable_type = variable_type,
      baseline_category = baseline_category
    )

    # Add names
    names(probs) <- data_columnnames[predict_vars]
    for(v in seq_along(probs)) {
      var_idx <- predict_vars[v]
      n_cats <- num_categories[var_idx] + 1
      colnames(probs[[v]]) <- paste0("cat_", 0:(n_cats - 1))
    }

  } else {
    # Use posterior samples
    pairwise_samples <- do.call(rbind, object$raw_samples$pairwise)
    main_samples <- do.call(rbind, object$raw_samples$main)

    total_draws <- nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws <- total_draws
    }
    ndraws <- min(ndraws, total_draws)

    draw_indices <- sample.int(total_draws, ndraws)

    # Collect probabilities from each draw
    all_probs <- vector("list", ndraws)

    for(i in seq_len(ndraws)) {
      idx <- draw_indices[i]

      # Reconstruct interaction matrix
      interactions <- matrix(0, nrow = num_variables, ncol = num_variables)
      interactions[lower.tri(interactions)] <- pairwise_samples[idx, ]
      interactions <- interactions + t(interactions)

      # Reconstruct threshold matrix
      thresholds <- reconstruct_thresholds(
        main_samples[idx, ],
        num_variables,
        num_categories,
        variable_type
      )

      all_probs[[i]] <- compute_conditional_probs(
        observations = newdata_recoded,
        predict_vars = predict_vars - 1L,
        interactions = interactions,
        thresholds = thresholds,
        no_categories = num_categories,
        variable_type = variable_type,
        baseline_category = baseline_category
      )
    }

    # Average over draws
    probs <- vector("list", length(predict_vars))
    probs_sd <- vector("list", length(predict_vars))
    names(probs) <- data_columnnames[predict_vars]
    names(probs_sd) <- data_columnnames[predict_vars]

    for(v in seq_along(predict_vars)) {
      # Stack probabilities from all draws: n x categories x ndraws
      var_probs <- lapply(all_probs, `[[`, v)
      prob_array <- array(unlist(var_probs),
                          dim = c(nrow(newdata), ncol(var_probs[[1]]), ndraws))

      probs[[v]] <- apply(prob_array, c(1, 2), mean)
      probs_sd[[v]] <- apply(prob_array, c(1, 2), sd)

      var_idx <- predict_vars[v]
      n_cats <- num_categories[var_idx] + 1
      colnames(probs[[v]]) <- paste0("cat_", 0:(n_cats - 1))
      colnames(probs_sd[[v]]) <- paste0("cat_", 0:(n_cats - 1))
    }

    attr(probs, "sd") <- probs_sd
  }

  if(type == "response") {
    # Return predicted categories (mode)
    pred_matrix <- sapply(probs, function(p) {
      apply(p, 1, which.max) - 1L  # Convert to 0-based category
    })
    if(is.vector(pred_matrix)) {
      pred_matrix <- matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) <- data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(probs)
}


# ==============================================================================
#   Helper Functions
# ==============================================================================

# Helper function to reconstruct threshold matrix from flat vector
reconstruct_thresholds <- function(main_vec, num_variables, num_categories, variable_type) {
  if(length(variable_type) == 1) {
    variable_type <- rep(variable_type, num_variables)
  }

  max_cats <- max(num_categories)
  thresholds <- matrix(NA, nrow = num_variables, ncol = max_cats)

  pos <- 1
  for(v in seq_len(num_variables)) {
    if(variable_type[v] != "blume-capel") {
      k <- num_categories[v]
      thresholds[v, 1:k] <- main_vec[pos:(pos + k - 1)]
      pos <- pos + k
    } else {
      thresholds[v, 1:2] <- main_vec[pos:(pos + 1)]
      pos <- pos + 2
    }
  }

  return(thresholds)
}


# Helper function to recode data for prediction
recode_data_for_prediction <- function(x, num_categories, is_ordinal) {
  x <- as.matrix(x)
  num_variables <- ncol(x)

  for(v in seq_len(num_variables)) {
    if(is_ordinal[v]) {
      # For ordinal variables, ensure values are in 0:num_categories[v]
      x[, v] <- as.integer(x[, v])
      if(min(x[, v], na.rm = TRUE) > 0) {
        # Shift to 0-based if necessary
        x[, v] <- x[, v] - min(x[, v], na.rm = TRUE)
      }
    }
  }

  return(x)
}
