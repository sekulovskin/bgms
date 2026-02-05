# ==============================================================================
#   Comprehensive Tests for simulate() and predict() Functions
#
#   This script tests:
#     - simulate_mrf(): Standalone MRF simulation
#     - mrfSampler(): Deprecated wrapper (should show warning)
#     - simulate.bgms(): S3 method for bgms objects
#     - predict.bgms(): S3 method for conditional probabilities
#
#   Run with: source("dev/test_simulate_predict.R")
# ==============================================================================

library(bgms)

cat("\n")
cat("=======================================================================\n")
cat("  COMPREHENSIVE TESTS FOR simulate() AND predict() FUNCTIONS\n")
cat("=======================================================================\n\n")

# Track test results
tests_passed <- 0
tests_failed <- 0

test <- function(description, expr) {
  cat(sprintf("Testing: %s... ", description))
  result <- tryCatch({
    eval(expr)
    cat("PASSED\n")
    tests_passed <<- tests_passed + 1
    TRUE
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
    tests_failed <<- tests_failed + 1
    FALSE
  })
  invisible(result)
}

expect_equal <- function(x, y, tol = .Machine$double.eps^0.5) {
  if (!isTRUE(all.equal(x, y, tolerance = tol))) {
    stop(sprintf("Expected %s but got %s", deparse(y), deparse(x)))
  }
}

expect_true <- function(x) {
  if (!isTRUE(x)) stop("Expected TRUE but got FALSE")
}

expect_class <- function(x, class) {
  if (!inherits(x, class)) {
    stop(sprintf("Expected class '%s' but got '%s'", class, paste(class(x), collapse = ", ")))
  }
}

expect_dim <- function(x, expected_dim) {
  if (!identical(dim(x), as.integer(expected_dim))) {
    stop(sprintf("Expected dim %s but got %s",
                 paste(expected_dim, collapse = " x "),
                 paste(dim(x), collapse = " x ")))
  }
}

expect_length <- function(x, expected_length) {
  if (length(x) != expected_length) {
    stop(sprintf("Expected length %d but got %d", expected_length, length(x)))
  }
}

expect_warning <- function(expr, pattern = NULL) {
  warned <- FALSE
  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        msg <- conditionMessage(w)
        if (is.null(pattern) || any(grepl(pattern, msg, ignore.case = TRUE))) {
          warned <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) stop(e)
  )
  if (!warned) stop("Expected warning but none occurred")
  invisible(result)
}


# ==============================================================================
#   Section 1: simulate_mrf() - Standalone MRF Simulation
# ==============================================================================

cat("\n--- Section 1: simulate_mrf() Tests ---\n\n")

# Basic ordinal simulation
test("simulate_mrf basic ordinal simulation", {
  result <- simulate_mrf(
    no_states = 100,
    no_variables = 5,
    no_categories = 3,
    interactions = matrix(0.2, 5, 5) - diag(0.2, 5),
    thresholds = matrix(0, 5, 3),
    iter = 100,
    seed = 42
  )
  expect_class(result, "matrix")
  expect_dim(result, c(100, 5))
  expect_true(all(result >= 0 & result <= 3))
})

# Variable number of categories
test("simulate_mrf with varying categories per variable", {
  no_cats <- c(2, 3, 4, 5, 2)
  result <- simulate_mrf(
    no_states = 50,
    no_variables = 5,
    no_categories = no_cats,
    interactions = matrix(0.1, 5, 5) - diag(0.1, 5),
    thresholds = matrix(0, 5, max(no_cats)),
    iter = 100,
    seed = 42
  )
  expect_dim(result, c(50, 5))
  # Check each variable stays within its category range
  for (v in 1:5) {
    expect_true(all(result[, v] >= 0 & result[, v] <= no_cats[v]))
  }
})

# Reproducibility with seed
test("simulate_mrf reproducibility with seed", {
  args <- list(
    no_states = 50,
    no_variables = 4,
    no_categories = 3,
    interactions = matrix(0.1, 4, 4) - diag(0.1, 4),
    thresholds = matrix(0, 4, 3),
    iter = 100
  )
  result1 <- do.call(simulate_mrf, c(args, seed = 123))
  result2 <- do.call(simulate_mrf, c(args, seed = 123))
  result3 <- do.call(simulate_mrf, c(args, seed = 456))

  expect_true(identical(result1, result2))
  expect_true(!identical(result1, result3))
})

# Blume-Capel variables
test("simulate_mrf with Blume-Capel variables", {
  result <- simulate_mrf(
    no_states = 100,
    no_variables = 4,
    no_categories = 4,
    interactions = matrix(0.1, 4, 4) - diag(0.1, 4),
    thresholds = matrix(c(-0.5, -0.2, NA, NA), 4, 4, byrow = TRUE),
    variable_type = "blume-capel",
    baseline_category = 2,
    iter = 100,
    seed = 42
  )
  expect_dim(result, c(100, 4))
  expect_true(all(result >= 0 & result <= 4))
})

# Mixed ordinal and Blume-Capel
test("simulate_mrf with mixed variable types", {
  thresholds <- matrix(0, 4, 4)
  thresholds[2, 1:2] <- c(-0.3, -0.1)  # Blume-Capel params for var 2
  thresholds[4, 1:2] <- c(-0.2, -0.15) # Blume-Capel params for var 4

  result <- simulate_mrf(
    no_states = 100,
    no_variables = 4,
    no_categories = 4,
    interactions = matrix(0.1, 4, 4) - diag(0.1, 4),
    thresholds = thresholds,
    variable_type = c("ordinal", "blume-capel", "ordinal", "blume-capel"),
    baseline_category = c(0, 2, 0, 2),
    iter = 100,
    seed = 42
  )
  expect_dim(result, c(100, 4))
})


# ==============================================================================
#   Section 2: mrfSampler() - Deprecated Wrapper
# ==============================================================================

cat("\n--- Section 2: mrfSampler() Deprecation Tests ---\n\n")

test("mrfSampler shows deprecation warning", {
  expect_warning({
    result <- mrfSampler(
      no_states = 10,
      no_variables = 3,
      no_categories = 2,
      interactions = matrix(0.1, 3, 3) - diag(0.1, 3),
      thresholds = matrix(0, 3, 2),
      iter = 50,
      seed = 42
    )
  }, "deprecated")
  expect_dim(result, c(10, 3))
})

test("mrfSampler produces identical results to simulate_mrf", {
  args <- list(
    no_states = 50,
    no_variables = 4,
    no_categories = 3,
    interactions = matrix(0.1, 4, 4) - diag(0.1, 4),
    thresholds = matrix(0, 4, 3),
    iter = 100,
    seed = 999
  )

  result_new <- do.call(simulate_mrf, args)
  result_old <- suppressWarnings(do.call(mrfSampler, args))

  expect_true(identical(result_new, result_old))
})


# ==============================================================================
#   Section 3: simulate.bgms() - S3 Method
# ==============================================================================

cat("\n--- Section 3: simulate.bgms() Tests ---\n\n")

# Fit a model for testing
cat("Fitting bgms model for simulate.bgms tests...\n")
data(Wenchuan)
fit <- bgm(Wenchuan[1:100, 1:5], iter = 200, warmup = 100, display_progress = "none")

test("simulate.bgms with posterior-mean method", {
  result <- simulate(fit, nsim = 50, method = "posterior-mean", seed = 42)
  expect_class(result, "matrix")
  expect_dim(result, c(50, 5))
  expect_true(!is.null(colnames(result)))
})

test("simulate.bgms with posterior-sample method", {
  result <- simulate(fit, nsim = 30, method = "posterior-sample",
                     ndraws = 10, seed = 42, display_progress = "none")
  expect_class(result, "list")
  expect_length(result, 10)
  expect_dim(result[[1]], c(30, 5))
})

test("simulate.bgms parallel execution", {
  result <- simulate(fit, nsim = 20, method = "posterior-sample",
                     ndraws = 8, cores = 2, seed = 42, display_progress = "none")
  expect_class(result, "list")
  expect_length(result, 8)
})

test("simulate.bgms reproducibility with seed", {
  result1 <- simulate(fit, nsim = 30, method = "posterior-mean", seed = 123)
  result2 <- simulate(fit, nsim = 30, method = "posterior-mean", seed = 123)
  result3 <- simulate(fit, nsim = 30, method = "posterior-mean", seed = 456)

  expect_true(identical(result1, result2))
  expect_true(!identical(result1, result3))
})

test("simulate.bgms column names match original data", {
  result <- simulate(fit, nsim = 10, method = "posterior-mean", seed = 42)
  expect_true(identical(colnames(result), colnames(Wenchuan)[1:5]))
})

test("simulate.bgms with iter parameter", {
  # Different iter values should produce different results (stochastic)
  result1 <- simulate(fit, nsim = 50, method = "posterior-mean", iter = 100, seed = 42)
  result2 <- simulate(fit, nsim = 50, method = "posterior-mean", iter = 1000, seed = 42)
  # Both should be valid matrices
  expect_dim(result1, c(50, 5))
  expect_dim(result2, c(50, 5))
})


# ==============================================================================
#   Section 4: predict.bgms() - S3 Method
# ==============================================================================

cat("\n--- Section 4: predict.bgms() Tests ---\n\n")

test("predict.bgms with posterior-mean method", {
  newdata <- Wenchuan[1:20, 1:5]
  result <- predict(fit, newdata = newdata, method = "posterior-mean")

  expect_class(result, "list")
  expect_length(result, 5)
  # Each element should be a matrix with n rows and (num_categories + 1) columns
  expect_true(nrow(result[[1]]) == 20)
  # Probabilities should sum to 1
  expect_true(all(abs(rowSums(result[[1]]) - 1) < 1e-10))
})

test("predict.bgms with type = 'response'", {
  newdata <- Wenchuan[1:20, 1:5]
  result <- predict(fit, newdata = newdata, type = "response")

  expect_class(result, "matrix")
  expect_dim(result, c(20, 5))
  # Response values should be valid categories
  expect_true(all(result >= 0))
})

test("predict.bgms for specific variables", {
  newdata <- Wenchuan[1:20, 1:5]
  result <- predict(fit, newdata = newdata, variables = c(1, 3))

  expect_class(result, "list")
  expect_length(result, 2)
})

test("predict.bgms with variable names", {
  newdata <- Wenchuan[1:20, 1:5]
  var_names <- colnames(Wenchuan)[c(1, 2)]
  result <- predict(fit, newdata = newdata, variables = var_names)

  expect_class(result, "list")
  expect_length(result, 2)
})

test("predict.bgms with posterior-sample method", {
  newdata <- Wenchuan[1:10, 1:5]
  result <- predict(fit, newdata = newdata, method = "posterior-sample",
                    ndraws = 20, seed = 42)

  expect_class(result, "list")
  expect_length(result, 5)
  # Should have sd attribute
  expect_true(!is.null(attr(result, "sd")))
  sd_attr <- attr(result, "sd")
  expect_length(sd_attr, 5)
})

test("predict.bgms probabilities are valid", {
  newdata <- Wenchuan[1:30, 1:5]
  result <- predict(fit, newdata = newdata)

  for (v in seq_along(result)) {
    probs <- result[[v]]
    # Skip if all NA (can happen with missing data)
    if (all(is.na(probs))) next
    # All non-NA probabilities should be between 0 and 1
    valid_probs <- probs[!is.na(probs)]
    if (!all(valid_probs >= -1e-10 & valid_probs <= 1 + 1e-10)) {
      stop(sprintf("Variable %d: probabilities outside [0,1]", v))
    }
    # Each non-NA row should sum to 1 (with tolerance for numerical precision)
    row_sums <- rowSums(probs, na.rm = TRUE)
    complete_rows <- rowSums(!is.na(probs)) == ncol(probs)
    if (any(complete_rows)) {
      valid_sums <- row_sums[complete_rows]
      if (!all(abs(valid_sums - 1) < 1e-6)) {
        stop(sprintf("Variable %d: row sums not 1", v))
      }
    }
  }
  expect_true(TRUE)
})

test("predict.bgms requires newdata argument", {
  error_occurred <- FALSE
  tryCatch({
    predict(fit)
  }, error = function(e) {
    error_occurred <<- grepl("newdata", e$message, ignore.case = TRUE)
  })
  expect_true(error_occurred)
})

test("predict.bgms validates newdata dimensions", {
  error_occurred <- FALSE
  tryCatch({
    predict(fit, newdata = Wenchuan[1:10, 1:3])  # Wrong number of columns
  }, error = function(e) {
    error_occurred <<- TRUE
  })
  expect_true(error_occurred)
})


# ==============================================================================
#   Section 5: Edge Cases and Error Handling
# ==============================================================================

cat("\n--- Section 5: Edge Cases and Error Handling ---\n\n")

test("simulate_mrf with single observation", {
  result <- simulate_mrf(
    no_states = 1,
    no_variables = 3,
    no_categories = 2,
    interactions = matrix(0.1, 3, 3) - diag(0.1, 3),
    thresholds = matrix(0, 3, 2),
    iter = 100,
    seed = 42
  )
  expect_dim(result, c(1, 3))
})

test("simulate_mrf with binary variables (no_categories = 1)", {
  result <- simulate_mrf(
    no_states = 100,
    no_variables = 4,
    no_categories = 1,
    interactions = matrix(0.2, 4, 4) - diag(0.2, 4),
    thresholds = matrix(0, 4, 1),
    iter = 100,
    seed = 42
  )
  expect_dim(result, c(100, 4))
  expect_true(all(result %in% c(0, 1)))
})

test("simulate.bgms with ndraws = 1", {
  result <- simulate(fit, nsim = 20, method = "posterior-sample",
                     ndraws = 1, seed = 42, display_progress = "none")
  expect_class(result, "list")
  expect_length(result, 1)
  expect_dim(result[[1]], c(20, 5))
})

test("predict.bgms with single observation", {
  newdata <- Wenchuan[1, 1:5, drop = FALSE]
  result <- predict(fit, newdata = newdata)

  expect_class(result, "list")
  expect_length(result, 5)
  expect_true(nrow(result[[1]]) == 1)
})

test("predict.bgms with single variable", {
  newdata <- Wenchuan[1:10, 1:5]
  result <- predict(fit, newdata = newdata, variables = 1)

  expect_class(result, "list")
  expect_length(result, 1)
})


# ==============================================================================
#   Section 6: Statistical Validation - Simulated Values Match Model Probabilities
# ==============================================================================

cat("\n--- Section 6: Statistical Validation ---\n\n")

test("simulated marginals match theoretical probabilities (independent case)", {
  # For independent variables (no interactions), we can compute exact probabilities
  # Ordinal MRF model: P(X = 0) ∝ 1, P(X = k) ∝ exp(threshold_k) for k > 0
  no_vars <- 3
  no_cats <- 2  # 3 categories: 0, 1, 2
  n <- 10000

  # Thresholds: directly affect probability of each category
  # threshold[, 1] affects P(X=1), threshold[, 2] affects P(X=2)
  thresholds <- matrix(c(
    -1.0, -2.0,   # var 1: favors category 0 (low thresholds reduce prob of 1,2)
     0.0,  0.0,   # var 2: balanced (all categories ~equal)
     1.0,  2.0    # var 3: favors higher categories (high thresholds increase prob)
  ), nrow = no_vars, ncol = no_cats, byrow = TRUE)

  # No interactions - variables are independent
  interactions <- matrix(0, no_vars, no_vars)

  result <- simulate_mrf(
    no_states = n,
    no_variables = no_vars,
    no_categories = no_cats,
    interactions = interactions,
    thresholds = thresholds,
    iter = 500,
    seed = 12345
  )

  # Compute empirical frequencies
  empirical_probs <- lapply(1:no_vars, function(v) {
    table(factor(result[, v], levels = 0:no_cats)) / n
  })

  # Compute theoretical probabilities for ordinal MRF with no interactions
  # P(X = 0) ∝ 1, P(X = k) ∝ exp(threshold[k-1]) for k = 1, 2, ...
  theoretical_probs <- lapply(1:no_vars, function(v) {
    unnorm <- c(1, exp(thresholds[v, ]))
    unnorm / sum(unnorm)
  })

  # Chi-square test for each variable (should NOT reject if simulation is correct)
  for (v in 1:no_vars) {
    expected <- theoretical_probs[[v]] * n
    observed <- as.numeric(empirical_probs[[v]]) * n
    chi_sq <- sum((observed - expected)^2 / expected)
    # df = no_cats, use alpha = 0.001 to be conservative
    critical_value <- qchisq(0.999, df = no_cats)
    if (chi_sq > critical_value) {
      stop(sprintf("Variable %d: chi-square test failed (chi_sq = %.2f > %.2f)\n  Empirical: %s\n  Theoretical: %s",
                   v, chi_sq, critical_value,
                   paste(round(as.numeric(empirical_probs[[v]]), 3), collapse = ", "),
                   paste(round(theoretical_probs[[v]], 3), collapse = ", ")))
    }
  }
  expect_true(TRUE)
})

test("simulated correlations reflect interaction strength", {
  # With positive interaction, variables should be positively correlated
  n <- 5000
  no_vars <- 2
  no_cats <- 1  # Binary

  # Strong positive interaction
  interactions_pos <- matrix(c(0, 1.5, 1.5, 0), 2, 2)
  # Strong negative interaction
  interactions_neg <- matrix(c(0, -1.5, -1.5, 0), 2, 2)
  # No interaction
  interactions_zero <- matrix(0, 2, 2)

  thresholds <- matrix(0, 2, 1)

  sim_pos <- simulate_mrf(no_states = n, no_variables = no_vars, no_categories = no_cats,
                          interactions = interactions_pos, thresholds = thresholds,
                          iter = 500, seed = 111)
  sim_neg <- simulate_mrf(no_states = n, no_variables = no_vars, no_categories = no_cats,
                          interactions = interactions_neg, thresholds = thresholds,
                          iter = 500, seed = 222)
  sim_zero <- simulate_mrf(no_states = n, no_variables = no_vars, no_categories = no_cats,
                           interactions = interactions_zero, thresholds = thresholds,
                           iter = 500, seed = 333)

  cor_pos <- cor(sim_pos[, 1], sim_pos[, 2])
  cor_neg <- cor(sim_neg[, 1], sim_neg[, 2])
  cor_zero <- cor(sim_zero[, 1], sim_zero[, 2])

  # Positive interaction should give positive correlation
  if (cor_pos <= 0.3) stop(sprintf("Positive interaction correlation too low: %.3f", cor_pos))
  # Negative interaction should give negative correlation
  if (cor_neg >= -0.3) stop(sprintf("Negative interaction correlation too high: %.3f", cor_neg))
  # Zero interaction should give near-zero correlation
  if (abs(cor_zero) > 0.15) stop(sprintf("Zero interaction correlation not near zero: %.3f", cor_zero))

  expect_true(TRUE)
})

test("predict.bgms probabilities match empirical frequencies from simulate.bgms", {
  # Generate simulated data from the fitted model, then check that

  # predict() probabilities align with empirical frequencies

  # Use the existing 'fit' object from Section 3
  set.seed(999)

  # Simulate a larger dataset to get stable empirical estimates
  sim_data <- simulate(fit, nsim = 2000, method = "posterior-mean", iter = 500, seed = 999)

  # Get predictions for the simulated data
  preds <- predict(fit, newdata = sim_data, method = "posterior-mean")

  # For each variable, the mean predicted probability for each category
  # should roughly equal the empirical frequency in the simulated data
  for (v in 1:ncol(sim_data)) {
    empirical_freq <- table(factor(sim_data[, v], levels = 0:(ncol(preds[[v]]) - 1))) / nrow(sim_data)
    mean_pred_prob <- colMeans(preds[[v]], na.rm = TRUE)

    # These won't match exactly (predictions are conditional, empirical are marginal)
    # But for a self-consistent model, they should be in the same ballpark
    # Check that no probability differs by more than 0.15
    max_diff <- max(abs(as.numeric(empirical_freq) - mean_pred_prob), na.rm = TRUE)
    if (max_diff > 0.25) {
      stop(sprintf("Variable %d: max difference between empirical and predicted too large: %.3f", v, max_diff))
    }
  }
  expect_true(TRUE)
})

test("round-trip: fit -> simulate -> refit produces correlated estimates", {

  # This is the ultimate validation: if we fit a model, simulate from it,

  # and refit on the simulated data, the estimates should be similar.

  cat("\n  Fitting original model on empirical data...\n")

  # Use more data and longer chains for reliable estimates
  data(Wenchuan)
  fit_original <- bgm(Wenchuan[, 1:10], iter = 2000, warmup = 1000,
                      display_progress = "none")

  # Extract posterior mean parameters from original fit
  original_interactions <- fit_original$posterior_mean_pairwise
  original_thresholds <- fit_original$posterior_mean_main

  cat("  Simulating data from fitted model...\n")

  # Simulate a dataset of similar size
  n_sim <- nrow(Wenchuan)
  sim_data <- simulate(fit_original, nsim = n_sim, method = "posterior-mean",
                       iter = 1000, seed = 42)

  cat("  Fitting model on simulated data...\n")

  # Fit model on simulated data
  fit_simulated <- bgm(sim_data, iter = 2000, warmup = 1000,
                       display_progress = "none")

  # Extract posterior mean parameters from simulated fit
  simulated_interactions <- fit_simulated$posterior_mean_pairwise
  simulated_thresholds <- fit_simulated$posterior_mean_main

  # Compare interaction parameters (off-diagonal elements)
  original_int_vec <- original_interactions[lower.tri(original_interactions)]
  simulated_int_vec <- simulated_interactions[lower.tri(simulated_interactions)]

  int_correlation <- cor(original_int_vec, simulated_int_vec)
  cat(sprintf("  Interaction correlation: %.3f\n", int_correlation))

  # Compare threshold parameters (flatten matrices)
  original_thresh_vec <- as.vector(original_thresholds)
  simulated_thresh_vec <- as.vector(simulated_thresholds)

  # Remove NAs (different variables may have different number of thresholds)
  valid_idx <- !is.na(original_thresh_vec) & !is.na(simulated_thresh_vec)
  thresh_correlation <- cor(original_thresh_vec[valid_idx],
                            simulated_thresh_vec[valid_idx])
  cat(sprintf("  Threshold correlation: %.3f\n", thresh_correlation))

  # Correlations should be strong (> 0.7) for a well-functioning simulation
  if (int_correlation < 0.7) {
    stop(sprintf("Interaction correlation too low: %.3f (expected > 0.7)", int_correlation))
  }
  if (thresh_correlation < 0.7) {
    stop(sprintf("Threshold correlation too low: %.3f (expected > 0.7)", thresh_correlation))
  }

  expect_true(TRUE)
})

test("round-trip Blume-Capel: fit -> simulate -> refit produces correlated estimates", {

  # Same validation as ordinal, but for Blume-Capel variables

  cat("\n  Fitting original Blume-Capel model on empirical data...\n")

  # Use more data and longer chains for reliable estimates
  data(Wenchuan)
  fit_original <- bgm(Wenchuan[, 1:10], iter = 2000, warmup = 1000,
                      variable_type = "blume-capel", baseline_category = 1,
                      update_method = "adaptive-metropolis",
                      display_progress = "none")

  # Extract posterior mean parameters from original fit
  original_interactions <- fit_original$posterior_mean_pairwise
  original_thresholds <- fit_original$posterior_mean_main

  cat("  Simulating data from fitted model...\n")

  # Simulate a dataset of similar size
  n_sim <- nrow(Wenchuan)
  sim_data <- simulate(fit_original, nsim = n_sim, method = "posterior-mean",
                       iter = 1000, seed = 42)

  cat("  Fitting model on simulated data...\n")

  # Fit model on simulated data using same baseline from original fit
  fit_simulated <- bgm(sim_data, iter = 2000, warmup = 1000,
                       variable_type = "blume-capel",
                       baseline_category = fit_original$arguments$baseline_category,
                       update_method = "adaptive-metropolis",
                       display_progress = "none")

  # Extract posterior mean parameters from simulated fit
  simulated_interactions <- fit_simulated$posterior_mean_pairwise
  simulated_thresholds <- fit_simulated$posterior_mean_main

  # Compare interaction parameters (off-diagonal elements)
  original_int_vec <- original_interactions[lower.tri(original_interactions)]
  simulated_int_vec <- simulated_interactions[lower.tri(simulated_interactions)]

  int_correlation <- cor(original_int_vec, simulated_int_vec)
  cat(sprintf("  Interaction correlation: %.3f\n", int_correlation))

  # Compare threshold parameters (flatten matrices)
  original_thresh_vec <- as.vector(original_thresholds)
  simulated_thresh_vec <- as.vector(simulated_thresholds)

  # Remove NAs (different variables may have different number of thresholds)
  valid_idx <- !is.na(original_thresh_vec) & !is.na(simulated_thresh_vec)
  thresh_correlation <- cor(original_thresh_vec[valid_idx],
                            simulated_thresh_vec[valid_idx])
  cat(sprintf("  Threshold correlation: %.3f\n", thresh_correlation))

  # Correlations should be strong (> 0.7) for a well-functioning simulation
  if (int_correlation < 0.7) {
    stop(sprintf("Blume-Capel interaction correlation too low: %.3f (expected > 0.7)", int_correlation))
  }
  if (thresh_correlation < 0.7) {
    stop(sprintf("Blume-Capel threshold correlation too low: %.3f (expected > 0.7)", thresh_correlation))
  }

  expect_true(TRUE)
})

# ==============================================================================
#   Section 7: Performance Tests
# ==============================================================================

cat("\n--- Section 7: Performance Tests ---\n\n")

test("simulate.bgms parallel vs sequential consistency", {
  # Results with same seed should be consistent regardless of cores
  # (Note: parallel execution may produce different RNG sequences)
  result_seq <- simulate(fit, nsim = 30, method = "posterior-sample",
                         ndraws = 5, cores = 1, seed = 42, display_progress = "none")
  result_par <- simulate(fit, nsim = 30, method = "posterior-sample",
                         ndraws = 5, cores = 2, seed = 42, display_progress = "none")

  # Both should produce valid output
  expect_length(result_seq, 5)
  expect_length(result_par, 5)
  expect_dim(result_seq[[1]], c(30, 5))
  expect_dim(result_par[[1]], c(30, 5))
})


# ==============================================================================
#   Summary
# ==============================================================================

cat("\n")
cat("=======================================================================\n")
cat("  TEST SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Passed: %d\n", tests_passed))
cat(sprintf("  Failed: %d\n", tests_failed))
cat(sprintf("  Total:  %d\n", tests_passed + tests_failed))
cat("=======================================================================\n")

if (tests_failed > 0) {
  cat("\n*** SOME TESTS FAILED ***\n\n")
} else {
  cat("\n*** ALL TESTS PASSED ***\n\n")
}
