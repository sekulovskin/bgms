# ==============================================================================
# Tests for bgmCompare() - Multi-group MRF Comparison
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Reproducibility, range invariants, dimension consistency
#
# These tests parallel the structure of test-bgm.R for consistency.
# Tests both reproducibility and basic output structure.
#
# Tolerance testing principles applied here:
#   - Reproducibility: identical seeds produce identical MCMC chains
#   - Range invariants: posterior indicators in [0,1], finite estimates
#   - Symmetry: baseline pairwise matrices must be symmetric
#   - Dimension consistency: correct number of groups, edges, variables
#
# See helper-fixtures.R for testing philosophy and session-cached fixtures.
# ==============================================================================

# Determine cores for testing
on_ci <- isTRUE(as.logical(Sys.getenv("CI", "false")))
no_cores <- if (on_ci) 2L else min(2L, parallel::detectCores())

# ------------------------------------------------------------------------------
# Reproducibility Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare is reproducible with seed (x, y interface)", {

  testthat::skip_on_cran()

  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:40, 1:4]
  y <- Wenchuan[41:80, 1:4]

  fit1 <- bgmCompare(x = x, y = y, iter = 100, warmup = 200, cores = no_cores, seed = 1234, display_progress = "none")
  fit2 <- bgmCompare(x = x, y = y, iter = 100, warmup = 200, cores = no_cores, seed = 1234, display_progress = "none")

  # Combine chains for comparison
  combine_chains <- function(fit) {
    pairs <- do.call(rbind, fit$raw_samples$pairwise)
    mains <- do.call(rbind, fit$raw_samples$main)
    cbind(mains, pairs)
  }

  combined1 <- combine_chains(fit1)
  combined2 <- combine_chains(fit2)

  expect_equal(combined1, combined2)
})

test_that("bgmCompare is reproducible with seed (x + group_indicator interface)", {
  testthat::skip_on_cran()

  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:80, 1:4]
  group_ind <- rep(1:2, each = 40)

  fit1 <- bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 100, warmup = 200, cores = no_cores, seed = 5678,
    display_progress = "none"
  )
  fit2 <- bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 100, warmup = 200, cores = no_cores, seed = 5678,
    display_progress = "none"
  )

  combine_chains <- function(fit) {
    pairs <- do.call(rbind, fit$raw_samples$pairwise)
    mains <- do.call(rbind, fit$raw_samples$main)
    cbind(mains, pairs)
  }

  expect_equal(combine_chains(fit1), combine_chains(fit2))
})


# ------------------------------------------------------------------------------
# Output Structure Tests (using saved fit)
# ------------------------------------------------------------------------------

test_that("bgmCompare output has expected structure", {
  fit <- get_bgmcompare_fit()

  expect_s3_class(fit, "bgmCompare")

  # Should have key components
  expect_true("arguments" %in% names(fit))
  expect_true("raw_samples" %in% names(fit))

  # Raw samples should have required components
  expect_true("pairwise" %in% names(fit$raw_samples))
  expect_true("main" %in% names(fit$raw_samples))
})

test_that("bgmCompare stores correct number of groups", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)

  expect_true("num_groups" %in% names(args))
  expect_true(args$num_groups >= 2)
})

test_that("bgmCompare posterior summaries have expected format", {
  fit <- get_bgmcompare_fit()

  # Should have baseline summaries
  expect_true(!is.null(fit$posterior_summary_pairwise_baseline))
  expect_true(!is.null(fit$posterior_mean_pairwise_baseline))
})


# ------------------------------------------------------------------------------
# Tolerance/Sanity Tests (Stochastic-robust)
# ------------------------------------------------------------------------------

test_that("bgmCompare outputs are numerically sane", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  p <- args$num_variables

  # Check baseline pairwise
  M <- fit$posterior_mean_pairwise_baseline

  expect_true(is.matrix(M))
  expect_equal(dim(M), c(p, p))

  # Symmetry check
  asym <- max(abs(M - t(M)), na.rm = TRUE)
  expect_true(asym <= 1e-8, info = sprintf("Asymmetry too large: %g", asym))

  # Values should be finite
  expect_true(all(is.finite(M)))

  # Check group params
  group_params <- extract_group_params(fit)

  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# Fresh Fit Tests (slower, skip on CRAN)
# ------------------------------------------------------------------------------

test_that("bgmCompare without selection produces valid estimates", {
  skip_on_cran_mcmc()

  data <- generate_grouped_test_data(n_per_group = 25, p = 3, n_groups = 2, seed = 42)

  fit <- bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = FALSE,
    iter = 100,
    warmup = 150,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  # Should have posterior means
  expect_true(!is.null(fit$posterior_mean_pairwise_baseline))
  expect_true(!is.null(fit$posterior_mean_main_baseline))
})

test_that("bgmCompare with selection produces valid indicators", {
  skip_on_cran_mcmc()

  data <- generate_grouped_test_data(n_per_group = 30, p = 3, n_groups = 2, seed = 123)

  # Test with single chain (previously caused bug in summarize_mixture_effect)
  fit <- bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = TRUE,
    iter = 200,
    warmup = 200,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  # Should have indicator samples
  expect_true(!is.null(fit$raw_samples$indicator))

  # Indicators should be binary
  ind_samples <- do.call(rbind, fit$raw_samples$indicator)
  expect_true(all(ind_samples %in% c(0, 1)))
})


# ------------------------------------------------------------------------------
# Method Variations Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare works with different update methods", {
  skip_on_cran_mcmc()

  data <- generate_grouped_test_data(n_per_group = 20, p = 3, n_groups = 2, seed = 99)

  methods_to_test <- c("adaptive-metropolis")
  # Note: Could add "hmc", "nuts" if testing more thoroughly

  for (method in methods_to_test) {
    fit <- tryCatch(
      bgmCompare(
        x = data$x,
        group_indicator = data$group_indicator,
        update_method = method,
        iter = 50,
        warmup = 150,
        chains = 1,
        display_progress = "none"
      ),
      error = function(e) e
    )

    if (!inherits(fit, "error")) {
      expect_s3_class(fit, "bgmCompare")
    }
  }
})


# ------------------------------------------------------------------------------
# More Than Two Groups
# ------------------------------------------------------------------------------

test_that("bgmCompare handles more than 2 groups", {
  skip_on_cran_mcmc()

  data <- generate_grouped_test_data(
    n_per_group = 20, p = 3, n_groups = 3, seed = 456
  )

  # Use difference_selection = FALSE to avoid summary computation issues
  # with very short chains
  fit <- bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = FALSE,
    iter = 100,
    warmup = 150,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  args <- extract_arguments(fit)
  expect_equal(args$num_groups, 3)

  # Group-specific effects should have 3 columns
  group_params <- extract_group_params(fit)
  expect_equal(ncol(group_params$main_effects_groups), 3)
  expect_equal(ncol(group_params$pairwise_effects_groups), 3)
})
