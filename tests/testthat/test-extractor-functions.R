# ==============================================================================
# Contract Tests for Extractor Functions (Parameterized)
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, symmetry checks, dimension consistency
#
# This file uses parameterized testing (specs + loop) to reduce code repetition.
# Each extractor is tested across multiple fixture types with shared assertions.
#
# IMPORTANT: Changes to extractor function output structure may break easybgm!
# ==============================================================================

# ------------------------------------------------------------------------------
# Fixture Specifications
# ------------------------------------------------------------------------------
# Define all fixtures to test against, with their properties

get_all_fixtures <- function() {
  list(
    list(
      label = "bgms_binary",
      get_fit = get_bgms_fit,
      type = "bgms",
      var_type = "binary"
    ),
    list(
      label = "bgms_ordinal",
      get_fit = get_bgms_fit_ordinal,
      type = "bgms",
      var_type = "ordinal"
    ),
    list(
      label = "bgms_blumecapel",
      get_fit = get_bgms_fit_blumecapel,
      type = "bgms",
      var_type = "blume-capel"
    ),
    list(
      label = "bgmCompare_binary",
      get_fit = get_bgmcompare_fit,
      type = "bgmCompare",
      var_type = "binary"
    ),
    list(
      label = "bgmCompare_ordinal",
      get_fit = get_bgmcompare_fit_ordinal,
      type = "bgmCompare",
      var_type = "ordinal"
    )
  )
}


# ------------------------------------------------------------------------------
# extract_arguments() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_arguments returns complete argument list for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    # Basic structure
    expect_true(is.list(args), info = paste(ctx, "should be list"))
    expect_true(length(args) > 0, info = ctx)
    
    # Essential fields present
    expect_true("num_variables" %in% names(args), info = paste(ctx, "missing num_variables"))
    expect_true("num_cases" %in% names(args), info = paste(ctx, "missing num_cases"))
    expect_true("data_columnnames" %in% names(args), info = paste(ctx, "missing data_columnnames"))
    
    # Values are sensible
    expect_true(args$num_variables >= 1, info = ctx)
    expect_true(args$num_cases >= 1, info = ctx)
    
    # Type-specific fields
    if (spec$type == "bgms") {
      expect_true(is.logical(args$edge_selection), info = paste(ctx, "edge_selection should be logical"))
    } else {
      expect_true(args$num_groups >= 2, info = paste(ctx, "bgmCompare should have >= 2 groups"))
    }
  }
})

test_that("extract_arguments errors on non-bgms objects", {
  expect_error(extract_arguments(list()), class = "error")
  expect_error(extract_arguments(data.frame()), class = "error")
})


# ------------------------------------------------------------------------------
# extract_pairwise_interactions() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_pairwise_interactions returns valid matrix for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    pairwise <- extract_pairwise_interactions(fit)
    
    # Structure checks
    expect_true(is.matrix(pairwise), info = paste(ctx, "should be matrix"))
    
    p <- args$num_variables
    expected_cols <- p * (p - 1) / 2
    expect_equal(ncol(pairwise), expected_cols, 
                 info = paste(ctx, "wrong number of edge columns"))
    
    # Values finite
    expect_true(all(is.finite(pairwise)), info = paste(ctx, "should have finite values"))
    
    # Has column names
    expect_true(!is.null(colnames(pairwise)), info = paste(ctx, "should have column names"))
  }
})


# ------------------------------------------------------------------------------
# extract_category_thresholds() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_category_thresholds returns valid output for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    thresholds <- extract_category_thresholds(fit)
    
    # Structure checks
    expect_true(is.matrix(thresholds), info = paste(ctx, "should be matrix"))
    
    # Values finite where not NA
    vals <- thresholds[!is.na(thresholds)]
    expect_true(all(is.finite(vals)), info = paste(ctx, "non-NA values should be finite"))
  }
})


# ------------------------------------------------------------------------------
# extract_indicators() and extract_posterior_inclusion_probabilities() Tests
# ------------------------------------------------------------------------------
# These only apply to fits with edge_selection = TRUE

test_that("extract_indicators returns binary matrix for edge-selection fits", {
  # Only test fixtures with edge selection
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    indicators <- extract_indicators(fit)
    
    # Structure
    expect_true(is.matrix(indicators), info = ctx)
    
    p <- args$num_variables
    expected_cols <- p * (p - 1) / 2
    expect_equal(ncol(indicators), expected_cols, info = paste(ctx, "wrong indicator columns"))
    
    # Binary values
    expect_true(all(indicators %in% c(0, 1)), 
                info = paste(ctx, "indicators should be 0 or 1"))
  }
})

test_that("extract_posterior_inclusion_probabilities returns symmetric PIP matrix", {
  # Only test fixtures with edge selection
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    pip <- extract_posterior_inclusion_probabilities(fit)
    p <- args$num_variables
    
    # Structure
    expect_true(is.matrix(pip), info = ctx)
    expect_equal(dim(pip), c(p, p), info = paste(ctx, "should be p x p"))
    
    # Symmetry
    expect_true(is_symmetric(pip), info = paste(ctx, "should be symmetric"))
    
    # Range [0, 1]
    expect_true(values_in_range(pip, 0, 1), info = paste(ctx, "PIPs should be in [0,1]"))
    
    # Diagonal is zero (no self-loops)
    expect_true(all(diag(pip) == 0), info = paste(ctx, "diagonal should be 0"))
    
    # Has variable names
    expect_equal(colnames(pip), args$data_columnnames, info = ctx)
  }
})

test_that("extract_indicators errors when edge_selection = FALSE", {
  skip_on_cran_mcmc()
  
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)
  
  expect_error(extract_indicators(fit), regexp = "edge_selection")
})


# ------------------------------------------------------------------------------
# extract_rhat() and extract_ess() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_rhat returns valid diagnostics for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    rhat <- extract_rhat(fit)
    
    expect_true(is.list(rhat), info = paste(ctx, "should be list"))
    
    if (spec$type == "bgms") {
      expect_true("pairwise" %in% names(rhat), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(rhat$pairwise), info = ctx)
      expect_true(all(is.na(rhat$pairwise) | rhat$pairwise > 0), 
                  info = paste(ctx, "R-hat should be positive"))
    } else {
      expect_true("pairwise_baseline" %in% names(rhat), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(rhat$pairwise_baseline), info = ctx)
      expect_true(all(is.na(rhat$pairwise_baseline) | rhat$pairwise_baseline > 0), info = ctx)
    }
  }
})

test_that("extract_ess returns valid diagnostics for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    ess <- extract_ess(fit)
    
    expect_true(is.list(ess), info = paste(ctx, "should be list"))
    
    if (spec$type == "bgms") {
      expect_true("pairwise" %in% names(ess), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(ess$pairwise), info = ctx)
      expect_true(all(is.na(ess$pairwise) | ess$pairwise > 0), 
                  info = paste(ctx, "ESS should be positive"))
    } else {
      expect_true("pairwise_baseline" %in% names(ess), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(ess$pairwise_baseline), info = ctx)
    }
  }
})

test_that("extract_rhat and extract_ess error on non-bgms objects", {
  expect_error(extract_rhat(list()), class = "error")
  expect_error(extract_rhat(data.frame()), class = "error")
  expect_error(extract_ess(list()), class = "error")
  expect_error(extract_ess(data.frame()), class = "error")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors() Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors returns prior specification", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$edge_selection)) {
    skip("Fit object does not have edge_selection = TRUE")
  }
  
  priors <- extract_indicator_priors(fit)
  
  expect_type(priors, "list")
  expect_true("type" %in% names(priors))
  
  valid_types <- c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
  expect_true(priors$type %in% valid_types)
  
  # Type-specific checks
  if (priors$type == "Bernoulli") {
    expect_true("prior_inclusion_probability" %in% names(priors))
    pip <- priors$prior_inclusion_probability
    expect_true(all(pip >= 0 & pip <= 1))
  }
  
  if (priors$type == "Beta-Bernoulli") {
    expect_true(all(c("alpha", "beta") %in% names(priors)))
    expect_true(priors$alpha > 0 && priors$beta > 0)
  }
})

test_that("extract_indicator_priors errors when no selection performed", {
  skip_on_cran_mcmc()
  
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)
  
  expect_error(extract_indicator_priors(fit), regexp = "selection")
})


# ------------------------------------------------------------------------------
# bgmCompare-specific Tests
# ------------------------------------------------------------------------------

test_that("extract_group_params returns group-level parameters", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  
  group_params <- extract_group_params(fit)
  
  expect_type(group_params, "list")
  expect_true("main_effects_groups" %in% names(group_params))
  expect_true("pairwise_effects_groups" %in% names(group_params))
  
  # Dimensions match number of groups
  n_groups <- args$num_groups
  expect_equal(ncol(group_params$main_effects_groups), n_groups)
  expect_equal(ncol(group_params$pairwise_effects_groups), n_groups)
  
  # Values finite
  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# Deprecated Function Tests
# ------------------------------------------------------------------------------

test_that("deprecated extractors still work with warning", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$edge_selection)) {
    skip("Fit object does not have edge_selection = TRUE")
  }
  
  expect_warning(result <- extract_edge_indicators(fit), regexp = "deprecated")
  expect_true(is.matrix(result))
  
  expect_warning(result <- extract_pairwise_thresholds(fit), regexp = "deprecated")
  expect_true(is.matrix(result))
})


# ------------------------------------------------------------------------------
# Cross-Function Consistency Tests
# ------------------------------------------------------------------------------

test_that("extractor outputs are dimensionally consistent", {
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    p <- args$num_variables
    n_edges <- p * (p - 1) / 2
    
    # All should agree on number of variables/edges
    pip <- extract_posterior_inclusion_probabilities(fit)
    expect_equal(nrow(pip), p, info = paste(ctx, "PIP rows"))
    
    indicators <- extract_indicators(fit)
    expect_equal(ncol(indicators), n_edges, info = paste(ctx, "indicator cols"))
    
    pairwise <- extract_pairwise_interactions(fit)
    expect_equal(ncol(pairwise), n_edges, info = paste(ctx, "pairwise cols"))
    
    thresholds <- extract_category_thresholds(fit)
    expect_equal(nrow(thresholds), p, info = paste(ctx, "threshold rows"))
  }
})


# ------------------------------------------------------------------------------
# Contract Tests for easybgm Integration
# ------------------------------------------------------------------------------

test_that("bgms fit contains all fields accessed by easybgm", {
  fixtures <- list(
    list(label = "bgms", get_fit = get_bgms_fit, type = "bgms"),
    list(label = "bgmCompare", get_fit = get_bgmcompare_fit, type = "bgmCompare")
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (spec$type == "bgms") {
      expect_true("posterior_summary_pairwise" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise), info = ctx)
      
      if (isTRUE(args$edge_selection)) {
        expect_true("posterior_summary_indicator" %in% names(fit), info = ctx)
        expect_true("n_eff" %in% names(fit$posterior_summary_indicator), info = ctx)
      }
    } else {
      expect_true("posterior_summary_pairwise_baseline" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
    }
  }
})
