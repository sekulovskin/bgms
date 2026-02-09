# ==============================================================================
# Core bgm() Function Tests
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Reproducibility, correlation with sufficient statistics
#
# These tests verify core bgm() functionality:
#   - Reproducibility: identical seeds produce identical MCMC chains
#   - Sanity: posterior means correlate with classical sufficient statistics
#
# See helper-fixtures.R for testing philosophy and session-cached fixtures.
# ==============================================================================

test_that("Posterior means correlate with sufficient precision statistics", {
  testthat::skip_on_cran()

  fit = bgm(Wenchuan, edge_selection = FALSE, iter = 100, warmup = 1000, seed = 1234, chains = 1, display_progress = "none")
  x = Wenchuan
  x = na.omit(x)
  alt = -solve(t(x) %*% x)
  alt = alt[lower.tri(alt)]
  posterior_means = colMeans(extract_pairwise_interactions(fit))
  testthat::expect_gte(cor(posterior_means, alt, method = "spearman"), .98)
})

on_ci <- isTRUE(as.logical(Sys.getenv("CI", "false")))
no_cores <- if(on_ci) 2L else min(4, parallel::detectCores())

test_that("bgm is reproducible", {
  testthat::skip_on_cran()

  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:50, 1:5]
  fit1 <- bgm(x = x, iter = 100, warmup = 1000, cores = no_cores, seed = 1234, display_progress = "none")
  fit2 <- bgm(x = x, iter = 100, warmup = 1000, cores = no_cores, seed = 1234, display_progress = "none")

  testthat::expect_equal(fit1$raw_samples, fit2$raw_samples)
})

test_that("bgmCompare is reproducible", {
  testthat::skip_on_cran()

  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:50, 1:5]
  y <- Wenchuan[1:50, c(1:4, 6)]
  fit1 <- bgmCompare(x = x, y = y, iter = 100, warmup = 1000, cores = no_cores, seed = 1234, display_progress = "none")
  fit2 <- bgmCompare(x = x, y = y, iter = 100, warmup = 1000, cores = no_cores, seed = 1234, display_progress = "none")

  combine_chains <- function(fit) {
    pairs <- fit$raw_samples$pairwise
    pair <- do.call(rbind, pairs)
    mains <- fit$raw_samples$main
    main <- do.call(rbind, mains)
    inds <- fit$raw_samples$indicator
    ind <- do.call(rbind, inds)
    return(cbind(main, pair, ind))
  }

  combined1 <- combine_chains(fit1)
  combined2 <- combine_chains(fit2)

  testthat::expect_equal(combined1, combined2)
})
