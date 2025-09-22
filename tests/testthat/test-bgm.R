test_that("inclusion probabilities correlate with posterior mode", {
  data("Wenchuan", package = "bgms")
  fit <- bgm(x = Wenchuan, iter = 1e2, burnin = 10)

  posterior_modes = extract_pairwise_interactions(fit)
  posterior_incl_probs = extract_posterior_inclusion_probabilities(fit)

  posterior_modes = posterior_modes[lower.tri(posterior_modes)]
  posterior_incl_probs = posterior_incl_probs[lower.tri(posterior_incl_probs)]

  testthat::expect_gte(cor(abs(posterior_modes), posterior_incl_probs, method = "spearman"), .9)

})

on_ci <-   isTRUE(as.logical(Sys.getenv("CI", "false")))
no_cores <- if (on_ci) 2L else min(4, parallel::detectCores())

test_that("bgm is reproducible", {
  data("Wenchuan", package = "bgms")
  x <-  Wenchuan[1:50, 1:5]
  fit1 <- bgm(x = x, iter = 100, burnin = 1000, cores = no_cores, seed = 1234)
  fit2 <- bgm(x = x, iter = 100, burnin = 1000, cores = no_cores, seed = 1234)

  testthat::expect_equal(fit1$raw_samples, fit2$raw_samples)
})

test_that("bgmCompare is reproducible", {
  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:50, 1:5]
  y <- Wenchuan[1:50, c(1:4, 6)]
  fit1 <- bgmCompare2(x = x, y = y, iter = 100, burnin = 1000, cores = no_cores, seed = 1234)
  fit2 <- bgmCompare2(x = x, y = y, iter = 100, burnin = 1000, cores = no_cores, seed = 1234)

  combine_chains <- function(lst) {
    # without abind
    element_names <- names(lst[[1]])

    out <- lapply(element_names, function(nm) {
      arrays <- lapply(lst, function(chain) chain[[nm]])
      dims <- dim(arrays[[1]])
      if (is.null(dims)) {
        # handle scalar/vector case (e.g. chain_id)
        return(unlist(arrays))
      } else {
        # bind along a new last dimension
        arr <- array(
          unlist(arrays),
          dim = c(dims, length(arrays))
        )
        return(arr)
      }
    })

    names(out) <- element_names
    out
  }

  combined1 <- combine_chains(samples1)
  combined2 <- combine_chains(samples2)

  testthat::expect_equal(combined1, combined2)
})
