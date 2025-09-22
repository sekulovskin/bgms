## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = identical(Sys.getenv("BUILD_VIGNETTE"), "true")
)

## ----setup--------------------------------------------------------------------
# library(bgms)

## ----usage2, eval=FALSE-------------------------------------------------------
# bgmCompare(x,
#            y,
#            difference_selection = TRUE,
#            main_difference_model = c("Collapse", "Constrain", "Free"),
#            variable_type = "ordinal",
#            reference_category,
#            pairwise_difference_scale = 1,
#            main_difference_scale = 1,
#            pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
#            main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
#            pairwise_difference_probability = 0.5,
#            main_difference_probability = 0.5,
#            pairwise_beta_bernoulli_alpha = 1,
#            pairwise_beta_bernoulli_beta = 1,
#            main_beta_bernoulli_alpha = 1,
#            main_beta_bernoulli_beta = 1,
#            interaction_scale = 2.5,
#            threshold_alpha = 0.5,
#            threshold_beta = 0.5,
#            iter = 1e4,
#            burnin = 1e3,
#            na.action = c("listwise", "impute"),
#            save = FALSE,
#            display_progress = TRUE)

