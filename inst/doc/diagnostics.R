## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----eval=FALSE---------------------------------------------------------------
# library(bgms)
# data = Wenchuan[, 1:5]
# fit = bgm(data, seed = 1234)

## ----include=FALSE------------------------------------------------------------
library(bgms)
data = Wenchuan[, 1:5]
fit <- readRDS(system.file("extdata", "fit_5items.rds", package = "bgms"))

## -----------------------------------------------------------------------------
summary(fit)$pairwise

## ----fig.width= 7, fig.height= 7----------------------------------------------
library(coda)

param_index = 1
chains = lapply(fit$raw_samples$pairwise, function(mat) mat[, param_index])
mcmc_obj = mcmc.list(lapply(chains, mcmc))

traceplot(mcmc_obj, col = c("firebrick", "steelblue", "darkgreen", "goldenrod"), 
          main = "Traceplot of pairwise[1]") 

## -----------------------------------------------------------------------------
coef(fit)$indicator

## -----------------------------------------------------------------------------
# Example for one edge
p = coef(fit)$indicator[1, 5]
BF_10 = p / (1 - p)
BF_10

## -----------------------------------------------------------------------------
1 / BF_10

