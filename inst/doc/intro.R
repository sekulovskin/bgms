## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(bgms)

# Analyse a subset of the Wenchuan dataset
?Wenchuan
data = Wenchuan[, 1:5]
head(data)

## ----include=FALSE------------------------------------------------------------
fit <- readRDS(system.file("extdata", "fit_5items.rds", package = "bgms"))

## ----eval=FALSE---------------------------------------------------------------
# fit = bgm(data, seed = 1234)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
coef(fit)

## ----fig.width= 7, fig.height= 7----------------------------------------------
library(qgraph)

median_probability_network = coef(fit)$pairwise
median_probability_network[coef(fit)$indicator < 0.5] = 0.0

qgraph(median_probability_network, 
       theme = "TeamFortress", 
       maximum = 1,
       fade = FALSE,
       color = c("#f0ae0e"), vsize = 10, repulsion = .9, 
       label.cex = 1, label.scale = "FALSE", 
       labels = colnames(data))

