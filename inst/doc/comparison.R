## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(bgms)

?Boredom
data_french = Boredom[Boredom$language == "fr", -1]
data_french = data_french[, 1:5]
data_english = Boredom[Boredom$language != "fr", -1]
data_english = data_english[, 1:5]

## ----include=FALSE------------------------------------------------------------
  fit <- readRDS(system.file("extdata", "fit_boredom.rds", package = "bgms"))

## ----eval = FALSE-------------------------------------------------------------
# fit = bgmCompare(x = data_french, y = data_english, seed = 1234)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
coef(fit)

## ----fig.width= 7, fig.height= 7----------------------------------------------
library(qgraph)

french_network = matrix(0, 5, 5)
french_network[lower.tri(french_network)] = coef(fit)$pairwise_effects_groups[, 1]
french_network = french_network + t(french_network)
colnames(french_network) = colnames(data_french)
rownames(french_network) = colnames(data_french)

qgraph(french_network,
       theme = "TeamFortress",
       maximum = 1,
       fade = FALSE,
       color = c("#f0ae0e"), vsize = 10, repulsion = .9,
       label.cex = 1, label.scale = "FALSE",
       labels = colnames(data_french))

