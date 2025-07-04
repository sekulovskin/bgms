#' Bayesian Estimation or Edge Selection for Markov Random Fields
#'
#' @description
#' The \code{bgm} function estimates the pseudoposterior distribution of threshold and
#' pairwise interaction parameters in a Markov Random Field (MRF) model for binary and/or
#' ordinal variables. Optionally, it performs Bayesian edge selection using spike-and-slab
#' priors to infer the network structure.
#'
#' @details
#' This function models the joint distribution of binary and ordinal variables
#' using a Markov Random Field, with support for edge selection through Bayesian
#' variable selection. Key components of the model are described in the sections below.
#'
#' @section Ordinal Variables:
#' The function supports two types of ordinal variables:
#'
#' \strong{Regular ordinal variables}:
#' Assign a threshold parameter to each response category except the lowest. The model imposes no additional constraints on the distribution of category responses.
#'
#' \strong{Blume-Capel ordinal variables}:
#' Assume a reference category (e.g., a “neutral” response) and score responses by distance from this reference. Thresholds are modeled quadratically:
#'
#' \deqn{\mu_{c} = \alpha \cdot c + \beta \cdot (c - r)^2}
#'
#' where:
#' \itemize{
#'   \item \eqn{\mu_{c}}: threshold for category \eqn{c}
#'   \item \eqn{\alpha}: linear trend across categories
#'   \item \eqn{\beta}: preference toward or away from the reference
#'    \itemize{
#'      \item If \eqn{\beta < 0}, the model favors responses near the reference category;
#'      \item if \eqn{\beta > 0}, it favors responses farther away (i.e., extremes).
#'    }
#'   \item \eqn{r}: reference category
#' }
#'
#'
#' @section Edge Selection:
#' When \code{edge_selection = TRUE}, the function performs Bayesian variable selection
#' on the pairwise interactions (edges) in the MRF using spike-and-slab priors.
#'
#' Supported priors for edge inclusion:
#' \itemize{
#'   \item \strong{Bernoulli}: Fixed inclusion probability across edges.
#'   \item \strong{Beta-Bernoulli}: Inclusion probability is assigned a Beta prior distribution.
#'   \item \strong{Stochastic Block Model}: Cluster-based edge priors with Beta, Dirichlet, and Poisson hyperpriors.
#' }
#'
#' All priors operate via binary indicator variables controlling the inclusion or exclusion of each edge.
#'
#' @section Prior Distributions:
#'
#' \itemize{
#'   \item \strong{Interaction parameters}: Modeled with a Cauchy slab prior.
#'   \item \strong{Threshold parameters}: Modeled using a beta-prime distribution.
#'   \item \strong{Edge indicators}: Use either a Bernoulli, Beta-Bernoulli, or SBM prior (as above).
#' }
#'
#' @section Gibbs Sampling:
#'
#' Parameters are estimated using a Metropolis-within-Gibbs sampling scheme.
#' When \code{edge_selection = TRUE}, the algorithm runs \code{2 * burnin} warmup iterations:
#' \itemize{
#'   \item First half without edge selection.
#'   \item Second half with edge selection enabled.
#' }
#' This warmup strategy improves stability of adaptive Metropolis-Hastings proposals and starting values.
#'
#'
#' @section Missing Data:
#'
#' If \code{na_action = "listwise"}, observations with missing values are removed.
#' If \code{na_action = "impute"}, missing values are imputed during MCMC.
#'

#' @param x A data frame or matrix with \code{n} rows and \code{p} columns containing binary and ordinal responses. Binary and ordinal variables are automatically recoded to non-negative integers (\code{0, 1, ..., m}). For regular ordinal variables, unobserved categories are collapsed; for Blume-Capel variables, all categories are retained.
#'
#' @param variable_type Character or character vector. Specifies the type of each variable in \code{x}. Allowed values: \code{"ordinal"} or \code{"blume-capel"}. Binary variables are automatically treated as \code{"ordinal"}. Default: \code{"ordinal"}.
#'
#' @param reference_category Integer or vector. Reference category used in Blume-Capel variables. Can be a single integer (applied to all) or a vector of length \code{p}. Required if at least one variable is of type \code{"blume-capel"}.
#'
#' @param iter Integer. Number of Gibbs sampling iterations. Default: \code{1e4}. For stable estimates, consider using at least \code{1e5}.
#'
#' @param burnin Integer. Number of burn-in iterations before saving samples. When \code{edge_selection = TRUE}, the function runs \code{2 * burnin} iterations: first half without edge selection, second half with edge selection. Default: \code{1e3}.
#'
#' @param interaction_scale Double. Scale of the Cauchy prior for pairwise interaction parameters. Default: \code{2.5}.
#'
#' @param threshold_alpha,threshold_beta Double. Shape parameters of the beta-prime prior for threshold parameters. Must be positive. If equal, the prior is symmetric. Defaults: \code{threshold_alpha = 0.5} and \code{threshold_beta = 0.5}.
#'
#' @param edge_selection Logical. Whether to perform Bayesian edge selection. If \code{FALSE}, the model estimates all edges. Default: \code{TRUE}.
#'
#' @param edge_prior Character. Specifies the prior for edge inclusion. Options:
#' \code{"Bernoulli"}, \code{"Beta-Bernoulli"}, or \code{"Stochastic-Block"}. Default: \code{"Bernoulli"}.
#'
#' @param inclusion_probability Numeric scalar or matrix. Prior inclusion probability of each edge (used with the Bernoulli prior). A single value applies to all edges; a matrix allows edge-specific probabilities. Default: \code{0.5}.
#'
#' @param beta_bernoulli_alpha,beta_bernoulli_beta Double. Shape parameters for the beta distribution in the Beta-Bernoulli prior. Must be positive. Defaults: \code{beta_bernoulli_alpha = 1} and \code{beta_bernoulli_beta = 1}.
#'
#' @param dirichlet_alpha Double. Concentration parameter of the Dirichlet prior on block assignments (used with the Stochastic Block Model). Default: \code{1}.
#'
#' @param lambda Double. Rate of the zero-truncated Poisson prior on the number of clusters in the Stochastic Block Model. Default: \code{1}.
#'
#' @param na_action Character. Specifies missing data handling. \code{"listwise"} deletes rows with missing values. \code{"impute"} imputes missing values during MCMC. Default: \code{"listwise"}.
#'
#' @param save Logical; \strong{deprecated}. Whether to return all sampled states from the Gibbs sampler. If \code{FALSE}, only posterior means are returned. Default: \code{FALSE}.
#'
#' @param display_progress Logical. Whether to show a progress bar during sampling. Default: \code{TRUE}.
#'
#' @return
#' A list of class \code{"bgms"} containing posterior summaries or sampled states, depending on the \code{save} option:
#'
#' \strong{If \code{save = FALSE}} (default), the list contains model-averaged posterior summaries:
#' \itemize{
#'   \item \code{indicator}: A \code{p × p} matrix of posterior inclusion probabilities for each edge.
#'   \item \code{interactions}: A \code{p × p} matrix of posterior means for pairwise interactions.
#'   \item \code{thresholds}: A \code{p × max(m)} matrix of posterior means for threshold parameters. For Blume-Capel variables, the first entry corresponds to the linear term and the second to the quadratic term.
#' }
#'
#' \strong{If \code{save = TRUE}}, the list also includes raw MCMC samples:
#' \itemize{
#'   \item \code{indicator}: A matrix with \code{iter} rows and \code{p × (p - 1) / 2} columns; sampled edge inclusion indicators.
#'   \item \code{interactions}: A matrix with \code{iter} rows and \code{p × (p - 1) / 2} columns; sampled interaction parameters.
#'   \item \code{thresholds}: A matrix with \code{iter} rows and \code{sum(m)} columns; sampled threshold parameters.
#' }
#'
#' \strong{If \code{edge_prior = "Stochastic-Block"}}, two additional components may be returned:
#' \itemize{
#'   \item \code{allocations}: A vector (or matrix if \code{save = TRUE}) with posterior cluster assignments for each node.
#'   \item \code{components}: A matrix of posterior probabilities for the number of clusters, based on sampled allocations.
#' }
#'
#' Column-wise averages of the sampled matrices yield posterior means, except for \code{allocations}, which should be summarized using \code{summarySBM()}.
#'
#' The returned list also includes some of the function call arguments, useful for post-processing.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Store user par() settings
#' op <- par(no.readonly = TRUE)
#'
#' # Run bgm on the Wenchuan dataset
#' # For reliable results, consider using at least 1e5 iterations
#' fit <- bgm(x = Wenchuan, iter = 1e4)
#'
#' #--- INCLUSION VS EDGE WEIGHT ----------------------------------------------
#' edge_weights <- fit$interactions[lower.tri(fit$interactions)]
#' incl_probs   <- fit$indicator[lower.tri(fit$indicator)]
#'
#' par(mar = c(5, 5, 1, 1) + 0.1, cex = 1.7)
#' plot(edge_weights, incl_probs,
#'      pch = 21, bg = "gray", cex = 1.3,
#'      ylim = c(0, 1), axes = FALSE,
#'      xlab = "", ylab = "")
#' abline(h = c(0, 0.5, 1), lty = 2, col = "gray")
#' axis(1); axis(2, las = 1)
#' mtext("Posterior Mean Edge Weight", side = 1, line = 3, cex = 1.7)
#' mtext("Posterior Inclusion Probability", side = 2, line = 3, cex = 1.7)
#'
#' #--- EVIDENCE PLOT ----------------------------------------------------------
#' prior_odds <- 1
#' post_odds <- incl_probs / (1 - incl_probs)
#' log_bf <- log(post_odds / prior_odds)
#' log_bf <- pmin(log_bf, 5)  # cap extreme values
#'
#' plot(edge_weights, log_bf,
#'      pch = 21, bg = "#bfbfbf", cex = 1.3,
#'      axes = FALSE, xlab = "", ylab = "",
#'      ylim = c(-5, 5.5), xlim = c(-0.5, 1.5))
#' axis(1); axis(2, las = 1)
#' abline(h = log(c(1/10, 10)), lwd = 2, col = "#bfbfbf")
#' text(1, log(1 / 10), "Evidence for Exclusion", pos = 1, cex = 1.1)
#' text(1, log(10),     "Evidence for Inclusion", pos = 3, cex = 1.1)
#' text(1, 0,           "Absence of Evidence", cex = 1.1)
#' mtext("Log-Inclusion Bayes Factor", side = 2, line = 3, cex = 1.7)
#' mtext("Posterior Mean Interactions", side = 1, line = 3.7, cex = 1.7)
#'
#' #--- MEDIAN PROBABILITY NETWORK --------------------------------------------
#' median_edges <- ifelse(incl_probs >= 0.5, edge_weights, 0)
#' n <- ncol(Wenchuan)
#' net <- matrix(0, n, n)
#' net[lower.tri(net)] <- median_edges
#' net <- net + t(net)
#' dimnames(net) <- list(colnames(Wenchuan), colnames(Wenchuan))
#'
#' par(cex = 1)
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   qgraph::qgraph(net,
#'     theme = "TeamFortress", maximum = 0.5, fade = FALSE,
#'     color = "#f0ae0e", vsize = 10, repulsion = 0.9,
#'     label.cex = 1.1, label.scale = FALSE,
#'     labels = colnames(Wenchuan)
#'   )
#' }
#'
#' # Restore user par() settings
#' par(op)
#' }
#'
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#'
#' @export
bgm = function(x,
               variable_type = "ordinal",
               reference_category,
               iter = 1e4,
               burnin = 5e2,
               interaction_scale = 2.5,
               threshold_alpha = 0.5,
               threshold_beta = 0.5,
               edge_selection = TRUE,
               edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
               inclusion_probability = 0.5,
               beta_bernoulli_alpha = 1,
               beta_bernoulli_beta = 1,
               dirichlet_alpha = 1,
               lambda = 1,
               na_action = c("listwise", "impute"),
               save = FALSE,
               display_progress = TRUE) {

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Check model input -----------------------------------------------------------
  model = check_model(x = x,
                      variable_type = variable_type,
                      reference_category = reference_category,
                      interaction_scale = interaction_scale,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      edge_selection = edge_selection,
                      edge_prior = edge_prior,
                      inclusion_probability = inclusion_probability,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      dirichlet_alpha = dirichlet_alpha,
                      lambda = lambda)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  variable_bool = model$variable_bool
  # ----------------------------------------------------------------------------

  reference_category = model$reference_category
  edge_selection = model$edge_selection
  edge_prior = model$edge_prior
  theta = model$theta

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > .Machine$double.eps)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(iter <= 0)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > .Machine$double.eps || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")
  if(burnin <= 0)
    stop("Parameter ``burnin'' needs to be a positive integer.")

  #Check na_action -------------------------------------------------------------
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if(inherits(na_action, what = "try-error"))
    stop(paste0("The na_action argument should equal listwise or impute, not ",
                na_action_input,
                "."))
  #Check save ------------------------------------------------------------------
  save_input = save
  save = as.logical(save)
  if(is.na(save))
    stop(paste0("The save argument should equal TRUE or FALSE, not ",
                save_input,
                "."))

  #Check display_progress ------------------------------------------------------
  display_progress = as.logical(display_progress)
  if(is.na(display_progress))
    stop("The display_progress argument should equal TRUE or FALSE.")

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x,
                       na_action = na_action,
                       variable_bool = variable_bool,
                       reference_category = reference_category)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na_impute = data$na_impute
  reference_category = data$reference_category

  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2
  no_thresholds = sum(no_categories)

  #Specify the variance of the (normal) proposal distribution ------------------
  proposal_sd = matrix(1,
                       nrow = no_variables,
                       ncol = no_variables)
  proposal_sd_blumecapel = matrix(1,
                                  nrow = no_variables,
                                  ncol = 2)

  # Starting value of model matrix ---------------------------------------------
  indicator = matrix(1,
                 nrow = no_variables,
                 ncol = no_variables)


  #Starting values of interactions and thresholds (posterior mode) -------------
  interactions = matrix(0, nrow = no_variables, ncol = no_variables)
  thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))

  #Precompute the number of observations per category for each variable --------
  n_cat_obs = matrix(0,
                     nrow = max(no_categories) + 1,
                     ncol = no_variables)
  for(variable in 1:no_variables) {
    for(category in 0:no_categories[variable]) {
      n_cat_obs[category + 1, variable] = sum(x[, variable] == category)
    }
  }

  #Precompute the sufficient statistics for the two Blume-Capel parameters -----
  sufficient_blume_capel = matrix(0, nrow = 2, ncol = no_variables)
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
    bc_vars = which(!variable_bool)
    for(i in bc_vars) {
      sufficient_blume_capel[1, i] = sum(x[, i])
      sufficient_blume_capel[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
    }
  }

  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0,
                 nrow = no_variables * (no_variables - 1) / 2,
                 ncol = 3)
  cntr = 0
  for(variable1 in 1:(no_variables - 1)) {
    for(variable2 in (variable1 + 1):no_variables) {
      cntr =  cntr + 1
      Index[cntr, 1] = cntr
      Index[cntr, 2] = variable1 - 1
      Index[cntr, 3] = variable2 - 1
    }
  }

  #Preparing the output --------------------------------------------------------
  arguments = list(
    no_variables = no_variables,
    no_cases = nrow(x),
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = theta,
    beta_bernoulli_alpha = beta_bernoulli_alpha ,
    beta_bernoulli_beta =  beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda,
    na_action = na_action,
    save = save,
    version = packageVersion("bgms")
  )

  #The Metropolis within Gibbs sampler -----------------------------------------
  out = gibbs_sampler(observations = x,
                      indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      no_categories  = no_categories,
                      interaction_scale = interaction_scale,
                      proposal_sd = proposal_sd,
                      proposal_sd_blumecapel = proposal_sd_blumecapel,
                      edge_prior = edge_prior,
                      theta = theta,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      dirichlet_alpha = dirichlet_alpha,
                      lambda = lambda,
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs,
                      sufficient_blume_capel = sufficient_blume_capel,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      na_impute = na_impute,
                      missing_index = missing_index,
                      variable_bool = variable_bool,
                      reference_category = reference_category,
                      save = save,
                      display_progress = display_progress,
                      edge_selection = edge_selection)

  if(save == FALSE) {
    if(edge_selection == TRUE) {
      indicator = out$indicator
    }
    interactions = out$interactions
    tresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames = paste0("variable ", 1:no_variables)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(indicator) = data_columnnames
        rownames(indicator) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    } else {
      data_columnnames <- colnames(x)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(indicator) = data_columnnames
        rownames(indicator) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    }

    if(any(variable_bool)) {
      colnames(thresholds) = paste0("category ", 1:max(no_categories))
    } else {
      thresholds = thresholds[, 1:2]
      colnames(thresholds) = c("linear", "quadratic")
    }

    arguments$data_columnnames = data_columnnames

    if(edge_selection == TRUE) {
      if(edge_prior == "Stochastic-Block"){
        output = list(
          indicator = indicator, interactions = interactions,
          thresholds = thresholds, allocations = out$allocations,
          arguments = arguments)
        class(output) = "bgms"
        summary_Sbm = summarySBM(output,
                                 internal_call = TRUE)

        output$components = summary_Sbm$components
        output$allocations = summary_Sbm$allocations
      } else {

        output = list(
          indicator = indicator, interactions = interactions,
          thresholds = thresholds, arguments = arguments)
        }

    } else {
      output = list(
        interactions = interactions, thresholds = thresholds,
        arguments = arguments)
    }

    class(output) = "bgms"
  } else {
    if(edge_selection == TRUE) {
      indicator = out$indicator
    }
    interactions = out$interactions
    thresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames <- 1:ncol(x)
    } else {
      data_columnnames <- colnames(x)
    }
    p <- ncol(x)
    names_bycol <- matrix(rep(data_columnnames, each = p), ncol = p)
    names_byrow <- matrix(rep(data_columnnames, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    names_vec <- names_comb[lower.tri(names_comb)]

    if(edge_selection == TRUE) {
      colnames(indicator) = names_vec
    }
    colnames(interactions) = names_vec
    names = character(length = sum(no_categories))
    cntr = 0
    for(variable in 1:no_variables) {
      for(category in 1:no_categories[variable]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",variable, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    if(edge_selection == TRUE) {
      dimnames(indicator) = list(Iter. = 1:iter, colnames(indicator))
    }
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))
    dimnames(thresholds) = list(Iter. = 1:iter, colnames(thresholds))

    arguments$data_columnnames = data_columnnames

    if(edge_selection == TRUE) {
      if(edge_prior == "Stochastic-Block"){
        output = list(indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      allocations = out$allocations,
                      arguments = arguments)
      } else {
        output = list(indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      arguments = arguments)
      }
    } else {
      output = list(interactions = interactions,
                    thresholds = thresholds,
                    arguments = arguments)
    }
    class(output) = "bgms"
  }
  return(output)
}