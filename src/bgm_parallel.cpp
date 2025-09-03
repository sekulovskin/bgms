// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "bgm_sampler.h"
#include <tbb/global_control.h>
#include <tbb/mutex.h>
#include "print_mutex.h"

using namespace Rcpp;
using namespace RcppParallel;

struct GibbsChainRunner : public Worker {
  const int num_chains;
  const arma::imat observations;
  const arma::ivec num_categories;
  const double interaction_scale;
  const std::string edge_prior;
  const arma::mat inclusion_probability;
  const double beta_bernoulli_alpha;
  const double beta_bernoulli_beta;
  const double dirichlet_alpha;
  const double lambda;
  const arma::imat interaction_index_matrix;
  const int iter;
  const int burnin;
  const arma::imat num_obs_categories;
  const arma::imat sufficient_blume_capel;
  const double threshold_alpha;
  const double threshold_beta;
  const bool na_impute;
  const arma::imat missing_index;
  const arma::uvec is_ordinal_variable;
  const arma::ivec reference_category;
  const bool edge_selection;
  const std::string update_method;
  const arma::imat pairwise_effect_indices;
  const double target_accept;
  const arma::imat sufficient_pairwise;
  const int hmc_num_leapfrogs;
  const int nuts_max_depth;
  const bool learn_mass_matrix;

  Rcpp::List& output;

  GibbsChainRunner(
    int num_chains,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const arma::imat& sufficient_pairwise,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    Rcpp::List& output
  ) :
    num_chains(num_chains),
    observations(observations), num_categories(num_categories),
    interaction_scale(interaction_scale), edge_prior(edge_prior),
    inclusion_probability(inclusion_probability),
    beta_bernoulli_alpha(beta_bernoulli_alpha), beta_bernoulli_beta(beta_bernoulli_beta),
    dirichlet_alpha(dirichlet_alpha), lambda(lambda),
    interaction_index_matrix(interaction_index_matrix),
    iter(iter), burnin(burnin),
    num_obs_categories(num_obs_categories), sufficient_blume_capel(sufficient_blume_capel),
    threshold_alpha(threshold_alpha), threshold_beta(threshold_beta),
    na_impute(na_impute), missing_index(missing_index),
    is_ordinal_variable(is_ordinal_variable), reference_category(reference_category),
    edge_selection(edge_selection), update_method(update_method),
    pairwise_effect_indices(pairwise_effect_indices), target_accept(target_accept),
    sufficient_pairwise(sufficient_pairwise),
    hmc_num_leapfrogs(hmc_num_leapfrogs), nuts_max_depth(nuts_max_depth),
    learn_mass_matrix(learn_mass_matrix), output(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      try {
        RNGScope scope;

        Rcpp::List result = run_gibbs_sampler_for_bgm(
          i + 1, // chain_id
          observations,
          num_categories,
          interaction_scale,
          edge_prior,
          inclusion_probability,
          beta_bernoulli_alpha,
          beta_bernoulli_beta,
          dirichlet_alpha,
          lambda,
          interaction_index_matrix,
          iter,
          burnin,
          num_obs_categories,
          sufficient_blume_capel,
          threshold_alpha,
          threshold_beta,
          na_impute,
          missing_index,
          is_ordinal_variable,
          reference_category,
          edge_selection,
          update_method,
          pairwise_effect_indices,
          target_accept,
          sufficient_pairwise,
          hmc_num_leapfrogs,
          nuts_max_depth,
          learn_mass_matrix
        );

        output[i] = result;

      } catch (std::exception& e) {
        output[i] = Rcpp::List::create(
          Rcpp::Named("error") = e.what(),
          Rcpp::Named("chain_id") = i + 1
        );
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List run_bgm_parallel(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const arma::imat& sufficient_pairwise,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const int num_chains,
    const int nThreads
) {
  Rcpp::List output(num_chains);

  GibbsChainRunner worker(
      num_chains, observations, num_categories, interaction_scale, edge_prior, inclusion_probability,
      beta_bernoulli_alpha, beta_bernoulli_beta, dirichlet_alpha, lambda,
      interaction_index_matrix, iter, burnin, num_obs_categories, sufficient_blume_capel,
      threshold_alpha, threshold_beta, na_impute, missing_index, is_ordinal_variable,
      reference_category, edge_selection, update_method, pairwise_effect_indices, target_accept,
      sufficient_pairwise, hmc_num_leapfrogs, nuts_max_depth, learn_mass_matrix, output
  );

  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, num_chains, worker);
  }

  return output;
}