// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "bgmCompare_sampler.h"
#include <tbb/global_control.h>
#include <vector>
#include <string>
#include "rng_utils.h"

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// Result struct
// -----------------------------------------------------------------------------
struct ChainResultCompare {
  bool error;
  std::string error_msg;
  int chain_id;
  SamplerOutput result;
};

// -----------------------------------------------------------------------------
// Worker struct
// -----------------------------------------------------------------------------
struct GibbsCompareChainRunner : public Worker {
  // inputs
  const arma::imat& observations;
  const int num_groups;
  const std::vector<arma::imat>& num_obs_categories_master;
  const std::vector<arma::imat>& sufficient_blume_capel_master;
  const std::vector<arma::mat>&  sufficient_pairwise_master;
  const arma::ivec& num_categories;
  const double main_alpha;
  const double main_beta;
  const double pairwise_scale;
  const double difference_scale;
  const double difference_selection_alpha;
  const double difference_selection_beta;
  const std::string& difference_prior;
  const int iter;
  const int burnin;
  const bool na_impute;
  const arma::imat& missing_data_indices;
  const arma::uvec& is_ordinal_variable;
  const arma::ivec& baseline_category;
  const bool difference_selection;
  const arma::imat& main_effect_indices;
  const arma::imat& pairwise_effect_indices;
  const double target_accept;
  const int nuts_max_depth;
  const bool learn_mass_matrix;
  const arma::mat& projection;
  const arma::ivec& group_membership;
  const arma::imat& group_indices;
  const arma::imat& interaction_index_matrix;
  const arma::mat& inclusion_probability_master;

  // RNG seeds
  const std::vector<SafeRNG>& chain_rngs;

  // output
  std::vector<ChainResultCompare>& results;

  GibbsCompareChainRunner(
    const arma::imat& observations,
    int num_groups,
    const std::vector<arma::imat>& num_obs_categories_master,
    const std::vector<arma::imat>& sufficient_blume_capel_master,
    const std::vector<arma::mat>&  sufficient_pairwise_master,
    const arma::ivec& num_categories,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    double difference_scale,
    double difference_selection_alpha,
    double difference_selection_beta,
    const std::string& difference_prior,
    int iter,
    int burnin,
    bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    const arma::mat& inclusion_probability_master,
    const std::vector<SafeRNG>& chain_rngs,
    std::vector<ChainResultCompare>& results
  ) :
    observations(observations),
    num_groups(num_groups),
    num_obs_categories_master(num_obs_categories_master),
    sufficient_blume_capel_master(sufficient_blume_capel_master),
    sufficient_pairwise_master(sufficient_pairwise_master),
    num_categories(num_categories),
    main_alpha(main_alpha),
    main_beta(main_beta),
    pairwise_scale(pairwise_scale),
    difference_scale(difference_scale),
    difference_selection_alpha(difference_selection_alpha),
    difference_selection_beta(difference_selection_beta),
    difference_prior(difference_prior),
    iter(iter),
    burnin(burnin),
    na_impute(na_impute),
    missing_data_indices(missing_data_indices),
    is_ordinal_variable(is_ordinal_variable),
    baseline_category(baseline_category),
    difference_selection(difference_selection),
    main_effect_indices(main_effect_indices),
    pairwise_effect_indices(pairwise_effect_indices),
    target_accept(target_accept),
    nuts_max_depth(nuts_max_depth),
    learn_mass_matrix(learn_mass_matrix),
    projection(projection),
    group_membership(group_membership),
    group_indices(group_indices),
    interaction_index_matrix(interaction_index_matrix),
    inclusion_probability_master(inclusion_probability_master),
    chain_rngs(chain_rngs),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      ChainResultCompare out;
      out.chain_id = static_cast<int>(i + 1);
      out.error = false;

      try {
        // per-chain RNG
        SafeRNG rng(chain_rngs[i]);

        // make per-chain copies
        std::vector<arma::imat> num_obs_categories = num_obs_categories_master;
        std::vector<arma::imat> sufficient_blume_capel = sufficient_blume_capel_master;
        std::vector<arma::mat>  sufficient_pairwise = sufficient_pairwise_master;
        arma::mat inclusion_probability = inclusion_probability_master;
        arma::imat observations_copy = observations;

        // run sampler (pure C++)
        SamplerOutput result = run_gibbs_sampler_for_bgmCompare(
          out.chain_id,
          observations_copy,
          num_groups,
          num_obs_categories,
          sufficient_blume_capel,
          sufficient_pairwise,
          num_categories,
          main_alpha,
          main_beta,
          pairwise_scale,
          difference_scale,
          difference_selection_alpha,
          difference_selection_beta,
          difference_prior,
          iter,
          burnin,
          na_impute,
          missing_data_indices,
          is_ordinal_variable,
          baseline_category,
          difference_selection,
          main_effect_indices,
          pairwise_effect_indices,
          target_accept,
          nuts_max_depth,
          learn_mass_matrix,
          projection,
          group_membership,
          group_indices,
          interaction_index_matrix,
          inclusion_probability,
          rng
        );

        out.result = result;

      } catch (std::exception& e) {
        out.error = true;
        out.error_msg = e.what();
      } catch (...) {
        out.error = true;
        out.error_msg = "Unknown error";
      }

      results[i] = out;
    }
  }
};

// -----------------------------------------------------------------------------
// Main entry point
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List run_bgmCompare_parallel(
    const arma::imat& observations,
    int num_groups,
    const std::vector<arma::imat>& num_obs_categories,
    const std::vector<arma::imat>& sufficient_blume_capel,
    const std::vector<arma::mat>&  sufficient_pairwise,
    const arma::ivec& num_categories,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    double difference_scale,
    double difference_selection_alpha,
    double difference_selection_beta,
    const std::string& difference_prior,
    int iter,
    int burnin,
    bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    const arma::mat& inclusion_probability,
    int num_chains,
    int nThreads,
    int seed
) {
  std::vector<ChainResultCompare> results(num_chains);

  // per-chain seeds
  std::vector<SafeRNG> chain_rngs(num_chains);
  for (int c = 0; c < num_chains; ++c) {
    chain_rngs[c] = SafeRNG(seed + c);
  }


  GibbsCompareChainRunner worker(
      observations, num_groups,
      num_obs_categories, sufficient_blume_capel, sufficient_pairwise,
      num_categories, main_alpha, main_beta, pairwise_scale, difference_scale,
      difference_selection_alpha, difference_selection_beta, difference_prior,
      iter, burnin, na_impute, missing_data_indices, is_ordinal_variable,
      baseline_category, difference_selection, main_effect_indices,
      pairwise_effect_indices, target_accept, nuts_max_depth, learn_mass_matrix,
      projection, group_membership, group_indices, interaction_index_matrix,
      inclusion_probability, chain_rngs, results
  );

  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, num_chains, worker);
  }

  // wrap results back into Rcpp::List
  Rcpp::List output(num_chains);
  for (int i = 0; i < num_chains; ++i) {
    if (results[i].error) {
      output[i] = Rcpp::List::create(
        Rcpp::Named("error") = results[i].error_msg,
        Rcpp::Named("chain_id") = results[i].chain_id
      );
    } else {
      const auto& r = results[i].result;
      Rcpp::List chain_out = Rcpp::List::create(
        Rcpp::Named("main_samples") = r.main_samples,
        Rcpp::Named("pairwise_samples") = r.pairwise_samples,
        Rcpp::Named("treedepth__") = r.treedepth_samples,
        Rcpp::Named("divergent__") = r.divergent_samples,
        Rcpp::Named("energy__") = r.energy_samples,
        Rcpp::Named("chain_id") = r.chain_id
      );
      if (r.has_indicator) {
        chain_out["indicator_samples"] = r.indicator_samples;
      }
      output[i] = chain_out;
    }
  }

  return output;
}
