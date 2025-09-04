// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>
#include <tbb/mutex.h>

#include "bgmCompare_sampler.h"
#include "print_mutex.h"

using namespace Rcpp;
using namespace RcppParallel;


// -----------------------------------------------------------------------------
// Worker: runs chains [begin, end) in parallel.
// -----------------------------------------------------------------------------
struct GibbsCompareChainRunner : public Worker {
  // immutable shared inputs (masters)
  const arma::imat observations_master;
  const int num_groups;

  // preconverted, thread-shared masters (read-only in worker)
  const std::vector<arma::imat>& num_obs_categories_cpp_master;
  const std::vector<arma::imat>& sufficient_blume_capel_cpp_master;
  const std::vector<arma::mat>&  sufficient_pairwise_cpp_master;

  // other immutable inputs
  const arma::ivec& num_categories;
  const double main_alpha;
  const double main_beta;
  const double pairwise_scale;
  const double difference_scale;
  const double difference_selection_alpha;
  const double difference_selection_beta;
  const std::string difference_prior;
  const int iter;
  const int burnin;
  const bool na_impute;
  const arma::imat& missing_data_indices;
  const arma::uvec& is_ordinal_variable;
  const arma::ivec& baseline_category;
  const bool difference_selection;
  const arma::imat main_effect_indices;
  const arma::imat pairwise_effect_indices;
  const double target_accept;
  const int nuts_max_depth;
  const bool learn_mass_matrix;
  const arma::mat projection;
  const arma::ivec& group_membership;
  const arma::imat& group_indices;
  const arma::imat& interaction_index_matrix;

  // base inclusion prob (will be deep-copied per chain)
  const arma::mat& inclusion_probability_master;

  // shared output (one slot per chain index)
  Rcpp::List& output;

  GibbsCompareChainRunner(
    arma::imat observations_master_,
    const int num_groups_,
    const std::vector<arma::imat>& num_obs_categories_cpp_master_,
    const std::vector<arma::imat>& sufficient_blume_capel_cpp_master_,
    const std::vector<arma::mat>&  sufficient_pairwise_cpp_master_,
    const arma::ivec& num_categories_,
    const double main_alpha_,
    const double main_beta_,
    const double pairwise_scale_,
    const double difference_scale_,
    const double difference_selection_alpha_,
    const double difference_selection_beta_,
    const std::string difference_prior_,
    const int iter_,
    const int burnin_,
    const bool na_impute_,
    const arma::imat& missing_data_indices_,
    const arma::uvec& is_ordinal_variable_,
    const arma::ivec& baseline_category_,
    const bool difference_selection_,
    const arma::imat main_effect_indices_,
    const arma::imat pairwise_effect_indices_,
    const double target_accept_,
    const int nuts_max_depth_,
    const bool learn_mass_matrix_,
    const arma::mat projection_,
    const arma::ivec& group_membership_,
    const arma::imat& group_indices_,
    const arma::imat& interaction_index_matrix_,
    const arma::mat& inclusion_probability_master_,
    Rcpp::List& output_
  )
    : observations_master(observations_master_), num_groups(num_groups_),
      num_obs_categories_cpp_master(num_obs_categories_cpp_master_),
      sufficient_blume_capel_cpp_master(sufficient_blume_capel_cpp_master_),
      sufficient_pairwise_cpp_master(sufficient_pairwise_cpp_master_),
      num_categories(num_categories_),
      main_alpha(main_alpha_), main_beta(main_beta_),
      pairwise_scale(pairwise_scale_), difference_scale(difference_scale_),
      difference_selection_alpha(difference_selection_alpha_),
      difference_selection_beta(difference_selection_beta_),
      difference_prior(difference_prior_),
      iter(iter_), burnin(burnin_), na_impute(na_impute_),
      missing_data_indices(missing_data_indices_),
      is_ordinal_variable(is_ordinal_variable_),
      baseline_category(baseline_category_),
      difference_selection(difference_selection_),
      main_effect_indices(main_effect_indices_),
      pairwise_effect_indices(pairwise_effect_indices_),
      target_accept(target_accept_), nuts_max_depth(nuts_max_depth_),
      learn_mass_matrix(learn_mass_matrix_),
      projection(projection_), group_membership(group_membership_),
      group_indices(group_indices_),
      interaction_index_matrix(interaction_index_matrix_),
      inclusion_probability_master(inclusion_probability_master_),
      output(output_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      try {
        RNGScope scope;  // per-thread RNG scope

        // ---- per-chain DEEP COPIES of potentially mutated data ----
        std::vector<arma::imat> num_obs_categories_cpp = num_obs_categories_cpp_master;   // deep copy
        std::vector<arma::imat> sufficient_blume_capel_cpp = sufficient_blume_capel_cpp_master; // deep copy
        std::vector<arma::mat>  sufficient_pairwise_cpp = sufficient_pairwise_cpp_master; // deep copy
        arma::mat inclusion_probability = inclusion_probability_master;                    // deep copy
        arma::imat observations = observations_master;                                     // deep copy

        // Wrap locals into Lists (no Rcpp::as in threads)
        Rcpp::List num_obs_categories(num_groups);
        Rcpp::List sufficient_blume_capel(num_groups);
        Rcpp::List sufficient_pairwise(num_groups);
        for (int g = 0; g < num_groups; ++g) {
          num_obs_categories[g]     = num_obs_categories_cpp[g];
          sufficient_blume_capel[g] = sufficient_blume_capel_cpp[g];
          sufficient_pairwise[g]    = sufficient_pairwise_cpp[g];
        }

        // Run existing sampler
        Rcpp::List res = run_gibbs_sampler_for_bgmCompare(
          static_cast<int>(i) + 1, // chain_id (1-based)
          observations,
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
          inclusion_probability
        );

        output[static_cast<int>(i)] = res;

      } catch (std::exception& e) {
        output[static_cast<int>(i)] = Rcpp::List::create(
          _["error"] = std::string("chain ") + std::to_string(i+1) + ": " + e.what()
        );
      } catch (...) {
        output[static_cast<int>(i)] = Rcpp::List::create(
          _["error"] = std::string("chain ") + std::to_string(i+1) + ": unknown C++ exception"
        );
      }
    }
  }
};

// -----------------------------------------------------------------------------
// Exported entry point
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List run_bgmCompare_parallel(
    arma::imat observations,
    const int num_groups,
    Rcpp::List num_obs_categories,      // per-group (R)
    Rcpp::List sufficient_blume_capel,  // per-group (R)
    Rcpp::List sufficient_pairwise,     // per-group (R)
    const arma::ivec& num_categories,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const double difference_scale,
    const double difference_selection_alpha,
    const double difference_selection_beta,
    const std::string difference_prior,
    const int iter,
    const int burnin,
    const bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,
    const arma::imat main_effect_indices,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    const arma::mat& inclusion_probability,
    const int num_chains,
    const int nThreads
) {
  // --- Preconvert Lists to Armadillo ON THE MAIN THREAD (no Rcpp::as in workers) ---
  std::vector<arma::imat> num_obs_categories_cpp(num_groups);
  std::vector<arma::imat> sufficient_blume_capel_cpp(num_groups);
  std::vector<arma::mat>  sufficient_pairwise_cpp(num_groups);
  for (int g = 0; g < num_groups; ++g) {
    num_obs_categories_cpp[g]     = Rcpp::as<arma::imat>(num_obs_categories[g]);
    sufficient_blume_capel_cpp[g] = Rcpp::as<arma::imat>(sufficient_blume_capel[g]);
    sufficient_pairwise_cpp[g]    = Rcpp::as<arma::mat>(sufficient_pairwise[g]);
  }

  // Output placeholder (one slot per chain)
  Rcpp::List out(num_chains);

  // Build worker
  GibbsCompareChainRunner worker(
      observations, num_groups,
      num_obs_categories_cpp,
      sufficient_blume_capel_cpp,
      sufficient_pairwise_cpp,
      num_categories, main_alpha, main_beta,
      pairwise_scale, difference_scale,
      difference_selection_alpha, difference_selection_beta, difference_prior,
      iter, burnin, na_impute,
      missing_data_indices, is_ordinal_variable, baseline_category,
      difference_selection, main_effect_indices, pairwise_effect_indices,
      target_accept, nuts_max_depth, learn_mass_matrix,
      projection, group_membership, group_indices, interaction_index_matrix,
      inclusion_probability, out
  );

  // Control parallelism (TBB)
  {
    tbb::global_control ctrl(tbb::global_control::max_allowed_parallelism,
                             std::max(1, nThreads));
    parallelFor(0, static_cast<std::size_t>(num_chains), worker);
  }

  return out;
}