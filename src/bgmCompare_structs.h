#pragma once
#include <RcppArmadillo.h>
#include <vector>

// Data for bgmCompare
struct DataBgmCompare {
  arma::imat observations;
  int num_groups;
  arma::ivec group_membership;
  arma::imat group_indices;
  arma::mat projection;
  std::vector<arma::imat> counts_per_category_group;
  std::vector<arma::imat> blume_capel_stats_group;
  std::vector<arma::mat> pairwise_stats_group;
  arma::imat missing_data_indices;
  arma::imat interaction_index_matrix;
};
