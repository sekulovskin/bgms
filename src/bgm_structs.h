#pragma once
#include <RcppArmadillo.h>

// Data for bgm
struct DataBgm {
  arma::imat observations;
  arma::imat counts_per_category;
  arma::imat blume_capel_stats;
  arma::imat pairwise_stats;
  arma::imat interaction_index_matrix;
  arma::imat missing_index;
};
