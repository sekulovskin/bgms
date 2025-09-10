#ifndef SAMPLEROUTPUT_H
#define SAMPLEROUTPUT_H

#include <RcppArmadillo.h>

// Plain C++ struct, no Rcpp types
struct SamplerOutput {
  arma::mat main_samples;
  arma::mat pairwise_samples;
  arma::imat indicator_samples;
  arma::ivec treedepth_samples;
  arma::ivec divergent_samples;
  arma::vec energy_samples;
  int chain_id;
  bool has_indicator;
};

#endif
