#ifndef MCMC_RWM_H
#define MCMC_RWM_H

#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc_utils.h"
using namespace Rcpp;

SamplerResult rwm_sampler(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post
);

#endif // MCMC_RWM.H