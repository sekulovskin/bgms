#pragma once

#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc_utils.h"
#include "rng_utils.h"

using namespace Rcpp;

SamplerResult rwm_sampler(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    dqrng::xoshiro256plus& rng
);