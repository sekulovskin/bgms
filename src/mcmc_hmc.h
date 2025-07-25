#ifndef MCMC_HMC_H
#define MCMC_HMC_H

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc_utils.h"



SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs
);

#endif // MCMC_HMC.H