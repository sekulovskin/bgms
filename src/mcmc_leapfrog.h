// mcmc_leapfrog.h
#ifndef MCMC_LEAPFROG_H
#define MCMC_LEAPFROG_H

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc_memoization.h"



/**
 * Function: leapfrog
 *
 * Performs num_leapfrogs standard leapfrog integration steps using the given gradient function.
 */
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs = 1
);



/**
 * Function: leapfrog_memo
 *
 * Performs a leapfrog step using a memoization wrapper to avoid redundant gradient evaluations.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo
);

#endif // MCMC_LEAPFROG_H
