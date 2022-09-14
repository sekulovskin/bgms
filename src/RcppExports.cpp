// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_omrf_gibbs
IntegerMatrix sample_omrf_gibbs(int no_persons, int no_nodes, IntegerVector no_categories, NumericMatrix Interactions, NumericMatrix Thresholds, int no_iterations);
RcppExport SEXP _Bmrf_sample_omrf_gibbs(SEXP no_personsSEXP, SEXP no_nodesSEXP, SEXP no_categoriesSEXP, SEXP InteractionsSEXP, SEXP ThresholdsSEXP, SEXP no_iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_persons(no_personsSEXP);
    Rcpp::traits::input_parameter< int >::type no_nodes(no_nodesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Interactions(InteractionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Thresholds(ThresholdsSEXP);
    Rcpp::traits::input_parameter< int >::type no_iterations(no_iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_omrf_gibbs(no_persons, no_nodes, no_categories, Interactions, Thresholds, no_iterations));
    return rcpp_result_gen;
END_RCPP
}
// em_gamma
NumericVector em_gamma(NumericMatrix interactions, NumericMatrix slab_var, double theta, double xi, int no_persons);
RcppExport SEXP _Bmrf_em_gamma(SEXP interactionsSEXP, SEXP slab_varSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP no_personsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type slab_var(slab_varSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type no_persons(no_personsSEXP);
    rcpp_result_gen = Rcpp::wrap(em_gamma(interactions, slab_var, theta, xi, no_persons));
    return rcpp_result_gen;
END_RCPP
}
// em_interaction_var
NumericVector em_interaction_var(NumericMatrix gamma, NumericMatrix slab_var, double theta, double xi, int no_persons);
RcppExport SEXP _Bmrf_em_interaction_var(SEXP gammaSEXP, SEXP slab_varSEXP, SEXP thetaSEXP, SEXP xiSEXP, SEXP no_personsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type slab_var(slab_varSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type no_persons(no_personsSEXP);
    rcpp_result_gen = Rcpp::wrap(em_interaction_var(gamma, slab_var, theta, xi, no_persons));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_samples
List gibbs_samples(IntegerMatrix observations, IntegerVector no_categories, String interaction_prior, double cauchy_scale, NumericMatrix unit_info, NumericMatrix proposal_sd, IntegerMatrix Index, int no_iterations, IntegerMatrix n_cat_obs, double a, double b, bool display_progress);
RcppExport SEXP _Bmrf_gibbs_samples(SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP interaction_priorSEXP, SEXP cauchy_scaleSEXP, SEXP unit_infoSEXP, SEXP proposal_sdSEXP, SEXP IndexSEXP, SEXP no_iterationsSEXP, SEXP n_cat_obsSEXP, SEXP aSEXP, SEXP bSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< String >::type interaction_prior(interaction_priorSEXP);
    Rcpp::traits::input_parameter< double >::type cauchy_scale(cauchy_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type unit_info(unit_infoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< int >::type no_iterations(no_iterationsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_cat_obs(n_cat_obsSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_samples(observations, no_categories, interaction_prior, cauchy_scale, unit_info, proposal_sd, Index, no_iterations, n_cat_obs, a, b, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_eap
List gibbs_eap(IntegerMatrix observations, IntegerVector no_categories, String interaction_prior, double cauchy_scale, NumericMatrix unit_info, NumericMatrix proposal_sd, IntegerMatrix Index, int no_iterations, IntegerMatrix n_cat_obs, double a, double b, bool display_progress);
RcppExport SEXP _Bmrf_gibbs_eap(SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP interaction_priorSEXP, SEXP cauchy_scaleSEXP, SEXP unit_infoSEXP, SEXP proposal_sdSEXP, SEXP IndexSEXP, SEXP no_iterationsSEXP, SEXP n_cat_obsSEXP, SEXP aSEXP, SEXP bSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< String >::type interaction_prior(interaction_priorSEXP);
    Rcpp::traits::input_parameter< double >::type cauchy_scale(cauchy_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type unit_info(unit_infoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< int >::type no_iterations(no_iterationsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_cat_obs(n_cat_obsSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_eap(observations, no_categories, interaction_prior, cauchy_scale, unit_info, proposal_sd, Index, no_iterations, n_cat_obs, a, b, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// gradient_interactions
NumericVector gradient_interactions(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, NumericMatrix interaction_var);
RcppExport SEXP _Bmrf_gradient_interactions(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP interaction_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interaction_var(interaction_varSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient_interactions(interactions, thresholds, observations, no_categories, interaction_var));
    return rcpp_result_gen;
END_RCPP
}
// gradient_interactions_cauchy
NumericVector gradient_interactions_cauchy(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, double cauchy_scale);
RcppExport SEXP _Bmrf_gradient_interactions_cauchy(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP cauchy_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< double >::type cauchy_scale(cauchy_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient_interactions_cauchy(interactions, thresholds, observations, no_categories, cauchy_scale));
    return rcpp_result_gen;
END_RCPP
}
// gradient_thresholds
NumericVector gradient_thresholds(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, NumericMatrix threshold_var);
RcppExport SEXP _Bmrf_gradient_thresholds(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP threshold_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold_var(threshold_varSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient_thresholds(interactions, thresholds, observations, no_categories, threshold_var));
    return rcpp_result_gen;
END_RCPP
}
// hessian_interactions
NumericMatrix hessian_interactions(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, NumericMatrix interaction_var);
RcppExport SEXP _Bmrf_hessian_interactions(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP interaction_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interaction_var(interaction_varSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_interactions(interactions, thresholds, observations, no_categories, interaction_var));
    return rcpp_result_gen;
END_RCPP
}
// hessian_interactions_cauchy
NumericMatrix hessian_interactions_cauchy(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, double cauchy_scale);
RcppExport SEXP _Bmrf_hessian_interactions_cauchy(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP cauchy_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< double >::type cauchy_scale(cauchy_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_interactions_cauchy(interactions, thresholds, observations, no_categories, cauchy_scale));
    return rcpp_result_gen;
END_RCPP
}
// hessian_thresholds
NumericMatrix hessian_thresholds(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, NumericMatrix threshold_var);
RcppExport SEXP _Bmrf_hessian_thresholds(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP threshold_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold_var(threshold_varSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_thresholds(interactions, thresholds, observations, no_categories, threshold_var));
    return rcpp_result_gen;
END_RCPP
}
// hessian_crossparameters
NumericMatrix hessian_crossparameters(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories);
RcppExport SEXP _Bmrf_hessian_crossparameters(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_crossparameters(interactions, thresholds, observations, no_categories));
    return rcpp_result_gen;
END_RCPP
}
// joint_log_density
double joint_log_density(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, NumericMatrix interaction_var, NumericMatrix threshold_var);
RcppExport SEXP _Bmrf_joint_log_density(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP interaction_varSEXP, SEXP threshold_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interaction_var(interaction_varSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold_var(threshold_varSEXP);
    rcpp_result_gen = Rcpp::wrap(joint_log_density(interactions, thresholds, observations, no_categories, interaction_var, threshold_var));
    return rcpp_result_gen;
END_RCPP
}
// joint_log_density_cauchy
double joint_log_density_cauchy(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, NumericMatrix threshold_var, double cauchy_scale, IntegerVector no_categories);
RcppExport SEXP _Bmrf_joint_log_density_cauchy(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP threshold_varSEXP, SEXP cauchy_scaleSEXP, SEXP no_categoriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold_var(threshold_varSEXP);
    Rcpp::traits::input_parameter< double >::type cauchy_scale(cauchy_scaleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    rcpp_result_gen = Rcpp::wrap(joint_log_density_cauchy(interactions, thresholds, observations, threshold_var, cauchy_scale, no_categories));
    return rcpp_result_gen;
END_RCPP
}
// emvs_joint_log_density
double emvs_joint_log_density(NumericMatrix interactions, NumericMatrix thresholds, IntegerMatrix observations, IntegerVector no_categories, double xi, NumericMatrix slab_var, NumericMatrix threshold_var, double theta, bool hierarchical, double alpha, double beta);
RcppExport SEXP _Bmrf_emvs_joint_log_density(SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP observationsSEXP, SEXP no_categoriesSEXP, SEXP xiSEXP, SEXP slab_varSEXP, SEXP threshold_varSEXP, SEXP thetaSEXP, SEXP hierarchicalSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type slab_var(slab_varSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold_var(threshold_varSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type hierarchical(hierarchicalSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(emvs_joint_log_density(interactions, thresholds, observations, no_categories, xi, slab_var, threshold_var, theta, hierarchical, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Bmrf_sample_omrf_gibbs", (DL_FUNC) &_Bmrf_sample_omrf_gibbs, 6},
    {"_Bmrf_em_gamma", (DL_FUNC) &_Bmrf_em_gamma, 5},
    {"_Bmrf_em_interaction_var", (DL_FUNC) &_Bmrf_em_interaction_var, 5},
    {"_Bmrf_gibbs_samples", (DL_FUNC) &_Bmrf_gibbs_samples, 12},
    {"_Bmrf_gibbs_eap", (DL_FUNC) &_Bmrf_gibbs_eap, 12},
    {"_Bmrf_gradient_interactions", (DL_FUNC) &_Bmrf_gradient_interactions, 5},
    {"_Bmrf_gradient_interactions_cauchy", (DL_FUNC) &_Bmrf_gradient_interactions_cauchy, 5},
    {"_Bmrf_gradient_thresholds", (DL_FUNC) &_Bmrf_gradient_thresholds, 5},
    {"_Bmrf_hessian_interactions", (DL_FUNC) &_Bmrf_hessian_interactions, 5},
    {"_Bmrf_hessian_interactions_cauchy", (DL_FUNC) &_Bmrf_hessian_interactions_cauchy, 5},
    {"_Bmrf_hessian_thresholds", (DL_FUNC) &_Bmrf_hessian_thresholds, 5},
    {"_Bmrf_hessian_crossparameters", (DL_FUNC) &_Bmrf_hessian_crossparameters, 4},
    {"_Bmrf_joint_log_density", (DL_FUNC) &_Bmrf_joint_log_density, 6},
    {"_Bmrf_joint_log_density_cauchy", (DL_FUNC) &_Bmrf_joint_log_density_cauchy, 6},
    {"_Bmrf_emvs_joint_log_density", (DL_FUNC) &_Bmrf_emvs_joint_log_density, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_Bmrf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}