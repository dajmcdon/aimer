// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// findThresholdAIMER
Rcpp::List findThresholdAIMER(arma::mat X, arma::colvec y, arma::colvec ncomps, arma::colvec nCovs, int nCovsMin, int nCovsMax, int nthresh, int kfold);
RcppExport SEXP aimer_findThresholdAIMER(SEXP XSEXP, SEXP ySEXP, SEXP ncompsSEXP, SEXP nCovsSEXP, SEXP nCovsMinSEXP, SEXP nCovsMaxSEXP, SEXP nthreshSEXP, SEXP kfoldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type ncomps(ncompsSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nCovs(nCovsSEXP);
    Rcpp::traits::input_parameter< int >::type nCovsMin(nCovsMinSEXP);
    Rcpp::traits::input_parameter< int >::type nCovsMax(nCovsMaxSEXP);
    Rcpp::traits::input_parameter< int >::type nthresh(nthreshSEXP);
    Rcpp::traits::input_parameter< int >::type kfold(kfoldSEXP);
    rcpp_result_gen = Rcpp::wrap(findThresholdAIMER(X, y, ncomps, nCovs, nCovsMin, nCovsMax, nthresh, kfold));
    return rcpp_result_gen;
END_RCPP
}
// AIMER
arma::colvec AIMER(arma::mat X, arma::colvec y, double t, double b, int d);
RcppExport SEXP aimer_AIMER(SEXP XSEXP, SEXP ySEXP, SEXP tSEXP, SEXP bSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(AIMER(X, y, t, b, d));
    return rcpp_result_gen;
END_RCPP
}
