// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lmcpp
List lmcpp(const arma::mat& X, const arma::colvec& y, const arma::colvec& w);
RcppExport SEXP _blblm_lmcpp(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type w(wSEXP);
    rcpp_result = Rcpp::wrap(lmcpp(X, y, w));
    return rcpp_result;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_blblm_lmcpp", (DL_FUNC) &_blblm_lmcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_blblm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
