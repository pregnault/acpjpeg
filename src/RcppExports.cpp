// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// moy
double moy(NumericVector x);
RcppExport SEXP _acpjpeg_moy(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(moy(x));
    return rcpp_result_gen;
END_RCPP
}
// ect
double ect(NumericVector x);
RcppExport SEXP _acpjpeg_ect(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ect(x));
    return rcpp_result_gen;
END_RCPP
}
// moy2
double moy2(NumericVector x);
RcppExport SEXP _acpjpeg_moy2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(moy2(x));
    return rcpp_result_gen;
END_RCPP
}
// ect2
double ect2(NumericVector x);
RcppExport SEXP _acpjpeg_ect2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ect2(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _acpjpeg_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_acpjpeg_moy", (DL_FUNC) &_acpjpeg_moy, 1},
    {"_acpjpeg_ect", (DL_FUNC) &_acpjpeg_ect, 1},
    {"_acpjpeg_moy2", (DL_FUNC) &_acpjpeg_moy2, 1},
    {"_acpjpeg_ect2", (DL_FUNC) &_acpjpeg_ect2, 1},
    {"_acpjpeg_rcpp_hello_world", (DL_FUNC) &_acpjpeg_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_acpjpeg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
