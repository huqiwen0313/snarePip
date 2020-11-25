// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// fast_multiply
Eigen::SparseMatrix<double> fast_multiply(Eigen::SparseMatrix<double> M1, Eigen::SparseMatrix<double> M2);
RcppExport SEXP _snarePip_fast_multiply(SEXP M1SEXP, SEXP M2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type M2(M2SEXP);
    rcpp_result_gen = Rcpp::wrap(fast_multiply(M1, M2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_snarePip_fast_multiply", (DL_FUNC) &_snarePip_fast_multiply, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_snarePip(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
