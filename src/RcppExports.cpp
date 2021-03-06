// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dgeShuffle
List dgeShuffle(int n_genes, int n_bcds, int num_nonzeros, NumericVector numis, NumericVector geneInd, NumericVector barcodeInd, int totalumi);
RcppExport SEXP _SiftCell_dgeShuffle(SEXP n_genesSEXP, SEXP n_bcdsSEXP, SEXP num_nonzerosSEXP, SEXP numisSEXP, SEXP geneIndSEXP, SEXP barcodeIndSEXP, SEXP totalumiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_genes(n_genesSEXP);
    Rcpp::traits::input_parameter< int >::type n_bcds(n_bcdsSEXP);
    Rcpp::traits::input_parameter< int >::type num_nonzeros(num_nonzerosSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type numis(numisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type geneInd(geneIndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type barcodeInd(barcodeIndSEXP);
    Rcpp::traits::input_parameter< int >::type totalumi(totalumiSEXP);
    rcpp_result_gen = Rcpp::wrap(dgeShuffle(n_genes, n_bcds, num_nonzeros, numis, geneInd, barcodeInd, totalumi));
    return rcpp_result_gen;
END_RCPP
}
// runDMM
NumericMatrix runDMM(NumericMatrix geneProfile, NumericMatrix mtx, NumericVector ub_p, NumericVector lb_p);
RcppExport SEXP _SiftCell_runDMM(SEXP geneProfileSEXP, SEXP mtxSEXP, SEXP ub_pSEXP, SEXP lb_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type geneProfile(geneProfileSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ub_p(ub_pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lb_p(lb_pSEXP);
    rcpp_result_gen = Rcpp::wrap(runDMM(geneProfile, mtx, ub_p, lb_p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SiftCell_dgeShuffle", (DL_FUNC) &_SiftCell_dgeShuffle, 7},
    {"_SiftCell_runDMM", (DL_FUNC) &_SiftCell_runDMM, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SiftCell(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
