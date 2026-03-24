#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> 
#include <R_ext/Rdynload.h>

extern SEXP super_mcmc_loop_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                                SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP knnx_index_predict_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"super_mcmc_loop_cpp", (DL_FUNC) &super_mcmc_loop_cpp, 21},
  {"knnx_index_predict_cpp", (DL_FUNC) &knnx_index_predict_cpp, 4},
  {NULL, NULL, 0}
};

void R_init_AddiVortes(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}