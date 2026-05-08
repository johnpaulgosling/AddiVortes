#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 1. Explicitly declare the C++ functions so the C compiler knows they exist */
extern SEXP super_mcmc_loop_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP knnx_index_predict_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP soft_predict_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);

/* 2. Map the R string names to the compiled C++ memory addresses and specify argument counts */
static const R_CallMethodDef CallEntries[] = {
  {"super_mcmc_loop_cpp", (DL_FUNC) &super_mcmc_loop_cpp, 37},
  {"knnx_index_predict_cpp", (DL_FUNC) &knnx_index_predict_cpp, 4},
  {"soft_predict_cpp", (DL_FUNC) &soft_predict_cpp, 5},
  {NULL, NULL, 0}
};

/* 3. Register the routines when the package loads */
void R_init_AddiVortes(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}