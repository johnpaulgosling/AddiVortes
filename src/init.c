#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP knnx_index_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP calculate_residuals_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP propose_tessellation_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"knnx_index_cpp",             (DL_FUNC) &knnx_index_cpp,             5},
  {"calculate_residuals_cpp",    (DL_FUNC) &calculate_residuals_cpp,    5},
  {"propose_tessellation_cpp",   (DL_FUNC) &propose_tessellation_cpp,   6},
  {NULL, NULL, 0}
};

void R_init_AddiVortes(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
