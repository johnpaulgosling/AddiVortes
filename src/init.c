#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP knnx_index_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP addi_vortes_mcmc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predict_sample_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"knnx_index_cpp",             (DL_FUNC) &knnx_index_cpp,             5},
  {"addi_vortes_mcmc_cpp",       (DL_FUNC) &addi_vortes_mcmc_cpp,      20},
  {"predict_sample_cpp",         (DL_FUNC) &predict_sample_cpp,         5},
  {NULL, NULL, 0}
};

void R_init_AddiVortes(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
