#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Forward declarations */
SEXP knnx_index_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP k_sexp);
SEXP calculate_residuals_cpp(SEXP R_j_sexp, SEXP indexes_sexp,
                             SEXP indexesStar_sexp, SEXP num_levels_old_sexp,
                             SEXP num_centres_new_sexp);
SEXP propose_tessellation_cpp(SEXP tess_j_sexp, SEXP dim_j_sexp, SEXP var_sexp, SEXP num_cov_sexp);

/* Registration table */
static const R_CallMethodDef CallEntries[] = {
  {"knnx_index_cpp",            (DL_FUNC) &knnx_index_cpp,            3},
  {"calculate_residuals_cpp",   (DL_FUNC) &calculate_residuals_cpp,   5},
  {"propose_tessellation_cpp",  (DL_FUNC) &propose_tessellation_cpp,  4},
  {NULL, NULL, 0}
};

void R_init_AddiVortes(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
