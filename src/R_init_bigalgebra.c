//
// 2021-04-10. Author F. Bertrand <frederic.bertrand@utt.fr>
//

#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "bigalgebra.h" // for NULL

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[]  = {
  {"_dgemm_wrapper", (DL_FUNC) &dgemm_wrapper, 17},
  {"_daxpy_wrapper", (DL_FUNC) &daxpy_wrapper, 5},
  {"_dpotrf_wrapper", (DL_FUNC) &dpotrf_wrapper, 6},
  {"_dcopy_wrapper", (DL_FUNC) &dcopy_wrapper, 7},
  {"_dscal_wrapper", (DL_FUNC) &dscal_wrapper, 5},
  {"_dgeqrf_wrapper", (DL_FUNC) &dgeqrf_wrapper, 11},
  {"_dgeev_wrapper", (DL_FUNC) &dgeev_wrapper, 20},
  {"_dgesdd_wrapper", (DL_FUNC) &dgesdd_wrapper, 18},
  {NULL, NULL, 0}
};

void R_init_bigalgebra(DllInfo * info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}


