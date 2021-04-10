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
  {NULL, NULL, 0}
};

void R_init_bigalgebra(DllInfo * info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}


