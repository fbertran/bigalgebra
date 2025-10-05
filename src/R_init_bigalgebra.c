//
// 2021-04-10. Author F. Bertrand <frederic.bertrand@lecnam.net>
//

#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "bigalgebra.h" // for NULL

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef callMethods[] = {
  {"_dgeev_lwork_query_wrapper", (DL_FUNC) &dgeev_lwork_query_wrapper, 3},
  {"_dgemm_wrapper",  (DL_FUNC) &dgemm_wrapper,  17},
  {"_dadd_wrapper",  (DL_FUNC) &dadd_wrapper,   7},
  {"_daxpy_wrapper",  (DL_FUNC) &daxpy_wrapper,   5},
  {"_dpotrf_wrapper", (DL_FUNC) &dpotrf_wrapper,  6},
  {"_dcopy_wrapper",  (DL_FUNC) &dcopy_wrapper,   7},
  {"_dswap_wrapper",  (DL_FUNC) &dswap_wrapper,   7},
  {"_dscal_wrapper",  (DL_FUNC) &dscal_wrapper,   5},
  {"_dset_wrapper",   (DL_FUNC) &dset_wrapper,    5},
  {"_dvcal_wrapper",  (DL_FUNC) &dvcal_wrapper,   9},
  {"_dsub_wrapper",   (DL_FUNC) &dsub_wrapper,    7},
  {"_ddot_wrapper",   (DL_FUNC) &ddot_wrapper,    7},
  {"_dqddot_wrapper", (DL_FUNC) &dqddot_wrapper,  7},
  {"_dhprod_wrapper", (DL_FUNC) &dhprod_wrapper, 10},
  {"_dxyz_wrapper",   (DL_FUNC) &dxyz_wrapper,    6},
  {"_dsum_wrapper",   (DL_FUNC) &dsum_wrapper,    4},
  {"_dasum_wrapper",  (DL_FUNC) &dasum_wrapper,   4},
  {"_dnrm2_wrapper",  (DL_FUNC) &dnrm2_wrapper,   4},
  {"_dprdct_wrapper", (DL_FUNC) &dprdct_wrapper,  4},
  {"_idmin_wrapper",  (DL_FUNC) &idmin_wrapper,   4},
  {"_idmax_wrapper",  (DL_FUNC) &idmax_wrapper,   4},
  {"_idamin_wrapper", (DL_FUNC) &idamin_wrapper,  4},
  {"_idamax_wrapper", (DL_FUNC) &idamax_wrapper,  4},
  {"_dsymm_wrapper",  (DL_FUNC) &dsymm_wrapper,  15},
  {"_dgeqrf_wrapper", (DL_FUNC) &dgeqrf_wrapper, 11},
  {"_dgeev_wrapper",  (DL_FUNC) &dgeev_wrapper,  20},
  {"_dgesdd_wrapper", (DL_FUNC) &dgesdd_wrapper, 18},
  {NULL, NULL, 0}
};

#ifdef __cplusplus
extern "C"
#endif
void R_init_bigalgebra(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}


