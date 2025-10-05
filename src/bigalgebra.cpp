/* Revised bigalgebra.cpp — rchk‑clean wrappers and C linkage
 *
 * Changes vs archive version:
 *  - make_double_ptr(): balanced PROTECT/UNPROTECT, early‑error safety
 *  - All .Call wrappers return SEXP (R_NilValue or an object) and use extern "C"
 *  - Added defensive NULL pointer checks after make_double_ptr()
 */

/*#define USE_FC_LEN_T*/
#define R_NO_REMAP  /* avoid macro remapping like length -> Rf_length that breaks C++ headers */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <cmath>

/* Portable helpers for hidden Fortran string lengths (FCONE). */
#ifndef FCLEN1
# ifdef FCONE
#  define FCLEN1 FCONE
#  define FCLEN2 FCONE FCONE
#  define FCLEN3 FCONE FCONE FCONE
# else
#  define FCLEN1
#  define FCLEN2
#  define FCLEN3
# endif
#endif

#include "bigalgebra.h"  /* must declare BigMatrix, index_type, LOGICAL_VALUE, etc. */
#include "bigmemory/BigMatrix.h"  /* needed for BigMatrix methods used below */


#ifdef REFBLAS
#include "refblas64longlong.h"
#define INT long long
#else
#define INT int
#endif

#ifndef NUMERIC_DATA
#define NUMERIC_DATA(x) REAL(x)
#endif
#ifndef DOUBLE_DATA
#define DOUBLE_DATA(x) REAL(x)
#endif
#ifndef CHARACTER_VALUE
#define CHARACTER_VALUE(x) CHAR(STRING_ELT((x), 0))
#endif
#ifndef LOGICAL_VALUE
#define LOGICAL_VALUE(x) (LOGICAL(x)[0])
#endif

/* ---------------------------------------------------------------------- */
/* Utility: returns double* for a BigMatrix or a regular R matrix.
 Ensures balanced PROTECT/UNPROTECT for rchk. */
extern "C" double* make_double_ptr(SEXP matrix, SEXP isBigMatrix)
{
  double* matrix_ptr = NULL;
  int nprot = 0;
  
  if ((Rboolean)LOGICAL_VALUE(isBigMatrix) == (Rboolean)TRUE) { /* BigMatrix case */
/* Cache "address" symbol (persistent once installed) */
static SEXP sym_address = NULL;
    if (sym_address == NULL)
      sym_address = Rf_install("address");
    
    SEXP address = R_do_slot(matrix, sym_address); /* replace GET_SLOT with R_do_slot to avoid dependency on Rdefines.h */
PROTECT(address); nprot++;

BigMatrix* pbm = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(address));
if (!pbm) {
  goto cleanup; /* matrix_ptr stays NULL */
}

/* Validate sub.big.matrix shape */
if (pbm->row_offset() > 0 && pbm->ncol() > 1) {
  UNPROTECT(nprot); /* must unprotect before Rf_error (long jump) */
Rf_error("sub.big.matrix objects cannot have row offset > 0 and number of columns > 1");
}

index_type offset = pbm->nrow() * pbm->col_offset();
matrix_ptr = reinterpret_cast<double*>(pbm->matrix()) + offset;

  } else { /* Regular R matrix case */
matrix_ptr = NUMERIC_DATA(matrix);
  }
  
  cleanup:
    if (nprot) UNPROTECT(nprot);
    return matrix_ptr;
}

/* ---------------------------------------------------------------------- */
/* .Call wrappers with C linkage                                               */
extern "C" {
  
  /* DGEMM: C := alpha * op(A) %*% op(B) + beta * C
   Args (17): TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC,
   A_isBM, B_isBM, C_isBM, C_offset */
  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                      SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
                      SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM,
                      SEXP C_offset)
  {
    const char *tA = Rf_translateChar(Rf_asChar(TRANSA));
    const char *tB = Rf_translateChar(Rf_asChar(TRANSB));
    
    INT MM    = (INT) Rf_asInteger(M);
    INT NN    = (INT) Rf_asInteger(N);
    INT KK    = (INT) Rf_asInteger(K);
    INT LDAA  = (INT) Rf_asInteger(LDA);
    INT LDBB  = (INT) Rf_asInteger(LDB);
    INT LDCC  = (INT) Rf_asInteger(LDC);
    double alpha = Rf_asReal(ALPHA);
    double beta  = Rf_asReal(BETA);
    
    double *pA = make_double_ptr(A, A_isBM);
    double *pB = make_double_ptr(B, B_isBM);
    double *pC = make_double_ptr(C, C_isBM);
    if (!pA) Rf_error("NULL pointer for A");
    if (!pB) Rf_error("NULL pointer for B");
    if (!pC) Rf_error("NULL pointer for C");
    
    /* Optional starting column offset in C (as seen in archive code) */
    INT j = (INT) Rf_asInteger(C_offset);
    pC += j;
    
#ifdef REFBLAS
    /* Standard Fortran interface without underscoring */
    int8_dgemm (tA, tB, &MM, &NN, &KK, &alpha, pA, &LDAA, pB, &LDBB,
                &beta, pC, &LDCC);
#else
    F77_CALL(dgemm) ((char*)tA, (char*)tB, &MM, &NN, &KK, &alpha, pA, &LDAA,
             pB, &LDBB, &beta, pC, &LDCC FCLEN2);
    
    
    
#endif
    return C;  // return the modified output matrix
  }
  
  /* DAXPY: Y := A * X + Y (returns modified Y for convenience) */
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM)
  {
    SEXP ans, Tr;
    double *pY;
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    
    INT incx = 1;
    INT incy = 1;
    INT NN = (INT) Rf_asInteger(N);
    
    PROTECT(ans = Y);
    PROTECT(Tr = Rf_allocVector(LGLSXP, 1));
    LOGICAL(Tr)[0] = 1; /* treat Y as BigMatrix in make_double_ptr */
  pY = make_double_ptr(Y, Tr);
  if (!pY) { UNPROTECT(2); Rf_error("NULL pointer for Y"); }
  
#ifdef REFBLAS
  /* Standard Fortran interface without underscoring */
  double alpha = Rf_asReal(A);
  int8_daxpy (NN, &alpha, pX, &incx, pY, &incy);
#else
  double alpha = Rf_asReal(A);
  F77_CALL(daxpy) (&NN, &alpha, pX, &incx, pY, &incy);
#endif
  
  UNPROTECT(2);
  return ans; /* return Y */
  }
  
  /* DPOTRF: Cholesky factorization (returns A) */
  SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM)
  {
    const char *UUPLO = Rf_translateChar(Rf_asChar(UPLO));
    INT NN   = (INT) Rf_asInteger(N);
    INT LLDA = (INT) Rf_asInteger(LDA);
    INT IINFO= (INT) Rf_asInteger(INFO);
    
    double *AA = make_double_ptr(A, A_isBM);
    if (!AA) Rf_error("NULL pointer for A");
    
#ifdef REFBLAS
    /* Standard Fortran interface without underscoring */
    int8_dpotrf (UUPLO, &NN, AA, &LLDA, &IINFO);
#else
    /* Fortran interface via F77_CALL */
    F77_CALL(dpotrf) ((char*)UUPLO, &NN, AA, &LLDA, &IINFO FCLEN1);
    
    
#endif
    return A;
  }
  
  /* DCOPY: Y := X (in‑place) */
  SEXP dcopy_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                      SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr (X, X_isBM);
    double *pY = make_double_ptr (Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);

#ifdef ACMLBLAS
    dcopy (NN, pX, INCXX, pY, INCYY);
#else
    F77_CALL(dcopy) (&NN, pX, &INCXX, pY, &INCYY);
#endif
    return R_NilValue;
  }

  /* DSWAP: swap X and Y in place */
  SEXP dswap_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                      SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);

    if (NN <= 0) {
      return R_NilValue;
    }

#ifdef ACMLBLAS
    dswap (NN, pX, INCXX, pY, INCYY);
#elif defined(REFBLAS)
    int8_dswap(&NN, pX, &INCXX, pY, &INCYY);
#else
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    for (INT i = 0; i < NN; ++i) {
      double tmp = pX[ix];
      pX[ix] = pY[iy];
      pY[iy] = tmp;
      ix += INCXX;
      iy += INCYY;
    }
#endif

    return R_NilValue;
  }

  /* DADD: Y := X + Y (in-place, alpha = 1) */
  SEXP dadd_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);
    double alpha = 1.0;

#ifdef REFBLAS
    int8_daxpy(&NN, &alpha, pX, &INCXX, pY, &INCYY);
#else
    F77_CALL(daxpy)(&NN, &alpha, pX, &INCXX, pY, &INCYY);
#endif

    return Y;
  }

  /* DSCAL: Y := alpha * Y (in-place) */
  SEXP dscal_wrapper (SEXP N, SEXP ALPHA, SEXP Y, SEXP INCY, SEXP Y_isBM)
  {
    double *pY = make_double_ptr (Y, Y_isBM);
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCYY = (INT) Rf_asInteger(INCY);
    double alpha = Rf_asReal(ALPHA);

#ifdef ACMLBLAS
    dscal (NN, alpha, pY, INCYY);
#else
    F77_CALL(dscal) (&NN, &alpha, pY, &INCYY);
#endif
    /* Return the modified object */
    return Y;
  }

  /* DSET: X := alpha (in-place) */
  SEXP dset_wrapper (SEXP N, SEXP ALPHA, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    double alpha = Rf_asReal(ALPHA);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    for (INT i = 0; i < NN; ++i) {
      pX[ix] = alpha;
      ix += INCXX;
    }

    return X;
  }

  /* DVCAL: Y := alpha * X + beta * Y */
  SEXP dvcal_wrapper (SEXP N, SEXP ALPHA, SEXP X, SEXP INCX,
                      SEXP BETA, SEXP Y, SEXP INCY,
                      SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);
    double alpha = Rf_asReal(ALPHA);
    double beta  = Rf_asReal(BETA);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    for (INT i = 0; i < NN; ++i) {
      pY[iy] = alpha * pX[ix] + beta * pY[iy];
      ix += INCXX;
      iy += INCYY;
    }

    return Y;
  }

  /* DSUB: Y := Y - X */
  SEXP dsub_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    for (INT i = 0; i < NN; ++i) {
      pY[iy] -= pX[ix];
      ix += INCXX;
      iy += INCYY;
    }

    return Y;
  }

  /* DDOT: dot product */
  SEXP ddot_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);

#ifdef REFBLAS
    double ans = 0.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    for (INT i = 0; i < NN; ++i) {
      ans += pX[ix] * pY[iy];
      ix += INCXX;
      iy += INCYY;
    }
#else
    double ans = F77_CALL(ddot)(&NN, pX, &INCXX, pY, &INCYY);
#endif

    return Rf_ScalarReal(ans);
  }

  /* DQDDOT: dot product with long double accumulation */
  SEXP dqddot_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                       SEXP X_isBM, SEXP Y_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);

    long double acc = 0.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    for (INT i = 0; i < NN; ++i) {
      acc += static_cast<long double>(pX[ix]) * static_cast<long double>(pY[iy]);
      ix += INCXX;
      iy += INCYY;
    }

    return Rf_ScalarReal(static_cast<double>(acc));
  }

  /* DHPROD: Z := X * Y (Hadamard product) */
  SEXP dhprod_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY, SEXP Z,
                       SEXP INCZ, SEXP X_isBM, SEXP Y_isBM, SEXP Z_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    double *pZ = make_double_ptr(Z, Z_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");
    if (!pZ) Rf_error("NULL pointer for Z");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT INCYY = (INT) Rf_asInteger(INCY);
    INT INCZZ = (INT) Rf_asInteger(INCZ);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    INT iy = (INCYY >= 0) ? 0 : (1 - NN) * INCYY;
    INT iz = (INCZZ >= 0) ? 0 : (1 - NN) * INCZZ;
    for (INT i = 0; i < NN; ++i) {
      pZ[iz] = pX[ix] * pY[iy];
      ix += INCXX;
      iy += INCYY;
      iz += INCZZ;
    }

    return Z;
  }

  /* DXYZ: Z := X x Y (cross product, length 3) */
  SEXP dxyz_wrapper (SEXP X, SEXP Y, SEXP Z,
                     SEXP X_isBM, SEXP Y_isBM, SEXP Z_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    double *pY = make_double_ptr(Y, Y_isBM);
    double *pZ = make_double_ptr(Z, Z_isBM);
    if (!pX) Rf_error("NULL pointer for X");
    if (!pY) Rf_error("NULL pointer for Y");
    if (!pZ) Rf_error("NULL pointer for Z");

    pZ[0] = pX[1] * pY[2] - pX[2] * pY[1];
    pZ[1] = pX[2] * pY[0] - pX[0] * pY[2];
    pZ[2] = pX[0] * pY[1] - pX[1] * pY[0];

    return Z;
  }

  /* DSUM: sum of elements */
  SEXP dsum_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

    long double acc = 0.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    for (INT i = 0; i < NN; ++i) {
      acc += static_cast<long double>(pX[ix]);
      ix += INCXX;
    }

    return Rf_ScalarReal(static_cast<double>(acc));
  }

  /* DASUM: sum of absolute values */
  SEXP dasum_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

#ifdef REFBLAS
    long double acc = 0.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    for (INT i = 0; i < NN; ++i) {
      acc += std::fabs(pX[ix]);
      ix += INCXX;
    }
    double ans = static_cast<double>(acc);
#else
    double ans = F77_CALL(dasum)(&NN, pX, &INCXX);
#endif

    return Rf_ScalarReal(ans);
  }

  /* DNRM2: Euclidean norm */
  SEXP dnrm2_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

#ifdef REFBLAS
    long double acc = 0.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    for (INT i = 0; i < NN; ++i) {
      long double val = static_cast<long double>(pX[ix]);
      acc += val * val;
      ix += INCXX;
    }
    double ans = std::sqrt(static_cast<double>(acc));
#else
    double ans = F77_CALL(dnrm2)(&NN, pX, &INCXX);
#endif

    return Rf_ScalarReal(ans);
  }

  /* DPRDCT: product of elements */
  SEXP dprdct_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

    long double prod = 1.0;
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    for (INT i = 0; i < NN; ++i) {
      prod *= static_cast<long double>(pX[ix]);
      ix += INCXX;
    }

    return Rf_ScalarReal(static_cast<double>(prod));
  }

  /* IDMIN: index of minimum element */
  SEXP idmin_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    double min_val = pX[ix];
    INT min_idx = 0;
    INT idx = ix;
    for (INT i = 1; i < NN; ++i) {
      idx += INCXX;
      double val = pX[idx];
      if (val < min_val) {
        min_val = val;
        min_idx = i;
      }
    }

    return Rf_ScalarInteger((int)(min_idx + 1));
  }

  /* IDMAX: index of maximum element */
  SEXP idmax_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);

    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    double max_val = pX[ix];
    INT max_idx = 0;
    INT idx = ix;
    for (INT i = 1; i < NN; ++i) {
      idx += INCXX;
      double val = pX[idx];
      if (val > max_val) {
        max_val = val;
        max_idx = i;
      }
    }

    return Rf_ScalarInteger((int)(max_idx + 1));
  }

  /* IDAMIN: index of minimum absolute value */
  SEXP idamin_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    double min_val = std::fabs(pX[ix]);
    INT min_idx = 0;
    INT idx = ix;
    for (INT i = 1; i < NN; ++i) {
      idx += INCXX;
      double val = std::fabs(pX[idx]);
      if (val < min_val) {
        min_val = val;
        min_idx = i;
      }
    }
    return Rf_ScalarInteger((int)(min_idx + 1));
  }

  /* IDAMAX: index of maximum absolute value */
  SEXP idamax_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM)
  {
    double *pX = make_double_ptr(X, X_isBM);
    if (!pX) Rf_error("NULL pointer for X");

    INT NN    = (INT) Rf_asInteger(N);
    INT INCXX = (INT) Rf_asInteger(INCX);
    INT ix = (INCXX >= 0) ? 0 : (1 - NN) * INCXX;
    double max_val = std::fabs(pX[ix]);
    INT max_idx = 0;
    INT idx = ix;
    for (INT i = 1; i < NN; ++i) {
      idx += INCXX;
      double val = std::fabs(pX[idx]);
      if (val > max_val) {
        max_val = val;
        max_idx = i;
      }
    }
    return Rf_ScalarInteger((int)(max_idx + 1));
  }

  /* DSYMM: symmetric matrix multiply */
  SEXP dsymm_wrapper (SEXP SIDE, SEXP UPLO, SEXP M, SEXP N, SEXP ALPHA,
                      SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
                      SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM)
  {
    const char *side = Rf_translateChar(Rf_asChar(SIDE));
    const char *uplo = Rf_translateChar(Rf_asChar(UPLO));

    double *pA = make_double_ptr(A, A_isBM);
    double *pB = make_double_ptr(B, B_isBM);
    double *pC = make_double_ptr(C, C_isBM);
    if (!pA) Rf_error("NULL pointer for A");
    if (!pB) Rf_error("NULL pointer for B");
    if (!pC) Rf_error("NULL pointer for C");

    INT MM    = (INT) Rf_asInteger(M);
    INT NN    = (INT) Rf_asInteger(N);
    INT LDAA  = (INT) Rf_asInteger(LDA);
    INT LDBB  = (INT) Rf_asInteger(LDB);
    INT LDCC  = (INT) Rf_asInteger(LDC);
    double alpha = Rf_asReal(ALPHA);
    double beta  = Rf_asReal(BETA);

#ifdef REFBLAS
    INT dimA = (*side == 'L' || *side == 'l') ? MM : NN;
    double *A_full = (double *) R_alloc((size_t)dimA * (size_t)dimA, sizeof(double));
    for (INT j = 0; j < dimA; ++j) {
      for (INT i = 0; i < dimA; ++i) {
        double val;
        if (*uplo == 'U' || *uplo == 'u') {
          val = (i <= j) ? pA[i + j * LDAA] : pA[j + i * LDAA];
        } else {
          val = (i >= j) ? pA[i + j * LDAA] : pA[j + i * LDAA];
        }
        A_full[i + j * dimA] = val;
      }
    }
    char trans = 'N';
    if (*side == 'L' || *side == 'l') {
      F77_CALL(dgemm)(&trans, &trans, &MM, &NN, &dimA, &alpha,
                      A_full, &dimA, pB, &LDBB, &beta, pC, &LDCC);
    } else {
      F77_CALL(dgemm)(&trans, &trans, &MM, &NN, &dimA, &alpha,
                      pB, &LDBB, A_full, &dimA, &beta, pC, &LDCC);
    }
#else
    F77_CALL(dsymm)((char*)side, (char*)uplo, &MM, &NN, &alpha,
                    pA, &LDAA, pB, &LDBB, &beta, pC, &LDCC FCLEN2);
#endif

    return C;
  }
  
  /* DGEQRF: QR factorization (in‑place) */
  SEXP dgeqrf_wrapper (SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP TAU, SEXP WORK,
                       SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM,
                       SEXP WORK_isBM)
  {
    double *pA   = make_double_ptr (A,   A_isBM);
    double *pTAU = make_double_ptr (TAU, TAU_isBM);
    double *pWORK= make_double_ptr (WORK,WORK_isBM);
    if (!pA)    Rf_error("NULL pointer for A");
    if (!pTAU)  Rf_error("NULL pointer for TAU");
    if (!pWORK) Rf_error("NULL pointer for WORK");
    
    INT MM    = (INT) Rf_asInteger(M);
    INT NN    = (INT) Rf_asInteger(N);
    INT LDAA  = (INT) Rf_asInteger(LDA);
    INT LWORKK= (INT) Rf_asInteger(LWORK);
    INT INFOO = (INT) Rf_asInteger(INFO);
    
#ifdef REFBLAS
    dgeqrf (MM, NN, pA, LDAA, pTAU, &INFOO); /* placeholder if REFBLAS provides it */
#else
  F77_CALL(dgeqrf) (&MM, &NN, pA, &LDAA, pTAU, pWORK, &LWORKK, &INFOO);
#endif
  return R_NilValue;
  }
  
  
  /* --- dgeev workspace query: returns optimal LWORK as an R integer --- */
  SEXP dgeev_lwork_query_wrapper(SEXP JOBVL, SEXP JOBVR, SEXP N)
  {
    const char *jvl = Rf_translateChar(Rf_asChar(JOBVL));
    const char *jvr = Rf_translateChar(Rf_asChar(JOBVR));
    int n    = Rf_asInteger(N);
    int lda  = (n > 0 ? n : 1);
    int ldvl = (*jvl == 'V' ? n : 1);
    int ldvr = (*jvr == 'V' ? n : 1);
    int lwork = -1;
    int info  = 0;
  
    /* Minimal dummy arrays; with LWORK=-1 dgeev will not use them */
    SEXP A   = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)lda * n));
    SEXP WR  = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP WI  = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP VL  = PROTECT(Rf_allocMatrix(REALSXP, ldvl, n));
    SEXP VR  = PROTECT(Rf_allocMatrix(REALSXP, ldvr, n));
    SEXP WORK= PROTECT(Rf_allocVector(REALSXP, 1));
  
    double *pA   = REAL(A);
    double *pWR  = REAL(WR);
    double *pWI  = REAL(WI);
    double *pVL  = REAL(VL);
    double *pVR  = REAL(VR);
    double *pWORK= REAL(WORK);
  
    F77_CALL(dgeev)((char*)jvl, (char*)jvr, &n, pA, &lda, pWR, pWI,
                                pVL, &ldvl, pVR, &ldvr, pWORK, &lwork, &info
                                FCONE FCONE);
                int opt = (info == 0) ? (int)pWORK[0] : info; /* negative info -> pass it back */
    UNPROTECT(6);
    return Rf_ScalarInteger(opt);
  }
  

  
  /* DGEEV: eigenvalues and (optionally) eigenvectors of a real matrix */
  SEXP dgeev_wrapper (SEXP JOBVL, SEXP JOBVR, SEXP N, SEXP A, SEXP LDA, SEXP WR,
                      SEXP WI, SEXP VL, SEXP LDVL, SEXP VR, SEXP LDVR, SEXP WORK,
                      SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP WR_isBM, SEXP WI_isBM,
                      SEXP VL_isBM, SEXP VR_isBM, SEXP WORK_isBM)
  {
    const char *jvl = Rf_translateChar(Rf_asChar(JOBVL));
    const char *jvr = Rf_translateChar(Rf_asChar(JOBVR));
    
    double *pA    = make_double_ptr (A,    A_isBM);
    double *pWR   = make_double_ptr (WR,   WR_isBM);
    double *pWI   = make_double_ptr (WI,   WI_isBM);
    double *pVL   = make_double_ptr (VL,   VL_isBM);
    double *pVR   = make_double_ptr (VR,   VR_isBM);
    double *pWORK = make_double_ptr (WORK, WORK_isBM);
    if (!pA)    Rf_error("NULL pointer for A");
    if (!pWR)   Rf_error("NULL pointer for WR");
    if (!pWI)   Rf_error("NULL pointer for WI");
    if (!pVL)   Rf_error("NULL pointer for VL");
    if (!pVR)   Rf_error("NULL pointer for VR");
    if (!pWORK) Rf_error("NULL pointer for WORK");
    
    INT NN     = (INT) Rf_asInteger(N);
    INT LDAA   = (INT) Rf_asInteger(LDA);
    INT LDVLL  = (INT) Rf_asInteger(LDVL);
    INT LDVRR  = (INT) Rf_asInteger(LDVR);
    INT LWORKK = (INT) Rf_asInteger(LWORK);
    INT INFOO  = (INT) Rf_asInteger(INFO);
    
#ifdef REFBLAS
    dgeev (jvl, jvr, NN, pA, LDAA, pWR, pWI, pVL, LDVLL, pVR, LDVRR, &INFOO);
#else
    F77_CALL(dgeev) ((char*)jvl, (char*)jvr, &NN, pA, &LDAA, pWR, pWI, pVL,
             &LDVLL, pVR, &LDVRR, pWORK, &LWORKK, &INFOO FCLEN2);
    
#endif
    return R_NilValue;
  }
  
  /* DGESDD: divide‑and‑conquer SVD */
  SEXP dgesdd_wrapper (SEXP JOBZ, SEXP M, SEXP N, SEXP A, SEXP LDA,
                       SEXP S, SEXP U, SEXP LDU, SEXP VT, SEXP LDVT, SEXP WORK,
                       SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP S_isBM,
                       SEXP U_isBM, SEXP VT_isBM, SEXP WORK_isBM)
  {
    const char *jz = Rf_translateChar(Rf_asChar(JOBZ));
    
    double *pA    = make_double_ptr (A,  A_isBM);
    double *pS    = make_double_ptr (S,  S_isBM);
    double *pU    = make_double_ptr (U,  U_isBM);
    double *pVT   = make_double_ptr (VT, VT_isBM);
    double *pWORK = make_double_ptr (WORK, WORK_isBM);
    if (!pA)    Rf_error("NULL pointer for A");
    if (!pS)    Rf_error("NULL pointer for S");
    if (!pU)    Rf_error("NULL pointer for U");
    if (!pVT)   Rf_error("NULL pointer for VT");
    if (!pWORK) Rf_error("NULL pointer for WORK");
    
    INT MM     = (INT) Rf_asInteger(M);
    INT NN     = (INT) Rf_asInteger(N);
    INT LDAA   = (INT) Rf_asInteger(LDA);
    INT LDUU   = (INT) Rf_asInteger(LDU);
    INT LDVTT  = (INT) Rf_asInteger(LDVT);
    INT LWORKK = (INT) Rf_asInteger(LWORK);
    INT INFOO  = (INT) Rf_asInteger(INFO);
    
#ifdef REFBLAS
    dgesdd (jz, MM, NN, pA, LDAA, pS, pU, LDUU, pVT, LDVTT, &INFOO);
#else
    INT piworkdim = 8 * MM; if (NN > MM) piworkdim = 8 * NN;
    INT *pIWORK = (INT *) R_alloc((size_t)piworkdim, sizeof(INT));
    F77_CALL(dgesdd) ((char*)jz, &MM, &NN, pA, &LDAA, pS, pU, &LDUU, pVT,
             &LDVTT, pWORK, &LWORKK, pIWORK, &INFOO FCLEN1);
    
#endif
    return R_NilValue;
  }
  
} /* extern "C" */

