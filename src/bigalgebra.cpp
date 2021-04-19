#include <string>
#include "bigmemory/BigMatrix.h"



#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


#ifdef REFBLAS
#include "refblas64longlong.h"
#define INT long long
#else
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#define INT int
#endif

#include "bigalgebra.h"



/* Pointer utility, returns a double pointer for either a BigMatrix or a
 * standard R matrix.
 */
double *
make_double_ptr (SEXP matrix, SEXP isBigMatrix)
{
  double *matrix_ptr;

  if (LOGICAL_VALUE (isBigMatrix) == (Rboolean) TRUE)   // Big Matrix
    {
      SEXP address = GET_SLOT (matrix, install ("address"));
      BigMatrix *pbm =
        reinterpret_cast < BigMatrix * >(R_ExternalPtrAddr (address));
      if (!pbm)
        return (NULL);

      // Check that have acceptable big.matrix
      if (pbm->row_offset () > 0 && pbm->ncol () > 1)
        {
          std::string errMsg =
            string ("sub.big.matrix objects cannoth have row ") +
            string
            ("offset greater than zero and number of columns greater than 1");
          Rf_error (errMsg.c_str ());
          return (NULL);
        }

      index_type offset = pbm->nrow () * pbm->col_offset ();
      matrix_ptr = reinterpret_cast < double *>(pbm->matrix ()) + offset;
    }
  else                          // Regular R Matrix
    {
      matrix_ptr = NUMERIC_DATA (matrix);
    }

  return (matrix_ptr);
};


/* Wrappers for miscellaneous BLAS and LAPACK routines. */
SEXP
dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
               SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
               SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM,
               SEXP C_offset)
{
  long j = *(DOUBLE_DATA (C_offset));
  double *pA = make_double_ptr (A, A_isBM);
  double *pB = make_double_ptr (B, B_isBM);
  double *pC;
  SEXP ans;
  INT MM = (INT) * (DOUBLE_DATA (M));
  INT NN = (INT) * (DOUBLE_DATA (N));
  INT KK = (INT) * (DOUBLE_DATA (K));
  INT LDAA = (INT) * (DOUBLE_DATA (LDA));
  INT LDBB = (INT) * (DOUBLE_DATA (LDB));
  INT LDCC = (INT) * (DOUBLE_DATA (LDC));
  if(LOGICAL_VALUE(C_isBM) == (Rboolean) TRUE)
  {
/* Return results in a big matrix */
    pC = make_double_ptr (C, C_isBM) + j;
    PROTECT(ans = C);
  } else {
/* Allocate an output R matrix and return results there
   XXX Add check for size of MM and NN XXX 
 */
    PROTECT(ans = allocMatrix(REALSXP, (int)MM, (int)NN));
    pC = NUMERIC_DATA(ans);
  }
#if REFBLAS
/* Standard Fortran interface without underscoring */
  int8_dgemm ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#else
/* Adaptative Fortran interface from F77_CALL */
  F77_CALL(dgemm) ((char *) CHARACTER_VALUE (TRANSA),
         (char *) CHARACTER_VALUE (TRANSB),
         &MM, &NN, &KK, NUMERIC_DATA (ALPHA), pA, &LDAA, pB,
         &LDBB, NUMERIC_DATA (BETA), pC, &LDCC);
#endif
  unprotect(1);
  return ans;
}



/* Compute A*X + Y for scalar a, vectors X and Y of length N.
 * Y must be a big.matrix, X can be an R vector or big.matrix.
 * The contents of Y are *replaced* by this routine and a reference
 * to Y is returned.
 */
SEXP
daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM)
{
  SEXP ans, Tr;
  double *pY;
  double *pA = DOUBLE_DATA(A);
  double *pX = make_double_ptr (X, X_isBM);
  INT incx = 1;
  INT incy = 1;
  INT NN = (INT) * (DOUBLE_DATA (N));
  PROTECT(ans = Y);
  PROTECT(Tr = allocVector(LGLSXP, 1));
  LOGICAL(Tr)[0] = 1;
  pY = make_double_ptr (Y, Tr);
#if REFBLAS
/* Standard Fortran interface without underscoring */
  int8_daxpy (&NN, pA, pX, &incx, pY, &incy);
#else
/* Adaptative Fortran interface from F77_CALL */
  F77_CALL(daxpy) (&NN, pA, pX, &incx, pY, &incy);
#endif
  unprotect(2);
  return ans;
}


SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM)
{
  SEXP ans;
  const char *_UPLO = CHAR(Rf_asChar(UPLO));
  INT _N = (INT)* (DOUBLE_DATA(N));
  double *_A = make_double_ptr(A, A_isBM);
  INT _LDA = (INT) *(DOUBLE_DATA(LDA));
  INT _INFO = (INT) *(DOUBLE_DATA(INFO));
#if REFBLAS
  /* Standard Fortran interface without underscoring */
  int8_dpotrf (_UPLO, &_N, _A, &_LDA, &_INFO);
#else
  /* Adaptative Fortran interface from F77_CALL */
  F77_CALL(dpotrf) (_UPLO, &_N, _A, &_LDA, &_INFO);
#endif
  PROTECT(ans = A);
  Rf_unprotect(1);
  return ans;
}





void dcopy_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY, SEXP X_isBM,
                    SEXP Y_isBM)
{
  double *pX = make_double_ptr (X, X_isBM);
  double *pY = make_double_ptr (Y, Y_isBM);
  INT NN = (INT)*(DOUBLE_DATA (N));
  INT INCXX = (INT)*(DOUBLE_DATA (INCX));
  INT INCYY = (INT)*(DOUBLE_DATA (INCY));
#ifdef ACMLBLAS
  dcopy (NN, pX, INCXX, pY, INCYY);
#else
  F77_CALL(dcopy) (&NN, pX, &INCXX, pY, &INCYY);
#endif
}

void dscal_wrapper (SEXP N, SEXP ALPHA, SEXP Y, SEXP INCY, SEXP Y_isBM)
{
  double *pY = make_double_ptr (Y, Y_isBM);
  INT NN = (INT)*(DOUBLE_DATA (N));
  INT INCYY = (INT)*(DOUBLE_DATA (INCY));
#ifdef ACMLBLAS
  dscal (NN, *(NUMERIC_DATA (ALPHA)), pY, INCYY);
#else
  F77_CALL(dscal) (&NN, NUMERIC_DATA (ALPHA), pY, &INCYY);
#endif
}

void dgeqrf_wrapper (SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP TAU, SEXP WORK,
                     SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM,
                     SEXP WORK_isBM)
{
  double *pA = make_double_ptr (A, A_isBM);
  double *pTAU = make_double_ptr (TAU, TAU_isBM);
  double *pWORK = make_double_ptr (WORK, WORK_isBM);
  INT MM = (INT)*(DOUBLE_DATA (M));
  INT NN = (INT)*(DOUBLE_DATA (N));
  INT LDAA = (INT)*(DOUBLE_DATA (LDA));
  INT LWORKK = (INT)*(DOUBLE_DATA (LWORK));
  INT INFOO = (INT)*(DOUBLE_DATA (INFO));
#ifdef ACMLBLAS
  dgeqrf (MM, NN, pA, LDAA, pTAU, &INFOO);
#else
  F77_CALL(dgeqrf) (&MM, &NN, pA, &LDAA, pTAU, pWORK, &LWORKK, &INFOO);
#endif
}

void dgeev_wrapper (SEXP JOBVL, SEXP JOBVR, SEXP N, SEXP A, SEXP LDA, SEXP WR,
                    SEXP WI, SEXP VL, SEXP LDVL, SEXP VR, SEXP LDVR, SEXP WORK,
                    SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP WR_isBM, SEXP WI_isBM,
                    SEXP VL_isBM, SEXP VR_isBM, SEXP WORK_isBM)
{
  double *pA = make_double_ptr (A, A_isBM);
  double *pWR = make_double_ptr (WR, WR_isBM);
  double *pWI = make_double_ptr (WI, WI_isBM);
  double *pVL = make_double_ptr (VL, VL_isBM);
  double *pVR = make_double_ptr (VR, VR_isBM);
  double *pWORK = make_double_ptr (WORK, WORK_isBM);
  INT NN = (INT)*(DOUBLE_DATA (N));
  INT LDAA = (INT)*(DOUBLE_DATA (LDA));
  INT LDVLL = (INT)*(DOUBLE_DATA (LDVL));
  INT LDVRR = (INT)*(DOUBLE_DATA (LDVR));
  INT LWORKK = (INT)*(DOUBLE_DATA (LWORK));
  INT INFOO = (INT)*(DOUBLE_DATA (INFO));
#ifdef ACMLBLAS
  dgeev (*((char *)CHARACTER_VALUE (JOBVL)),
         *((char *)CHARACTER_VALUE (JOBVR)),
         NN, pA, LDAA, pWR, pWI, pVL, LDVLL, pVR, LDVRR, &INFOO);
#else
  F77_CALL(dgeev) ((char *)CHARACTER_VALUE (JOBVL), (char *)CHARACTER_VALUE (JOBVR),
         &NN, pA, &LDAA, pWR, pWI, pVL, &LDVLL, pVR, &LDVRR, pWORK,
         &LWORKK, &INFOO);
#endif
}

void dgesdd_wrapper (SEXP JOBZ, SEXP M, SEXP N, SEXP A, SEXP LDA, 
                     SEXP S, SEXP U,
                     SEXP LDU, SEXP VT, SEXP LDVT, SEXP WORK, SEXP LWORK,
                     SEXP INFO, SEXP A_isBM, SEXP S_isBM, SEXP U_isBM,
                     SEXP VT_isBM, SEXP WORK_isBM)
{
  INT *pIWORK;
  INT piworkdim;
  double *pA = make_double_ptr (A, A_isBM);
  double *pS = make_double_ptr (S, S_isBM);
  double *pU = make_double_ptr (U, U_isBM);
  double *pVT = make_double_ptr (VT, VT_isBM);
  double *pWORK = make_double_ptr (WORK, WORK_isBM);
  INT MM = (INT)*(DOUBLE_DATA (N));
  INT NN = (INT)*(DOUBLE_DATA (N));
  INT LDAA = (INT)*(DOUBLE_DATA (LDA));
  INT LDUU = (INT)*(DOUBLE_DATA (LDU));
  INT LDVTT = (INT)*(DOUBLE_DATA (LDVT));
  INT LWORKK = (INT)*(DOUBLE_DATA (LWORK));
  INT INFOO = (INT)*(DOUBLE_DATA (INFO));
  
  piworkdim = 8*MM;
  if(NN>MM) piworkdim = 8*NN;
#ifdef ACMLBLAS
  dgesdd (*((char *)CHARACTER_VALUE (JOBZ)), MM, NN, pA,
          LDAA, pS, pU, LDUU, pVT, LDVTT, &INFOO);
#else
  pIWORK = (INT *)malloc(piworkdim*sizeof(INT));
  F77_CALL(dgesdd) ((char *)CHARACTER_VALUE (JOBZ), &MM, &NN, pA,
          &LDAA, pS, pU, &LDUU, pVT,
          &LDVTT, pWORK, &LWORKK, pIWORK, &INFOO);
          free(pIWORK);
#endif
}

