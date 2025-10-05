
#ifdef __cplusplus
extern "C"
{
#endif
  
  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);
  
  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                       SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB,
                       SEXP BETA, SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM,
                       SEXP C_isBM, SEXP C_offset);
  SEXP dadd_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM);
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM);
  SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM);
  
  SEXP dcopy_wrapper(SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                            SEXP X_isBM, SEXP Y_isBM);
  SEXP dscal_wrapper (SEXP N, SEXP ALPHA, SEXP Y, SEXP INCY, SEXP Y_isBM);
  SEXP dgeqrf_wrapper (SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP TAU, SEXP WORK,
                       SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP TAU_isBM,
                       SEXP WORK_isBM);
  SEXP dgeev_lwork_query_wrapper(SEXP JOBVL, SEXP JOBVR, SEXP N);
  SEXP dgeev_wrapper (SEXP JOBVL, SEXP JOBVR, SEXP N, SEXP A, SEXP LDA, SEXP WR,
                      SEXP WI, SEXP VL, SEXP LDVL, SEXP VR, SEXP LDVR, SEXP WORK,
                      SEXP LWORK, SEXP INFO, SEXP A_isBM, SEXP WR_isBM, SEXP WI_isBM,
                      SEXP VL_isBM, SEXP VR_isBM, SEXP WORK_isBM);
  SEXP dgesdd_wrapper (SEXP JOBZ, SEXP M, SEXP N, SEXP A, SEXP LDA, 
                       SEXP S, SEXP U,
                       SEXP LDU, SEXP VT, SEXP LDVT, SEXP WORK, SEXP LWORK,
                       SEXP INFO, SEXP A_isBM, SEXP S_isBM, SEXP U_isBM,
                       SEXP VT_isBM, SEXP WORK_isBM);
    
#ifdef __cplusplus
}
#endif

