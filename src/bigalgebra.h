
#ifdef __cplusplus
extern "C"
{
#endif
  
  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);
  
  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                       SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB,
                       SEXP BETA, SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM,
                       SEXP C_isBM, SEXP C_offset);
  SEXP dswap_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                      SEXP X_isBM, SEXP Y_isBM);
  SEXP dadd_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM);
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM);
  SEXP dpotrf_wrapper(SEXP UPLO, SEXP N, SEXP A, SEXP LDA, SEXP INFO, SEXP A_isBM);

  SEXP dcopy_wrapper(SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                            SEXP X_isBM, SEXP Y_isBM);
  SEXP dscal_wrapper (SEXP N, SEXP ALPHA, SEXP Y, SEXP INCY, SEXP Y_isBM);
  SEXP dset_wrapper (SEXP N, SEXP ALPHA, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP dvcal_wrapper (SEXP N, SEXP ALPHA, SEXP X, SEXP INCX,
                      SEXP BETA, SEXP Y, SEXP INCY,
                      SEXP X_isBM, SEXP Y_isBM);
  SEXP dsub_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM);
  SEXP dsqrt_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP ddot_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                     SEXP X_isBM, SEXP Y_isBM);
  SEXP dqddot_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY,
                       SEXP X_isBM, SEXP Y_isBM);
  SEXP dhprod_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP Y, SEXP INCY, SEXP Z,
                       SEXP INCZ, SEXP X_isBM, SEXP Y_isBM, SEXP Z_isBM);
  SEXP dxyz_wrapper (SEXP X, SEXP Y, SEXP Z,
                     SEXP X_isBM, SEXP Y_isBM, SEXP Z_isBM);
  SEXP dsum_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP dasum_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP dnrm2_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP dprdct_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP idmin_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP idmax_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP idamin_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP idamax_wrapper (SEXP N, SEXP X, SEXP INCX, SEXP X_isBM);
  SEXP dsymm_wrapper (SEXP SIDE, SEXP UPLO, SEXP M, SEXP N, SEXP ALPHA,
                      SEXP A, SEXP LDA, SEXP B, SEXP LDB, SEXP BETA,
                      SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM, SEXP C_isBM);
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

