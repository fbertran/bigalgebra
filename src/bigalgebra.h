
#ifdef __cplusplus
extern "C"
{
#endif
  
  double *make_double_ptr (SEXP matrix, SEXP isBigMatrix);
  
  SEXP dgemm_wrapper (SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
                       SEXP ALPHA, SEXP A, SEXP LDA, SEXP B, SEXP LDB,
                       SEXP BETA, SEXP C, SEXP LDC, SEXP A_isBM, SEXP B_isBM,
                       SEXP C_isBM, SEXP C_offset);
  SEXP daxpy_wrapper (SEXP N, SEXP A, SEXP X, SEXP Y, SEXP X_isBM);
  
#ifdef __cplusplus
}
#endif

