
is_transposed = function( tcode )
{
  if ( sum(tcode == c('n', 'N')) > 0 )
    return(FALSE)
  if ( sum(tcode == c('T', 't', 'C', 'c')) > 0 )
    return(TRUE)
  stop("Invalid transpose code given") 
}

check_matrix = function(A, classes=c('big.matrix', 'matrix'), 
  types='double')
{
  if (!any( class(A) %in% classes))
  {
    stop("A is not the correct class type")
  }
  if (!any(typeof(A) == types))
  {
    stop("The matrix type is not correct")
  }
  return( ifelse( inherits(A, 'big.matrix'), TRUE, FALSE ) )
}

# Do I need a function to add a scalar to each element of a matrix?

#' @title Copy a vector.
#'
#' @description Copy double precision DX to double precision DY.
#' For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
#' where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
#' defined in a similar way using INCY.
#' @param N number of elements in input vector(s)
#' @param X double precision vector with N elements
#' @param INCX storage spacing between elements of DX
#' @param Y double precision vector with N elements
#' @param INCY storage spacing between elements of DY
#'
#' @return  DY copy of vector DX (unchanged if N .LE. 0)
#' @references C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. Krogh, Basic linear algebra subprograms for Fortran usage, Algorithm No. 539, Transactions on Mathematical Software 5, 3 (September 1979), pp. 308-323.
#' @export
#'
#' @examples
#' set.seed(4669)
#' A = big.matrix(3, 2, type="double", init=1, dimnames=list(NULL, 
#' c("alpha", "beta")), shared=FALSE)
#' B = big.matrix(3, 2, type="double", init=0, dimnames=list(NULL, 
#' c("alpha", "beta")), shared=FALSE)
#' 
#' dcopy(X=A,Y=B)
#' A[,]-B[,]
#' 
#' # The big.matrix file backings will be deleted when garbage collected.
#' rm(A,B)
#' gc()
#'
# Copy a matrix
# Y := X
dcopy = function(N=NULL, X, INCX=1, Y, INCY=1)
{
  X.is.bm = check_matrix(X)
  Y.is.bm = check_matrix(Y)
  if (is.null(N))
  {
    N = as.double(nrow(X))*as.double(ncol(X))
  }
  .Call(`_dcopy_wrapper`, as.integer(N), X, as.integer(INCX),
        Y, as.integer(INCY), X.is.bm, Y.is.bm)
  return(0)
}

#' @title Swap two double-precision vectors.
#'
#' @description Exchange the elements of two double precision vectors in place.
#' For I = 0 to N-1, swap DX(LX + I * INCX) with DY(LY + I * INCY) where LX and
#' LY depend on the increment signs. When an optimized BLAS is available the
#' implementation dispatches to \code{DSWAP}; otherwise a portable C loop
#' performs the exchange.
#' 
#' @details
#' When an optimized BLAS is available the implementation delegates to the
#' Fortran \code{DSWAP} routine. Otherwise a portable C fallback performs the
#' exchange directly while respecting the supplied vector increments.
#'
#' @param N Number of elements in the input vectors. Defaults to the full
#'   length of \code{X} if \code{NULL}.
#' @param X Double precision vector or matrix providing the first data block.
#' @param INCX Storage spacing between elements of \code{X}.
#' @param Y Double precision vector or matrix providing the second data block.
#' @param INCY Storage spacing between elements of \code{Y}.
#'
#' @return Invisibly returns \code{NULL}; both \code{X} and \code{Y} are
#'   modified in place.
#'   
#' @seealso [dcopy()], [daxpy()] and [dscal()].
#' 
#' @export
#'
#' @examples
#' set.seed(4670)
#' X <- matrix(runif(6), 3, 2)
#' Y <- matrix(runif(6), 3, 2)
#' X_original <- X
#' Y_original <- Y
#' dswap(X = X, Y = Y)
#' all.equal(X, Y_original)
#' all.equal(Y, X_original)
dswap <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L) {
  X.is.bm <- check_matrix(X)
  Y.is.bm <- check_matrix(Y)
  if (is.null(N)) {
    N <- as.integer(nrow(X) * ncol(X))
  }
  N <- as.integer(N)
  INCX <- as.integer(INCX)
  INCY <- as.integer(INCY)

  .Call(`_dswap_wrapper`, N, X, INCX, Y, INCY, X.is.bm, Y.is.bm)
  invisible(NULL)
}

#' @title Add two double-precision vectors.
#'
#' @description Compute double precision DY := DX + DY for N elements, applying the
#' specified storage increments. This routine mirrors the BLAS DAXPY operation with
#' a unit scaling factor.
#' 
#' @details
#' The implementation delegates to the BLAS \code{DAXPY} routine with a unit scaling
#' factor, making it equivalent to \code{daxpy(1.0, X, Y)} while exposing an interface
#' consistent with other low-level wrappers such as \code{dcopy} and \code{dscal}.
#' 
#' @param N number of elements in the input vectors. Defaults to the length of \code{X} if \code{NULL}.
#' @param X double precision vector or matrix providing the addend.
#' @param INCX storage spacing between elements of \code{X}.
#' @param Y double precision vector or matrix that is updated in place.
#' @param INCY storage spacing between elements of \code{Y}.
#'
#' @return The modified object \code{Y} containing the element-wise sums.
#' @export
#' 
#' \seealso{[daxpy()], [dcopy()] and [dscal()].}
#' 
#' @examples
#' set.seed(4669)
#' X <- matrix(runif(6), 3, 2)
#' Y <- matrix(runif(6), 3, 2)
#' dadd(X = X, Y = Y)
#' all.equal(Y, X + Y)
#'
dadd <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L) {
  X.is.bm <- check_matrix(X)
  Y.is.bm <- check_matrix(Y)
  if (is.null(N)) {
    N <- as.integer(nrow(X) * ncol(X))
  }
  N <- as.integer(N)
  INCX <- as.integer(INCX)
  INCY <- as.integer(INCY)

  .Call(`_dadd_wrapper`, N, X, INCX, Y, INCY, X.is.bm, Y.is.bm)
}

#' @title Scales a vector by a constant.
#'
#' @param N an integer. Number of elements in input vector(s)
#' @param ALPHA a real number. The scalar alpha
#' @param Y a big matrix to scale by ALPHA
#' @param INCY an integer. Storage spacing between elements of Y.
#'
#' @return Update Y.
#' @export
#'
#' @examples
#' set.seed(4669)
#' A = big.matrix(3, 2, type="double", init=1, dimnames=list(NULL, 
#' c("alpha", "beta")), shared=FALSE)
#' dscal(ALPHA=2,Y=A)
#' A[,]
#' 
#' # The big.matrix file backings will be deleted when garbage collected.
#' rm(A)
#' gc()
#' 
dscal <- function(N = NULL, ALPHA, Y, INCY = 1L) {
  Y.is.bm <- check_matrix(Y)
  #if (!is.double(Y)) storage.mode(Y) <- "double"
  
  N    <- if (is.null(N)) as.integer(nrow(Y) * ncol(Y)) else as.integer(N)
  INCY <- as.integer(INCY)
  ALPHA <- as.numeric(ALPHA)
  return(.Call(`_dscal_wrapper`, as.integer(N), ALPHA, Y, 
               as.integer(INCY), Y.is.bm))
}


#' @title QR factorization
#'
#' @description DGEQRF computes a QR factorization of a real M-by-N matrix A: A = Q * R.
#' @param M an integer. The number of rows of the matrix A.  M >= 0.
#' @param N an integer. The number of columns of the matrix A.  N >= 0.
#' @param A the M-by-N big matrix A.
#' @param LDA an integer. The leading dimension of the array A.  LDA >= max(1,M).
#' @param TAU a min(M,N) matrix. The scalar factors of the elementary reflectors.
#' @param WORK a (MAX(1,LWORK)) matrix. On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
#' @param LWORK an integer. The dimension of th array WORK.
#'
#' @return M-by-N big matrix A. The elements on and above the diagonal of the 
#' array contain the min(M,N)-by-N upper trapezoidal matrix R (R is upper 
#' triangular if m >= n); the elements below the diagonal, with the array TAU, 
#' represent the orthogonal matrix Q as a product of min(m,n) elementary 
#' reflectors.
#' @export 
#'
#' @examples
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#' h9 <- hilbert(9); h9
#' qr(h9)$rank           #--> only 7
#' qrh9 <- qr(h9, tol = 1e-10)
#' qrh9$rank 
#' C <- as.big.matrix(h9)
#' dgeqrf(A=C)
#' 
#' # The big.matrix file backings will be deleted when garbage collected.
#' rm(C)
#' gc()
#' 
# QR factorization
# return 0 if successful, -i if ith argument has illegal value
dgeqrf=function(M=NULL, N=NULL, A, LDA=NULL, TAU=NULL, WORK=NULL,
  LWORK=NULL)
{
  A.is.bm = check_matrix(A)
  if (is.null(M))
  {
    M = nrow(A)
  }
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  if (is.null(TAU))
  {
    TAU = as.matrix(rep(0.0, min(M,N)))
  }
  if (is.null(LWORK))
  {
    LWORK = max(1, N)
  }
  if (is.null(WORK))
  { 
    WORK = as.matrix(rep(0.0, max(1, LWORK)))
  }
  TAU.is.bm = check_matrix(TAU)
  WORK.is.bm = check_matrix(WORK)
  INFO = 0
  .Call(`_dgeqrf_wrapper`, as.integer(M), as.integer(N), A, as.integer(LDA), 
    TAU, WORK, as.integer(LWORK), as.integer(INFO), A.is.bm, TAU.is.bm, 
    WORK.is.bm)
  return(INFO) 
}




#' @title Cholesky factorization 
#'
#' @description DPOTRF computes the Cholesky factorization of a real symmetric positive definite matrix A.
#' 
#' The factorization has the form
#' \describe{
#'   \item{A =}{ U**T * U,  if UPLO = 'U', or}
#'   \item{A =}{ L  * L**T,  if UPLO = 'L',}
#' }
#' where U is an upper triangular matrix and L is lower triangular.
#' 
#' This is the block version of the algorithm, calling Level 3 BLAS.
#' @param UPLO a character. 
#' \describe{
#'   \item{'U':}{ Upper triangle of A is stored;}
#'   \item{'L':}{ Lower triangle of A is stored.}
#' }
#' @param N an integer. The order of the matrix A.  N >= 0.
#' @param A a big.matrix, dimension (LDA,N).
#' @param LDA an integer. Dimension of the array A.  LDA >= max(1,N).
#'
#' @return updates the big matrix A with the result, INFO is an integer 
#' \describe{
#'   \item{= 0:}{ successful exit}
#'   \item{< 0:}{ if INFO = -i, the i-th argument had an illegal value}
#'   \item{> 0:}{ if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed.}
#' }
#' Terms laying out of the computed triangle should be discarded.
#' @export
#'
#' @examples
#' set.seed(4669)
#' A = matrix(rnorm(16),4)
#' B = as.big.matrix(A %*% t(A))
#' C = A %*% t(A)
#' chol(C)
#' dpotrf(UPLO='U', N=4, A=B, LDA=4)
#' D <- B[,]
#' D[lower.tri(D)]<-0
#' D
#' D-chol(C)
#' t(D)%*%D-C
#' 
#' #' # The big.matrix file backings will be deleted when garbage collected.
#' rm(A,B,C,D)
#' gc()
# return 0 if successful, <0 if -i-th argument is invalid, > 0 if leading minor
# is not positive definite
dpotrf=function(UPLO='U', N=NULL, A, LDA=NULL)
{
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  A.is.bm = check_matrix(A)
  INFO = 0
  .Call(`_dpotrf_wrapper`, as.character(UPLO), as.integer(N), A, as.integer(LDA),
        as.integer(INFO), A.is.bm)
  return(INFO)
}

#' @title DGEEV computes eigenvalues and eigenvectors.
#'
#' @description DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices.
#' 
#' DGEEV computes for an N-by-N real nonsymmetric matrix A, the eigenvalues and, optionally, the left and/or right eigenvectors.
#' The right eigenvector v(j) of A satisfies A * v(j) = lambda(j) * v(j) where lambda(j) is its eigenvalue.
#' The left eigenvector u(j) of A satisfies u(j)**H * A = lambda(j) * u(j)**H where u(j)**H denotes the conjugate-transpose of u(j).
#' 
#' The computed eigenvectors are normalized to have Euclidean norm equal to 1 and largest component real.
#' 
#' @param JOBVL a character.
#' \describe{
#'   \item{= 'N':}{left eigenvectors of A are not computed;}
#'   \item{= 'V':}{left eigenvectors of A are computed.}
#' }
#' @param JOBVR a character.
#' \describe{
#'   \item{= 'N':}{right eigenvectors of A are not computed;}
#'   \item{= 'V':}{right eigenvectors of A are computed.}
#' }
#' @param N an integer. The order of the matrix A. N >= 0.
#' @param A a matrix of dimension (LDA,N), the N-by-N matrix A.
#' @param LDA an integer. The leading dimension of the matrix A. LDA >= max(1,N).
#' @param WR a vector of dimension (N). WR contain the real part of the computed eigenvalues. Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
#' @param WI a vector of dimension (N). WI contain the imaginary part of the computed eigenvalues. Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
#' @param VL a matrx of dimension (LDVL,N)
#' \describe{
#'   \item{If}{ JOBVL = 'V', the left eigenvectors u(j) are stored one
#'   after another in the columns of VL, in the same order
#'   as their eigenvalues.}
#'   \item{If}{ JOBVL = 'N', VL is not referenced.}
#'   \item{If}{ the j-th eigenvalue is real, then u(j) = VL(:,j),
#'   the j-th column of VL.}
#'   \item{If}{ the j-th and (j+1)-st eigenvalues form a complex
#'   conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
#'   u(j+1) = VL(:,j) - i*VL(:,j+1).}
#' }
#' @param LDVL an integer. The leading dimension of the array VL. LDVL >= 1; if JOBVL = 'V', LDVL >= N.
#' @param VR  a matrix of dimension (LDVR,N). 
#' \describe{
#'   \item{If}{ JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.}
#'   \item{If}{ JOBVR = 'N', VR is not referenced.}
#'   \item{If}{ the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.}
#'   \item{If}{ the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).}
#'  }
#' @param LDVR an integer. The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
#' @param WORK a matrix of dimension (MAX(1,LWORK))
#' @param LWORK an integer. The dimension of the array WORK.LWORK >= max(1,3*N), and if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N. For good performance, LWORK must generally be larger. If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
#'
#' @return WR, WI, VR, VL and Work. On exit, A has been overwritten. 
#' @export
#'
#' @examples
#' set.seed(4669)
#' A = matrix(rnorm(16),4)
#' WR= matrix(0,nrow=4,ncol=1)
#' WI= matrix(0,nrow=4,ncol=1)
#' VL = matrix(0,ncol=4,nrow=4)
#' eigen(A)
#' dgeev(A=A,WR=WR,WI=WI,VL=VL)
#' VL
#' WR
#' WI
#' 
#' rm(A,WR,WI,VL)
#' 
#' A = as.big.matrix(matrix(rnorm(16),4))
#' WR= matrix(0,nrow=4,ncol=1)
#' WI= matrix(0,nrow=4,ncol=1)
#' VL = as.big.matrix(matrix(0,ncol=4,nrow=4))
#' eigen(A[,])
#' dgeev(A=A,WR=WR,WI=WI,VL=VL)
#' VL[,]
#' WR[,]
#' WI[,]
#' 
# The big.matrix file backings will be deleted when garbage collected.
#' rm(A,WR,WI,VL)
#' gc()
#'
# General eigenvalue
# return 0 if successful, <0 i-th argument has illegal value, >0 QR 
# algorithm failed.
# for now, VL and VR have to be matrices but they could be NULL

dgeev=function(JOBVL=NULL, JOBVR=NULL, N=NULL, A, LDA=NULL, WR, WI, VL=NULL, 
               LDVL=NULL,  VR=NULL, LDVR=NULL, WORK=NULL, LWORK=NULL)
{
  if (is.null(N))
  {
    N = ncol(A)
  }
  if (is.null(LDA))
  {
    LDA = nrow(A)
  }
  
  # Derive JOB flags from object presence if not explicitly set
  if (is.null(JOBVL)) JOBVL <- if (is.null(VL)) 'N' else 'V'
  if (is.null(JOBVR)) JOBVR <- if (is.null(VR)) 'N' else 'V'
  # Leading dimensions: 1 when that side is not requested
  if (is.null(LDVL)) LDVL <- if (JOBVL == 'V') nrow(VL) else 1L
  if (is.null(LDVR)) LDVR <- if (JOBVR == 'V') nrow(VR) else 1L
  
  # If caller didn't supply WORK/LWORK, query LAPACK for the optimal size.
  if (is.null(LWORK) || is.null(WORK)) {
    opt <- .Call(`_dgeev_lwork_query_wrapper`,
                 as.character(JOBVL), as.character(JOBVR), as.integer(N))
    if (opt < 0L) stop(sprintf("dgeev workspace query failed: INFO=%d", opt))
    LWORK <- max(1L, as.integer(opt))
    WORK  <- as.matrix(numeric(LWORK))
  }
  
  # Take car of the case where someone doesn't want to get the 
  # eigen vectors and passed NULL.
  if (is.null(VL))
  {
    VL = matrix(0.0, nrow=1, ncol=1)
  }
  if (is.null(VR))
  {
    VR = matrix(0.0, nrow=1, ncol=1)
  }
  INFO = 0
  A.is.bm = check_matrix(A)
  WR.is.bm = check_matrix(WR)
  WI.is.bm = check_matrix(WI)
  VL.is.bm = check_matrix(VL)
  VR.is.bm = check_matrix(VR)
  WORK.is.bm = check_matrix(WORK)
  INFO=0
  .Call(`_dgeev_wrapper`, as.character(JOBVL), as.character(JOBVR), as.integer(N), A, as.integer(LDA),
        WR, WI, VL, as.integer(LDVL), VR, as.integer(LDVR),
        WORK, as.integer(LWORK), as.integer(INFO), A.is.bm, WR.is.bm, WI.is.bm,
        VL.is.bm, VR.is.bm, WORK.is.bm)
  return(INFO)
}


#' @title DGESDD computes the singular value decomposition (SVD) of a real matrix.
#'
#' @description DGESDD computes the singular value decomposition (SVD) of a real M-by-N matrix A, optionally computing the left and right singular vectors.  If singular vectors are desired, it uses a divide-and-conquer algorithm.
#' 
#' The SVD is written
#' 
#' A = U * SIGMA * transpose(V)
#' 
#' where SIGMA is an M-by-N matrix which is zero except for its min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA are the singular values of A; they are real and non-negative, and are returned in descending order.  The first min(m,n) columns of U and V are the left and right singular vectors of A.
#' 
#' Note that the routine returns VT = V**T, not V.
#' 
#' @param JOBZ a character. Specifies options for computing all or part of the matrix U:
#' \describe{
#'   \item{= 'A':}{ all M columns of U and all N rows of V**T are returned in the arrays U and VT;}
#'   \item{= 'S':}{ the first min(M,N) columns of U and the first min(M,N) rows of V**T are returned in the arrays U and VT;}
#'   \item{= 'O':}{ If M >= N, the first N columns of U are overwritten on the array A and all rows of V**T are returned in the array VT; otherwise, all columns of U are returned in the array U and the first M rows of V**T are overwritten in the array A;}
#'   \item{= 'N':}{ no columns of U or rows of V**T are computed.}
#' }
#' @param M an integer. The number of rows of the input matrix A. M >= 0.
#' @param N an integer. The number of columns of the input matrix A. N >= 0.
#' @param A the M-by-N matrix A.
#' @param LDA an integer. The leading dimension of the matrix A.  LDA >= max(1,M).
#' @param S  a matrix of dimension (min(M,N)). The singular values of A, sorted so that S(i) >= S(i+1).
#' @param U  U is a matrx of dimension (LDU,UCOL)
#' \describe{
#'   \item{UCOL = M if}{ JOBZ = 'A' or JOBZ = 'O' and M < N; UCOL = min(M,N) if JOBZ = 'S'.}
#'   \item{If}{ JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M orthogonal matrix U;}
#'   \item{if}{ JOBZ = 'S', U contains the first min(M,N) columns of U (the left singular vectors, stored columnwise);}
#'   \item{if}{ JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.}
#' }
#' @param LDU an integer. The leading dimension of the matrix U.  LDU >= 1; if JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
#' @param VT VT is matrix of dimension (LDVT,N)
#' \describe{
#'   \item{If}{ JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the N-by-N orthogonal matrix V**T;}
#'   \item{if}{ JOBZ = 'S', VT contains the first min(M,N) rows of V**T (the right singular vectors, stored rowwise);}
#'   \item{if}{ JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.}
#'  }
#' @param LDVT an integer. The leading dimension of the matrix VT.  LDVT >= 1; if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N; if JOBZ = 'S', LDVT >= min(M,N).
#' @param WORK  a matrix of dimension (MAX(1,LWORK))
#' @param LWORK an integer. The dimension of the array WORK. LWORK >= 1.
#' If LWORK = -1, a workspace query is assumed.  The optimal
#' size for the WORK array is calculated and stored in WORK(1),
#' and no other work except argument checking is performed.
#' 
#' Let mx = max(M,N) and mn = min(M,N).
#' \describe{
#'   \item{If}{ JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ).}
#'   \item{If}{ JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ).}
#'   \item{If}{ JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn.}
#'   \item{If}{ JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx.}
#'  }
#'  These are not tight minimums in all cases; see comments inside code.
#'  For good performance, LWORK should generally be larger;
#'  a query is recommended.
#'
#' @return IWORK an integer matrix dimension of (8*min(M,N))
#' A is updated.
#' \describe{
#'   \item{if}{ JOBZ = 'O',  A is overwritten with the first N columns of U (the left singular vectors, stored columnwise) if M >= N; A is overwritten with the first M rows of V**T (the right singular vectors, stored rowwise) otherwise.}
#'   \item{if}{ JOBZ .ne. 'O', the contents of A are destroyed.}
#' }
#' INFO an integer
#' \describe{
#'   \item{= 0:}{ successful exit.}
#'   \item{< 0:}{ if INFO = -i, the i-th argument had an illegal value.}
#'   \item{> 0:}{ DBDSDC did not converge, updating process failed.}
#' }
#' @export
#'
#' @examples
#' set.seed(4669)
#' A = matrix(rnorm(12),4,3)
#' S = matrix(0,nrow=3,ncol=1)
#' U = matrix(0,nrow=4,ncol=4)
#' VT = matrix(0,ncol=3,nrow=3)
#' dgesdd(A=A,S=S,U=U,VT=VT)
#' S
#' U
#' VT
#' 
#' rm(A,S,U,VT)
#' 
#' A = as.big.matrix(matrix(rnorm(12),4,3))
#' S = as.big.matrix(matrix(0,nrow=3,ncol=1))
#' U = as.big.matrix(matrix(0,nrow=4,ncol=4))
#' VT = as.big.matrix(matrix(0,ncol=3,nrow=3))
#' dgesdd(A=A,S=S,U=U,VT=VT)
#' S[,]
#' U[,]
#' VT[,]
#' 
#' rm(A,S,U,VT)
#' gc()
#' 
# Singular value decomposition (SVD) 
# Returns: = 0 if successful
#          < 0 if INFO = -i had an illegal value
#          > 0 if DBDSDC did not converge
dgesdd = function( JOBZ='A', M=NULL, N=NULL, A, LDA=NULL, S, U, LDU=NULL, 
  VT, LDVT=NULL, WORK=NULL, LWORK=NULL)
{
  A.is.bm = check_matrix(A)
  S.is.bm = check_matrix(S)
  U.is.bm = check_matrix(U)
  VT.is.bm = check_matrix(VT)
  if (is.null(M))
  {
    M=nrow(A)
  }
  if (is.null(N))
  {
    N=ncol(A)
  }
  if (is.null(LDA))
  {
    LDA=nrow(A)
  }
  if (is.null(LDU))
  {
    LDU=nrow(U)
  }
  if (is.null(LDVT))
  {
    LDVT=nrow(VT)
  }
  if (is.null(LWORK) && is.null(WORK))
  {
    if (JOBZ=='N')
    {
      LWORK = 3 * min(M, N) + max( max(M, N), 7*min(M, N) )
    }
    else if (JOBZ=='O')
    {
      LWORK = 3 * (min(M, N))^2 + max( max(M, N), 5 * (min(M, N))^2
        + 4 * min(M, N) )
    }
    else if (JOBZ == 'S' || JOBZ == 'A')
    {
      LWORK = 3 * (min(M, N))^2 + max( max(M, N), 4 * (min(M, N))^2
        + 4 * min(M, N) )
    }
    else
    {
      stop("Invalid JOBZ argument specified")
    }
    WORK = as.matrix(rep(0.0, max(1, LWORK) ))
  }
  if (is.null(LWORK))
  {
    LWORK = length(WORK)
  }
  if (is.null(WORK))
  {
    WORK = as.matrix(rep(0.0, max(1, LWORK)))
  }
  WORK.is.bm = check_matrix(WORK)
  INFO = 0
  .Call(`_dgesdd_wrapper`, as.character(JOBZ), as.integer(M), as.integer(N), A, 
        as.integer(LDA), S, U, as.integer(LDU), VT, as.integer(LDVT), WORK, 
        as.integer(LWORK), as.integer(INFO), A.is.bm, S.is.bm, U.is.bm, VT.is.bm, 
    WORK.is.bm)
  return(INFO)
}
