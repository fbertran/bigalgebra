# XXX TODO
#
# Each routine below may return either a big matrix or an R matrix
# depending on how they are called (or a set of such matrices for
# the decomposition methods). In each case the returned value is
# allocated by the routine.
#
# When big matrices are returned, we need to register a finalizer
# with the big matrix that removes its allocated files when the
# garbage collector is called on the big matrix, to clean up.
# Add this!
#


#' Matrix Multiply
#'
#' @description This is function provides dgemm functionality, which DGEMM 
#' performs one of the matrix-matrix operations.
#' C := ALPHA * op(A) * op(B) + BETA * C.
#' 
#' @param TRANSA a character. TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:
#' \describe{
#'   \item{TRANSA =}{ 'N' or 'n',  op( A ) = A.}
#'   \item{TRANSA =}{ 'T' or 't',  op( A ) = A**T.}
#'   \item{TRANSA =}{ 'C' or 'c',  op( A ) = A**T.}
#' }
#' @param TRANSB a character. TRANSB specifies the form of op( B ) to be used in the matrix multiplication as follows:
#' #' \describe{
#'   \item{TRANSA =}{ 'N' or 'n',  op( B ) = B.}
#'   \item{TRANSA =}{ 'T' or 't',  op( B ) = B**T.}
#'   \item{TRANSA =}{ 'C' or 'c',  op( B ) = B**T.}
#' }
#' @param M an integer. M specifies the number of rows  of the  matrix op( A )  and of the  matrix  C.  M  must  be at least  zero.
#' @param N an integer. N specifies the number of columns  of the  matrix op( B )  and of the  matrix  C.  N  must  be at least  zero.
#' @param K an integer. K specifies the number of columns  of the  matrix op( A )  and the number of rows of the  matrix  op( B ).  K  must  be at least  zero.
#' @param ALPHA a real number. Specifies the scalar alpha.
#' @param A a matrix of dimension (LDA, ka), where ka is k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. Before entry with  TRANSA = 'N' or 'n',  the leading  m by k part of the array  A  must contain the matrix  A,  otherwise the leading  k by m  part of the array  A  must contain  the matrix A.
#' @param LDA an integer. 
#' @param B a matrix of dimension ( LDB, kb ), where kb is n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. Before entry with  TRANSB = 'N' or 'n',  the leading  k by n part of the array  B  must contain the matrix  B,  otherwise the leading  n by k  part of the array  B  must contain  the matrix B.
#' @param LDB an integer. 
#' @param BETA a real number. Specifies the scalar beta
#' @param C a matrix of dimension ( LDC, N ). Before entry, the leading  m by n  part of the array  C must contain the matrix  C,  except when  beta  is zero, in which case C need not be set on entry. On exit, the array  C  is overwritten by the  m by n  matrix ( alpha*op( A )*op( B ) + beta*C ).
#' @param LDC an integer. 
#' @param COFF offset for C.
#'
#' @return Update C with the result.
#' @export
#'
#' @examples
#' require(bigmemory)
#' A = as.big.matrix(matrix(1, nrow=3, ncol=2))
#' B <- big.matrix(2, 3, type="double", init=-1,
#'                 dimnames=list(NULL, c("alpha", "beta")), shared=FALSE)
#' C = big.matrix(3, 3, type="double", init=1,
#'                dimnames=list(NULL, c("alpha", "beta", "gamma")), shared=FALSE)  
#' 2*A[,]%*%B[,]+0.5*C[,]
#' E = dgemm(ALPHA=2.0, A=A, B=B, BETA=0.5, C=C)
#' E[,] # Same result
#' 
#' # The big.matrix file backings will be deleted when garbage collected.
#' rm(A,B,C,E)
#' gc()
# Matrix Multiply
# C := ALPHA * op(A) * op(B) + BETA * C
# This is function provides dgemm functionality.
dgemm = function(TRANSA='N', TRANSB='N', M=NULL, N=NULL, K=NULL,
  ALPHA=1, A, LDA=NULL, B, LDB=NULL, BETA=0, C, LDC=NULL, COFF=0) 
{
  A.is.bm = check_matrix(A)
  B.is.bm = check_matrix(B)
# The matrices look OK.  Now, if they haven't been specified, let's
# specify some reasonable dimension information.
  if ( is.null(M) )
  {
    M = ifelse ( is_transposed(TRANSA), ncol(A), nrow(A) )
  }
  if ( is.null(N) ) 
  {
    N = ifelse ( is_transposed(TRANSB), nrow(B), ncol(B) )
  }
  if ( is.null(K) )
  {
    K = ifelse ( is_transposed(TRANSA), nrow(A), ncol(A) )
  }
  if ( is.null(LDA) ) LDA = ifelse (is_transposed(TRANSA), K, M)
  if ( is.null(LDB) ) LDB = ifelse (is_transposed(TRANSB), N, K) 
  if ( is.null(LDC) ) LDC = M

# Default to big matrix output
  if(missing(C)) C = anon_matrix(M, N)
  C.is.bm = "big.matrix" %in% class(C)

  .Call(`_dgemm_wrapper`, as.character(TRANSA), as.character(TRANSB),
    as.double(M), as.double(N), as.double(K), as.double(ALPHA), A, 
    as.double(LDA), B, as.double(LDB),
    as.double(BETA), C, as.double(LDC), as.logical(A.is.bm), 
    as.logical(B.is.bm), as.logical(C.is.bm), COFF)
}

#' @title BLAS daxpy functionality
#'
#' @description This function implements the function Y := A * X + Y where X and Y may be either native double-precision valued R matrices or numeric vectors, or double-precision valued \code{\link[bigmemory]{big.matrix}} objects, and A is a scalar.
#' @param A Optional numeric scalar value to scale the matrix \code{X} by, with a default value of 1.
#' @param X Requried to be either a native \R \code{\link{matrix}} or numeric vector, or a \code{\link[bigmemory]{big.matrix}} object
#' @param Y Optional native \R \code{\link{matrix}} or numeric vector, or a \code{\link[bigmemory]{big.matrix}} object
#'
#' @details At least one of either \code{X} or \code{Y} must be a \code{big.matrix}. All values must be of type \code{double} (the only type presently supported by the bigalgebra package).
#' 
#' This function is rarely necessary to use directly since the bigalgebra package defines standard arithmetic operations and scalar multiplication. It is more efficient to use \code{daxpy} directly when both scaling and matrix addition are required, in which case both operations are performed in one step.
#' 
#' @return The output value depends on the classes of input values \code{X} and \code{Y} and on the value of the global option \code{bigalgebra.mixed_arithmetic_returns_R_matrix}.
#' 
#' If \code{X} and \code{Y} are both big matrices, or \code{Y} is missing, \code{options("bigalgebra.mixed_arithmetic_returns_R_matrix")} is \code{FALSE}, then a \code{big.matrix} is returned. The returned \code{big.matrix} is backed by a temporary file mapping that will be deleted when the returned result is garbage collected by R (see the examples).
#' 
#' Otherwise, a standard R matrix is returned. The dimensional shape of the output is taken from \code{X}. If input \code{X} is dimensionless (that is, lacks a dimension attribute), then the output is a column vector.
#' 
#' @references \url{https://www.netlib.org/blas/daxpy.f}
#' @author Michael J. Kane
#' @seealso \code{\link[bigmemory]{bigmemory}}
#' @export
#'
#' @examples
#' require(bigmemory)
#'A = matrix(1, nrow=3, ncol=2)
#'B <- big.matrix(3, 2, type="double", init=0,
#'                dimnames=list(NULL, c("alpha", "beta")), shared=FALSE)
#'C = B + B   # C is a new big matrix
#'D = A + B   # D defaults to a regular R matrix, to change this, set the option:
#'# options(bigalgebra.mixed_arithmetic_returns_R_matrix=FALSE)
#'E = daxpy(A=1.0, X=B, Y=B)  # Same kind of result as C
#'print(C[])
#'print(D)
#'print(E[])
#'
#'# The C and E big.matrix file backings will be deleted when garbage collected:
#'# (We enable debugging to see this explicitly)
#'options(bigalgebra.DEBUG=TRUE)
#'rm(C,E)
#'gc()
#' 
# Vector addition and scaling
# Y := A * X  + Y
# A is a scalar double
# X is either a big.matrix or regular R matrix.
# Y is an optional matrix or vector of the same length as X.
# Returns a new matrix or big matrix with the same dimensions as X. If
# X is a dimension-less R vector, returns a column. Returned value type
# depends on the arguments and the value of the option
# options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]].
daxpy = function(A=1, X, Y)
{
  mixed = FALSE
  X.is.bm = check_matrix(X,classes=c('big.matrix','matrix','vector','numeric'))
# default to a column big matrix output
  M = length(X)
  L = M
  N = 1L
  D = dim(X)
  if(!is.null(D) && length(D)==2)
  {
    M = D[1]
    N = D[2]
  }
  Z = anon_matrix(M,N,val=0.0)
  if(!missing(Y))
  {
# Check conformity of Y and duplicate
    if(length(Y)!=length(X)) stop("Lengths of X and Y must match")
    mixed = (X.is.bm != check_matrix(Y,classes=c('big.matrix','matrix','vector','numeric')))
    Z[] = Y[]
  }
  ans = .Call(`_daxpy_wrapper`, as.double(L), as.double(A), X, Z, X.is.bm)
  if(mixed) return(ans[])
  ans
}
