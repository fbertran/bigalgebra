# DGESDD computes the singular value decomposition (SVD) of a real matrix.

DGESDD computes the singular value decomposition (SVD) of a real M-by-N
matrix A, optionally computing the left and right singular vectors. If
singular vectors are desired, it uses a divide-and-conquer algorithm.

The SVD is written

A = U \* SIGMA \* transpose(V)

where SIGMA is an M-by-N matrix which is zero except for its min(m,n)
diagonal elements, U is an M-by-M orthogonal matrix, and V is an N-by-N
orthogonal matrix. The diagonal elements of SIGMA are the singular
values of A; they are real and non-negative, and are returned in
descending order. The first min(m,n) columns of U and V are the left and
right singular vectors of A.

Note that the routine returns VT = V\*\*T, not V.

## Usage

``` r
dgesdd(
  JOBZ = "A",
  M = NULL,
  N = NULL,
  A,
  LDA = NULL,
  S,
  U,
  LDU = NULL,
  VT,
  LDVT = NULL,
  WORK = NULL,
  LWORK = NULL
)
```

## Arguments

- JOBZ:

  a character. Specifies options for computing all or part of the matrix
  U:

  = 'A':

  :   all M columns of U and all N rows of V\*\*T are returned in the
      arrays U and VT;

  = 'S':

  :   the first min(M,N) columns of U and the first min(M,N) rows of
      V\*\*T are returned in the arrays U and VT;

  = 'O':

  :   If M \>= N, the first N columns of U are overwritten on the array
      A and all rows of V\*\*T are returned in the array VT; otherwise,
      all columns of U are returned in the array U and the first M rows
      of V\*\*T are overwritten in the array A;

  = 'N':

  :   no columns of U or rows of V\*\*T are computed.

- M:

  an integer. The number of rows of the input matrix A. M \>= 0.

- N:

  an integer. The number of columns of the input matrix A. N \>= 0.

- A:

  the M-by-N matrix A.

- LDA:

  an integer. The leading dimension of the matrix A. LDA \>= max(1,M).

- S:

  a matrix of dimension (min(M,N)). The singular values of A, sorted so
  that S(i) \>= S(i+1).

- U:

  U is a matrx of dimension (LDU,UCOL)

  UCOL = M if

  :   JOBZ = 'A' or JOBZ = 'O' and M \< N; UCOL = min(M,N) if JOBZ =
      'S'.

  If

  :   JOBZ = 'A' or JOBZ = 'O' and M \< N, U contains the M-by-M
      orthogonal matrix U;

  if

  :   JOBZ = 'S', U contains the first min(M,N) columns of U (the left
      singular vectors, stored columnwise);

  if

  :   JOBZ = 'O' and M \>= N, or JOBZ = 'N', U is not referenced.

- LDU:

  an integer. The leading dimension of the matrix U. LDU \>= 1; if JOBZ
  = 'S' or 'A' or JOBZ = 'O' and M \< N, LDU \>= M.

- VT:

  VT is matrix of dimension (LDVT,N)

  If

  :   JOBZ = 'A' or JOBZ = 'O' and M \>= N, VT contains the N-by-N
      orthogonal matrix V\*\*T;

  if

  :   JOBZ = 'S', VT contains the first min(M,N) rows of V\*\*T (the
      right singular vectors, stored rowwise);

  if

  :   JOBZ = 'O' and M \< N, or JOBZ = 'N', VT is not referenced.

- LDVT:

  an integer. The leading dimension of the matrix VT. LDVT \>= 1; if
  JOBZ = 'A' or JOBZ = 'O' and M \>= N, LDVT \>= N; if JOBZ = 'S', LDVT
  \>= min(M,N).

- WORK:

  a matrix of dimension (MAX(1,LWORK))

- LWORK:

  an integer. The dimension of the array WORK. LWORK \>= 1. If LWORK =
  -1, a workspace query is assumed. The optimal size for the WORK array
  is calculated and stored in WORK(1), and no other work except argument
  checking is performed.

  Let mx = max(M,N) and mn = min(M,N).

  If

  :   JOBZ = 'N', LWORK \>= 3\*mn + max( mx, 7\*mn ).

  If

  :   JOBZ = 'O', LWORK \>= 3\*mn + max( mx, 5\*mn\*mn + 4\*mn ).

  If

  :   JOBZ = 'S', LWORK \>= 4\*mn\*mn + 7\*mn.

  If

  :   JOBZ = 'A', LWORK \>= 4\*mn\*mn + 6\*mn + mx.

  These are not tight minimums in all cases; see comments inside code.
  For good performance, LWORK should generally be larger; a query is
  recommended.

## Value

IWORK an integer matrix dimension of (8\*min(M,N)) A is updated.

- if:

  JOBZ = 'O', A is overwritten with the first N columns of U (the left
  singular vectors, stored columnwise) if M \>= N; A is overwritten with
  the first M rows of V\*\*T (the right singular vectors, stored
  rowwise) otherwise.

- if:

  JOBZ .ne. 'O', the contents of A are destroyed.

INFO an integer

- = 0::

  successful exit.

- \< 0::

  if INFO = -i, the i-th argument had an illegal value.

- \> 0::

  DBDSDC did not converge, updating process failed.

## Examples

``` r
set.seed(4669)
A = matrix(rnorm(12),4,3)
S = matrix(0,nrow=3,ncol=1)
U = matrix(0,nrow=4,ncol=4)
VT = matrix(0,ncol=3,nrow=3)
dgesdd(A=A,S=S,U=U,VT=VT)
#> [1] 0
S
#>           [,1]
#> [1,] 2.5683690
#> [2,] 1.9865047
#> [3,] 0.1726739
U
#>            [,1]        [,2]       [,3]       [,4]
#> [1,]  0.1711751  0.09221058 0.78855772  0.5834150
#> [2,]  0.8605882 -0.33631329 0.12096045 -0.3628360
#> [3,]  0.2092122  0.91796420 0.08793762 -0.3253290
#> [4,] -0.4316449 -0.18902993 0.59650002 -0.6497215
VT
#>            [,1]        [,2]        [,3]
#> [1,]  0.2496213 -0.04136007  0.96745984
#> [2,]  0.1202120 -0.99003537 -0.07334196
#> [3,] -0.9608529 -0.13460796  0.24216198

rm(A,S,U,VT)

A = as.big.matrix(matrix(rnorm(12),4,3))
S = as.big.matrix(matrix(0,nrow=3,ncol=1))
U = as.big.matrix(matrix(0,nrow=4,ncol=4))
VT = as.big.matrix(matrix(0,ncol=3,nrow=3))
dgesdd(A=A,S=S,U=U,VT=VT)
#> [1] 0
S[,]
#> [1] 2.6819327 1.2859400 0.8128895
U[,]
#>            [,1]       [,2]       [,3]        [,4]
#> [1,] -0.1658242 -0.5205485 -0.3243924 -0.77220537
#> [2,]  0.4356896 -0.6461869 -0.3770659  0.50043820
#> [3,] -0.8648680 -0.1078115 -0.3028043  0.38560286
#> [4,] -0.1862263 -0.5475842  0.8129578  0.06760851
VT[,]
#>           [,1]       [,2]       [,3]
#> [1,] 0.5814564  0.0745909  0.8101510
#> [2,] 0.4372976  0.8110586 -0.3885289
#> [3,] 0.6860607 -0.5801897 -0.4389768

rm(A,S,U,VT)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1150951 61.5    2286010 122.1         NA  2286010 122.1
#> Vcells 2060435 15.8    8388608  64.0      65536  6482741  49.5
```
