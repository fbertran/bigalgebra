# DGEEV computes eigenvalues and eigenvectors.

DGEEV computes the eigenvalues and, optionally, the left and/or right
eigenvectors for GE matrices.

DGEEV computes for an N-by-N real nonsymmetric matrix A, the eigenvalues
and, optionally, the left and/or right eigenvectors. The right
eigenvector v(j) of A satisfies A \* v(j) = lambda(j) \* v(j) where
lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
u(j)\*\*H \* A = lambda(j) \* u(j)\*\*H where u(j)\*\*H denotes the
conjugate-transpose of u(j).

The computed eigenvectors are normalized to have Euclidean norm equal to
1 and largest component real.

## Usage

``` r
dgeev(
  JOBVL = NULL,
  JOBVR = NULL,
  N = NULL,
  A,
  LDA = NULL,
  WR,
  WI,
  VL = NULL,
  LDVL = NULL,
  VR = NULL,
  LDVR = NULL,
  WORK = NULL,
  LWORK = NULL
)
```

## Arguments

- JOBVL:

  a character.

  = 'N':

  :   left eigenvectors of A are not computed;

  = 'V':

  :   left eigenvectors of A are computed.

- JOBVR:

  a character.

  = 'N':

  :   right eigenvectors of A are not computed;

  = 'V':

  :   right eigenvectors of A are computed.

- N:

  an integer. The order of the matrix A. N \>= 0.

- A:

  a matrix of dimension (LDA,N), the N-by-N matrix A.

- LDA:

  an integer. The leading dimension of the matrix A. LDA \>= max(1,N).

- WR:

  a vector of dimension (N). WR contain the real part of the computed
  eigenvalues. Complex conjugate pairs of eigenvalues appear
  consecutively with the eigenvalue having the positive imaginary part
  first.

- WI:

  a vector of dimension (N). WI contain the imaginary part of the
  computed eigenvalues. Complex conjugate pairs of eigenvalues appear
  consecutively with the eigenvalue having the positive imaginary part
  first.

- VL:

  a matrx of dimension (LDVL,N)

  If

  :   JOBVL = 'V', the left eigenvectors u(j) are stored one after
      another in the columns of VL, in the same order as their
      eigenvalues.

  If

  :   JOBVL = 'N', VL is not referenced.

  If

  :   the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column
      of VL.

  If

  :   the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
      then u(j) = VL(:,j) + i\*VL(:,j+1) and u(j+1) = VL(:,j) -
      i\*VL(:,j+1).

- LDVL:

  an integer. The leading dimension of the array VL. LDVL \>= 1; if
  JOBVL = 'V', LDVL \>= N.

- VR:

  a matrix of dimension (LDVR,N).

  If

  :   JOBVR = 'V', the right eigenvectors v(j) are stored one after
      another in the columns of VR, in the same order as their
      eigenvalues.

  If

  :   JOBVR = 'N', VR is not referenced.

  If

  :   the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column
      of VR.

  If

  :   the j-th and (j+1)-st eigenvalues form a complex conjugate pair,
      then v(j) = VR(:,j) + i\*VR(:,j+1) and v(j+1) = VR(:,j) -
      i\*VR(:,j+1).

- LDVR:

  an integer. The leading dimension of the array VR. LDVR \>= 1; if
  JOBVR = 'V', LDVR \>= N.

- WORK:

  a matrix of dimension (MAX(1,LWORK))

- LWORK:

  an integer. The dimension of the array WORK.LWORK \>= max(1,3\*N), and
  if JOBVL = 'V' or JOBVR = 'V', LWORK \>= 4\*N. For good performance,
  LWORK must generally be larger. If LWORK = -1, then a workspace query
  is assumed; the routine only calculates the optimal size of the WORK
  array, returns this value as the first entry of the WORK array, and no
  error message related to LWORK is issued by XERBLA.

## Value

WR, WI, VR, VL and Work. On exit, A has been overwritten.

## Examples

``` r
set.seed(4669)
A = matrix(rnorm(16),4)
WR= matrix(0,nrow=4,ncol=1)
WI= matrix(0,nrow=4,ncol=1)
VL = matrix(0,ncol=4,nrow=4)
eigen(A)
#> eigen() decomposition
#> $values
#> [1]  0.85730911+1.690442i  0.85730911-1.690442i -0.86120947+0.000000i
#> [4] -0.04044732+0.000000i
#> 
#> $vectors
#>                         [,1]                    [,2]          [,3]
#> [1,] -0.20131416-0.08051662i -0.20131416+0.08051662i -0.4396218+0i
#> [2,] -0.73251338+0.00000000i -0.73251338+0.00000000i  0.4547709+0i
#> [3,] -0.06644069-0.55880008i -0.06644069+0.55880008i -0.1693250+0i
#> [4,]  0.22714677+0.21942261i  0.22714677-0.21942261i -0.7558076+0i
#>               [,4]
#> [1,] -0.9125418+0i
#> [2,] -0.2656930+0i
#> [3,]  0.2527695+0i
#> [4,]  0.1810590+0i
#> 
dgeev(A=A,WR=WR,WI=WI,VL=VL)
#> [1] 0
VL
#>              [,1]       [,2]       [,3]       [,4]
#> [1,] -0.099694024  0.1149983 -0.2126753 -0.7057922
#> [2,]  0.657326284  0.0000000 -0.1918649  0.3331385
#> [3,]  0.007357843  0.5515939 -0.3234797  0.3139870
#> [4,]  0.451853982 -0.1904644 -0.9018414  0.5406369
WR
#>             [,1]
#> [1,]  0.85730911
#> [2,]  0.85730911
#> [3,] -0.86120947
#> [4,] -0.04044732
WI
#>           [,1]
#> [1,]  1.690442
#> [2,] -1.690442
#> [3,]  0.000000
#> [4,]  0.000000

rm(A,WR,WI,VL)

A = as.big.matrix(matrix(rnorm(16),4))
WR= matrix(0,nrow=4,ncol=1)
WI= matrix(0,nrow=4,ncol=1)
VL = as.big.matrix(matrix(0,ncol=4,nrow=4))
eigen(A[,])
#> eigen() decomposition
#> $values
#> [1] -0.7802559+0.6873447i -0.7802559-0.6873447i  0.6930662+0.4810716i
#> [4]  0.6930662-0.4810716i
#> 
#> $vectors
#>                       [,1]                  [,2]                   [,3]
#> [1,]  0.1471182+0.4334710i  0.1471182-0.4334710i -0.2986975+0.07375248i
#> [2,]  0.2638302+0.1607608i  0.2638302-0.1607608i -0.4446219-0.18640256i
#> [3,] -0.4094915+0.0733955i -0.4094915-0.0733955i  0.7778787+0.00000000i
#> [4,] -0.7224532+0.0000000i -0.7224532+0.0000000i  0.2526100-0.06323534i
#>                        [,4]
#> [1,] -0.2986975-0.07375248i
#> [2,] -0.4446219+0.18640256i
#> [3,]  0.7778787+0.00000000i
#> [4,]  0.2526100+0.06323534i
#> 
dgeev(A=A,WR=WR,WI=WI,VL=VL)
#> [1] 0
VL[,]
#>             [,1]       [,2]        [,3]        [,4]
#> [1,] -0.63149291  0.0000000  0.35751797  0.04531232
#> [2,] -0.22007417  0.1890068 -0.81587447  0.00000000
#> [3,] -0.33976651  0.2889619 -0.32444761 -0.26761276
#> [4,] -0.08779556 -0.5571465 -0.04124336  0.16091204
WR[,]
#> [1] -0.7802559 -0.7802559  0.6930662  0.6930662
WI[,]
#> [1]  0.6873447 -0.6873447  0.4810716 -0.4810716

rm(A,WR,WI,VL)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1014132 54.2    1995290 106.6         NA  1995290 106.6
#> Vcells 1821032 13.9    8388608  64.0      65536  2724101  20.8
```
