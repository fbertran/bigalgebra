# Matrix Multiply

This is function provides dgemm functionality, which DGEMM performs one
of the matrix-matrix operations. C := ALPHA \* op(A) \* op(B) + BETA \*
C.

## Usage

``` r
dgemm(
  TRANSA = "N",
  TRANSB = "N",
  M = NULL,
  N = NULL,
  K = NULL,
  ALPHA = 1,
  A,
  LDA = NULL,
  B,
  LDB = NULL,
  BETA = 0,
  C,
  LDC = NULL,
  COFF = 0
)
```

## Arguments

- TRANSA:

  a character. TRANSA specifies the form of op( A ) to be used in the
  matrix multiplication as follows:

  TRANSA =

  :   'N' or 'n', op( A ) = A.

  TRANSA =

  :   'T' or 't', op( A ) = A\*\*T.

  TRANSA =

  :   'C' or 'c', op( A ) = A\*\*T.

- TRANSB:

  a character. TRANSB specifies the form of op( B ) to be used in the
  matrix multiplication as follows: \#'

  TRANSA =

  :   'N' or 'n', op( B ) = B.

  TRANSA =

  :   'T' or 't', op( B ) = B\*\*T.

  TRANSA =

  :   'C' or 'c', op( B ) = B\*\*T.

- M:

  an integer. M specifies the number of rows of the matrix op( A ) and
  of the matrix C. M must be at least zero.

- N:

  an integer. N specifies the number of columns of the matrix op( B )
  and of the matrix C. N must be at least zero.

- K:

  an integer. K specifies the number of columns of the matrix op( A )
  and the number of rows of the matrix op( B ). K must be at least zero.

- ALPHA:

  a real number. Specifies the scalar alpha.

- A:

  a matrix of dimension (LDA, ka), where ka is k when TRANSA = 'N' or
  'n', and is m otherwise. Before entry with TRANSA = 'N' or 'n', the
  leading m by k part of the array A must contain the matrix A,
  otherwise the leading k by m part of the array A must contain the
  matrix A.

- LDA:

  an integer.

- B:

  a matrix of dimension ( LDB, kb ), where kb is n when TRANSB = 'N' or
  'n', and is k otherwise. Before entry with TRANSB = 'N' or 'n', the
  leading k by n part of the array B must contain the matrix B,
  otherwise the leading n by k part of the array B must contain the
  matrix B.

- LDB:

  an integer.

- BETA:

  a real number. Specifies the scalar beta

- C:

  a matrix of dimension ( LDC, N ). Before entry, the leading m by n
  part of the array C must contain the matrix C, except when beta is
  zero, in which case C need not be set on entry. On exit, the array C
  is overwritten by the m by n matrix ( alpha\*op( A )\*op( B ) +
  beta\*C ).

- LDC:

  an integer.

- COFF:

  offset for C.

## Value

Update C with the result.

## Examples

``` r
require(bigmemory)
A = as.big.matrix(matrix(1, nrow=3, ncol=2))
B <- big.matrix(2, 3, type="double", init=-1,
                dimnames=list(NULL, c("alpha", "beta")), shared=FALSE)
C = big.matrix(3, 3, type="double", init=1,
               dimnames=list(NULL, c("alpha", "beta", "gamma")), shared=FALSE)  
2*A[,]%*%B[,]+0.5*C[,]
#>      alpha beta gamma
#> [1,]  -3.5 -3.5  -3.5
#> [2,]  -3.5 -3.5  -3.5
#> [3,]  -3.5 -3.5  -3.5
E = dgemm(ALPHA=2.0, A=A, B=B, BETA=0.5, C=C)
E[,] # Same result
#>      alpha beta gamma
#> [1,]  -3.5 -3.5  -3.5
#> [2,]  -3.5 -3.5  -3.5
#> [3,]  -3.5 -3.5  -3.5

# The big.matrix file backings will be deleted when garbage collected.
rm(A,B,C,E)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1149438 61.4    2286010 122.1         NA  2286010 122.1
#> Vcells 2056964 15.7    8388608  64.0      65536  6482741  49.5
```
