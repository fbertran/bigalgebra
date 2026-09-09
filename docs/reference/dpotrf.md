# Cholesky factorization

DPOTRF computes the Cholesky factorization of a real symmetric positive
definite matrix A.

The factorization has the form

- A =:

  U\*\*T \* U, if UPLO = 'U', or

- A =:

  L \* L\*\*T, if UPLO = 'L',

where U is an upper triangular matrix and L is lower triangular.

This is the block version of the algorithm, calling Level 3 BLAS.

## Usage

``` r
dpotrf(UPLO = "U", N = NULL, A, LDA = NULL)
```

## Arguments

- UPLO:

  a character.

  'U':

  :   Upper triangle of A is stored;

  'L':

  :   Lower triangle of A is stored.

- N:

  an integer. The order of the matrix A. N \>= 0.

- A:

  a big.matrix, dimension (LDA,N).

- LDA:

  an integer. Dimension of the array A. LDA \>= max(1,N).

## Value

updates the big matrix A with the result, INFO is an integer

- = 0::

  successful exit

- \< 0::

  if INFO = -i, the i-th argument had an illegal value

- \> 0::

  if INFO = i, the leading minor of order i is not positive definite,
  and the factorization could not be completed.

Terms laying out of the computed triangle should be discarded.

## Examples

``` r
set.seed(4669)
A = matrix(rnorm(16),4)
B = as.big.matrix(A %*% t(A))
C = A %*% t(A)
chol(C)
#>           [,1]     [,2]       [,3]        [,4]
#> [1,] 0.8840437 0.876387  1.9545734 -0.49320694
#> [2,] 0.0000000 2.139012 -0.8930167 -0.83237775
#> [3,] 0.0000000 0.000000  1.2194774 -0.68016374
#> [4,] 0.0000000 0.000000  0.0000000  0.05426794
dpotrf(UPLO='U', N=4, A=B, LDA=4)
#> [1] 0
D <- B[,]
D[lower.tri(D)]<-0
D
#>           [,1]     [,2]       [,3]        [,4]
#> [1,] 0.8840437 0.876387  1.9545734 -0.49320694
#> [2,] 0.0000000 2.139012 -0.8930167 -0.83237775
#> [3,] 0.0000000 0.000000  1.2194774 -0.68016374
#> [4,] 0.0000000 0.000000  0.0000000  0.05426794
D-chol(C)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0
t(D)%*%D-C
#>               [,1]         [,2]         [,3] [,4]
#> [1,] -1.110223e-16 0.000000e+00 0.000000e+00    0
#> [2,]  0.000000e+00 0.000000e+00 8.326673e-17    0
#> [3,]  0.000000e+00 8.326673e-17 0.000000e+00    0
#> [4,]  0.000000e+00 0.000000e+00 0.000000e+00    0

#' # The big.matrix file backings will be deleted when garbage collected.
rm(A,B,C,D)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1016550 54.3    1995290 106.6         NA  1995290 106.6
#> Vcells 1826359 14.0    8388608  64.0      65536  2724101  20.8
```
