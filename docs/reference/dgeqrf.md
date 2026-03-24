# QR factorization

DGEQRF computes a QR factorization of a real M-by-N matrix A: A = Q \*
R.

## Usage

``` r
dgeqrf(
  M = NULL,
  N = NULL,
  A,
  LDA = NULL,
  TAU = NULL,
  WORK = NULL,
  LWORK = NULL
)
```

## Arguments

- M:

  an integer. The number of rows of the matrix A. M \>= 0.

- N:

  an integer. The number of columns of the matrix A. N \>= 0.

- A:

  the M-by-N big matrix A.

- LDA:

  an integer. The leading dimension of the array A. LDA \>= max(1,M).

- TAU:

  a min(M,N) matrix. The scalar factors of the elementary reflectors.

- WORK:

  a (MAX(1,LWORK)) matrix. On exit, if INFO = 0, WORK(1) returns the
  optimal LWORK.

- LWORK:

  an integer. The dimension of th array WORK.

## Value

M-by-N big matrix A. The elements on and above the diagonal of the array
contain the min(M,N)-by-N upper trapezoidal matrix R (R is upper
triangular if m \>= n); the elements below the diagonal, with the array
TAU, represent the orthogonal matrix Q as a product of min(m,n)
elementary reflectors.

## Examples

``` r
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
h9 <- hilbert(9); h9
#>            [,1]      [,2]       [,3]       [,4]       [,5]       [,6]
#>  [1,] 1.0000000 0.5000000 0.33333333 0.25000000 0.20000000 0.16666667
#>  [2,] 0.5000000 0.3333333 0.25000000 0.20000000 0.16666667 0.14285714
#>  [3,] 0.3333333 0.2500000 0.20000000 0.16666667 0.14285714 0.12500000
#>  [4,] 0.2500000 0.2000000 0.16666667 0.14285714 0.12500000 0.11111111
#>  [5,] 0.2000000 0.1666667 0.14285714 0.12500000 0.11111111 0.10000000
#>  [6,] 0.1666667 0.1428571 0.12500000 0.11111111 0.10000000 0.09090909
#>  [7,] 0.1428571 0.1250000 0.11111111 0.10000000 0.09090909 0.08333333
#>  [8,] 0.1250000 0.1111111 0.10000000 0.09090909 0.08333333 0.07692308
#>  [9,] 0.1111111 0.1000000 0.09090909 0.08333333 0.07692308 0.07142857
#>             [,7]       [,8]       [,9]
#>  [1,] 0.14285714 0.12500000 0.11111111
#>  [2,] 0.12500000 0.11111111 0.10000000
#>  [3,] 0.11111111 0.10000000 0.09090909
#>  [4,] 0.10000000 0.09090909 0.08333333
#>  [5,] 0.09090909 0.08333333 0.07692308
#>  [6,] 0.08333333 0.07692308 0.07142857
#>  [7,] 0.07692308 0.07142857 0.06666667
#>  [8,] 0.07142857 0.06666667 0.06250000
#>  [9,] 0.06666667 0.06250000 0.05882353
qr(h9)$rank           #--> only 7
#> [1] 7
qrh9 <- qr(h9, tol = 1e-10)
qrh9$rank 
#> [1] 9
C <- as.big.matrix(h9)
dgeqrf(A=C)
#> [1] 0

# The big.matrix file backings will be deleted when garbage collected.
rm(C)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1150291 61.5    2286010 122.1         NA  2286010 122.1
#> Vcells 2059149 15.8    8388608  64.0      65536  6482741  49.5
```
