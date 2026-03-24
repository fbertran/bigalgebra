# Swap two double-precision vectors.

Exchange the elements of two double precision vectors in place. For I =
0 to N-1, swap DX(LX + I \* INCX) with DY(LY + I \* INCY) where LX and
LY depend on the increment signs. When an optimized BLAS is available
the implementation dispatches to `DSWAP`; otherwise a portable C loop
performs the exchange.

## Usage

``` r
dswap(N = NULL, X, INCX = 1L, Y, INCY = 1L)
```

## Arguments

- N:

  Number of elements in the input vectors. Defaults to the full length
  of `X` if `NULL`.

- X:

  Double precision vector or matrix providing the first data block.

- INCX:

  Storage spacing between elements of `X`.

- Y:

  Double precision vector or matrix providing the second data block.

- INCY:

  Storage spacing between elements of `Y`.

## Value

Invisibly returns `NULL`; both `X` and `Y` are modified in place.

## Details

When an optimized BLAS is available the implementation delegates to the
Fortran `DSWAP` routine. Otherwise a portable C fallback performs the
exchange directly while respecting the supplied vector increments.

## See also

\[dcopy()\], \[daxpy()\] and \[dscal()\].

## Examples

``` r
set.seed(4670)
X <- matrix(runif(6), 3, 2)
Y <- matrix(runif(6), 3, 2)
X_original <- X
Y_original <- Y
dswap(X = X, Y = Y)
all.equal(X, Y_original)
#> [1] "Mean relative difference: 0.274648"
all.equal(Y, X_original)
#> [1] "Mean relative difference: 0.2364376"
```
