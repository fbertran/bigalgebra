# Add two double-precision vectors.

Compute double precision DY := DX + DY for N elements, applying the
specified storage increments. This routine mirrors the BLAS DAXPY
operation with a unit scaling factor.

## Usage

``` r
dadd(N = NULL, X, INCX = 1L, Y, INCY = 1L)
```

## Arguments

- N:

  number of elements in the input vectors. Defaults to the length of `X`
  if `NULL`.

- X:

  double precision vector or matrix providing the addend.

- INCX:

  storage spacing between elements of `X`.

- Y:

  double precision vector or matrix that is updated in place.

- INCY:

  storage spacing between elements of `Y`.

## Value

The modified object `Y` containing the element-wise sums.

## Details

The implementation delegates to the BLAS `DAXPY` routine with a unit
scaling factor, making it equivalent to `daxpy(1.0, X, Y)` while
exposing an interface consistent with other low-level wrappers such as
`dcopy` and `dscal`.

## Examples

``` r
set.seed(4669)
X <- matrix(runif(6), 3, 2)
Y <- matrix(runif(6), 3, 2)
dadd(X = X, Y = Y)
#>           [,1]     [,2]
#> [1,] 0.8373056 1.273518
#> [2,] 0.9437275 1.347313
#> [3,] 1.0879022 1.192258
all.equal(Y, X + Y)
#> [1] "Mean relative difference: 0.5023438"
```
