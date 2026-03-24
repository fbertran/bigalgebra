# BLAS daxpy functionality

This function implements the function Y := A \* X + Y where X and Y may
be either native double-precision valued R matrices or numeric vectors,
or double-precision valued
[`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects, and A is a scalar.

## Usage

``` r
daxpy(A = 1, X, Y)
```

## Arguments

- A:

  Optional numeric scalar value to scale the matrix `X` by, with a
  default value of 1.

- X:

  Requried to be either a native R
  [`matrix`](https://rdrr.io/r/base/matrix.html) or numeric vector, or a
  [`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  object

- Y:

  Optional native R [`matrix`](https://rdrr.io/r/base/matrix.html) or
  numeric vector, or a
  [`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  object

## Value

The output value depends on the classes of input values `X` and `Y` and
on the value of the global option
`bigalgebra.mixed_arithmetic_returns_R_matrix`.

If `X` and `Y` are both big matrices, or `Y` is missing,
`options("bigalgebra.mixed_arithmetic_returns_R_matrix")` is `FALSE`,
then a `big.matrix` is returned. The returned `big.matrix` is backed by
a temporary file mapping that will be deleted when the returned result
is garbage collected by R (see the examples).

Otherwise, a standard R matrix is returned. The dimensional shape of the
output is taken from `X`. If input `X` is dimensionless (that is, lacks
a dimension attribute), then the output is a column vector.

## Details

At least one of either `X` or `Y` must be a `big.matrix`. All values
must be of type `double` (the only type presently supported by the
bigalgebra package).

This function is rarely necessary to use directly since the bigalgebra
package defines standard arithmetic operations and scalar
multiplication. It is more efficient to use `daxpy` directly when both
scaling and matrix addition are required, in which case both operations
are performed in one step.

## References

<https://www.netlib.org/blas/daxpy.f>

## See also

[`bigmemory`](https://rdrr.io/pkg/bigmemory/man/bigmemory-package.html)

## Author

Michael J. Kane

## Examples

``` r
require(bigmemory)
A = matrix(1, nrow=3, ncol=2)
B <- big.matrix(3, 2, type="double", init=0,
                dimnames=list(NULL, c("alpha", "beta")), shared=FALSE)
C = B + B   # C is a new big matrix
D = A + B   # D defaults to a regular R matrix, to change this, set the option:
# options(bigalgebra.mixed_arithmetic_returns_R_matrix=FALSE)
E = daxpy(A=1.0, X=B, Y=B)  # Same kind of result as C
print(C[])
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    0    0
print(D)
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    1
#> [3,]    1    1
print(E[])
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    0    0

# The C and E big.matrix file backings will be deleted when garbage collected:
# (We enable debugging to see this explicitly)
options(bigalgebra.DEBUG=TRUE)
rm(C,E)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1147699 61.3    2286010 122.1         NA  2286010 122.1
#> Vcells 2053999 15.7    8388608  64.0      65536  6482741  49.5
```
