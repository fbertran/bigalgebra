# High-accuracy dot product

Forms the dot product using long-double accumulation to mitigate
rounding error.

## Usage

``` r
dqddot(N = NULL, X, INCX = 1L, Y, INCY = 1L)
```

## Arguments

- N:

  Optional integer giving the number of elements. Defaults to
  \`length(X)\`.

- X:

  Double-precision vector, matrix or \[\`bigmemory::big.matrix\`\]
  input.

- INCX:

  Integer stride for traversing \`X\`.

- Y:

  Double-precision object updated in place.

- INCY:

  Integer stride for traversing \`Y\`.

## Value

Numeric scalar equal to the dot product of \`X\` and \`Y\`.

## Examples

``` r
dqddot(X = 1:3, Y = c(2, 4, 6))
#> [1] 28
```
