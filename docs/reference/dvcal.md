# Generalised AXPBY operation

Computes the linear combination \\Y := \alpha X + \beta Y\\ in place.

## Usage

``` r
dvcal(N = NULL, ALPHA = 1, X, INCX = 1L, BETA = 1, Y, INCY = 1L)
```

## Arguments

- N:

  Optional integer giving the number of elements. Defaults to
  \`length(X)\`.

- ALPHA:

  Numeric scalar multiplying \`X\`.

- X:

  Double-precision vector, matrix or \[\`bigmemory::big.matrix\`\]
  input.

- INCX:

  Integer stride for traversing \`X\`.

- BETA:

  Numeric scalar multiplying \`Y\`.

- Y:

  Double-precision object updated in place.

- INCY:

  Integer stride for traversing \`Y\`.

## Value

Invisibly returns the modified \`Y\`.

## Examples

``` r
x <- 1:5
y <- rep(2, 5)
dvcal(ALPHA = 2, X = x, BETA = -1, Y = y)
y
#> [1] 0 2 4 6 8
```
