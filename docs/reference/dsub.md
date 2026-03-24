# In-place vector subtraction

Updates \`Y\` by subtracting \`X\`, i.e. \\Y := Y - X\\.

## Usage

``` r
dsub(N = NULL, X, INCX = 1L, Y, INCY = 1L)
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

Invisibly returns the modified \`Y\`.

## Examples

``` r
x <- 1:4
y <- rep(10, 4)
dsub(X = x, Y = y)
y
#> [1] 9 8 7 6
```
