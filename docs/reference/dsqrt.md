# Element-wise square root

Applies the square root to each entry of \`X\` in place, supporting both
base R and \[\`bigmemory::big.matrix\`\] inputs.

## Usage

``` r
dsqrt(N = NULL, X, INCX = 1L)
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

## Value

Invisibly returns the modified object \`X\`.

## Examples

``` r
vals <- matrix(c(1, 4, 9, 16), 2)
dsqrt(X = vals)
vals
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```
