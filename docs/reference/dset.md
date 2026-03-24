# Fill a vector or matrix with a constant value

Fill a vector or matrix with a constant value

## Usage

``` r
dset(N = NULL, ALPHA, X, INCX = 1L)
```

## Arguments

- N:

  Optional integer specifying the number of elements to modify. Defaults
  to the length of \`X\`.

- ALPHA:

  Numeric scalar used to populate \`X\`.

- X:

  Double-precision vector, matrix or \[\`bigmemory::big.matrix\`\] to be
  filled in place.

- INCX:

  Integer stride between successive elements of \`X\`.

## Value

Invisibly returns \`X\` after modification.

## Examples

``` r
x <- matrix(0, 2, 3)
dset(ALPHA = 5, X = x)
x
#>      [,1] [,2] [,3]
#> [1,]    5    5    5
#> [2,]    5    5    5
```
