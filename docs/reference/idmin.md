# Index of the minimum element

Index of the minimum element

## Usage

``` r
idmin(N = NULL, X, INCX = 1L)
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

Integer index (1-based) of the smallest entry in \`X\`.
