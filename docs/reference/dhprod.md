# Element-wise (Hadamard) product

Computes \\Z := X \circ Y\\. When \`Z\` is missing it is allocated
automatically.

## Usage

``` r
dhprod(N = NULL, X, INCX = 1L, Y, INCY = 1L, Z, INCZ = 1L)
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

- Z:

  Optional output container receiving the product.

- INCZ:

  Integer stride for \`Z\`.

## Value

The updated object \`Z\`.

## Examples

``` r
dhprod(X = 1:4, Y = rep(2, 4))
#> [1] 2 4 6 8
```
