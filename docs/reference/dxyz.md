# Three-dimensional cross product

Three-dimensional cross product

## Usage

``` r
dxyz(X, Y, Z)
```

## Arguments

- X:

  Numeric vector of length three, matrix with three rows, or big.matrix.

- Y:

  Numeric object matching the shape of \`X\`.

- Z:

  Optional output container.

## Value

The updated object \`Z\` containing the cross product.

## Examples

``` r
dxyz(X = c(1, 0, 0), Y = c(0, 1, 0))
#> [1] 0 0 1
```
