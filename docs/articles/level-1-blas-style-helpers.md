# Level 1 BLAS-Style Helpers

## Overview

Level 1 BLAS routines operate on vectors and provide the building blocks
for more elaborate linear algebra workflows. The `bigalgebra` package
exposes a collection of convenience wrappers that mirror these classic
operations while extending them to work seamlessly with
\[[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)\]
objects. This vignette introduces each helper, illustrating how to use
them with both base R vectors and matrices backed by disk.

## Filling and combining vectors

The helpers
[`dset()`](https://fbertran.github.io/bigalgebra/reference/dset.md),
[`dvcal()`](https://fbertran.github.io/bigalgebra/reference/dvcal.md)
and [`dsub()`](https://fbertran.github.io/bigalgebra/reference/dsub.md)
all operate in place, making it easy to reuse pre-allocated storage when
working with large data sets.

``` r

# Fill an existing matrix with a scalar
X <- matrix(0, 2, 3)
dset(ALPHA = 3.5, X = X)
X
#>      [,1] [,2] [,3]
#> [1,]  3.5  3.5  3.5
#> [2,]  3.5  3.5  3.5

# Form alpha * X + beta * Y in place
dx <- as.numeric(1:4)
dy <- rep(2, 4)
dvcal(ALPHA = 2, X = dx, BETA = -1, Y = dy)
dy
#> [1] 0 2 4 6

# Subtract X from Y element-wise
dy2 <- rep(10, 4)
dsub(X = dx, Y = dy2)
dy2
#> [1] 9 8 7 6
```

When either argument is a `big.matrix`, the wrapper automatically keeps
the data on disk while still updating it in place:

``` r

library(bigmemory)
# Allocate a file-backed big.matrix
dir.create(tmp <- tempfile())
X_bm <- filebacked.big.matrix(4, 1, type = "double",
                              backingfile = "vec.bin",
                              backingpath = tmp,
                              descriptorfile = "vec.desc")
X_bm[] <- 0:3

dvcal(ALPHA = -1, X = X_bm, BETA = 2, Y = X_bm)
X_bm[]
#> [1] 0 1 2 3
```

## Dot products and element-wise operations

[`ddot()`](https://fbertran.github.io/bigalgebra/reference/ddot.md)
evaluates the standard dot product. For problems where numerical
stability is paramount,
[`dqddot()`](https://fbertran.github.io/bigalgebra/reference/dqddot.md)
accumulates in extended precision. The helper
[`dhprod()`](https://fbertran.github.io/bigalgebra/reference/dhprod.md)
forms the Hadamard (element-wise) product and returns the populated
output container. When you need to transform each entry of a vector or
matrix in place,
[`dsqrt()`](https://fbertran.github.io/bigalgebra/reference/dsqrt.md)
applies the square root while preserving the original storage mode.

``` r

# Classic dot product and its extended-precision counterpart
v1 <- as.numeric(1:5)
v2 <- seq(2, 10, by = 2)
list(ddot = ddot(X = v1, Y = v2), dqddot = dqddot(X = v1, Y = v2))
#> $ddot
#> [1] 110
#> 
#> $dqddot
#> [1] 110

# Hadamard product stored in a new matrix
A <- matrix(as.numeric(1:4), 2, 2)
B <- matrix(rep(2, 4), 2, 2)
Z <- dhprod(X = A, Y = B)
Z
#>      [,1] [,2]
#> [1,]    2    6
#> [2,]    4    8

# Element-wise square root performed in place
sqrt_vals <- matrix(c(1, 4, 9, 16), 2)
dsqrt(X = sqrt_vals)
sqrt_vals
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```

For three-dimensional vectors the cross product helper
[`dxyz()`](https://fbertran.github.io/bigalgebra/reference/dxyz.md)
offers a specialised convenience routine:

``` r

ux <- c(1, 0, 0)
uy <- c(0, 1, 0)
dxyz(X = ux, Y = uy)
#> [1] 0 0 1
```

## Reduction helpers

The reduction helpers mirror the Level 1 BLAS routines for computing
sums, absolute sums, Euclidean norms and cumulative products. They
return numeric scalars and work with strided access patterns via the
`INCX` parameter when needed.

``` r

vals <- c(-1, 2, -3, 4)
list(
  sum = dsum(X = vals),
  abs_sum = dasum(X = vals),
  euclidean_norm = dnrm2(X = vals),
  product = dprdct(X = vals)
)
#> $sum
#> [1] 2
#> 
#> $abs_sum
#> [1] 10
#> 
#> $euclidean_norm
#> [1] 5.477226
#> 
#> $product
#> [1] 24
```

The `id*` family returns one-based indices of extrema based on value or
absolute value criteria.

``` r

idx_vals <- c(-2, 5, -7, 3)
list(
  min_index = idmin(X = idx_vals),
  max_index = idmax(X = idx_vals),
  min_abs_index = idamin(X = idx_vals),
  max_abs_index = idamax(X = idx_vals)
)
#> $min_index
#> [1] 3
#> 
#> $max_index
#> [1] 2
#> 
#> $min_abs_index
#> [1] 1
#> 
#> $max_abs_index
#> [1] 3
```

Level 1 helpers can be combined with higher-level matrix routines to
implement more complex algorithms while retaining the performance
benefits of the underlying compiled implementations.

``` r

unlink(file.path(tmp, "vec.bin"))
unlink(file.path(tmp, "vec.desc"))
unlink(tmp, recursive = TRUE)
```
