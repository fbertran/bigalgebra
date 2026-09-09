# Matrix Wrapper Helpers

## Overview

Beyond its Level 1 vector helpers, `bigalgebra` includes wrappers around
the Level 3 BLAS routines and related utilities for manipulating
matrices stored in memory or as
\[[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)\]
objects. This vignette demonstrates the core workflows for the symmetric
matrix-matrix multiply
[`dsymm()`](https://fbertran.github.io/bigalgebra/reference/dsymm.md),
the general matrix multiply
[`dgemm()`](https://fbertran.github.io/bigalgebra/reference/dgemm.md)
and the AXPY-style updater
[`daxpy()`](https://fbertran.github.io/bigalgebra/reference/daxpy.md).

## Symmetric matrix products with `dsymm()`

[`dsymm()`](https://fbertran.github.io/bigalgebra/reference/dsymm.md)
mirrors the BLAS symmetric matrix-matrix multiplication kernel. It can
multiply from the left or right and accepts optional leading-dimension
arguments when working with non-standard strides.

``` r

A_sym <- matrix(c(2, 1, 1, 3), 2, 2)
B_rhs <- diag(2)
C_out <- matrix(0, 2, 2)
dsymm(A = A_sym, B = B_rhs, C = C_out, SIDE = "L", UPLO = "U")
C_out
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    3
```

When either input is a `big.matrix`, the wrapper performs the
computation in place without materialising an intermediate R matrix.

``` r

library(bigmemory)
dir.create(tmp_sym <- tempfile())
A_bm <- filebacked.big.matrix(2, 2, type = "double",
                              backingpath = tmp_sym,
                              backingfile = "A.bin",
                              descriptorfile = "A.desc")
A_bm[,] <- A_sym
B_bm <- filebacked.big.matrix(2, 2, type = "double",
                              backingpath = tmp_sym,
                              backingfile = "B.bin",
                              descriptorfile = "B.desc")
B_bm[,] <- B_rhs
C_bm <- filebacked.big.matrix(2, 2, type = "double",
                              backingpath = tmp_sym,
                              backingfile = "C.bin",
                              descriptorfile = "C.desc")

dsymm(A = A_bm, B = B_bm, C = C_bm)
C_bm[]
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    3
```

## General matrix multiplication with `dgemm()`

[`dgemm()`](https://fbertran.github.io/bigalgebra/reference/dgemm.md)
exposes the workhorse Level 3 BLAS routine that computes
`C := alpha * op(A) %*% op(B) + beta * C`. All matrices can be ordinary
R matrices or `big.matrix` objects, and any missing output will be
created automatically. When you construct native R matrices, coerce
integer literals to doubles (for example with
[`as.numeric()`](https://rdrr.io/r/base/numeric.html)) so they satisfy
the package’s double-precision requirement.

``` r

A <- matrix(as.numeric(1:6), nrow = 2)
B <- matrix(seq(2, 12, by = 2), nrow = 3)
C <- matrix(0, nrow = 2, ncol = 4)
dgemm(TRANSA = "N", TRANSB = "N", A = A, B = B, C = C, BETA = 0)
C
#>      [,1] [,2] [,3] [,4]
#> [1,]   44   98    0    0
#> [2,]   56  128    0    0
```

Transposed operands are also supported via `TRANSA` and `TRANSB`:

``` r

C_t <- matrix(0, nrow = 3, ncol = 3)
dgemm(TRANSA = "T", TRANSB = "N", A = A, B = B, C = C_t, BETA = 0)
C_t
#>      [,1] [,2] [,3]
#> [1,]   10   22    0
#> [2,]   22   50    0
#> [3,]   34   78    0
```

## Updating matrices with `daxpy()`

[`daxpy()`](https://fbertran.github.io/bigalgebra/reference/daxpy.md)
provides a convenient wrapper for the AXPY operation
`Y := alpha * X + Y`. The helper respects the package option
`bigalgebra.mixed_arithmetic_returns_R_matrix` to decide when to return
a native matrix versus a `big.matrix`.

``` r

X <- matrix(1, nrow = 2, ncol = 2)
Y <- matrix(c(0, 1, 2, 3), nrow = 2)
daxpy(A = 0.5, X = X, Y = Y)
#> An object of class "big.matrix"
#> Slot "address":
#> <pointer: 0x155fcf7c0>
```

When both operands are `big.matrix` objects, the result stays on disk:

``` r

dir.create(tmp_axpy <- tempfile())
X_bm <- filebacked.big.matrix(2, 2, type = "double",
                              backingpath = tmp_axpy,
                              backingfile = "X.bin",
                              descriptorfile = "X.desc")
X_bm[,] <- 1
daxpy(A = 3, X = X_bm)
#> An object of class "big.matrix"
#> Slot "address":
#> <pointer: 0x155fdb6b0>
X_bm[]
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    1
```

``` r

unlink(file.path(tmp_sym, c("A.bin", "A.desc", "B.bin", "B.desc", "C.bin", "C.desc")))
unlink(tmp_sym, recursive = TRUE)
unlink(file.path(tmp_axpy, c("X.bin", "X.desc")))
unlink(tmp_axpy, recursive = TRUE)
```
