# LAPACK Decompositions with bigalgebra

## Overview

Beyond BLAS-level helpers, `bigalgebra` ships several wrappers around
[LAPACK](https://www.netlib.org/lapack/) routines so that factorisations
can run against in-memory `matrix` objects or file-backed
\[[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)\]
containers without changing workflows. This vignette highlights the four
decompositions currently provided:

- [`dgeqrf()`](https://fbertran.github.io/bigalgebra/reference/dgeqrf.md)
  — QR factorisation driven by Householder reflectors.
- [`dpotrf()`](https://fbertran.github.io/bigalgebra/reference/dpotrf.md)
  — Cholesky factorisation of symmetric positive definite matrices.
- [`dgeev()`](https://fbertran.github.io/bigalgebra/reference/dgeev.md)
  — real eigenvalues and (optionally) eigenvectors of general matrices.
- [`dgesdd()`](https://fbertran.github.io/bigalgebra/reference/dgesdd.md)
  — divide-and-conquer singular value decomposition (SVD).

Each example illustrates how to prepare workspace arguments, inspect the
results, and clean up temporary `big.matrix` files when necessary.

## QR factorisation with `dgeqrf()`

The
[`dgeqrf()`](https://fbertran.github.io/bigalgebra/reference/dgeqrf.md)
helper overwrites its input with the R factor in the upper triangle
while storing Householder reflector coefficients in a user-supplied
`TAU` vector. Supplying both `TAU` and `WORK` explicitly makes it easy
to inspect the outputs afterwards.

``` r
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
H <- hilbert(4)
H_big <- as.big.matrix(H)
TAU <- matrix(0, nrow = min(nrow(H), ncol(H)))
WORK <- matrix(0, nrow = max(1, ncol(H)))
dgeqrf(A = H_big, TAU = TAU, WORK = WORK)
#> [1] 0

# Extract the R factor from the overwritten big.matrix
R_big <- H_big[,]
R_big[lower.tri(R_big)] <- 0
R_big
#>           [,1]       [,2]         [,3]          [,4]
#> [1,] -1.193152 -0.6704931 -0.474932601 -0.3698354709
#> [2,]  0.000000 -0.1185333 -0.125655095 -0.1175419928
#> [3,]  0.000000  0.0000000 -0.006221774 -0.0095660929
#> [4,]  0.000000  0.0000000  0.000000000  0.0001879049

# Compare against base R's QR decomposition
all.equal(R_big, qr.R(qr(H)))
#> [1] TRUE
TAU
#>          [,1]
#> [1,] 1.838116
#> [2,] 1.560868
#> [3,] 1.413383
#> [4,] 0.000000
```

When the input is file-backed, the updates are written directly to disk
and the factorisation scales to matrices that exceed available RAM.

``` r
tmp <- tempfile()
H_fb <- filebacked.big.matrix(nrow(H), ncol(H), init = H,
                              backingfile = basename(tmp),
                              backingpath = dirname(tmp))
#> Warning in filebacked.big.matrix(nrow(H), ncol(H), init = H, backingfile =
#> basename(tmp), : No descriptor file given, it will be named
#> fileb50174135d2.desc
dgeqrf(A = H_fb)
#> [1] 0
H_fb[,][upper.tri(H_fb[,], diag = TRUE)]
#>  [1] -2.000000e+00 -2.000000e+00 -9.614813e-17 -2.000000e+00 -9.614813e-17
#>  [6] -1.739150e-33 -2.000000e+00 -9.614813e-17 -1.739150e-33  4.056603e-49
rm(H_fb); gc()
#>           used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
#> Ncells  816477 43.7    1543642 82.5         NA  1543642 82.5
#> Vcells 1481039 11.3    8388608 64.0      65536  2969954 22.7
```

## Cholesky factorisation with `dpotrf()`

[`dpotrf()`](https://fbertran.github.io/bigalgebra/reference/dpotrf.md)
computes an in-place Cholesky factorisation. The helper returns `0` when
the matrix is positive definite and leaves the result in the selected
triangle (upper by default).

``` r
set.seed(42)
A <- matrix(rnorm(9), 3)
SPD <- crossprod(A)  # symmetric positive definite
SPD_big <- as.big.matrix(SPD)
info <- dpotrf(UPLO = "U", A = SPD_big)
info
#> [1] 0

U <- SPD_big[,]
U[lower.tri(U)] <- 0
all.equal(U, chol(SPD))
#> [1] TRUE
```

If the input matrix is not positive definite, the return value indicates
the leading minor that failed so you can diagnose numerical issues
before proceeding with downstream solves.

## Eigenvalues and eigenvectors via `dgeev()`

[`dgeev()`](https://fbertran.github.io/bigalgebra/reference/dgeev.md)
wraps LAPACK’s general eigen solver and accepts optional storage for
left and right eigenvectors. By default the helper queries LAPACK for an
optimal workspace size, making small examples straightforward.

``` r
set.seed(123)
M <- matrix(rnorm(16), 4)
WR <- matrix(0, nrow = ncol(M))
WI <- matrix(0, nrow = ncol(M))
VL <- matrix(0, nrow = nrow(M), ncol = ncol(M))
dgeev(A = M, WR = WR, WI = WI, VL = VL)
#> [1] 0

# Compare eigenvalues with base R
complex_eigs <- WR[, 1] + 1i * WI[, 1]
all.equal(sort(complex_eigs), sort(eigen(M)$values))
#> [1] TRUE
```

For large matrices stored on disk you can pass `big.matrix` containers
for `A` and (optionally) `VL`/`VR`. The wrapper automatically converts
between big storage and native `matrix` arguments as required.

## Singular value decomposition with `dgesdd()`

The divide-and-conquer SVD routine
[`dgesdd()`](https://fbertran.github.io/bigalgebra/reference/dgesdd.md)
requires workspace for the singular values and (optionally) the left and
right singular vectors. Providing `big.matrix` containers lets you
persist decompositions to disk without copying through R’s heap.

``` r
set.seed(101)
X <- matrix(rnorm(12), 4)
S <- matrix(0, nrow = min(dim(X)))
U <- matrix(0, nrow = nrow(X), ncol = nrow(X))
VT <- matrix(0, nrow = ncol(X), ncol = ncol(X))
dgesdd(A = X, S = S, U = U, VT = VT)
#> [1] 0

svd_base <- svd(X)
all.equal(drop(S), svd_base$d)
#> [1] "Mean relative difference: 0.3311388"
all.equal(U, svd_base$u)
#> [1] "Attributes: < Component \"dim\": Mean relative difference: 0.25 >"
#> [2] "Numeric: lengths (16, 12) differ"
all.equal(VT, t(svd_base$v))
#> [1] "Mean relative difference: 1.484757"
```

Because
[`dgesdd()`](https://fbertran.github.io/bigalgebra/reference/dgesdd.md)
accepts both in-memory and file-backed matrices, it enables scalable
dimensionality reduction pipelines directly on datasets that are larger
than available RAM.
