<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigalgebra <img src="man/figures/logo.png" align="right" width="200"/>

# Arithmetic routines for native R matrices and big.matrix objects
## Frédéric Bertrand, Michael J. Kane, Bryan Lewis, John W. Emerson

<https://doi.org/10.32614/CRAN.package.bigalgebra>

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/doi-10.32614/CRAN.package.bigalgebra-blue.svg)](https://doi.org/10.32614/CRAN.package.bigalgebra)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/fbertran/bigalgebra/branch/master/graph/badge.svg)](https://app.codecov.io/gh/fbertran/bigalgebra?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/bigalgebra)](https://cran.r-project.org/package=bigalgebra)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/bigalgebra)](https://cran.r-project.org/package=bigalgebra)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/bigalgebra?style=social)](https://github.com/fbertran/bigalgebra)
[![R-CMD-check](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`bigalgebra` provides fast linear algebra primitives that operate seamlessly on
base `matrix` objects and [`bigmemory::big.matrix`] containers. The package
wraps BLAS and LAPACK routines with R-friendly helpers so that vector updates,
matrix products, and classic decompositions work the same way in memory or on
disk.

## Package highlights

* **big.matrix workflows** – Guidance on creating, sharing, and cleaning up
  file-backed matrices is collected in the
  [Working with big.matrix Objects vignette](https://fbertran.github.io/bigalgebra/articles/big-matrix-workflows.html).
* **Vector kernels** – Level 1 BLAS-style helpers such as `dset()`, `dsub()` and
  `ddot()` extend familiar vector algebra to `big.matrix` inputs. See the
  [Level 1 BLAS-Style Helpers vignette](https://fbertran.github.io/bigalgebra/articles/level-1-blas-style-helpers.html).
* **Matrix products** – Wrappers including `dgemm()` and `dsymm()` expose Level 3
  BLAS routines for dense matrix multiplication with optional file-backed
  outputs. Explore the
  [Matrix Wrapper Helpers vignette](https://fbertran.github.io/bigalgebra/articles/matrix-wrapper-helpers.html).
* **LAPACK decompositions** – QR, Cholesky, eigenvalue, and SVD helpers (`dgeqrf()`,
  `dpotrf()`, `dgeev()`, `dgesdd()`) bring advanced factorisations to large
  datasets. Walk through the
  [LAPACK Decompositions vignette](https://fbertran.github.io/bigalgebra/articles/lapack-decompositions.html).

## Package options

The package defines a number of global options that begin with `bigalgebra`:

Option  Default value
* `bigalgebra.temp_pattern` with default `matrix_`
* `bigalgebra.tempdir` with default `tempdir`
* `bigalgebra.mixed_arithmetic_returns_R_matrix` with default `TRUE`
* `bigalgebra.DEBUG` with default `FALSE`

The `bigalgebra.tempdir` option must be a function that returns a temporary
directory path used to store big matrix results of BLAS and LAPACK operations.
The default value is simply the base R `tempdir()` function.

The `bigalgebra.temp_pattern` option is a name prefix for file names of generated
big matrix objects output as a result of BLAS and LAPACK operations.

The `bigalgebra.mixed_arithmetic_returns_R_matrix` option determines whether
arithmetic operations involving an R matrix or vector and a `big.matrix`
matrix or vector return a big matrix (when the option is `FALSE`), or return a
normal R matrix (`TRUE`).

## BLAS and LAPACK backends

The package is built, by default, with R's native BLAS libraries, which use
32-bit signed integer indexing. The default build is limited to vectors of at
most 2^31 − 1 entries and matrices with at most 2^31 − 1 rows and 2^31 − 1
columns (note that standard R matrices are limited to 2^31 − 1 total entries).

The package includes a reference BLAS implementation that supports 64-bit
integer indexing, relaxing the limitation on vector lengths and matrix row and
column limits. Installation of this package with the 64-bit reference BLAS
implementation may be performed from the command-line install:

`REFBLAS=1 R CMD INSTALL bigalgebra`

where `bigalgebra` is the source package (for example, `bigalgebra_0.9.0.tar.gz`).

The package may also be built with user-supplied external BLAS and LAPACK
libraries, in either 32- or 64-bit varieties. This is an advanced topic that
requires additional `Makevars` modification, and may include adjustment of the
low-level calling syntax depending on the library used.

Feel free to contact us for help installing and running the package.

This website, the unit tests, some C code fixes and improvements as well as
these examples were created by F. Bertrand.

Maintainer: Frédéric Bertrand <frederic.bertrand@lecnam.net>.

## Installation

You can install the released version of bigalgebra from [CRAN](https://CRAN.R-project.org) with:


``` r
install.packages("bigalgebra")
```

You can install the development version of bigalgebra from [GitHub](https://github.com) with:


``` r
devtools::install_github("fbertran/bigalgebra")
```

## Quick tour of the functionality

The snippets below mirror the worked examples in the vignettes and show how the
helpers behave with in-memory and file-backed matrices.

### Level 1 BLAS helpers


``` r
library(bigmemory)
library(bigalgebra)

x <- bigmemory::big.matrix(5, 1, init = 0)
dset(ALPHA = 3, X = x)

y <- bigmemory::big.matrix(5, 1, init = 1)
dvcal(ALPHA = 0.5, X = x, BETA = 2, Y = y)
y[]
#> [1] 3.5 3.5 3.5 3.5 3.5
```

### Matrix products with `dgemm()`


``` r
A <- bigmemory::big.matrix(5, 4, init = 1)
B <- bigmemory::big.matrix(4, 4, init = 2)
C <- bigmemory::big.matrix(5, 4, init = 0)

dgemm(A = A, B = B, C = C, ALPHA = 1, BETA = 0)
C[]
#>      [,1] [,2] [,3] [,4]
#> [1,]    8    8    8    8
#> [2,]    8    8    8    8
#> [3,]    8    8    8    8
#> [4,]    8    8    8    8
#> [5,]    8    8    8    8
```

### LAPACK decompositions


``` r
set.seed(1)
M <- matrix(rnorm(9), 3)
SPD <- crossprod(M)
SPD_big <- as.big.matrix(SPD)
dpotrf(A = SPD_big)
#> [1] 0
chol_factor <- SPD_big[,]
chol_factor[lower.tri(chol_factor)] <- 0
chol_factor
#>          [,1]       [,2]       [,3]
#> [1,] 1.060398 -0.2388263 -0.6138286
#> [2,] 0.000000  1.8082109  0.2222424
#> [3,] 0.000000  0.0000000  0.8294922
```

### File-backed `big.matrix` workflows


``` r
tmpdir <- tempdir()
file_big <- filebacked.big.matrix(3, 3, init = diag(3),
                                  backingpath = tmpdir,
                                  backingfile = "example.bin")
#> Warning in filebacked.big.matrix(3, 3, init = diag(3), backingpath = tmpdir, : No
#> descriptor file given, it will be named example.bin.desc
#> Error in filebacked.big.matrix(3, 3, init = diag(3), backingpath = tmpdir, : Backing file already exists! Either remove or specify
#>            different backing file
file_big[1, 3] <- 5
#> Error: object 'file_big' not found
file_big[]
#> Error: object 'file_big' not found
rm(file_big)
#> Warning in rm(file_big): object 'file_big' not found
gc()
#>           used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
#> Ncells 1017419 54.4    1698357 90.8         NA  1698357 90.8
#> Vcells 2413994 18.5    8388608 64.0      65536  3432096 26.2
```

## Available vignettes

The full vignette set expands on the topics above and demonstrates how the
routines interact:

1. **Working with big.matrix Objects** – Managing shared memory, file backing,
   and clean-up for large datasets.
2. **Level 1 BLAS-Style Helpers** – Filling vectors, Hadamard products, and
   reductions on disk-backed data.
3. **Matrix Wrapper Helpers** – Symmetric and general matrix products, including
   strategies for chaining operations.
4. **LAPACK Decompositions with bigalgebra** – QR, Cholesky, eigenvalue, and SVD
   workflows for both in-memory and file-backed matrices.
