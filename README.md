<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigalgebra <img src="man/figures/logo.png" align="right" width="200"/>

# Arithmetic routines for native R matrices and big.matrix objects
## Frédéric Bertrand, Michael J. Kane, Bryan Lewis, John W. Emerson


<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/fbertran/bigalgebra/workflows/R-CMD-check/badge.svg)](https://github.com/fbertran/bigalgebra/actions)
[![Codecov test coverage](https://codecov.io/gh/fbertran/bigalgebra/branch/master/graph/badge.svg)](https://app.codecov.io/gh/fbertran/bigalgebra?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/bigalgebra)](https://cran.r-project.org/package=bigalgebra)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/bigalgebra)](https://cran.r-project.org/package=bigalgebra)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/bigalgebra?style=social)](https://github.com/fbertran/bigalgebra)
[![DOI](https://zenodo.org/badge/353292865.svg)](https://zenodo.org/badge/latestdoi/353292865)

<!-- badges: end -->

This package provides arithmetic functions for native `R` matrices and `bigmemory::big.matrix` objects  as well as functions for QR factorization, Cholesky factorization, General eigenvalue, and Singular value decomposition (SVD). A method matrix multiplication and an arithmetic method -for matrix addition, matrix difference- allows for mixed type operation -a matrix class object and a big.matrix class object- and pure type operation for two big.matrix class objects.


The package defines a number of global options that begin with `bigalgebra`.

They include:

Option  Default value
* `bigalgebra.temp_pattern` with default `matrix_`
* `bigalgebra.tempdir` with default `tempdir`
* `bigalgebra.mixed_arithmetic_returns_R_matrix` with default `TRUE`
* `bigalgebra.DEBUG` with default `FALSE`

The `bigalgebra.tempdir` option must be a function that returns
a temporary directory path used to big matrix results of BLAS and
LAPACK operations. The deault value is simply the default R `tempdir`
function.

The `bigalgebra.temp_pattern` is a name prefix for file names of generated
big matrix objects output as a result of BLAS and LAPACK operations.

The `bigalgebra.mixed_arithmetic_returns_R_matrix` option determines
whether arithmetic operations involving an R matrix or vector and a big.matrix
matrix or vector return a big matrix (when the option is `FALSE`), or
return a normal R matrix (`TRUE`).

The package is built, by default, with `R`'s native BLAS libraries, which use
32-bit signed integer indexing. The default build is limited to vectors of at
most 2^31 - 1 entries and matrices with at most 2^31 - 1 rows and 2^31 - 1
columns (note that standard R matrices are limtied to 2^31 - 1 total entries).

The package includes a reference BLAS implementation that supports 64-bit
integer indexing, relaxing the limitation on vector lengths and matrix
row and column limits. Installation of this package with the 64-bit reference
BLAS implementation may be performed from the command-line install:

`REFBLAS=1 R CMD INSTALL bigalgebra`

where `bigalgebra` is the source package (for example, `bigalgebra_0.9.0.tar.gz`).

The package may also be build with user-supplied external BLAS and LAPACK
libraries, in either 32- or 64-bit varieties. This is an advanced topic
that requires additional Makevars modification, and may include adjustment
of the low-level calling syntax depending on the library used.

Feel free to contact us for help installing and running the package.


This website and these examples were created by F. Bertrand.

Maintainer: Frédéric Bertrand <frederic.bertrand@utt.fr>.


## Installation

You can install the released version of bigalgebra from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("bigalgebra")
```

You can install the development version of bigalgebra from [github](https://github.com) with:


```r
devtools::install_github("fbertran/bigalgebra")
```


## Examples


```r
library("bigmemory")
A <- bigmemory::big.matrix(5,4,init = 1)
B <- bigmemory::big.matrix(4,4,init = 2)

C <- A %*% B       # Returns a new big.matrix object
#> Error in A %*% B: nécessite des arguments numériques/complexes matrice/vecteur
D <- A[] %*% B[]   # Compute the same thing in R

print(C - D)       # Compare the results (subtraction of an R matrix from a
#> Error in h(simpleError(msg, call)): erreur d'ï¿½valuation de l'argument 'x' lors de la sï¿½lection d'une mï¿½thode pour la fonction 'print' : argument non numérique pour un opérateur binaire
                   # big.matrix)

# The next example illustrates mixing R and big.matrix objects. It returns by
# default (see # options("bigalgebra.mixed_arithmetic_returns_R_matrix")
D <- matrix(rnorm(16),4)
E <- A %*% D
#> Error in A %*% D: nécessite des arguments numériques/complexes matrice/vecteur
```

