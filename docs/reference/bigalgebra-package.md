# Arithmetic routines for native R matrices and big.matrix objects.

This package provides arithmetic functions for native R matrices and
[`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects.

## Details

This package provides arithmetic functions for native R matrices and
[`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects.

The package defines a number of global options that begin with
`bigalgebra`. They include:

|                                                |               |
|------------------------------------------------|---------------|
| Option                                         | Default value |
| `bigalgebra.temp_pattern`                      | "matrix\_"    |
| `bigalgebra.tempdir`                           | `tempdir`     |
| `bigalgebra.mixed_arithmetic_returns_R_matrix` | TRUE          |
| `bigalgebra.DEBUG`                             | FALSE         |

The `bigalgebra.tempdir` option must be a function that returns a
temporary directory path used to big matrix results of BLAS and LAPACK
operations. The deault value is simply the default R `tempdir` function.

The `bigalgebra.temp_pattern` is a name prefix for file names of
generated big matrix objects output as a result of BLAS and LAPACK
operations.

The `bigalgebra.mixed_arithmetic_returns_R_matrix` option determines
whether arithmetic operations involving an R matrix or vector and a
big.matrix matrix or vector return a big matrix (when the option is
`FALSE`), or return a normal R matrix (`TRUE`).

The package is built, by default, with R's native BLAS libraries, which
use 32-bit signed integer indexing. The default build is limited to
vectors of at most 2\*\*31 - 1 entries and matrices with at most
2\*\*31 - 1 rows and 2\*\*31 - 1 columns (note that standard R matrices
are limtied to 2\*\*31 - 1 total entries).

The package includes a reference BLAS implementation that supports
64-bit integer indexing, relaxing the limitation on vector lengths and
matrix row and column limits. Installation of this package with the
64-bit reference BLAS implementation may be performed from the
command-line install:

`REFBLAS=1 R CMD INSTALL bigalgebra`

where "`bigalgebra`" is the source package (for example,
`bigalgebra_0.8.4.tar.gz`).

The package may also be build with user-supplied external BLAS and
LAPACK libraries, in either 32- or 64-bit varieties. This is an advanced
topic that requires additional Makevars modification, and may include
adjustment of the low-level calling syntax depending on the library
used.

Feel free to contact us for help installing and running the package.

## Author

Frédéric Bertrand, Michael J. Kane, Bryan Lewis, John W. Emerson

Maintainer: Frédéric Bertrand \<frederic.bertrand@lecnam.net\>

## References

<https://www.netlib.org/blas/>

<https://www.netlib.org/lapack/>

## See also

[`bigmemory`](https://rdrr.io/pkg/bigmemory/man/bigmemory-package.html),
[`big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)

## Examples

``` r

# Testing the development of the user-friendly operators:
# if you have any problems, please email us!  - Jay & Mike 4/29/2010

library("bigmemory")
A <- big.matrix(5,4, type="double", init=0,
                dimnames=list(NULL, c("alpha", "beta")))
B <- big.matrix(4,4, type="double", init=0,
                dimnames=list(NULL, c("alpha", "beta")))

C <- A 
D <- A[] 

print(C - D)       # Compare the results (subtraction of an R matrix from a
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0
#> [5,]    0    0    0    0
                   # big.matrix)

# The next example illustrates mixing R and big.matrix objects. It returns by
# default (see # options("bigalgebra.mixed_arithmetic_returns_R_matrix")
D <- matrix(rnorm(16),4)
E <- A 
```
