# Working with big.matrix Objects

## Overview

`bigalgebra` is designed to interoperate with the `bigmemory` ecosystem.
This vignette demonstrates how to create in-memory and file-backed
`big.matrix` objects, interact with them via the package’s wrappers, and
manage the underlying resources safely.

## Creating in-memory `big.matrix` objects

In-memory matrices behave much like ordinary R matrices but reside in
shared memory, allowing multiple R sessions to access the same data.

``` r
X <- big.matrix(3, 3, type = "double", init = 0)
X[,] <- matrix(1:9, nrow = 3)
X[]
#>      [,1] [,2] [,3]
#> [1,]    1    4    7
#> [2,]    2    5    8
#> [3,]    3    6    9
```

Once created, the objects can be passed directly to Level 1 helpers:

``` r
dvcal(ALPHA = 2, X = X, BETA = -1, Y = X)
X[]
#>      [,1] [,2] [,3]
#> [1,]    1    4    7
#> [2,]    2    5    8
#> [3,]    3    6    9
```

## Working with file-backed matrices

File-backed matrices persist their contents on disk, making them
suitable for data sets that exceed available RAM.

``` r
dir.create(tmp_fb <- tempfile())
Y <- filebacked.big.matrix(4, 2, type = "double",
                           backingpath = tmp_fb,
                           backingfile = "fb.bin",
                           descriptorfile = "fb.desc",
                           init = 0)
Y[,] <- matrix(runif(8), nrow = 4)
Y[]
#>            [,1]        [,2]
#> [1,] 0.08075014 0.007399441
#> [2,] 0.83433304 0.466393497
#> [3,] 0.60076089 0.497777389
#> [4,] 0.15720844 0.289767245
```

These objects participate in higher-level operations without being
loaded into memory.

``` r
Z <- filebacked.big.matrix(4, 2, type = "double",
                           backingpath = tmp_fb,
                           backingfile = "res.bin",
                           descriptorfile = "res.desc",
                           init = 0)
dvcal(ALPHA = 1.5, X = Y, BETA = 0, Y = Z)
Z[]
#>           [,1]       [,2]
#> [1,] 0.1211252 0.01109916
#> [2,] 1.2514996 0.69959025
#> [3,] 0.9011413 0.74666608
#> [4,] 0.2358127 0.43465087
```

## Sharing matrices between sessions

The descriptor file records the metadata needed to reopen a file-backed
matrix in a new R session. The
[`attach.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/attach.big.matrix.html)
helper reconstructs the object:

``` r
Y_desc <- dget(file.path(tmp_fb, "fb.desc"))
Y_again <- attach.big.matrix(Y_desc)
identical(Y[,], Y_again[,])
#> [1] TRUE
```

Any operations performed via `bigalgebra` update the shared backing
file, allowing all attached references to observe the change.

``` r
dsub(X = Z, Y = Y_again)
Y_again[]
#>             [,1]         [,2]
#> [1,] -0.04037507 -0.003699721
#> [2,] -0.41716652 -0.233196749
#> [3,] -0.30038044 -0.248888694
#> [4,] -0.07860422 -0.144883622
```

## Cleaning up backing files

File-backed matrices allocate resources on disk. Deleting the backing
and descriptor files once they are no longer needed helps keep the
workspace tidy.

``` r
unlink(file.path(tmp_fb, c("fb.bin", "fb.desc", "res.bin", "res.desc")))
unlink(tmp_fb, recursive = TRUE)
```
