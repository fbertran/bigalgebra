# Copy a vector.

Copy double precision DX to double precision DY. For I = 0 to N-1, copy
DX(LX+I\*INCX) to DY(LY+I\*INCY), where LX = 1 if INCX .GE. 0, else LX =
1+(1-N)\*INCX, and LY is defined in a similar way using INCY.

## Usage

``` r
dcopy(N = NULL, X, INCX = 1, Y, INCY = 1)
```

## Arguments

- N:

  number of elements in input vector(s)

- X:

  double precision vector with N elements

- INCX:

  storage spacing between elements of DX

- Y:

  double precision vector with N elements

- INCY:

  storage spacing between elements of DY

## Value

DY copy of vector DX (unchanged if N .LE. 0)

## References

C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. Krogh, Basic linear
algebra subprograms for Fortran usage, Algorithm No. 539, Transactions
on Mathematical Software 5, 3 (September 1979), pp. 308-323.

## Examples

``` r
set.seed(4669)
A = big.matrix(3, 2, type="double", init=1, dimnames=list(NULL, 
c("alpha", "beta")), shared=FALSE)
B = big.matrix(3, 2, type="double", init=0, dimnames=list(NULL, 
c("alpha", "beta")), shared=FALSE)

dcopy(X=A,Y=B)
#> [1] 0
A[,]-B[,]
#>      alpha beta
#> [1,]     0    0
#> [2,]     0    0
#> [3,]     0    0

# The big.matrix file backings will be deleted when garbage collected.
rm(A,B)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1147654 61.3    2286010 122.1         NA  2286010 122.1
#> Vcells 2052770 15.7    8388608  64.0      65536  6482741  49.5
```
