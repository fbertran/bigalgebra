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
#> Ncells 1012511 54.1    1995290 106.6         NA  1995290 106.6
#> Vcells 1817079 13.9    8388608  64.0      65536  2724101  20.8
```
