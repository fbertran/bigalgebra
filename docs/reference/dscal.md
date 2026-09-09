# Scales a vector by a constant.

Scales a vector by a constant.

## Usage

``` r
dscal(N = NULL, ALPHA, Y, INCY = 1L)
```

## Arguments

- N:

  an integer. Number of elements in input vector(s)

- ALPHA:

  a real number. The scalar alpha

- Y:

  a big matrix to scale by ALPHA

- INCY:

  an integer. Storage spacing between elements of Y.

## Value

Update Y.

## Examples

``` r
set.seed(4669)
A = big.matrix(3, 2, type="double", init=1, dimnames=list(NULL, 
c("alpha", "beta")), shared=FALSE)
dscal(ALPHA=2,Y=A)
#> An object of class "big.matrix"
#> Slot "address":
#> <pointer: 0x1078049d0>
#> 
A[,]
#>      alpha beta
#> [1,]     2    2
#> [2,]     2    2
#> [3,]     2    2

# The big.matrix file backings will be deleted when garbage collected.
rm(A)
gc()
#>           used (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
#> Ncells 1017315 54.4    1995290 106.6         NA  1995290 106.6
#> Vcells 1828205 14.0    8388608  64.0      65536  2724101  20.8
```
