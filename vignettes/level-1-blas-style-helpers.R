## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(bigalgebra)

## -----------------------------------------------------------------------------
# Fill an existing matrix with a scalar
X <- matrix(0, 2, 3)
dset(ALPHA = 3.5, X = X)
X

# Form alpha * X + beta * Y in place
dx <- as.numeric(1:4)
dy <- rep(2, 4)
dvcal(ALPHA = 2, X = dx, BETA = -1, Y = dy)
dy

# Subtract X from Y element-wise
dy2 <- rep(10, 4)
dsub(X = dx, Y = dy2)
dy2

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
# Classic dot product and its extended-precision counterpart
v1 <- as.numeric(1:5)
v2 <- seq(2, 10, by = 2)
list(ddot = ddot(X = v1, Y = v2), dqddot = dqddot(X = v1, Y = v2))

# Hadamard product stored in a new matrix
A <- matrix(as.numeric(1:4), 2, 2)
B <- matrix(rep(2, 4), 2, 2)
Z <- dhprod(X = A, Y = B)
Z

# Element-wise square root performed in place
sqrt_vals <- matrix(c(1, 4, 9, 16), 2)
dsqrt(X = sqrt_vals)
sqrt_vals

## -----------------------------------------------------------------------------
ux <- c(1, 0, 0)
uy <- c(0, 1, 0)
dxyz(X = ux, Y = uy)

## -----------------------------------------------------------------------------
vals <- c(-1, 2, -3, 4)
list(
  sum = dsum(X = vals),
  abs_sum = dasum(X = vals),
  euclidean_norm = dnrm2(X = vals),
  product = dprdct(X = vals)
)

## -----------------------------------------------------------------------------
idx_vals <- c(-2, 5, -7, 3)
list(
  min_index = idmin(X = idx_vals),
  max_index = idmax(X = idx_vals),
  min_abs_index = idamin(X = idx_vals),
  max_abs_index = idamax(X = idx_vals)
)

## -----------------------------------------------------------------------------
unlink(file.path(tmp, "vec.bin"))
unlink(file.path(tmp, "vec.desc"))
unlink(tmp, recursive = TRUE)

