## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(bigalgebra)

## -----------------------------------------------------------------------------
A_sym <- matrix(c(2, 1, 1, 3), 2, 2)
B_rhs <- diag(2)
C_out <- matrix(0, 2, 2)
dsymm(A = A_sym, B = B_rhs, C = C_out, SIDE = "L", UPLO = "U")
C_out

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
A <- matrix(as.numeric(1:6), nrow = 2)
B <- matrix(seq(2, 12, by = 2), nrow = 3)
C <- matrix(0, nrow = 2, ncol = 4)
dgemm(TRANSA = "N", TRANSB = "N", A = A, B = B, C = C, BETA = 0)
C

## -----------------------------------------------------------------------------
C_t <- matrix(0, nrow = 3, ncol = 3)
dgemm(TRANSA = "T", TRANSB = "N", A = A, B = B, C = C_t, BETA = 0)
C_t

## -----------------------------------------------------------------------------
X <- matrix(1, nrow = 2, ncol = 2)
Y <- matrix(c(0, 1, 2, 3), nrow = 2)
daxpy(A = 0.5, X = X, Y = Y)

## -----------------------------------------------------------------------------
dir.create(tmp_axpy <- tempfile())
X_bm <- filebacked.big.matrix(2, 2, type = "double",
                              backingpath = tmp_axpy,
                              backingfile = "X.bin",
                              descriptorfile = "X.desc")
X_bm[,] <- 1
daxpy(A = 3, X = X_bm)
X_bm[]

## -----------------------------------------------------------------------------
unlink(file.path(tmp_sym, c("A.bin", "A.desc", "B.bin", "B.desc", "C.bin", "C.desc")))
unlink(tmp_sym, recursive = TRUE)
unlink(file.path(tmp_axpy, c("X.bin", "X.desc")))
unlink(tmp_axpy, recursive = TRUE)

