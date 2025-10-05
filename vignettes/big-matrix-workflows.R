## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(bigalgebra)
library(bigmemory)

## -----------------------------------------------------------------------------
X <- big.matrix(3, 3, type = "double", init = 0)
X[,] <- matrix(1:9, nrow = 3)
X[]

## -----------------------------------------------------------------------------
dvcal(ALPHA = 2, X = X, BETA = -1, Y = X)
X[]

## -----------------------------------------------------------------------------
dir.create(tmp_fb <- tempfile())
Y <- filebacked.big.matrix(4, 2, type = "double",
                           backingpath = tmp_fb,
                           backingfile = "fb.bin",
                           descriptorfile = "fb.desc",
                           init = 0)
Y[,] <- matrix(runif(8), nrow = 4)
Y[]

## -----------------------------------------------------------------------------
Z <- filebacked.big.matrix(4, 2, type = "double",
                           backingpath = tmp_fb,
                           backingfile = "res.bin",
                           descriptorfile = "res.desc",
                           init = 0)
dvcal(ALPHA = 1.5, X = Y, BETA = 0, Y = Z)
Z[]

## -----------------------------------------------------------------------------
Y_desc <- dget(file.path(tmp_fb, "fb.desc"))
Y_again <- attach.big.matrix(Y_desc)
identical(Y[,], Y_again[,])

## -----------------------------------------------------------------------------
dsub(X = Z, Y = Y_again)
Y_again[]

## -----------------------------------------------------------------------------
unlink(file.path(tmp_fb, c("fb.bin", "fb.desc", "res.bin", "res.desc")))
unlink(tmp_fb, recursive = TRUE)

