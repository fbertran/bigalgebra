## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(bigalgebra)
library(bigmemory)

## -----------------------------------------------------------------------------
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
H <- hilbert(4)
H_big <- as.big.matrix(H)
TAU <- matrix(0, nrow = min(nrow(H), ncol(H)))
WORK <- matrix(0, nrow = max(1, ncol(H)))
dgeqrf(A = H_big, TAU = TAU, WORK = WORK)

# Extract the R factor from the overwritten big.matrix
R_big <- H_big[,]
R_big[lower.tri(R_big)] <- 0
R_big

# Compare against base R's QR decomposition
all.equal(R_big, qr.R(qr(H)))
TAU

## -----------------------------------------------------------------------------
tmp <- tempfile()
H_fb <- filebacked.big.matrix(nrow(H), ncol(H), init = H,
                              backingfile = basename(tmp),
                              backingpath = dirname(tmp))
dgeqrf(A = H_fb)
H_fb[,][upper.tri(H_fb[,], diag = TRUE)]
rm(H_fb); gc()

## -----------------------------------------------------------------------------
set.seed(42)
A <- matrix(rnorm(9), 3)
SPD <- crossprod(A)  # symmetric positive definite
SPD_big <- as.big.matrix(SPD)
info <- dpotrf(UPLO = "U", A = SPD_big)
info

U <- SPD_big[,]
U[lower.tri(U)] <- 0
all.equal(U, chol(SPD))

## -----------------------------------------------------------------------------
set.seed(123)
M <- matrix(rnorm(16), 4)
WR <- matrix(0, nrow = ncol(M))
WI <- matrix(0, nrow = ncol(M))
VL <- matrix(0, nrow = nrow(M), ncol = ncol(M))
dgeev(A = M, WR = WR, WI = WI, VL = VL)

# Compare eigenvalues with base R
complex_eigs <- WR[, 1] + 1i * WI[, 1]
all.equal(sort(complex_eigs), sort(eigen(M)$values))

## -----------------------------------------------------------------------------
set.seed(101)
X <- matrix(rnorm(12), 4)
S <- matrix(0, nrow = min(dim(X)))
U <- matrix(0, nrow = nrow(X), ncol = nrow(X))
VT <- matrix(0, nrow = ncol(X), ncol = ncol(X))
dgesdd(A = X, S = S, U = U, VT = VT)

svd_base <- svd(X)
all.equal(drop(S), svd_base$d)
all.equal(U, svd_base$u)
all.equal(VT, t(svd_base$v))

