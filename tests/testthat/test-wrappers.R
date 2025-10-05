# test-wrappers.R — unit tests for BLAS/LAPACK wrappers via .Call
# All comments in English.

library(testthat)
library(bigalgebra)

testthat::local_edition(3)

# Helper: compare matrices/vectors with numeric tolerance
expect_equal_num <- function(x, y, tol = 1e-10) {
  expect_equal(as.numeric(x), as.numeric(y), tolerance = tol)
}

# --- DCOPY ---------------------------------------------------------------
test_that("dcopy copies matrix contents", {
  set.seed(101)
  X <- matrix(runif(12), 3, 4)
  Y <- matrix(0, 3, 4)
  # R front-end allocates N if NULL
  expect_silent(dcopy(X = X, Y = Y))
  expect_equal(Y, X, tolerance = 1e-12)
})

# --- DADD ---------------------------------------------------------------
test_that("dadd adds matrix contents in place", {
  set.seed(1011)
  X <- matrix(runif(10), 5, 2)
  Y <- matrix(runif(10), 5, 2)
  Y_ref <- X + Y
  expect_silent(dadd(X = X, Y = Y))
  expect_equal(Y, Y_ref, tolerance = 1e-12)
})

# --- DSWAP ---------------------------------------------------------------
test_that("dswap exchanges matrix contents", {
  set.seed(1012)
  X <- matrix(runif(8), 4, 2)
  Y <- matrix(runif(8), 4, 2)
  X_ref <- matrix(X, nrow(X), ncol(X))
  Y_ref <- matrix(Y, nrow(Y), ncol(Y))
  expect_silent(dswap(X = X, Y = Y))
  expect_equal(X, Y_ref, tolerance = 1e-12)
  expect_equal(Y, X_ref, tolerance = 1e-12)
})

# --- DSCAL ---------------------------------------------------------------
test_that("dscal scales a matrix in place", {
  set.seed(102)
  Y <- matrix(runif(9), 3, 3)
  Y0 <- matrix(Y, nrow(Y), ncol(Y))  # deep copy to avoid shared memory
  alpha <- 2.5
  expect_silent(dscal(ALPHA = alpha, Y = Y))
  expect_equal(Y, alpha * Y0, tolerance = 1e-12)
})

# --- DSET ---------------------------------------------------------------
test_that("dset fills the target with a constant", {
  X <- matrix(0, 2, 3)
  expect_invisible(dset(ALPHA = 3.14, X = X))
  expect_equal(X, matrix(3.14, 2, 3))
})

# --- DVCAL --------------------------------------------------------------
test_that("dvcal computes alpha * X + beta * Y", {
  X <- matrix(as.numeric(1:6), 3, 2)
  Y <- matrix(as.numeric(6:1), 3, 2)
  expect_invisible(dvcal(ALPHA = 2, X = X, BETA = -1, Y = Y))
  expect_equal(Y, 2 * X - matrix(6:1, 3, 2))
})

# --- DSUB ---------------------------------------------------------------
test_that("dsub subtracts matrices in place", {
  X <- matrix(1, 4, 2)
  Y <- matrix(5, 4, 2)
  expect_invisible(dsub(X = X, Y = Y))
  expect_equal(Y, matrix(4, 4, 2))
})

# --- DDOT / DQDDOT ------------------------------------------------------
test_that("ddot matches base crossprod", {
  X <- as.numeric(1:5)
  Y <- seq(2, 10, by = 2)
  expect_equal(ddot(X = X, Y = Y), sum(X * Y))
})

test_that("dqddot matches ddot", {
  set.seed(42)
  X <- runif(10)
  Y <- runif(10)
  expect_equal(dqddot(X = X, Y = Y), ddot(X = X, Y = Y), tolerance = 1e-15)
})

# --- DHPROD -------------------------------------------------------------
test_that("dhprod returns element-wise products", {
  X <- matrix(as.numeric(1:6), 3, 2)
  Y <- matrix(rep(2, 6), 3, 2)
  Z <- dhprod(X = X, Y = Y)
  expect_equal(Z, 2 * X)
})

# --- DXYZ ---------------------------------------------------------------
test_that("dxyz computes cross product", {
  X <- c(1, 0, 0)
  Y <- c(0, 1, 0)
  Z <- dxyz(X = X, Y = Y)
  expect_equal(Z, c(0, 0, 1))
})

# --- DSUM / DASUM / DNRM2 / DPRDCT --------------------------------------
test_that("scalar summary helpers match base R", {
  X <- c(-1, 2, -3, 4)
  expect_equal(dsum(X = X), sum(X))
  expect_equal(dasum(X = X), sum(abs(X)))
  expect_equal(dnrm2(X = X), sqrt(sum(X^2)))
  expect_equal(dprdct(X = X), prod(X))
})

# --- IDMIN / IDMAX / IDAMIN / IDAMAX -----------------------------------
test_that("index helpers return 1-based extrema", {
  X <- c(-2, 5, -7, 3)
  expect_equal(idmin(X = X), which.min(X))
  expect_equal(idmax(X = X), which.max(X))
  expect_equal(idamin(X = X), which.min(abs(X)))
  expect_equal(idamax(X = X), which.max(abs(X)))
})

# --- DSYMM --------------------------------------------------------------
test_that("dsymm matches base tcrossprod when SIDE = 'L'", {
  skip_on_cran()
  set.seed(321)
  A <- crossprod(matrix(rnorm(9), 3, 3))
  B <- matrix(rnorm(9), 3, 3)
  C <- matrix(0, 3, 3)
  expect_silent(dsymm(A = A, B = B, C = C, SIDE = "L", UPLO = "U"))
  expect_equal(C, A %*% B, tolerance = 1e-12)
})

# --- DGEMM ---------------------------------------------------------------
test_that("dgemm matches base %*% for N,N", {
  skip_on_cran()
  set.seed(103)
  A <- matrix(rnorm(6), 2, 3)
  B <- matrix(rnorm(9), 3, 3)
  C <- matrix(0,nrow(A),ncol(B))
  C_pkg <- dgemm(A = A, B = B, C = C, TRANSA = "N", TRANSB = "N", ALPHA = 1, BETA = 0)
  C_ref <- A %*% B
  expect_equal(C_pkg, C_ref, tolerance = 1e-12)
})

# --- DPOTRF --------------------------------------------------------------
test_that("dpotrf produces a valid Cholesky factor (upper)", {
  skip_on_cran()
  set.seed(104)
  A <- crossprod(matrix(rnorm(16), 4, 4)) + diag(1e-8, 4)
  A0 <- matrix(A, nrow(A), ncol(A))  # deep copy to avoid shared memory
  info <- dpotrf(UPLO = "U", N = nrow(A), A = A, LDA = nrow(A))
  expect_identical(info, 0)
  R <- matrix(A, nrow(A), ncol(A)); R[lower.tri(R)] <- 0
  expect_equal(t(R) %*% R, A0, tolerance = 1e-8)
})

# --- DGEQRF --------------------------------------------------------------
test_that("dgeqrf writes R in upper triangle compatible with base qr", {
  skip_on_cran()
  set.seed(105)
  A <- matrix(rnorm(15), 5, 3)
  A0 <- matrix(A, nrow(A), ncol(A))  # deep copy to avoid shared memory
  info <- dgeqrf(M = nrow(A), N = ncol(A), A = A, LDA = nrow(A))
  expect_identical(info, 0)
  R_mod <- matrix(A[1:ncol(A), 1:ncol(A)], ncol(A), ncol(A)); R_mod[lower.tri(R_mod)] <- 0
  R_ref <- qr.R(qr(A0))
  expect_equal(R_mod, R_ref, tolerance = 1e-8)
})

# --- DGEEV ---------------------------------------------------------------
test_that("dgeev returns correct eigenvalues (values only)", {
  skip_on_cran()
  A <- matrix(c(2,1,0,3), 2, 2)  # eigenvalues 2 and 3
  WR <- matrix(0,nrow=2,ncol=1); WI <- matrix(0,nrow=2,ncol=1)
  VL <- matrix(0,2,2); VR <- matrix(0,2,2)
  LWORK <- 8; WORK <- matrix(0,nrow=max(1,LWORK),ncol=max(1,LWORK)); INFO <- 0
  info <- dgeev(JOBVL = "N", JOBVR = "N", N = 2, A = A, LDA = 2,
                WR = WR, WI = WI, VL = VL, LDVL = 2, VR = VR, LDVR = 2,
                WORK = WORK, LWORK = LWORK)
  expect_identical(info, 0)
  vals <- sort(WR + 1i * WI)
  expect_equal(vals, as.complex(sort(eigen(A)$values)), tolerance = 1e-10)
})

# --- DGESDD --------------------------------------------------------------
test_that("dgesdd returns correct singular values (JOBZ = 'N')", {
  skip_on_cran()
  set.seed(106)
  A <- matrix(rnorm(12), 4, 3)
  A0 <- matrix(A, nrow(A), ncol(A))  # deep copy: dgesdd overwrites A
  S <- matrix(0,nrow=min(dim(A)),ncol=1)
  U <- matrix(0, nrow(A), nrow(A))
  VT <- matrix(0, ncol(A), ncol(A))
  # Heuristic workspace; the R front-end can allocate if needed
#  WORK <- matrix(0,nrow=max(1, 5 * min(200) + 7 * max(200)),ncol=1);
#  LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "N", M = nrow(A), N = ncol(A), A = A, LDA = nrow(A),
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  #,
   #              WORK = WORK, LWORK = LWORK)
  expect_identical(info, 0)
  S_ref <- svd(A0, nu = 0, nv = 0)$d  # compare to original A
  expect_equal(sort(S), sort(S_ref), tolerance = 1e-10)
})



















# --- BIGMEMORY INTEGRATION TESTS ----------------------------------------

test_that("dcopy works with big.matrix", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(201)
  Xr <- matrix(runif(12), 3, 4)
  Xbm <- big.matrix(nrow = 3, ncol = 4, type = "double", init = NA_real_)
  Ybm <- big.matrix(nrow = 3, ncol = 4, type = "double", init = 0)
  Xbm[,] <- Xr
  Ybm[,] <- 0
  expect_silent(dcopy(X = Xbm, Y = Ybm))
  expect_equal(as.matrix(Ybm[]), Xr, tolerance = 1e-12)
  rm(Xbm,Ybm)
  gc()
})

test_that("dadd works with big.matrix", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(2011)
  Xr <- matrix(runif(10), 5, 2)
  Yr <- matrix(runif(10), 5, 2)
  Xbm <- big.matrix(nrow = 5, ncol = 2, type = "double")
  Ybm <- big.matrix(nrow = 5, ncol = 2, type = "double")
  Xbm[,] <- Xr
  Ybm[,] <- Yr
  expect_silent(dadd(X = Xbm, Y = Ybm))
  expect_equal(as.matrix(Ybm[]), Xr + Yr, tolerance = 1e-12)
  rm(Xbm, Ybm)
  gc()
})

test_that("dswap works with big.matrix", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(2012)
  Xr <- matrix(runif(8), 4, 2)
  Yr <- matrix(runif(8), 4, 2)
  Xbm <- big.matrix(nrow = 4, ncol = 2, type = "double")
  Ybm <- big.matrix(nrow = 4, ncol = 2, type = "double")
  Xbm[,] <- Xr
  Ybm[,] <- Yr
  X_ref <- matrix(Xr, nrow(Xr), ncol(Xr))
  Y_ref <- matrix(Yr, nrow(Yr), ncol(Yr))
  expect_silent(dswap(X = Xbm, Y = Ybm))
  expect_equal(as.matrix(Xbm[]), Y_ref, tolerance = 1e-12)
  expect_equal(as.matrix(Ybm[]), X_ref, tolerance = 1e-12)
  rm(Xbm, Ybm)
  gc()
})

test_that("dscal scales big.matrix in place", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(202)
  Yr <- matrix(runif(9), 3, 3)
  Ybm <- big.matrix(nrow = 3, ncol = 3, type = "double", init = NA_real_)
  Ybm[,] <- Yr
  alpha <- 1.7
  expect_silent(dscal(ALPHA = alpha, Y = Ybm))
  expect_equal(as.matrix(Ybm[]), alpha * Yr, tolerance = 1e-12)
  rm(Ybm)
  gc()
})

test_that("dgemm with big.matrix matches base %*%", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(203)
  A <- matrix(rnorm(6), 2, 3)
  B <- matrix(rnorm(9), 3, 3)
  C_ref <- A %*% B
  Abm <- big.matrix(nrow = 2, ncol = 3, type = "double")
  Bbm <- big.matrix(nrow = 3, ncol = 3, type = "double")
  Cbm <- big.matrix(nrow = 2, ncol = 3, type = "double", init = 0)
  Abm[,] <- A; Bbm[,] <- B; Cbm[,] <- 0
  expect_silent(dgemm(A = Abm, B = Bbm, C = Cbm, TRANSA = "N", TRANSB = "N", ALPHA = 1, BETA = 0))
  expect_equal(as.matrix(Cbm[]), C_ref, tolerance = 1e-12)
  rm(Abm,Bbm,Cbm)
  gc()
})


test_that("dpotrf on big.matrix produces a valid Cholesky factor (upper)", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(301)
  Ar <- crossprod(matrix(rnorm(25), 5, 5)) + diag(1e-8, 5)
  Abm <- big.matrix(nrow = 5, ncol = 5, type = "double")
  Abm[,] <- Ar
  info <- dpotrf(UPLO = "U", N = nrow(Ar), A = Abm, LDA = nrow(Ar))
  expect_identical(info, 0)
  R <- as.matrix(Abm[]); R[lower.tri(R)] <- 0
  expect_equal(t(R) %*% R, Ar, tolerance = 1e-8)
  rm(Abm)
  gc()
})

test_that("dgeqrf on big.matrix matches base qr R factor", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(302)
  Ar <- matrix(rnorm(20), 5, 4)
  Abm <- big.matrix(nrow = nrow(Ar), ncol = ncol(Ar), type = "double")
  Abm[,] <- Ar
  info <- dgeqrf(M = nrow(Ar), N = ncol(Ar), A = Abm, LDA = nrow(Ar))
  expect_identical(info, 0)
  R_mod <- as.matrix(Abm[]);
  R_mod <- matrix(R_mod[1:ncol(R_mod), 1:ncol(R_mod)], ncol(R_mod), ncol(R_mod))
  R_mod[lower.tri(R_mod)] <- 0
  R_ref <- qr.R(qr(Ar))
  expect_equal(R_mod, R_ref, tolerance = 1e-8)
  rm(Abm)
  gc()
})

test_that("dgeev on big.matrix returns correct eigenvalues (values only)", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  Ar <- matrix(c(2,1,0,3), 2, 2)
  Abm <- big.matrix(nrow = 2, ncol = 2, type = "double")
  Abm[,] <- Ar
  WR <- matrix(0,nrow=2,ncol=1); WI <- matrix(0,nrow=2,ncol=1)
  VL <- matrix(0,2,2); VR <- matrix(0,2,2)
  LWORK <- 8; WORK <- matrix(0,nrow=max(1,LWORK),ncol=max(1,LWORK)); INFO <- 0
  info <- dgeev(JOBVL = "N", JOBVR = "N", N = 2, A = Abm, LDA = 2,
                WR = WR, WI = WI, VL = VL, LDVL = 2, VR = VR, LDVR = 2,
                WORK = WORK, LWORK = LWORK)
  expect_identical(info, 0)
  vals <- sort(WR + 1i * WI)
  expect_equal(vals, as.complex(sort(eigen(Ar)$values)), tolerance = 1e-10)
  rm(Abm)
  gc()
})

test_that("dgesdd on big.matrix returns correct singular values (JOBZ = 'N')", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(303)
  Ar <- matrix(rnorm(20), 5, 4)
  A0 <- matrix(Ar, nrow(Ar), ncol(Ar)) # deep copy: dgesdd overwrites A
  Abm <- big.matrix(nrow = nrow(Ar), ncol = ncol(Ar), type = "double")
  Abm[,] <- Ar
  S <- matrix(0,nrow=min(dim(Ar)),ncol=1)
  U <- matrix(0, nrow(Ar), nrow(Ar))
  VT <- matrix(0, ncol(Ar), ncol(Ar))
  # Heuristic workspace; the R front-end can allocate if needed
  #WORK <- numeric(max(1, 5 * min(dim(Ar)) + 7 * max(dim(Ar))))
  #LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "N", M = nrow(Ar), N = ncol(Ar), A = Abm, LDA = nrow(Ar),
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  expect_identical(info, 0)
  expect_equal(sort(S), sort(svd(A0, nu = 0, nv = 0)$d), tolerance = 1e-10)
  rm(Abm)
  gc()
})




# --- VARIANTS ------------------------------------------------------------

# DGEMM transpose variants (non-big)
test_that("dgemm matches base %*% for transpose variants", {
  skip_on_cran()
  set.seed(401)
  A <- matrix(rnorm(8), 4, 2)  # 4x2
  B <- matrix(rnorm(12), 2, 6) # 2x6
  # (N,N)
  C <- matrix(0, 4, 6)
  dgemm(A = A, B = B, C = C, TRANSA = "N", TRANSB = "N", ALPHA = 1, BETA = 0)
  expect_equal(C, A %*% B, tolerance = 1e-12)
  # (T,N)
  At <- t(A)
  C <- matrix(0, nrow(t(At)), ncol(B))
  dgemm(A = At, B = B, C = C, TRANSA = "T", TRANSB = "N", ALPHA = 1, BETA = 0)
  expect_equal(C, t(At) %*% B, tolerance = 1e-12)
  # (N,T)
  Bt <- t(B)
  C <- matrix(0, nrow(A), ncol(t(Bt)))
  dgemm(A = A, B = Bt, C = C, TRANSA = "N", TRANSB = "T", ALPHA = 1, BETA = 0)
  expect_equal(C, A %*% t(Bt), tolerance = 1e-12)
  # (T,T)
  C <- matrix(0, nrow(t(At)), ncol(t(Bt)))
  dgemm(A = At, B = Bt, C = C, TRANSA = "T", TRANSB = "T", ALPHA = 1, BETA = 0)
  expect_equal(C, t(At) %*% t(Bt), tolerance = 1e-12)
})

# DPOTRF lower-triangle variant (non-big)
test_that("dpotrf(L) produces valid lower Cholesky", {
  skip_on_cran()
  set.seed(402)
  A <- crossprod(matrix(rnorm(25), 5, 5)) + diag(1e-8, 5)
  A0 <- matrix(A,nrow(A),ncol(A))
  info <- dpotrf(UPLO = "L", N = nrow(A), A = A, LDA = nrow(A))
  expect_identical(info, 0)
  A[upper.tri(A)] <- 0
  expect_equal(A %*% t(A), A0, tolerance = 1e-8)
})

# DGEEV with eigenvectors (values + right eigenvectors) — symmetric case so real
test_that("dgeev returns correct eigenvalues and right eigenvectors (symmetric)", {
  skip_on_cran()
  set.seed(403)
  A <- crossprod(matrix(rnorm(16), 4, 4))  # symmetric positive semidefinite
  A0 <- matrix(A,nrow(A),ncol(A))
  WR <- matrix(0,nrow=nrow(A),ncol=1); WI <- matrix(0,nrow=nrow(A),ncol=1)
  VL <- matrix(0,nrow(A),nrow(A)); VR <- matrix(0,nrow(A),nrow(A))
  LWORK <- 20; WORK <- matrix(0,nrow=max(1,LWORK),ncol=max(1,LWORK)); INFO <- 0
  info <- dgeev(JOBVL = "N", JOBVR = "V", N = nrow(A), A = A, LDA = nrow(A),
                WR = WR, WI = WI, VL = VL, LDVL = nrow(VL), VR = VR, LDVR = nrow(VR),
                WORK = WORK, LWORK = LWORK)
  expect_identical(info, 0)
  eig <- eigen(A, symmetric = TRUE)
  # Eigenvalues match (up to ordering)
  expect_equal(sort(WR), sort(eig$values), tolerance = 1e-8)
  # Eigenvectors satisfy A %*% v ~= lambda * v (check first few)
  ord <- order(WR, decreasing = TRUE)
  for (j in ord) {
    v <- VR[, j]
    lam <- WR[j]
    expect_equal(as.vector(A0 %*% v), as.vector(lam * v), tolerance = 1e-6)
  }
})

# DGESDD JOBZ='S': thin U/VT with reconstruction (non-big)
test_that("dgesdd(JOBZ='S') returns thin factors that reconstruct A", {
  skip_on_cran()
  set.seed(404)
  A <- matrix(rnorm(20), 5, 4)  # m > n
  A0 <- matrix(A,nrow(A),ncol(A))
  m <- nrow(A); n <- ncol(A); r <- min(m, n)
  S <- matrix(0, r, 1)
  U <- matrix(0, m, r)
  VT <- matrix(0, r, n)
  #WORK <- numeric(max(1, 5 * r + 7 * max(m, n)))
  #LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "S", M = m, N = n, A = A, LDA = m,
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  expect_identical(info, 0)
  Arec <- U %*% diag(as.vector(S), nrow = r, ncol = r) %*% VT
  expect_equal(Arec, A0, tolerance = 1e-8)
})


# DGESDD JOBZ='A': full U/VT with reconstruction (non-big)
test_that("dgesdd(JOBZ='A') returns full factors that reconstruct A", {
  skip_on_cran()
  set.seed(405)
  A <- matrix(rnorm(20), 5, 4)
  A0 <- matrix(A,nrow(A),ncol(A))
  m <- nrow(A); n <- ncol(A); r <- min(m, n)
  S <- matrix(0, r, 1)
  U <- matrix(0, m, m)
  VT <- matrix(0, n, n)
  #WORK <- numeric(max(1, 5 * r + 7 * max(m, n)))
  #LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "A", M = m, N = n, A = A, LDA = m,
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  expect_identical(info, 0)
  # Reconstruct using the first r singular vectors
  Arec <- U[, seq_len(r), drop = FALSE] %*% diag(as.vector(S), nrow = r, ncol = r) %*%
    VT[seq_len(r), , drop = FALSE]
  expect_equal(Arec, A0, tolerance = 1e-8)
})

# BIGMATRIX variants ------------------------------------------------------

test_that("dgeev with big.matrix returns real eigenpairs for symmetric A", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(406)
  A <- crossprod(matrix(rnorm(16), 4, 4))
  Abm <- big.matrix(nrow = 4, ncol = 4, type = "double")
  Abm[,] <- A
  WR <- matrix(0, 4, 1); WI <- matrix(0, 4, 1)
  VL <- matrix(0, 4, 4); VR <- matrix(0, 4, 4)
  WORK <- matrix(0,64,1); LWORK <- length(WORK); INFO <- 0
  info <- dgeev(JOBVL = "N", JOBVR = "V", N = 4, A = Abm, LDA = 4,
                WR = WR, WI = WI, VL = VL, LDVL = 4, VR = VR, LDVR = 4,
                WORK = WORK, LWORK = LWORK)
  expect_identical(info, 0)
  eig <- eigen(A, symmetric = TRUE)
  expect_equal(sort(WR), sort(eig$values), tolerance = 1e-8)
  # Check eigen equation A v = lambda v for two columns
  for (j in 1:2) {
    v <- VR[, j]
    lam <- WR[j]
    expect_equal(as.vector(A %*% v), as.vector(lam * v), tolerance = 1e-6)
  }
  rm(Abm)
  gc()
}
)

test_that("dgesdd(JOBZ='S') on big.matrix reconstructs A", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(407)
  A <- matrix(rnorm(20), 5, 4)
  A0 <- matrix(A,nrow(A),ncol(A))
  Abm <- big.matrix(nrow = 5, ncol = 4, type = "double"); Abm[,] <- A
  m <- 5; n <- 4; r <- 4
  S <- matrix(0,r,1)
  U <- matrix(0, m, r)
  VT <- matrix(0, r, n)
  #WORK <- matrix(0,max(1, 5 * r + 7 * max(m, n)),1)
  #LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "S", M = m, N = n, A = Abm, LDA = m,
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  expect_identical(info, 0)
  Arec <- U %*% diag(as.vector(S), nrow = r, ncol = r) %*% VT
  expect_equal(Arec, A0, tolerance = 1e-8)
  rm(Abm)
  gc()
})

test_that("dgesdd(JOBZ='A') on big.matrix reconstructs A", {
  skip_on_cran()
  skip_if_not_installed("bigmemory")
  library(bigmemory)
  set.seed(408)
  A <- matrix(rnorm(20), 5, 4)
  A0 <- matrix(A,nrow(A),ncol(A))
  Abm <- big.matrix(nrow = 5, ncol = 4, type = "double"); Abm[,] <- A
  m <- 5; n <- 4; r <- 4
  S <- matrix(0, r, 1)
  U <- matrix(0, m, m)
  VT <- matrix(0, n, n)
  #WORK <- numeric(max(1, 5 * r + 7 * max(m, n)))
  #LWORK <- length(WORK); INFO <- 0
  info <- dgesdd(JOBZ = "A", M = m, N = n, A = Abm, LDA = m,
                 S = S, U = U, LDU = nrow(U), VT = VT, LDVT = nrow(VT))
  expect_identical(info, 0)
  Arec <- U[, seq_len(r), drop = FALSE] %*% diag(as.vector(S), nrow = r, ncol = r) %*%
    VT[seq_len(r), , drop = FALSE]
  expect_equal(Arec, A0, tolerance = 1e-8)
  rm(Abm)
  gc()
})

# --- DGEEV workspace query & auto-allocation tests ----------------------

test_that("dgeev workspace query returns >= 4N when vectors requested", {
  skip_on_cran()
  n <- 5L
  # requires the C symbol `_dgeev_lwork_query_wrapper` to be registered
  optV <- .Call(`_dgeev_lwork_query_wrapper`, "V", "N", as.integer(n))
  expect_type(optV, "integer")
  expect_gte(as.integer(optV), 4L * n)
  optN <- .Call(`_dgeev_lwork_query_wrapper`, "N", "N", as.integer(n))
  expect_type(optN, "integer")
  expect_gte(as.integer(optN), 3L * n)
})

test_that("dgeev auto-allocates WORK when VL requested (non-big)", {
  skip_on_cran()
  set.seed(4669)
  A <- matrix(rnorm(16), 4, 4)
  WR <- matrix(0, 4, 1); WI <- matrix(0, 4, 1)
  VL <- matrix(0, 4, 4)
  info <- dgeev(A = A, WR = WR, WI = WI, VL = VL) # no WORK/LWORK provided
  expect_identical(info, 0)
  expect_equal(sort(WR + 1i*WI), sort(eigen(A)$values), tolerance = 1e-8)
})

test_that("dgeev auto-allocates WORK when VR requested (non-big)", {
  skip_on_cran()
  set.seed(4670)
  A <- matrix(rnorm(25), 5, 5)
  WR <- matrix(0, 5, 1); WI <- matrix(0, 5, 1)
  VR <- matrix(0, 5, 5)
  info <- dgeev(A = A, WR = WR, WI = WI, VR = VR) # vectors on right
  expect_identical(info, 0)
  expect_equal(sort(WR + 1i*WI), sort(eigen(A)$values), tolerance = 1e-8)
})

test_that("dgeev auto-allocates WORK when neither VL or VR requested (non-big)", {
  skip_on_cran()
  set.seed(4669)
  A <- matrix(rnorm(16), 4, 4)
  WR <- matrix(0, 4, 1); WI <- matrix(0, 4, 1)
  VL <- matrix(0, 4, 4)
  info <- dgeev(A = A, WR = WR, WI = WI) # no WORK/LWORK provided
  expect_identical(info, 0)
  expect_equal(sort(WR + 1i*WI), sort(eigen(A)$values), tolerance = 1e-8)
})