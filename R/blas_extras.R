# Internal helper to ensure base R vectors and matrices use double storage.
ensure_double <- function(x) {
  if (inherits(x, "big.matrix")) {
    return(x)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  x
}

#' Fill a vector or matrix with a constant value
#' 
#' @param N Optional integer specifying the number of elements to modify. Defaults to the length of `X`.
#' @param ALPHA Numeric scalar used to populate `X`.
#' @param X Double-precision vector, matrix or [`bigmemory::big.matrix`] to be filled in place.
#' @param INCX Integer stride between successive elements of `X`.
#'
#' @return Invisibly returns `X` after modification.
#' @export
#'
#' @examples
#' x <- matrix(0, 2, 3)
#' dset(ALPHA = 5, X = x)
#' x
#'
dset <- function(N = NULL, ALPHA, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dset_wrapper`, as.integer(N), as.numeric(ALPHA), X,
        as.integer(INCX), X.is.bm)
  invisible(X)
}

#' Generalised AXPBY operation
#'
#' @description
#' Computes the linear combination \eqn{Y := \alpha X + \beta Y} in place.
#'
#' @param N Optional integer giving the number of elements. Defaults to `length(X)`.
#' @param ALPHA Numeric scalar multiplying `X`.
#' @param X Double-precision vector, matrix or [`bigmemory::big.matrix`] input.
#' @param INCX Integer stride for traversing `X`.
#' @param BETA Numeric scalar multiplying `Y`.
#' @param Y Double-precision object updated in place.
#' @param INCY Integer stride for traversing `Y`.
#'
#' @return Invisibly returns the modified `Y`.
#' @export
#'
#' @examples
#' x <- 1:5
#' y <- rep(2, 5)
#' dvcal(ALPHA = 2, X = x, BETA = -1, Y = y)
#' y
#'
dvcal <- function(N = NULL, ALPHA = 1, X, INCX = 1L, BETA = 1, Y, INCY = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dvcal_wrapper`, as.integer(N), as.numeric(ALPHA), X,
        as.integer(INCX), as.numeric(BETA), Y, as.integer(INCY),
        X.is.bm, Y.is.bm)
  invisible(Y)
}

#' In-place vector subtraction
#'
#' @description
#' Updates `Y` by subtracting `X`, i.e. \eqn{Y := Y - X}.
#'
#' @inheritParams dvcal
#'
#' @return Invisibly returns the modified `Y`.
#' @export
#'
#' @examples
#' x <- 1:4
#' y <- rep(10, 4)
#' dsub(X = x, Y = y)
#' y
#'
dsub <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dsub_wrapper`, as.integer(N), X, as.integer(INCX), Y,
        as.integer(INCY), X.is.bm, Y.is.bm)
  invisible(Y)
}

#' Dot product of two vectors
#'
#' @inheritParams dsub
#'
#' @return Numeric scalar containing the dot product.
#' @export
#'
#' @examples
#' ddot(X = 1:3, Y = c(2, 4, 6))
#'
ddot <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_ddot_wrapper`, as.integer(N), X, as.integer(INCX), Y,
        as.integer(INCY), X.is.bm, Y.is.bm)
}

#' High-accuracy dot product
#'
#' @description
#' Forms the dot product using long-double accumulation to mitigate rounding error.
#'
#' @inheritParams ddot
#'
#' @return Numeric scalar equal to the dot product of `X` and `Y`.
#' @export
#'
#' @examples
#' dqddot(X = 1:3, Y = c(2, 4, 6))
#'
dqddot <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dqddot_wrapper`, as.integer(N), X, as.integer(INCX), Y,
        as.integer(INCY), X.is.bm, Y.is.bm)
}

#' Element-wise (Hadamard) product
#'
#' @description
#' Computes \eqn{Z := X \circ Y}. When `Z` is missing it is allocated automatically.
#'
#' @inheritParams dsub
#' @param Z Optional output container receiving the product.
#' @param INCZ Integer stride for `Z`.
#'
#' @return The updated object `Z`.
#' @export
#'
#' @examples
#' dhprod(X = 1:4, Y = rep(2, 4))
#'
dhprod <- function(N = NULL, X, INCX = 1L, Y, INCY = 1L,
                   Z, INCZ = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  if (missing(Z)) {
    dims <- if (inherits(X, "big.matrix")) c(nrow(X), ncol(X)) else dim(X)
    if (X.is.bm || Y.is.bm) {
      if (is.null(dims)) {
        dims <- c(N, 1L)
      }
      if (length(dims) == 1L) {
        dims <- c(dims, 1L)
      }
      Z <- anon_matrix(dims[1], dims[2], val = 0)
      Z.is.bm <- TRUE
    } else {
      if (is.null(dims)) {
        Z <- numeric(N)
      } else {
        Z <- matrix(0, dims[1], ifelse(length(dims) == 1L, 1L, dims[2]))
      }
      Z.is.bm <- FALSE
    }
  } else {
    Z <- ensure_double(Z)
    Z.is.bm <- check_matrix(Z, classes = classes)
  }
  .Call(`_dhprod_wrapper`, as.integer(N), X, as.integer(INCX), Y,
        as.integer(INCY), Z, as.integer(INCZ), X.is.bm, Y.is.bm, Z.is.bm)
  Z
}

#' Three-dimensional cross product
#'
#' @param X Numeric vector of length three, matrix with three rows, or big.matrix.
#' @param Y Numeric object matching the shape of `X`.
#' @param Z Optional output container.
#'
#' @return The updated object `Z` containing the cross product.
#' @export
#'
#' @examples
#' dxyz(X = c(1, 0, 0), Y = c(0, 1, 0))
#'
dxyz <- function(X, Y, Z) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  Y <- ensure_double(Y)
  X.is.bm <- check_matrix(X, classes = classes)
  Y.is.bm <- check_matrix(Y, classes = classes)
  if (missing(Z)) {
    if (X.is.bm || Y.is.bm) {
      Z <- anon_matrix(3, 1, val = 0)
      Z.is.bm <- TRUE
    } else {
      Z <- numeric(3)
      Z.is.bm <- FALSE
    }
  } else {
    Z <- ensure_double(Z)
    Z.is.bm <- check_matrix(Z, classes = classes)
  }
  .Call(`_dxyz_wrapper`, X, Y, Z, X.is.bm, Y.is.bm, Z.is.bm)
  Z
}

#' Sum of elements
#'
#' @inheritParams ddot
#'
#' @return Numeric scalar giving the sum of elements of `X`.
#' @export
#'
dsum <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dsum_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Sum of absolute values
#'
#' @inheritParams dsum
#'
#' @return Numeric scalar.
#' @export
#'
dasum <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dasum_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Euclidean norm (2-norm)
#'
#' @inheritParams dsum
#'
#' @return Numeric scalar containing the Euclidean norm.
#' @export
#'
dnrm2 <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dnrm2_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Product of vector elements
#'
#' @inheritParams dsum
#'
#' @return Numeric scalar equal to the product of elements of `X`.
#' @export
#'
dprdct <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_dprdct_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Index of the minimum element
#'
#' @inheritParams dsum
#'
#' @return Integer index (1-based) of the smallest entry in `X`.
#' @export
#'
idmin <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_idmin_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Index of the maximum element
#'
#' @inheritParams idmin
#'
#' @return Integer index (1-based).
#' @export
#'
idmax <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_idmax_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Index of the minimum absolute value
#'
#' @inheritParams idmin
#'
#' @return Integer index (1-based).
#' @export
#'
idamin <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_idamin_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Index of the maximum absolute value
#'
#' @inheritParams idmin
#'
#' @return Integer index (1-based).
#' @export
#'
idamax <- function(N = NULL, X, INCX = 1L) {
  classes <- c("big.matrix", "matrix", "vector", "numeric")
  X <- ensure_double(X)
  X.is.bm <- check_matrix(X, classes = classes)
  if (is.null(N)) {
    N <- length(X)
  }
  .Call(`_idamax_wrapper`, as.integer(N), X, as.integer(INCX), X.is.bm)
}

#' Symmetric matrix-matrix multiplication
#'
#' @description
#' Computes \eqn{C := \alpha \operatorname{op}(A) B + \beta C} when `A` is symmetric.
#'
#' @param SIDE Character specifying whether `A` multiplies from the left (`"L"`) or right (`"R"`).
#' @param UPLO Character indicating whether `A` stores the upper (`"U"`) or lower (`"L"`) triangle.
#' @param M,N Optional integers for the output dimensions.
#' @param ALPHA,BETA Numeric scalars.
#' @param A Symmetric matrix or big.matrix.
#' @param LDA,LDB,LDC Leading dimensions.
#' @param B Input matrix.
#' @param C Optional output container updated in place.
#'
#' @return Invisibly returns `C`.
#' @export
#'
#' @examples
#' A <- matrix(c(2, 1, 1, 3), 2, 2)
#' B <- diag(2)
#' C <- matrix(0, 2, 2)
#' dsymm(A = A, B = B, C = C)
#' C
#'
dsymm <- function(SIDE = "L", UPLO = "U", M = NULL, N = NULL,
                  ALPHA = 1, A, LDA = NULL, B, LDB = NULL,
                  BETA = 0, C, LDC = NULL) {
  A <- ensure_double(A)
  B <- ensure_double(B)
  A.is.bm <- check_matrix(A)
  B.is.bm <- check_matrix(B)

  side <- toupper(SIDE)
  if (!side %in% c("L", "R")) {
    stop("SIDE must be 'L' or 'R'.")
  }

  if (is.null(M)) {
    M <- nrow(B)
  }
  if (is.null(N)) {
    N <- ncol(B)
  }

  if (is.null(LDA)) {
    LDA <- nrow(A)
  }
  if (is.null(LDB)) {
    LDB <- nrow(B)
  }

  if (missing(C)) {
    C <- anon_matrix(M, N, val = 0)
  }
  C <- ensure_double(C)
  C.is.bm <- check_matrix(C)
  if (is.null(LDC)) {
    LDC <- nrow(C)
  }

  .Call(`_dsymm_wrapper`, as.character(side), as.character(UPLO),
        as.integer(M), as.integer(N), as.numeric(ALPHA), A, as.integer(LDA),
        B, as.integer(LDB), as.numeric(BETA), C, as.integer(LDC),
        A.is.bm, B.is.bm, C.is.bm)
  invisible(C)
}

