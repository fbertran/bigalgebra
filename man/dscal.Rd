% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigalgebra.R
\name{dscal}
\alias{dscal}
\title{Scales a vector by a constant.}
\usage{
dscal(N = NULL, ALPHA, Y, INCY = 1)
}
\arguments{
\item{N}{an integer. Number of elements in input vector(s)}

\item{ALPHA}{a real number. The scalar alpha}

\item{Y}{a big matrix to scale by ALPHA}

\item{INCY}{an integer. Storage spacing between elements of Y.}
}
\value{
Update Y.
}
\description{
Scales a vector by a constant.
}
\examples{
\dontrun{
set.seed(4669)
A = big.matrix(3, 2, type="double", init=1, dimnames=list(NULL, 
c("alpha", "beta")), shared=FALSE)
dscal(ALPHA=2,Y=A)
A[,]

# The big.matrix file backings will be deleted when garbage collected.
rm(A)
gc()
}
}
