# Changelog

## bigalgebra 3.1.1

- Fixed compilation with recent Clang and libc++ by including the
  standard headers required by `bigmemory/BigMatrix.h`.

## bigalgebra 3.1.0

CRAN release: 2026-03-24

- Dealt with upcoming Rf_error issue with Rcpp

## bigalgebra 3.0.0

CRAN release: 2025-10-05

- Added equivalent of the BLAS FORTRAN functions: DSET, DVCAL, DADD,
  DSUB, DSWAP, DDOT, DQDDOT, DHPROD, DXYZ, DSUM, DASUM, DNRM2, DPRDCT,
  IDMIN, IDMAX, IDAMIN, IDAMAX, DSYMM. More precisely, the wrappers that
  call into BLAS are: dadd(), ddot(), dasum(), dnrm2() and dsymm(). The
  other routines are implemented with explicit C++ loops.
- Added vignettes to the package
- Updated README
- Created a dsqrt() helper wired through the C++ backend to apply
  element-wise square roots to double vectors, matrices, and big.matrix
  objects in place.

## bigalgebra 2.0.2

CRAN release: 2025-09-08

- Maintainer email update
- Code fix to get rid of random errors for gcc 1.5 or windows based.

## bigalgebra 2.0.0

CRAN release: 2025-08-20

- Code cleaning and improvement.
- Bug fixes and updates to adapt to new CRAN checks.
- Added tests

## bigalgebra 1.1.2

CRAN release: 2024-09-23

- Update to adapt to evolution in CRAN checks.

## bigalgebra 1.1.1

CRAN release: 2023-12-09

- Update to adapt to evolution in CRAN checks.

## bigalgebra 1.1

- Update FORTRAN calls as requested by CRAN for R4.2.

## bigalgebra 1.0.2

- As recommended, updated the link to the JSS.

## bigalgebra 1.0.1

CRAN release: 2021-05-12

- Fixes for CRAN checks with Fedora.

## bigalgebra 1.0.0

CRAN release: 2021-04-18

- Github actions, package logo, pkgdown site and readme.

## bigalgebra 0.9.0

- Fixes for CRAN check by F. Bertrand.
