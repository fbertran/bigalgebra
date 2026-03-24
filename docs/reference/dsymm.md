# Symmetric matrix-matrix multiplication

Computes \\C := \alpha \operatorname{op}(A) B + \beta C\\ when \`A\` is
symmetric.

## Usage

``` r
dsymm(
  SIDE = "L",
  UPLO = "U",
  M = NULL,
  N = NULL,
  ALPHA = 1,
  A,
  LDA = NULL,
  B,
  LDB = NULL,
  BETA = 0,
  C,
  LDC = NULL
)
```

## Arguments

- SIDE:

  Character specifying whether \`A\` multiplies from the left (\`"L"\`)
  or right (\`"R"\`).

- UPLO:

  Character indicating whether \`A\` stores the upper (\`"U"\`) or lower
  (\`"L"\`) triangle.

- M, N:

  Optional integers for the output dimensions.

- ALPHA, BETA:

  Numeric scalars.

- A:

  Symmetric matrix or big.matrix.

- LDA, LDB, LDC:

  Leading dimensions.

- B:

  Input matrix.

- C:

  Optional output container updated in place.

## Value

Invisibly returns \`C\`.

## Examples

``` r
A <- matrix(c(2, 1, 1, 3), 2, 2)
B <- diag(2)
C <- matrix(0, 2, 2)
dsymm(A = A, B = B, C = C)
C
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    1    3
```
