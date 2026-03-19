#' Kernel energy balancing weights via closed-form solution
#'
#' Computes balancing weights that minimize a kernelized energy distance between
#' the weighted treated and control distributions and the overall sample. The
#' weights are obtained via a closed-form solution to a linear system derived
#' from the kernel energy distance objective.
#'
#' @param trt A binary (0/1) integer or numeric vector indicating treatment
#'   assignment (\code{1} = treated, \code{0} = control).
#' @param kern A symmetric \eqn{n \times n} kernel matrix, where \eqn{n} is the
#'   length of \code{trt}. Typically a random forest proximity kernel from
#'   \code{\link{forest_kernel}}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{weights}{A numeric vector of length \eqn{n} containing the balancing
#'     weights. Treated weights sum to \eqn{n_1} and control weights sum to
#'     \eqn{n_0}.}
#' }
#'
#' @details
#' The kernel energy distance between two distributions \eqn{P} and \eqn{Q}
#' measured via kernel \eqn{K} has a quadratic form. By setting up the
#' optimization to balance the treated and control groups towards the overall
#' sample, subject to group-wise sum constraints, the solution reduces to a
#' linear system involving the kernel matrix.
#'
#' Specifically, the method solves:
#' \deqn{w = K^{-1} z,}
#' where \eqn{z} is adjusted to satisfy the constraints
#' \eqn{\sum_{i: A_i=1} w_i = n_1} and \eqn{\sum_{i: A_i=0} w_i = n_0}.
#'
#' @references
#' De, B. and Huling, J. (2024). Kernel Energy Balancing with Random Forest
#' Similarity. \emph{arXiv preprint arXiv:2512.18069}.
#'
#' @examples
#' \donttest{
#' library(grf)
#' n <- 200
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' A <- rbinom(n, 1, plogis(X[, 1]))
#' Y <- X[, 1] + rnorm(n)
#'
#' forest <- multi_regression_forest(X, cbind(A, Y), num.trees = 500)
#' K <- forest_kernel(forest)
#' bal <- kernel_balance(A, K)
#'
#' # Weighted ATE estimate
#' w <- bal$weights
#' ate <- weighted.mean(Y[A == 1], w[A == 1]) -
#'        weighted.mean(Y[A == 0], w[A == 0])
#' }
#'
#' @importFrom Matrix Matrix
#' @export
kernel_balance <- function(trt, kern) {
  n  <- length(trt)
  n1 <- sum(trt == 1)
  n0 <- n - n1

  if (n1 == 0 || n0 == 0) {
    stop("Treatment vector must contain both treated (1) and control (0) units.")
  }
  if (nrow(kern) != n || ncol(kern) != n) {
    stop("Kernel matrix dimensions must match the length of 'trt'.")
  }

  # Ensure kern is a regular matrix for the element-wise operations below.
  # The final linear system is solved with a sparse representation of K.
  kern <- as.matrix(kern)

  # Quadratic form matrix for within-group kernel sums
  K <- trt * t(trt * t(kern)) / n1^2 +
       (1 - trt) * t((1 - trt) * t(kern)) / n0^2

  # Constraint matrix: column 1 selects treated, column 2 selects control
  A <- matrix(0, nrow = n, ncol = 2)
  A[, 1] <- trt
  A[, 2] <- 1 - trt

  # Right-hand side: cross-group kernel sums (target)
  b <- as.vector(rowSums(trt * kern)) / (n1 * n) +
       as.vector(rowSums((1 - trt) * kern)) / (n0 * n)

  # Solve via sparse Cholesky
  K_sparse <- Matrix::Matrix(K, sparse = TRUE)
  f1 <- solve(K_sparse, trt)
  f2 <- solve(K_sparse, (1 - trt))
  M1 <- cbind(f1, f2)
  X  <- t(A) %*% M1
  f  <- c(n1, n0)
  YY <- t(A) %*% solve(K_sparse, b) - f
  z  <- b - A %*% solve(X) %*% YY
  w  <- drop(solve(K_sparse, z))

  list(weights = w)
}
