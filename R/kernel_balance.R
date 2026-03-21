#' Kernel energy balancing weights via closed-form solution
#'
#' Computes balancing weights that minimize a kernelized energy distance between
#' the weighted treated and control distributions and the overall sample. The
#' weights are obtained via a closed-form solution to a linear system derived
#' from the kernel energy distance objective.
#'
#' @param trt A binary (0/1) integer or numeric vector indicating treatment
#'   assignment (\code{1} = treated, \code{0} = control).
#' @param kern A symmetric \eqn{n \times n} kernel matrix (dense or sparse), or
#'   \code{NULL} if \code{Z} is provided.
#' @param Z Optional sparse indicator matrix from
#'   \code{\link{leaf_node_kernel_Z}} such that \eqn{K = Z Z^\top / B}. When
#'   supplied, the solver can avoid forming the full kernel matrix. If both
#'   \code{kern} and \code{Z} are given, \code{Z} takes priority when the CG
#'   solver is selected.
#' @param num.trees Number of trees \eqn{B}. Required when \code{Z} is
#'   provided.
#' @param solver Which linear solver to use. \code{"auto"} (default) selects
#'   \code{"direct"} for \eqn{n \le 5000} and \code{"cg"} for
#'   \eqn{n > 5000}. \code{"direct"} uses sparse Cholesky on the treated and
#'   control sub-blocks of the kernel. \code{"cg"} uses conjugate gradient
#'   iterations with the factored \eqn{Z} representation, avoiding formation of
#'   any kernel matrix.
#' @param tol Convergence tolerance for the CG solver. Default is \code{1e-8}.
#'   Ignored when \code{solver = "direct"}.
#' @param maxiter Maximum CG iterations. Default is 1000.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{weights}{A numeric vector of length \eqn{n} containing the balancing
#'     weights. Treated weights sum to \eqn{n_1} and control weights sum to
#'     \eqn{n_0}.}
#'   \item{solver}{The solver that was used (\code{"direct"} or \code{"cg"}).}
#' }
#'
#' @details
#' The modified kernel \eqn{K_q} used in the optimization is block-diagonal:
#' the treated--control cross-blocks are zero because
#' \eqn{K_q(i,j) = 0} whenever \eqn{A_i \neq A_j}. Both solvers exploit this
#' structure by working on the treated and control blocks independently.
#'
#' The \strong{direct} solver extracts the sub-blocks \eqn{K_{tt}} and
#' \eqn{K_{cc}} and solves via sparse Cholesky. This gives exact solutions
#' but requires forming (at least sub-blocks of) the kernel matrix.
#'
#' The \strong{CG} solver uses the factored representation \eqn{K = Z Z^\top / B}
#' to perform matrix--vector products without forming any kernel matrix,
#' via \eqn{K v = Z (Z^\top v) / B}. This is much faster and more
#' memory-efficient at large \eqn{n} (e.g., \eqn{n > 5000}). The CG iterates
#' converge to the exact solution; the default tolerance of \code{5e-11} yields
#' weight vectors that agree with the direct solution to several decimal places.
#'
#' @references
#' De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for causal
#' effect estimation for high dimensional data.
#' \emph{arXiv preprint arXiv:2512.18069}.
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
#' @importFrom methods as
#' @importFrom Matrix Matrix rowSums forceSymmetric sparseMatrix colSums
#' @importFrom Matrix crossprod tcrossprod
#' @importMethodsFrom Matrix solve
#' @export
kernel_balance <- function(trt, kern = NULL, Z = NULL, num.trees = NULL,
                           solver = c("auto", "direct", "cg"),
                           tol = 1e-8, maxiter = 2000L) {
  solver <- match.arg(solver)

  trt <- as.double(trt)
  n  <- length(trt)
  n1 <- sum(trt)
  n0 <- n - n1

  if (n1 == 0 || n0 == 0) {
    stop("Treatment vector must contain both treated (1) and control (0) units.")
  }
  if (is.null(kern) && is.null(Z)) {
    stop("Either 'kern' or 'Z' must be provided.")
  }
  if (!is.null(Z) && is.null(num.trees)) {
    stop("'num.trees' is required when 'Z' is provided.")
  }

  # Choose solver adaptively
  if (solver == "auto") {
    solver <- if (n > 5000 && !is.null(Z)) "cg" else "direct"
  }

  idx_t <- which(trt == 1)
  idx_c <- which(trt == 0)
  ones_t <- rep(1, n1)
  ones_c <- rep(1, n0)

  if (solver == "cg") {
    # ------------------------------------------------------------------
    # CG solver: uses Z factor, never forms the kernel matrix
    # ------------------------------------------------------------------
    if (is.null(Z)) {
      stop("CG solver requires the Z matrix. Supply it via the 'Z' argument.")
    }
    B <- num.trees

    # b vector via Z: rowSums(K) = Z %*% colSums(Z) / B
    rs <- as.numeric(Z %*% Matrix::colSums(Z)) / B
    b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)

    Z_t <- Z[idx_t, ]
    Z_c <- Z[idx_c, ]

    # CG helpers: solve Z_g Z_g^T x = B * rhs  for group g (C++ solver)
    Z_t_csc <- as(Z_t, "dgCMatrix")
    Z_c_csc <- as(Z_c, "dgCMatrix")
    cg_t <- function(rhs) as.numeric(cg_solve_cpp(Z_t_csc, B * rhs, tol, maxiter))
    cg_c <- function(rhs) as.numeric(cg_solve_cpp(Z_c_csc, B * rhs, tol, maxiter))

    # Treated block: only 2 solves needed.
    # The 3rd RHS z_t = b_t - c*1 is a linear combination, so
    # solve(K, z_t) = solve(K, b_t) - c * solve(K, 1).
    s1 <- cg_t(ones_t);   sb <- cg_t(b[idx_t])
    X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
    w_t <- n1^2 * (sb - (YY1 / X11) * s1)

    # Control block: same 2-solve optimization.
    s1 <- cg_c(ones_c);   sb <- cg_c(b[idx_c])
    X22 <- n0^2 * sum(s1); YY2 <- n0^2 * sum(sb) - n0
    w_c <- n0^2 * (sb - (YY2 / X22) * s1)

  } else {
    # ------------------------------------------------------------------
    # Direct solver: sparse Cholesky on sub-blocks of K
    # ------------------------------------------------------------------
    if (is.null(kern)) {
      # Build kernel from Z if not provided
      kern <- Matrix::tcrossprod(Z) / num.trees
    }
    if (nrow(kern) != n || ncol(kern) != n) {
      stop("Kernel matrix dimensions must match the length of 'trt'.")
    }

    # b vector: rowSums(trt * K) = trt * rowSums(K)
    rs <- as.numeric(Matrix::rowSums(kern))
    b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)

    # Sub-blocks of the block-diagonal modified kernel
    K_tt <- kern[idx_t, idx_t]
    K_cc <- kern[idx_c, idx_c]

    # Treated block: 2 solves (Cholesky cached after first).
    # The 3rd RHS is a linear combination, so we recombine solutions.
    s1 <- as.numeric(solve(K_tt, ones_t))
    sb <- as.numeric(solve(K_tt, b[idx_t]))
    X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
    w_t <- n1^2 * (sb - (YY1 / X11) * s1)

    # Control block: same 2-solve optimization.
    s1 <- as.numeric(solve(K_cc, ones_c))
    sb <- as.numeric(solve(K_cc, b[idx_c]))
    X22 <- n0^2 * sum(s1); YY2 <- n0^2 * sum(sb) - n0
    w_c <- n0^2 * (sb - (YY2 / X22) * s1)
  }

  # Reassemble
  w <- numeric(n)
  w[idx_t] <- w_t
  w[idx_c] <- w_c

  list(weights = w, solver = solver)
}



# Note: CG solver is now implemented in C++ (src/cg_solve.cpp) as cg_solve_cpp().
# The C++ implementation uses Eigen sparse mat-vecs with a compiled loop,
# giving ~2x speedup over the equivalent R loop by eliminating interpreter
# overhead, temporary vector allocations, and garbage collection.
