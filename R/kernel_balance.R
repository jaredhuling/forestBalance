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
#'   \code{NULL} if \code{Z} is provided. Required for \code{solver = "direct"}
#'   (if not provided but \code{Z} is available, the kernel is formed
#'   automatically, though this is \eqn{O(n^2)} and may be slow for large
#'   \eqn{n}).
#' @param Z Optional sparse indicator matrix from
#'   \code{\link{leaf_node_kernel_Z}} such that \eqn{K = Z Z^\top / B}. When
#'   supplied, the iterative solvers (\code{"cg"}, \code{"bj"}) can perform
#'   matrix-free products without forming the full kernel. Required for
#'   \code{solver = "cg"} and \code{solver = "bj"}.
#' @param leaf_matrix Optional integer matrix of leaf node assignments
#'   (observations x trees), as returned by \code{\link{get_leaf_node_matrix}}.
#'   Required for \code{solver = "bj"} (Block Jacobi preconditioner uses
#'   tree 1's leaf partition). If \code{NULL} and \code{solver = "bj"}, falls
#'   back to \code{"cg"} with a warning.
#' @param num.trees Number of trees \eqn{B}. Required when \code{Z} is
#'   provided.
#' @param solver Which linear solver to use. \code{"auto"} (default) selects
#'   the best available solver based on the inputs: \code{"cg"} when \code{Z}
#'   is available and \eqn{n > 5000}, or \code{"direct"} otherwise. See
#'   Details for solver requirements.
#' @param tol Convergence tolerance for iterative solvers. Default is
#'   \code{1e-8}.
#' @param maxiter Maximum iterations for iterative solvers. Default is 2000.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{weights}{A numeric vector of length \eqn{n} containing the balancing
#'     weights. Treated weights sum to \eqn{n_1} and control weights sum to
#'     \eqn{n_0}.}
#'   \item{solver}{The solver that was used.}
#' }
#'
#' @details
#' The modified kernel \eqn{K_q} used in the optimization is block-diagonal:
#' the treated--control cross-blocks are zero because
#' \eqn{K_q(i,j) = 0} whenever \eqn{A_i \neq A_j}. All solvers exploit this
#' structure by working on the treated and control blocks independently.
#'
#' \strong{Solver requirements:}
#' \tabular{lll}{
#'   Solver \tab Required inputs \tab Optional inputs \cr
#'   \code{"direct"} \tab \code{kern} (or \code{Z} + \code{num.trees}) \tab \cr
#'   \code{"cg"} \tab \code{Z} + \code{num.trees} \tab \cr
#'   \code{"bj"} \tab \code{Z} + \code{num.trees} + \code{leaf_matrix}
#'     \tab (falls back to \code{"cg"} if \code{leaf_matrix} is missing)
#' }
#'
#' The \strong{direct} solver extracts sub-blocks of the kernel and solves via
#' sparse Cholesky. If only \code{Z} is provided, the kernel is formed as
#' \eqn{K = Z Z^\top / B}, which requires \eqn{O(n^2)} time and memory.
#'
#' The \strong{CG} solver uses the factored representation \eqn{K = Z Z^\top / B}
#' to perform matrix--vector products without forming any kernel matrix.
#'
#' The \strong{Block Jacobi} solver (\code{"bj"}) uses the first tree's leaf
#' partition (from \code{leaf_matrix}) to define a block-diagonal
#' preconditioner for CG. Each leaf block is a small dense system that is
#' cheap to factor.
#'
#' Only 2 linear solves per block are needed (not 3) because the third
#' right-hand side is a linear combination of the first two.
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
kernel_balance <- function(trt, kern = NULL, Z = NULL, leaf_matrix = NULL,
                           num.trees = NULL,
                           solver = c("auto", "direct", "cg", "bj"),
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

  # Choose solver adaptively based on available inputs.
  if (solver == "auto") {
    if (!is.null(Z) && n > 5000) {
      solver <- "cg"
    } else {
      solver <- "direct"
    }
  }

  # Validate solver/input compatibility
  if (solver %in% c("cg", "bj") && is.null(Z)) {
    stop("solver = \"", solver, "\" requires the 'Z' matrix ",
         "(sparse indicator from leaf_node_kernel_Z). ",
         "Use solver = \"direct\" with 'kern', or provide 'Z'.")
  }
  if (solver == "bj" && is.null(leaf_matrix)) {
    warning("solver = \"bj\" requires 'leaf_matrix' for the block Jacobi ",
            "preconditioner. Falling back to solver = \"cg\".")
    solver <- "cg"
  }
  if (solver == "direct" && is.null(kern) && is.null(Z)) {
    stop("solver = \"direct\" requires either 'kern' or 'Z' + 'num.trees'.")
  }

  idx_t <- which(trt == 1)
  idx_c <- which(trt == 0)
  ones_t <- rep(1, n1)
  ones_c <- rep(1, n0)

  if (solver %in% c("cg", "bj")) {
    # ------------------------------------------------------------------
    # Iterative solver path (CG or Block Jacobi preconditioned CG)
    # ------------------------------------------------------------------
    B <- num.trees

    # b vector via Z: rowSums(K) = Z %*% colSums(Z) / B
    rs <- as.numeric(Z %*% Matrix::colSums(Z)) / B
    b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)

    Z_t <- Z[idx_t, ]
    Z_c <- Z[idx_c, ]

    if (solver == "bj") {
      # Block Jacobi: build preconditioner from tree 1's leaf partition
      solve_t <- .bj_pcg_solver(Z_t, leaf_matrix[idx_t, ], B, tol, maxiter)
      solve_c <- .bj_pcg_solver(Z_c, leaf_matrix[idx_c, ], B, tol, maxiter)
    } else {
      # Plain CG (Rcpp)
      Z_t_csc <- as(Z_t, "dgCMatrix")
      Z_c_csc <- as(Z_c, "dgCMatrix")
      solve_t <- function(rhs) as.numeric(cg_solve_cpp(Z_t_csc, B * rhs, tol, maxiter))
      solve_c <- function(rhs) as.numeric(cg_solve_cpp(Z_c_csc, B * rhs, tol, maxiter))
    }

    # Treated block: 2 solves, 3rd by linear combination
    s1 <- solve_t(ones_t);   sb <- solve_t(b[idx_t])
    X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
    w_t <- n1^2 * (sb - (YY1 / X11) * s1)

    # Control block
    s1 <- solve_c(ones_c);   sb <- solve_c(b[idx_c])
    X22 <- n0^2 * sum(s1); YY2 <- n0^2 * sum(sb) - n0
    w_c <- n0^2 * (sb - (YY2 / X22) * s1)

  } else {
    # ------------------------------------------------------------------
    # Direct solver: sparse Cholesky on sub-blocks of K
    # ------------------------------------------------------------------
    if (is.null(kern)) {
      kern <- Matrix::tcrossprod(Z) / num.trees
    }
    if (nrow(kern) != n || ncol(kern) != n) {
      stop("Kernel matrix dimensions must match the length of 'trt'.")
    }

    rs <- as.numeric(Matrix::rowSums(kern))
    b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)

    K_tt <- kern[idx_t, idx_t]
    K_cc <- kern[idx_c, idx_c]

    # Treated block: 2 solves, 3rd by linear combination
    s1 <- as.numeric(solve(K_tt, ones_t))
    sb <- as.numeric(solve(K_tt, b[idx_t]))
    X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
    w_t <- n1^2 * (sb - (YY1 / X11) * s1)

    # Control block
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


# Block Jacobi preconditioned CG solver.
# Returns a function(rhs) that solves Z_g Z_g^T x = B * rhs using
# tree 1's leaf partition as a block-diagonal preconditioner.
# @noRd
.bj_pcg_solver <- function(Z_g, lm_g, B, tol, maxiter) {
  Z_g_csc <- as(Z_g, "dgCMatrix")
  ng <- nrow(Z_g)

  # Build block-diagonal preconditioner from tree 1's leaves.
  # If any block is singular, fall back to identity (no preconditioning).
  leaves <- lm_g[, 1]
  groups <- split(seq_len(ng), leaves)
  block_solvers <- lapply(groups, function(idx) {
    K_block <- Matrix::tcrossprod(Z_g[idx, , drop = FALSE]) / B
    tryCatch(
      { ch <- chol(as.matrix(K_block)); function(v) backsolve(ch, forwardsolve(t(ch), v)) },
      error = function(e) function(v) v  # identity fallback
    )
  })

  precondition <- function(v) {
    result <- numeric(ng)
    for (g in seq_along(groups)) {
      result[groups[[g]]] <- block_solvers[[g]](v[groups[[g]]])
    }
    result
  }

  # Return a solve function
  function(rhs) {
    rhs_scaled <- B * rhs
    Kv <- function(v) as.numeric(Z_g_csc %*% Matrix::crossprod(Z_g_csc, v))

    x <- numeric(ng)
    r <- rhs_scaled - Kv(x)
    z <- precondition(r)
    p <- z
    rz <- sum(r * z)
    rhs_norm <- sqrt(sum(rhs_scaled^2))

    for (i in seq_len(maxiter)) {
      Ap <- Kv(p)
      pAp <- sum(p * Ap)
      if (pAp <= 0 || !is.finite(pAp)) break
      alpha <- rz / pAp
      x <- x + alpha * p
      r <- r - alpha * Ap
      rnorm2 <- sum(r * r)
      if (!is.finite(rnorm2) || sqrt(rnorm2) / rhs_norm < tol) break
      z <- precondition(r)
      rz_new <- sum(r * z)
      if (!is.finite(rz_new) || rz_new == 0) break
      p <- z + (rz_new / rz) * p
      rz <- rz_new
    }
    x
  }
}


# Note: Plain CG solver is implemented in C++ (src/cg_solve.cpp) as cg_solve_cpp().
