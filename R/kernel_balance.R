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
#' @param estimand Target estimand: \code{"ATE"} (default), \code{"ATT"}, or \code{"ATC"}.
#' @param solver Which linear solver to use. \code{"auto"} (default) selects
#'   the best available solver based on the inputs: \code{"cg"} when \code{Z}
#'   is available and \eqn{n > 5000}, or \code{"direct"} otherwise. See
#'   Details for solver requirements.
#' @param lambda Nonnegative ridge penalty on the weights. The group-\eqn{a} problem
#'   solved is \eqn{\tfrac12 \| \sum_{i \in S_a} w_i K(X_i,\cdot) - n_a P_n \|^2 + \tfrac{\lambda}{2}\|w_a\|^2}
#'   subject to \eqn{\sum_{i \in S_a} w_i = n_a}, whose KKT system is
#'   \eqn{(K_{aa} + \lambda I) w_a = b_a + \gamma_a 1}. Default \code{0} (no ridge).
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
#'   \item{diagnostics}{Residuals and convergence checks for every linear solve, arm-total errors, and projected KKT residuals. For positive ridge and a positive semidefinite kernel, the reported weight error bound includes the arm-total correction. Failed accuracy checks raise an error carrying diagnostics.}
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
#' # --- Direct solver (using the kernel matrix) ---
#' forest <- multi_regression_forest(X, cbind(A, Y), num.trees = 500)
#' K <- forest_kernel(forest)
#' bal_direct <- kernel_balance(A, kern = K, solver = "direct")
#'
#' # --- CG solver (using the Z matrix, avoids forming K) ---
#' # Step 1: extract leaf node assignments (n x B matrix)
#' leaf_mat <- get_leaf_node_matrix(forest, X)
#'
#' # Step 2: build sparse indicator matrix Z such that K = Z Z' / B
#' Z <- leaf_node_kernel_Z(leaf_mat)
#'
#' # Step 3: solve with CG (matrix-free, no kernel formed)
#' bal_cg <- kernel_balance(A, Z = Z, num.trees = 500, solver = "cg")
#'
#' # Both solvers give the same weights
#' max(abs(bal_direct$weights - bal_cg$weights))
#'
#' # Weighted ATE estimate
#' w <- bal_cg$weights
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
                           estimand = c("ATE", "ATT", "ATC"),
                           solver = c("auto", "direct", "cg", "bj"),
                           lambda = 0,
                           tol = 1e-8, maxiter = 2000L) {
  solver <- match.arg(solver)
  estimand <- match.arg(estimand)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("'lambda' must be a single nonnegative number.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("'tol' must be finite and positive.")
  if (length(maxiter) != 1L || !is.finite(maxiter) || maxiter < 1 || maxiter != floor(maxiter))
    stop("'maxiter' must be a positive integer.")
  if (any(!is.finite(trt)) || !all(trt %in% 0:1)) stop("Treatment vector must be binary (0/1).")
  if (!is.null(Z) && (nrow(Z) != length(trt) || !.matrix_all_finite(Z)))
    stop("Invalid Z dimensions or nonfinite entries.")
  if (!is.null(kern) && !.matrix_all_finite(kern)) stop("Nonfinite kernel entries.")
  if (!is.null(num.trees) && (length(num.trees) != 1L || !is.finite(num.trees) || num.trees <= 0))
    stop("'num.trees' must be finite and positive.")
  diagnostics <- list(solves = list(), arms = list())
  fail <- function(message) stop(structure(list(message = message, call = NULL,
    diagnostics = diagnostics), class = c("forest_balance_solver_error", "error", "condition")))
  # Recompute actual residuals with the original operator. Sparse paths stay sparse.
  audit_solve <- function(res, rhs, arm) {
    it <- attr(res, "iters")
    x <- as.numeric(res)
    ix <- if (arm == "treated") which(trt == 1) else which(trt == 0)
    op <- function(x) {
      if (solver %in% c("cg", "bj")) {
        za <- Z[ix, , drop = FALSE]
        as.numeric(za %*% Matrix::crossprod(za, x)) / num.trees + lambda * x
      } else as.numeric(kern[ix, ix, drop = FALSE] %*% x) + lambda * x
    }
    residual <- sqrt(sum((op(x) - rhs)^2))
    # CG stops on an absolute B-scaled residual; BJ uses a relative residual.
    target <- if (solver == "cg") tol / num.trees else tol * sqrt(sum(rhs^2))
    op_bound <- if (solver %in% c("cg", "bj")) {
      za_abs <- abs(Z[ix, , drop = FALSE])
      max(as.numeric(za_abs %*% Matrix::colSums(za_abs))) / num.trees + lambda
    } else max(Matrix::rowSums(abs(kern[ix, ix, drop = FALSE]))) + lambda
    roundoff <- 100 * .Machine$double.eps * (sqrt(sum(rhs^2)) + op_bound * sqrt(sum(x^2)))
    limit <- 10 * target + roundoff
    ok <- all(is.finite(x)) && is.finite(residual) && residual <= limit
    diagnostics$solves[[length(diagnostics$solves) + 1L]] <<- list(
      arm = arm, rhs = if (all(rhs == 1)) "ones" else "target",
      iterations = if (is.null(it)) NA_integer_ else it,
      residual_l2 = residual, residual_limit = limit, converged = ok)
    if (!ok) fail(paste("Linear solve did not converge for", arm, "arm."))
    note_iters(res)
  }
  iters_env <- new.env()
  iters_env$max <- NA_integer_
  note_iters <- function(res) {
    it <- attr(res, "iters")
    if (!is.null(it)) iters_env$max <- max(iters_env$max, it, na.rm = TRUE)
    as.numeric(res)
  }

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

    Z_t <- Z[idx_t, , drop = FALSE]
    Z_c <- Z[idx_c, , drop = FALSE]

    # Compute b vector based on estimand
    if (estimand == "ATE") {
      rs <- as.numeric(Z %*% Matrix::colSums(Z)) / B
      b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)
      b_t <- b[idx_t]
      b_c <- b[idx_c]
    } else if (estimand == "ATT") {
      # Target is treated distribution; only control block needs solving
      cs_zt <- Matrix::colSums(Z_t)
      rs_to_t <- as.numeric(Z_c %*% cs_zt) / B
      b_c <- rs_to_t / (n0 * n1)
    } else {
      # ATC: target is control distribution; only treated block needs solving
      cs_zc <- Matrix::colSums(Z_c)
      rs_to_c <- as.numeric(Z_t %*% cs_zc) / B
      b_t <- rs_to_c / (n1 * n0)
    }

    # Build solvers only for blocks that need them
    if (estimand != "ATC") {
      # Need control solver
      if (solver == "bj") {
        solve_c0 <- .bj_pcg_solver(Z_c, leaf_matrix[idx_c, ], B, tol, maxiter, lambda)
        solve_c <- function(rhs) audit_solve(solve_c0(rhs), rhs, "control")
      } else {
        Z_c_csc <- as(Z_c, "dgCMatrix")
        solve_c <- function(rhs) audit_solve(cg_solve_cpp(Z_c_csc, B * rhs, tol, maxiter, B * lambda), rhs, "control")
      }
    }
    if (estimand != "ATT") {
      # Need treated solver
      if (solver == "bj") {
        solve_t0 <- .bj_pcg_solver(Z_t, leaf_matrix[idx_t, ], B, tol, maxiter, lambda)
        solve_t <- function(rhs) audit_solve(solve_t0(rhs), rhs, "treated")
      } else {
        Z_t_csc <- as(Z_t, "dgCMatrix")
        solve_t <- function(rhs) audit_solve(cg_solve_cpp(Z_t_csc, B * rhs, tol, maxiter, B * lambda), rhs, "treated")
      }
    }

    # Solve each block
    if (estimand == "ATT") {
      w_t <- ones_t
    } else {
      s1 <- solve_t(ones_t);   sb <- solve_t(b_t)
      X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
      w_t <- n1^2 * (sb - (YY1 / X11) * s1)
    }

    if (estimand == "ATC") {
      w_c <- ones_c
    } else {
      s1 <- solve_c(ones_c);   sb <- solve_c(b_c)
      X22 <- n0^2 * sum(s1); YY2 <- n0^2 * sum(sb) - n0
      w_c <- n0^2 * (sb - (YY2 / X22) * s1)
    }

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

    # Compute b vector based on estimand
    if (estimand == "ATE") {
      rs <- as.numeric(Matrix::rowSums(kern))
      b  <- trt * rs / (n1 * n) + (1 - trt) * rs / (n0 * n)
      b_t <- b[idx_t]
      b_c <- b[idx_c]
    } else if (estimand == "ATT") {
      b_c <- as.numeric(kern[idx_c, idx_t, drop = FALSE] %*% ones_t) / (n0 * n1)
    } else {
      b_t <- as.numeric(kern[idx_t, idx_c, drop = FALSE] %*% ones_c) / (n1 * n0)
    }

    # Solve each block
    if (estimand == "ATT") {
      w_t <- ones_t
    } else {
      K_tt <- kern[idx_t, idx_t, drop = FALSE]
      if (lambda > 0) K_tt <- K_tt + lambda * Matrix::Diagonal(n1)
      s1 <- audit_solve(solve(K_tt, ones_t), ones_t, "treated")
      sb <- audit_solve(solve(K_tt, b_t), b_t, "treated")
      X11 <- n1^2 * sum(s1); YY1 <- n1^2 * sum(sb) - n1
      w_t <- n1^2 * (sb - (YY1 / X11) * s1)
    }

    if (estimand == "ATC") {
      w_c <- ones_c
    } else {
      K_cc <- kern[idx_c, idx_c, drop = FALSE]
      if (lambda > 0) K_cc <- K_cc + lambda * Matrix::Diagonal(n0)
      s1 <- audit_solve(solve(K_cc, ones_c), ones_c, "control")
      sb <- audit_solve(solve(K_cc, b_c), b_c, "control")
      X22 <- n0^2 * sum(s1); YY2 <- n0^2 * sum(sb) - n0
      w_c <- n0^2 * (sb - (YY2 / X22) * s1)
    }
  }

  # Projected KKT residual also checks amplification of the small target solve.
  for (a in c("treated", "control")) {
    ix <- if (a == "treated") idx_t else idx_c
    wa <- if (a == "treated") w_t else w_c
    na <- length(ix)
    fixed <- (estimand == "ATT" && a == "treated") || (estimand == "ATC" && a == "control")
    mass_error <- abs(sum(wa) - na)
    if (!all(is.finite(wa)) || mass_error > 1e-7 * na) fail("Invalid arm weights or arm total.")
    if (!fixed) {
      ba <- if (a == "treated") b_t else b_c
      feasible_w <- wa + (na - sum(wa)) / na
      if (solver %in% c("cg", "bj")) {
        za <- Z[ix, , drop = FALSE]
        grad <- as.numeric(za %*% Matrix::crossprod(za, feasible_w))/num.trees + lambda * feasible_w - na^2 * ba
      } else grad <- as.numeric(kern[ix, ix, drop = FALSE] %*% feasible_w) + lambda * feasible_w - na^2 * ba
      projected <- sqrt(sum((grad - mean(grad))^2))
      bound <- if (lambda > 0) projected / lambda + mass_error / sqrt(na) else NA_real_
      # For PSD kernels and positive ridge this bounds distance to the constrained
      # optimum, after the negligible recorded arm-total error is accounted for.
      error_limit <- 1e-5 * (1 + sqrt(sum(wa^2)))

    } else { projected <- 0; bound <- 0; error_limit <- 0 }
    diagnostics$arms[[a]] <- list(mass_error = mass_error,
      projected_residual_l2 = projected, weight_error_bound = bound,
      weight_error_limit = error_limit)
    if (!is.finite(projected) || (lambda > 0 && bound > error_limit))
      fail("Final weight KKT accuracy check failed.")
  }
  # Reassemble
  w <- numeric(n)
  w[idx_t] <- w_t
  w[idx_c] <- w_c

  list(weights = w, solver = solver, estimand = estimand, lambda = lambda,
       iters = iters_env$max, diagnostics = diagnostics)
}


# Block Jacobi preconditioned CG solver.
# Returns a function(rhs) that solves (Z_g Z_g^T + B * lambda * I) x = B * rhs using
# tree 1's leaf partition (plus the ridge) as a block-diagonal preconditioner.
# The result carries the iteration count as attribute "iters".
# @noRd
.bj_pcg_solver <- function(Z_g, lm_g, B, tol, maxiter, lambda = 0) {
  Z_g_csc <- as(Z_g, "dgCMatrix")
  ng <- nrow(Z_g)

  # Build block-diagonal preconditioner from tree 1's leaves.
  # If any block is singular, fall back to identity (no preconditioning).
  leaves <- lm_g[, 1]
  groups <- split(seq_len(ng), leaves)
  block_solvers <- lapply(groups, function(idx) {
    K_block <- Matrix::tcrossprod(Z_g[idx, , drop = FALSE]) / B
    if (lambda > 0) K_block <- K_block + lambda * Matrix::Diagonal(length(idx))
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
    Kv <- function(v) as.numeric(Z_g_csc %*% Matrix::crossprod(Z_g_csc, v)) + B * lambda * v

    x <- numeric(ng)
    r <- rhs_scaled - Kv(x)
    z <- precondition(r)
    p <- z
    rz <- sum(r * z)
    rhs_norm <- sqrt(sum(rhs_scaled^2))
    iters <- 0L

    for (i in seq_len(maxiter)) {
      iters <- i
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
    attr(x, "iters") <- iters
    x
  }
}


# Note: Plain CG solver is implemented in C++ (src/cg_solve.cpp) as cg_solve_cpp().

# Avoid densifying sparse matrices when implicit zeros are necessarily finite.
.matrix_all_finite <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(x)) all(is.finite(methods::slot(x, "x"))) else TRUE
  } else all(is.finite(x))
}
