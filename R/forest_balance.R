#' Estimate ATE using forest-based kernel energy balancing
#'
#' Fits a multivariate random forest that jointly models the relationship
#' between covariates, treatment, and outcome, computes a random forest
#' proximity kernel, and then uses kernel energy balancing to produce weights
#' for estimating the average treatment effect (ATE). By default, K-fold
#' cross-fitting is used to avoid overfitting bias from estimating the kernel
#' on the same data used for treatment effect estimation.
#'
#' @param X A numeric matrix or data frame of covariates (\eqn{n \times p}).
#' @param A A binary (0/1) vector of treatment assignments.
#' @param Y A numeric vector of outcomes.
#' @param num.trees Number of trees to grow in the forest. Default is 1000.
#' @param min.node.size Minimum number of observations per leaf node. If
#'   \code{NULL} (default), an adaptive heuristic is used:
#'   \code{max(20, min(floor(n/200) + p, floor(n/50)))}. This scales the leaf
#'   size with both the sample size and the number of covariates, which
#'   empirically yields low bias. See Details.
#' @param cross.fitting Logical; if \code{TRUE} (default), use K-fold
#'   cross-fitting to construct the kernel from held-out data, reducing
#'   overfitting bias. If \code{FALSE}, the kernel is estimated on the full
#'   sample.
#' @param num.folds Number of cross-fitting folds. Default is 2. Only used when
#'   \code{cross.fitting = TRUE}.
#' @param scale.outcomes If \code{TRUE} (default), the joint outcome matrix
#'   \code{cbind(A, Y)} is column-standardized before fitting the forest. This
#'   ensures that treatment and outcome contribute equally to the splits.
#' @param solver Which linear solver to use for the balancing weights.
#'   \code{"auto"} (default) selects \code{"direct"} for small fold sizes
#'   and \code{"cg"} for large fold sizes. See \code{\link{kernel_balance}} for
#'   details.
#' @param tol Convergence tolerance for the CG solver. Default is \code{5e-11}.
#' @param ... Additional arguments passed to
#'   \code{\link[grf]{multi_regression_forest}}.
#'
#' @return An object of class \code{"forest_balance"} (a list) with the
#'   following elements:
#' \describe{
#'   \item{ate}{The estimated average treatment effect. When cross-fitting is
#'     used, this is the average of per-fold Hajek estimates (DML1).}
#'   \item{weights}{The balancing weight vector (length \eqn{n}). When
#'     cross-fitting is used, these are the concatenated per-fold weights.}
#'   \item{kernel}{The \eqn{n \times n} forest proximity kernel (sparse matrix),
#'     or \code{NULL} when cross-fitting or the CG solver is used.}
#'   \item{forest}{The trained forest object. When cross-fitting is used, this
#'     is the last fold's forest.}
#'   \item{X, A, Y}{The input data.}
#'   \item{n, n1, n0}{Total, treated, and control sample sizes.}
#'   \item{solver}{The solver that was used (\code{"direct"} or \code{"cg"}).}
#'   \item{crossfit}{Logical indicating whether cross-fitting was used.}
#'   \item{num.folds}{Number of folds (if cross-fitting was used).}
#'   \item{fold_ates}{Per-fold ATE estimates (if cross-fitting was used).}
#'   \item{fold_ids}{Fold assignments (if cross-fitting was used).}
#' }
#' The object has \code{print} and \code{summary} methods. Use
#' \code{\link{summary.forest_balance}} for covariate balance diagnostics.
#'
#' @details
#' The method proceeds in three steps:
#' \enumerate{
#'   \item A \code{\link[grf]{multi_regression_forest}} is fit on covariates
#'     \code{X} with a bivariate response \code{(A, Y)}. This jointly models
#'     the relationship between covariates, treatment assignment, and outcome.
#'   \item The forest's leaf co-membership structure defines a proximity kernel:
#'     \eqn{K(i,j)} is the proportion of trees where \eqn{i} and \eqn{j} share
#'     a leaf. Because the forest splits on both \eqn{A} and \eqn{Y}, this
#'     kernel captures confounding structure.
#'   \item \code{\link{kernel_balance}} computes balancing weights via the
#'     closed-form kernel energy distance solution. The ATE is then estimated
#'     using the Hajek (ratio) estimator with these weights.
#' }
#'
#' \strong{Cross-fitting} (default): For each fold \eqn{k}, the forest is
#' trained on all data \emph{except} fold \eqn{k}, and the kernel for fold
#' \eqn{k} is built from that held-out forest's leaf predictions. This breaks
#' the dependence between the kernel and the outcomes, reducing overfitting
#' bias. The final ATE is the average of the per-fold Hajek estimates (DML1).
#'
#' \strong{Adaptive leaf size}: The default \code{min.node.size} is set
#' adaptively via \code{max(20, min(floor(n/200) + p, floor(n/50)))}. Larger
#' leaves produce smoother kernels that generalize better, while the cap at
#' \code{n/50} prevents kernel degeneracy. This heuristic was calibrated
#' empirically to minimize RMSE across a range of sample sizes and dimensions.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#' Newey, W. and Robins, J. (2018). Double/debiased machine learning for
#' treatment and structural parameters. \emph{The Econometrics Journal},
#' 21(1), C1--C68.
#'
#' De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for causal
#' effect estimation for high dimensional data.
#' \emph{arXiv preprint arXiv:2512.18069}.
#'
#' @examples
#' \donttest{
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
#' Y <- X[, 1] + rnorm(n)  # true ATE = 0
#'
#' # Default: cross-fitting with adaptive leaf size
#' result <- forest_balance(X, A, Y)
#' result
#'
#' # Without cross-fitting
#' result_nocf <- forest_balance(X, A, Y, cross.fitting = FALSE)
#' }
#'
#' @importFrom grf multi_regression_forest
#' @importFrom stats weighted.mean
#' @export
forest_balance <- function(X, A, Y,
                           num.trees = 1000,
                           min.node.size = NULL,
                           cross.fitting = TRUE,
                           num.folds = 2,
                           scale.outcomes = TRUE,
                           solver = c("auto", "direct", "cg"),
                           tol = 5e-11,
                           ...) {
  solver <- match.arg(solver)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (length(A) != n || length(Y) != n) {
    stop("X, A, and Y must have the same number of observations.")
  }
  if (!all(A %in% c(0, 1))) {
    stop("Treatment vector A must be binary (0/1).")
  }

  # Adaptive min.node.size heuristic
  if (is.null(min.node.size)) {
    min.node.size <- max(20L, min(floor(n / 200) + p, floor(n / 50)))
  }

  if (cross.fitting) {
    # ------------------------------------------------------------------
    # Cross-fitting path
    # ------------------------------------------------------------------
    if (num.folds < 2 || num.folds > n) {
      stop("num.folds must be between 2 and n.")
    }

    fold_ids  <- sample(rep(seq_len(num.folds), length.out = n))
    fold_ates <- numeric(num.folds)
    weights   <- numeric(n)
    last_forest <- NULL
    last_solver <- NULL

    for (k in seq_len(num.folds)) {
      idx_k    <- which(fold_ids == k)
      idx_notk <- which(fold_ids != k)
      n_k <- length(idx_k)

      A_test <- A[idx_k]; Y_test <- Y[idx_k]

      # Skip if a treatment group is empty in this fold
      if (sum(A_test == 1) == 0 || sum(A_test == 0) == 0) {
        fold_ates[k] <- NA
        next
      }

      # Train forest on held-out folds
      response_train <- cbind(A[idx_notk], Y[idx_notk])
      if (scale.outcomes) response_train <- scale(response_train)

      forest_k <- grf::multi_regression_forest(
        X[idx_notk, , drop = FALSE], Y = response_train,
        num.trees = num.trees, min.node.size = min.node.size, ...
      )

      # Predict leaf nodes for held-in fold
      leaf_mat_k <- get_leaf_node_matrix(forest_k, newdata = X[idx_k, , drop = FALSE])

      # Build kernel / Z and solve.
      # CG is preferred when the fold is large or the kernel would be dense
      # (large min.node.size relative to fold size creates denser kernels).
      use_cg <- (solver == "cg") ||
                (solver == "auto" && (n_k > 5000 || min.node.size > n_k / 20))
      if (use_cg) {
        Z_k <- leaf_node_kernel_Z(leaf_mat_k)
        bal_k <- kernel_balance(trt = A_test, Z = Z_k, num.trees = num.trees,
                                solver = "cg", tol = tol)
      } else {
        K_k <- leaf_node_kernel(leaf_mat_k)
        bal_k <- kernel_balance(trt = A_test, kern = K_k, solver = "direct")
      }

      w_k <- bal_k$weights
      fold_ates[k] <- weighted.mean(Y_test[A_test == 1], w_k[A_test == 1]) -
                       weighted.mean(Y_test[A_test == 0], w_k[A_test == 0])
      weights[idx_k] <- w_k
      last_forest <- forest_k
      last_solver <- bal_k$solver
    }

    ate <- mean(fold_ates, na.rm = TRUE)

    out <- list(
      ate       = ate,
      weights   = weights,
      fold_ates = fold_ates,
      fold_ids  = fold_ids,
      kernel    = NULL,
      forest    = last_forest,
      X         = X,
      A         = A,
      Y         = Y,
      n         = n,
      n1        = as.integer(sum(A == 1)),
      n0        = as.integer(sum(A == 0)),
      solver    = last_solver,
      crossfit  = TRUE,
      num.folds = num.folds
    )

  } else {
    # ------------------------------------------------------------------
    # No cross-fitting path
    # ------------------------------------------------------------------
    response <- cbind(A, Y)
    if (scale.outcomes) response <- scale(response)

    forest <- grf::multi_regression_forest(
      X, Y = response,
      num.trees = num.trees, min.node.size = min.node.size, ...
    )

    leaf_mat <- get_leaf_node_matrix(forest, newdata = X)
    use_cg <- (solver == "cg") ||
              (solver == "auto" && (n > 5000 || min.node.size > n / 20))

    if (use_cg) {
      Z <- leaf_node_kernel_Z(leaf_mat)
      bal <- kernel_balance(trt = A, Z = Z, num.trees = num.trees,
                            solver = "cg", tol = tol)
      K <- NULL
    } else {
      K <- leaf_node_kernel(leaf_mat)
      bal <- kernel_balance(trt = A, kern = K, solver = "direct")
    }
    w <- bal$weights

    ate <- weighted.mean(Y[A == 1], w[A == 1]) -
           weighted.mean(Y[A == 0], w[A == 0])

    out <- list(
      ate       = ate,
      weights   = w,
      kernel    = K,
      forest    = forest,
      X         = X,
      A         = A,
      Y         = Y,
      n         = n,
      n1        = as.integer(sum(A == 1)),
      n0        = as.integer(sum(A == 0)),
      solver    = bal$solver,
      crossfit  = FALSE
    )
  }

  class(out) <- "forest_balance"
  out
}
