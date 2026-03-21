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
#' @param augmented Logical; if \code{TRUE}, use an augmented (doubly-robust)
#'   estimator that combines the kernel energy balancing weights with
#'   group-specific outcome regression models. This reduces bias when either the
#'   kernel or the outcome models are correctly specified. Default is
#'   \code{FALSE}. See Details.
#' @param mu.hat Optional list with components \code{mu1} and \code{mu0}, each
#'   a numeric vector of length \eqn{n}, containing user-supplied predictions of
#'   \eqn{E[Y \mid X, A=1]} and \eqn{E[Y \mid X, A=0]}. When provided, these
#'   are used instead of fitting internal outcome models. If \code{NULL}
#'   (default) and \code{augmented = TRUE}, two
#'   \code{\link[grf]{regression_forest}} models are fit automatically (one on
#'   treated, one on control). When supplying \code{mu.hat} with
#'   \code{cross.fitting = TRUE}, the user is responsible for ensuring the
#'   predictions were cross-fitted externally.
#' @param scale.outcomes If \code{TRUE} (default), the joint outcome matrix
#'   \code{cbind(A, Y)} is column-standardized before fitting the forest. This
#'   ensures that treatment and outcome contribute equally to the splits.
#' @param solver Which linear solver to use for the balancing weights.
#'   \code{"auto"} (default) selects \code{"direct"} for small fold sizes
#'   and \code{"cg"} for large fold sizes. See \code{\link{kernel_balance}} for
#'   details.
#' @param tol Convergence tolerance for the CG solver. Default is \code{5e-11}.
#' @param parallel Logical or integer. If \code{FALSE} (default), folds are
#'   processed sequentially. If \code{TRUE}, folds are processed in parallel
#'   using all available cores via \code{\link[parallel]{mclapply}}. An integer
#'   value specifies the exact number of cores. Only used when
#'   \code{cross.fitting = TRUE}. Note: parallel processing is not supported on
#'   Windows.
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
#'   \item{mu1.hat}{Predictions of \eqn{E[Y|X, A=1]} (length \eqn{n}), or
#'     \code{NULL} if \code{augmented = FALSE}.}
#'   \item{mu0.hat}{Predictions of \eqn{E[Y|X, A=0]} (length \eqn{n}), or
#'     \code{NULL} if \code{augmented = FALSE}.}
#'   \item{kernel}{The \eqn{n \times n} forest proximity kernel (sparse matrix),
#'     or \code{NULL} when cross-fitting or the CG solver is used.}
#'   \item{forest}{The trained forest object. When cross-fitting is used, this
#'     is the last fold's forest.}
#'   \item{X, A, Y}{The input data.}
#'   \item{n, n1, n0}{Total, treated, and control sample sizes.}
#'   \item{solver}{The solver that was used (\code{"direct"} or \code{"cg"}).}
#'   \item{crossfit}{Logical indicating whether cross-fitting was used.}
#'   \item{augmented}{Logical indicating whether augmentation was used.}
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
#' \strong{Augmented estimator}: When \code{augmented = TRUE}, two
#' group-specific outcome models \eqn{\hat\mu_1(X) = E[Y|X, A=1]} and
#' \eqn{\hat\mu_0(X) = E[Y|X, A=0]} are fit, and the ATE is estimated via
#' the doubly-robust formula:
#' \deqn{\hat\tau = \frac{1}{n}\sum_i [\hat\mu_1(X_i) - \hat\mu_0(X_i)]
#'   + \frac{\sum w_i A_i (Y_i - \hat\mu_1(X_i))}{\sum w_i A_i}
#'   - \frac{\sum w_i (1-A_i)(Y_i - \hat\mu_0(X_i))}{\sum w_i (1-A_i)}.}
#' The first term is the regression-based estimate of the ATE; the remaining
#' terms are weighted bias corrections. This is consistent if either the kernel
#' (balancing weights) or the outcome models are correctly specified. When
#' combined with cross-fitting, the outcome models are automatically
#' cross-fitted in lockstep with the kernel.
#'
#' \strong{Adaptive leaf size}: The default \code{min.node.size} is set
#' adaptively via \code{max(20, min(floor(n/200) + p, floor(n/50)))}. Larger
#' leaves produce smoother kernels that generalize better, while the cap at
#' \code{n/50} prevents kernel degeneracy.
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
#' # Augmented (doubly-robust) estimator
#' result_aug <- forest_balance(X, A, Y, augmented = TRUE)
#'
#' # Without cross-fitting
#' result_nocf <- forest_balance(X, A, Y, cross.fitting = FALSE)
#' }
#'
#' @importFrom grf multi_regression_forest regression_forest
#' @importFrom stats predict weighted.mean
#' @export
forest_balance <- function(X, A, Y,
                           num.trees = 1000,
                           min.node.size = NULL,
                           cross.fitting = TRUE,
                           num.folds = 2,
                           augmented = FALSE,
                           mu.hat = NULL,
                           scale.outcomes = TRUE,
                           solver = c("auto", "direct", "cg", "bj"),
                           tol = 1e-8,
                           parallel = FALSE,
                           ...) {
  solver <- match.arg(solver)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  .validate_inputs(X, A, Y, mu.hat, cross.fitting, num.folds, n)
  if (!is.null(mu.hat)) augmented <- TRUE

  if (is.null(min.node.size)) {
    min.node.size <- .adaptive_min_node_size(n, p)
  }

  if (cross.fitting) {
    result <- .fit_crossfitted(X, A, Y, num.trees, min.node.size, num.folds,
                               augmented, mu.hat, scale.outcomes, solver,
                               tol, parallel, ...)
  } else {
    result <- .fit_full_sample(X, A, Y, num.trees, min.node.size,
                               augmented, mu.hat, scale.outcomes, solver,
                               tol, ...)
  }

  out <- c(result, list(
    X         = X,
    A         = A,
    Y         = Y,
    n         = n,
    n1        = as.integer(sum(A == 1)),
    n0        = as.integer(sum(A == 0)),
    crossfit  = cross.fitting,
    augmented = augmented
  ))
  class(out) <- "forest_balance"
  out
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' Validate inputs to forest_balance
#' @noRd
.validate_inputs <- function(X, A, Y, mu.hat, cross.fitting, num.folds, n) {
  if (length(A) != n || length(Y) != n) {
    stop("X, A, and Y must have the same number of observations.")
  }
  if (!all(A %in% c(0, 1))) {
    stop("Treatment vector A must be binary (0/1).")
  }
  if (!is.null(mu.hat)) {
    if (!is.list(mu.hat) || is.null(mu.hat$mu1) || is.null(mu.hat$mu0)) {
      stop("mu.hat must be a list with components 'mu1' and 'mu0'.")
    }
    if (length(mu.hat$mu1) != n || length(mu.hat$mu0) != n) {
      stop("mu.hat$mu1 and mu.hat$mu0 must each have length n.")
    }
    if (cross.fitting) {
      message("Note: user-supplied mu.hat used with cross.fitting = TRUE. ",
              "Ensure predictions were cross-fitted externally.")
    }
  }
  if (cross.fitting && (num.folds < 2 || num.folds > n)) {
    stop("num.folds must be between 2 and n.")
  }
}


#' Compute adaptive min.node.size
#' @noRd
.adaptive_min_node_size <- function(n, p) {
  max(20L, min(floor(n / 200) + p, floor(n / 50)))
}


#' Choose solver based on sample size and kernel density
#' @noRd
.choose_solver <- function(solver, n_obs, min.node.size) {
  if (solver != "auto") return(solver)
  if (n_obs > 5000) "bj"
  else if (min.node.size > n_obs / 20) "cg"
  else "direct"
}


#' Train joint forest and compute balancing weights for a set of observations
#' @return List with components: weights, forest, solver, kernel (or NULL)
#' @noRd
.fit_kernel_and_balance <- function(X_train, A_train, Y_train,
                                    X_pred, A_pred,
                                    num.trees, min.node.size,
                                    scale.outcomes, solver, tol, ...) {
  # Train joint forest
  response <- cbind(A_train, Y_train)
  if (scale.outcomes) response <- scale(response)

  forest <- grf::multi_regression_forest(
    X_train, Y = response,
    num.trees = num.trees, min.node.size = min.node.size, ...
  )

  # Extract leaf nodes and build kernel / Z
  leaf_mat <- get_leaf_node_matrix(forest, newdata = X_pred)
  n_pred <- nrow(X_pred)
  eff_solver <- .choose_solver(solver, n_pred, min.node.size)

  if (eff_solver %in% c("bj", "cg")) {
    Z <- leaf_node_kernel_Z(leaf_mat)
    bal <- kernel_balance(trt = A_pred, Z = Z, leaf_matrix = leaf_mat,
                          num.trees = num.trees,
                          solver = eff_solver, tol = tol)
    K <- NULL
  } else {
    K <- leaf_node_kernel(leaf_mat)
    bal <- kernel_balance(trt = A_pred, kern = K, solver = "direct")
    Z <- NULL
  }

  list(weights = bal$weights, forest = forest,
       solver = bal$solver, kernel = K)
}


#' Fit group-specific outcome models and predict for target observations
#' @return List with mu1 and mu0 prediction vectors
#' @noRd
.fit_outcome_models <- function(X_train, A_train, Y_train,
                                X_pred, num.trees) {
  idx_t <- which(A_train == 1)
  idx_c <- which(A_train == 0)

  mu1_forest <- grf::regression_forest(
    X_train[idx_t, , drop = FALSE], Y_train[idx_t], num.trees = num.trees
  )
  mu0_forest <- grf::regression_forest(
    X_train[idx_c, , drop = FALSE], Y_train[idx_c], num.trees = num.trees
  )

  list(
    mu1 = predict(mu1_forest, newdata = X_pred)$predictions,
    mu0 = predict(mu0_forest, newdata = X_pred)$predictions
  )
}


#' Compute ATE from weights, outcomes, and optional outcome model predictions
#' @noRd
.compute_ate <- function(Y, A, w, augmented, mu1 = NULL, mu0 = NULL) {
  if (augmented) {
    # Doubly-robust: regression term + weighted bias correction
    reg_term <- mean(mu1 - mu0)
    corr_1 <- weighted.mean(Y[A == 1] - mu1[A == 1], w[A == 1])
    corr_0 <- weighted.mean(Y[A == 0] - mu0[A == 0], w[A == 0])
    reg_term + corr_1 - corr_0
  } else {
    # Hajek weighted estimator
    weighted.mean(Y[A == 1], w[A == 1]) -
      weighted.mean(Y[A == 0], w[A == 0])
  }
}


#' Process a single cross-fitting fold
#' @return List with fold-level results (ate, weights, mu1, mu0, etc.)
#' @noRd
.fit_one_fold <- function(k, fold_ids, X, A, Y, num.trees, min.node.size,
                          augmented, mu.hat, scale.outcomes, solver, tol, ...) {
  idx_k    <- which(fold_ids == k)
  idx_notk <- which(fold_ids != k)

  A_k <- A[idx_k]; Y_k <- Y[idx_k]

  # Skip if a treatment group is empty in this fold
  if (sum(A_k == 1) == 0 || sum(A_k == 0) == 0) {
    return(list(idx = idx_k, ate = NA, weights = rep(NA, length(idx_k)),
                mu1 = NULL, mu0 = NULL, forest = NULL, solver = NULL))
  }

  # Fit kernel and compute weights
  kb <- .fit_kernel_and_balance(
    X_train = X[idx_notk, , drop = FALSE],
    A_train = A[idx_notk], Y_train = Y[idx_notk],
    X_pred = X[idx_k, , drop = FALSE], A_pred = A_k,
    num.trees = num.trees, min.node.size = min.node.size,
    scale.outcomes = scale.outcomes, solver = solver, tol = tol, ...
  )

  # Outcome model predictions for augmentation
  mu1_k <- NULL; mu0_k <- NULL
  if (augmented) {
    if (is.null(mu.hat)) {
      mu_k <- .fit_outcome_models(
        X_train = X[idx_notk, , drop = FALSE],
        A_train = A[idx_notk], Y_train = Y[idx_notk],
        X_pred = X[idx_k, , drop = FALSE], num.trees = num.trees
      )
      mu1_k <- mu_k$mu1
      mu0_k <- mu_k$mu0
    } else {
      mu1_k <- mu.hat$mu1[idx_k]
      mu0_k <- mu.hat$mu0[idx_k]
    }
  }

  ate_k <- .compute_ate(Y_k, A_k, kb$weights, augmented, mu1_k, mu0_k)

  list(idx = idx_k, ate = ate_k, weights = kb$weights,
       mu1 = mu1_k, mu0 = mu0_k,
       forest = kb$forest, solver = kb$solver)
}


#' Cross-fitted estimation path
#' @noRd
.fit_crossfitted <- function(X, A, Y, num.trees, min.node.size, num.folds,
                             augmented, mu.hat, scale.outcomes, solver,
                             tol, parallel, ...) {
  n <- nrow(X)
  fold_ids <- sample(rep(seq_len(num.folds), length.out = n))

  # Determine number of cores
  if (isTRUE(parallel)) {
    ncores <- parallel::detectCores()
  } else if (is.numeric(parallel) && parallel > 1) {
    ncores <- as.integer(parallel)
  } else {
    ncores <- 1L
  }

  # Run folds sequentially or in parallel
  fold_args <- list(fold_ids = fold_ids, X = X, A = A, Y = Y,
                    num.trees = num.trees, min.node.size = min.node.size,
                    augmented = augmented, mu.hat = mu.hat,
                    scale.outcomes = scale.outcomes, solver = solver,
                    tol = tol, ...)

  run_fold <- function(k) {
    do.call(.fit_one_fold, c(list(k = k), fold_args))
  }

  if (ncores > 1L) {
    fold_results <- parallel::mclapply(seq_len(num.folds), run_fold,
                                       mc.cores = ncores)
  } else {
    fold_results <- lapply(seq_len(num.folds), run_fold)
  }

  # Assemble results from fold outputs
  fold_ates <- numeric(num.folds)
  weights   <- numeric(n)
  mu1_hat   <- if (augmented) numeric(n) else NULL
  mu0_hat   <- if (augmented) numeric(n) else NULL
  last_forest <- NULL
  last_solver <- NULL

  for (k in seq_len(num.folds)) {
    res_k <- fold_results[[k]]
    fold_ates[k] <- res_k$ate
    weights[res_k$idx] <- res_k$weights
    if (augmented) {
      mu1_hat[res_k$idx] <- res_k$mu1
      mu0_hat[res_k$idx] <- res_k$mu0
    }
    if (!is.null(res_k$forest)) {
      last_forest <- res_k$forest
      last_solver <- res_k$solver
    }
  }

  # Assign user-supplied mu.hat if provided
  if (augmented && !is.null(mu.hat)) {
    mu1_hat <- mu.hat$mu1
    mu0_hat <- mu.hat$mu0
  }

  list(
    ate       = mean(fold_ates, na.rm = TRUE),
    weights   = weights,
    mu1.hat   = mu1_hat,
    mu0.hat   = mu0_hat,
    fold_ates = fold_ates,
    fold_ids  = fold_ids,
    kernel    = NULL,
    forest    = last_forest,
    solver    = last_solver,
    num.folds = num.folds
  )
}


#' Full-sample (no cross-fitting) estimation path
#' @noRd
.fit_full_sample <- function(X, A, Y, num.trees, min.node.size,
                             augmented, mu.hat, scale.outcomes, solver,
                             tol, ...) {
  n <- nrow(X)

  # Fit kernel and compute weights
  kb <- .fit_kernel_and_balance(
    X_train = X, A_train = A, Y_train = Y,
    X_pred = X, A_pred = A,
    num.trees = num.trees, min.node.size = min.node.size,
    scale.outcomes = scale.outcomes, solver = solver, tol = tol, ...
  )

  # Outcome model predictions for augmentation
  mu1_hat <- NULL; mu0_hat <- NULL
  if (augmented) {
    if (is.null(mu.hat)) {
      mu_preds <- .fit_outcome_models(X, A, Y, X, num.trees)
      mu1_hat <- mu_preds$mu1
      mu0_hat <- mu_preds$mu0
    } else {
      mu1_hat <- mu.hat$mu1
      mu0_hat <- mu.hat$mu0
    }
  }

  ate <- .compute_ate(Y, A, kb$weights, augmented, mu1_hat, mu0_hat)

  list(
    ate     = ate,
    weights = kb$weights,
    mu1.hat = mu1_hat,
    mu0.hat = mu0_hat,
    kernel  = kb$kernel,
    forest  = kb$forest,
    solver  = kb$solver
  )
}
