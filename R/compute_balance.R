#' Compute covariate balance diagnostics for a set of weights
#'
#' Computes standardized mean differences (SMD), effective sample sizes (ESS),
#' and optionally the weighted energy distance for a given set of balancing
#' weights. Can also assess balance on user-supplied nonlinear transformations
#' of the covariates.
#'
#' @param X An \eqn{n \times p} numeric covariate matrix.
#' @param trt A binary (0/1) vector of treatment assignments of length \eqn{n}.
#' @param weights A numeric weight vector of length \eqn{n}. Treated weights
#'   should sum to \eqn{n_1} and control weights to \eqn{n_0} (as returned by
#'   \code{\link{kernel_balance}} or \code{\link{forest_balance}}).
#' @param X.trans An optional matrix of nonlinear or transformed covariates
#'   (\eqn{n \times q}) on which to additionally assess balance (e.g.,
#'   interactions, squared terms). If \code{NULL} (default), only linear
#'   covariate balance is reported.
#' @param energy.dist Logical; if \code{TRUE}, compute the weighted energy
#'   distance between the treated and control groups. This requires computing an
#'   \eqn{n \times n} distance matrix and is only feasible for moderate \eqn{n}
#'   (automatically skipped when \eqn{n > 5000}). Default is \code{TRUE}.
#'
#' @return An object of class \code{"forest_balance_diag"} containing:
#' \describe{
#'   \item{smd}{Named vector of |SMD| for each covariate.}
#'   \item{max_smd}{Maximum |SMD| across covariates.}
#'   \item{smd_trans}{Named vector of |SMD| for transformed covariates (if
#'     \code{X.trans} was supplied), otherwise \code{NULL}.}
#'   \item{max_smd_trans}{Maximum |SMD| for transformed covariates, or
#'     \code{NA}.}
#'   \item{energy_dist}{Weighted energy distance, or \code{NA} if not
#'     computed.}
#'   \item{ess_treated}{Effective sample size for the treated group as a
#'     fraction of \eqn{n_1}.}
#'   \item{ess_control}{Effective sample size for the control group as a
#'     fraction of \eqn{n_0}.}
#'   \item{n}{Total sample size.}
#'   \item{n1}{Number of treated units.}
#'   \item{n0}{Number of control units.}
#' }
#'
#' @details
#' The standardized mean difference for covariate \eqn{j} is defined as
#' \deqn{\text{SMD}_j = \frac{|\bar X_{j,1}^w - \bar X_{j,0}^w|}{s_j},}
#' where \eqn{\bar X_{j,a}^w} is the weighted mean of covariate \eqn{j} in
#' group \eqn{a} and \eqn{s_j} is the pooled (unweighted) standard deviation.
#'
#' The effective sample size for a group is
#' \deqn{\text{ESS} = \frac{(\sum_i w_i)^2}{\sum_i w_i^2},}
#' reported as a fraction of the group size.
#'
#' The weighted energy distance is
#' \deqn{E = 2 \sum_{i,j} p_i q_j \|X_i - X_j\|
#'         - \sum_{i,j} p_i p_j \|X_i - X_j\|
#'         - \sum_{i,j} q_i q_j \|X_i - X_j\|,}
#' where \eqn{p} and \eqn{q} are the normalized treated and control weights.
#'
#' @examples
#' \donttest{
#' n <- 500; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
#' Y <- X[, 1] + rnorm(n)
#'
#' fit <- forest_balance(X, A, Y)
#' bal <- compute_balance(X, A, fit$weights)
#' bal
#'
#' # With nonlinear features
#' X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
#' colnames(X.nl) <- c("X1^2", "X1*X2")
#' bal2 <- compute_balance(X, A, fit$weights, X.trans = X.nl)
#' bal2
#' }
#'
#' @importFrom stats dist sd quantile
#' @export
compute_balance <- function(X, trt, weights, X.trans = NULL,
                            energy.dist = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  n1 <- as.numeric(sum(trt == 1))
  n0 <- as.numeric(n - n1)
  idx_t <- which(trt == 1L)
  idx_c <- which(trt == 0L)

  # Normalize weights within each group so they sum to group size
  wt_norm <- weights[idx_t] / mean(weights[idx_t])
  wc_norm <- weights[idx_c] / mean(weights[idx_c])

  # --- Linear covariate SMD ---
  pooled_sd <- apply(X, 2L, sd)
  mean_t <- colSums(X[idx_t, , drop = FALSE] * wt_norm) / sum(wt_norm)
  mean_c <- colSums(X[idx_c, , drop = FALSE] * wc_norm) / sum(wc_norm)
  smd_vec <- abs(mean_t - mean_c) / pmax(pooled_sd, 1e-10)

  if (!is.null(colnames(X))) {
    names(smd_vec) <- colnames(X)
  } else {
    names(smd_vec) <- paste0("X", seq_len(p))
  }
  max_smd <- max(smd_vec)

  # --- Transformed covariate SMD (optional) ---
  smd_trans_vec <- NULL
  max_smd_trans <- NA_real_
  if (!is.null(X.trans)) {
    X.trans <- as.matrix(X.trans)
    if (nrow(X.trans) != n) {
      stop("X.trans must have the same number of rows as X.")
    }
    sd_trans <- apply(X.trans, 2L, sd)
    mt_trans <- colSums(X.trans[idx_t, , drop = FALSE] * wt_norm) / sum(wt_norm)
    mc_trans <- colSums(X.trans[idx_c, , drop = FALSE] * wc_norm) / sum(wc_norm)
    smd_trans_vec <- abs(mt_trans - mc_trans) / pmax(sd_trans, 1e-10)
    if (!is.null(colnames(X.trans))) {
      names(smd_trans_vec) <- colnames(X.trans)
    } else {
      names(smd_trans_vec) <- paste0("T", seq_len(ncol(X.trans)))
    }
    max_smd_trans <- max(smd_trans_vec)
  }

  # --- Effective sample sizes ---
  ess_t <- (sum(weights[idx_t]))^2 / sum(weights[idx_t]^2) / n1
  ess_c <- (sum(weights[idx_c]))^2 / sum(weights[idx_c]^2) / n0

  # --- Energy distance (optional, expensive) ---
  ed <- NA_real_
  if (energy.dist && n <= 5000L) {
    wtp <- weights[idx_t] / sum(weights[idx_t])
    wcp <- weights[idx_c] / sum(weights[idx_c])
    D   <- as.matrix(dist(X))
    ed  <- 2 * sum(outer(wtp, wcp) * D[idx_t, idx_c]) -
               sum(outer(wtp, wtp) * D[idx_t, idx_t]) -
               sum(outer(wcp, wcp) * D[idx_c, idx_c])
  }

  out <- list(
    smd           = smd_vec,
    max_smd       = max_smd,
    smd_trans     = smd_trans_vec,
    max_smd_trans = max_smd_trans,
    energy_dist   = ed,
    ess_treated   = ess_t,
    ess_control   = ess_c,
    n             = n,
    n1            = as.integer(n1),
    n0            = as.integer(n0)
  )
  class(out) <- "forest_balance_diag"
  out
}


#' @rdname compute_balance
#' @param x A \code{forest_balance_diag} object.
#' @param threshold SMD threshold for flagging imbalanced covariates. Default is
#'   0.1, a standard threshold in the causal inference literature.
#' @param ... Ignored.
#' @export
print.forest_balance_diag <- function(x, threshold = 0.1, ...) {
  cat("Covariate Balance Diagnostics\n")
  cat(sprintf("  n = %s  (n_treated = %d, n_control = %d)\n",
              format(x$n, big.mark = ","), x$n1, x$n0))
  cat(strrep("-", 60), "\n")

  # ESS
  cat(sprintf("  ESS (treated): %5.1f%%\n", x$ess_treated * 100))
  cat(sprintf("  ESS (control): %5.1f%%\n", x$ess_control * 100))

  # Energy distance
  if (!is.na(x$energy_dist)) {
    cat(sprintf("  Energy distance: %.4f\n", x$energy_dist))
  }
  cat("\n")

  # Helper for printing an SMD block
  print_smd_block <- function(smd_vec, label) {
    n_feat  <- length(smd_vec)
    n_imbal <- sum(smd_vec > threshold)
    qq <- quantile(smd_vec, probs = c(0.5, 0.75, 0.9, 1.0))
    cat(sprintf("  %s (%d features):\n", label, n_feat))
    cat(sprintf("    median = %.4f   Q75 = %.4f   Q90 = %.4f   max = %.4f\n",
                qq[1], qq[2], qq[3], qq[4]))
    if (n_imbal == 0L) {
      cat(sprintf("    All |SMD| below %.2f\n", threshold))
    } else {
      cat(sprintf("    %d/%d above %.2f:\n", n_imbal, n_feat, threshold))
      imbal_idx <- order(smd_vec, decreasing = TRUE)
      imbal_idx <- imbal_idx[smd_vec[imbal_idx] > threshold]
      for (i in imbal_idx) {
        cat(sprintf("      %-25s %.4f\n", names(smd_vec)[i], smd_vec[i]))
      }
    }
  }

  print_smd_block(x$smd, "|SMD| for covariates")

  if (!is.null(x$smd_trans)) {
    cat("\n")
    print_smd_block(x$smd_trans, "|SMD| for transformed covariates")
  }

  cat(strrep("-", 60), "\n")
  invisible(x)
}
