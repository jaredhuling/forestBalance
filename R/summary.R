#' Print a forest_balance object
#'
#' Displays a concise summary of the forest balance fit, including the ATE
#' estimate, sample sizes, effective sample sizes, and a brief covariate
#' balance overview.
#'
#' @param x A \code{forest_balance} object returned by
#'   \code{\link{forest_balance}}.
#' @param ... Ignored.
#'
#' @export
print.forest_balance <- function(x, ...) {
  n  <- x$n
  n1 <- x$n1
  n0 <- x$n0
  w  <- x$weights
  is_cf <- isTRUE(x$crossfit)

  if (is_cf) {
    cat("Forest Kernel Energy Balancing (cross-fitted)\n")
  } else {
    cat("Forest Kernel Energy Balancing\n")
  }
  cat(strrep("-", 50), "\n")
  cat(sprintf("  n = %s  (n_treated = %d, n_control = %d)\n",
              format(n, big.mark = ","), n1, n0))
  cat(sprintf("  Trees: %d\n", x$forest[["_num_trees"]]))
  if (is_cf) {
    cat(sprintf("  Cross-fitting: %d folds\n", x$num.folds))
  }
  if (isTRUE(x$augmented)) {
    cat("  Augmented: yes (doubly-robust)\n")
  }
  cat(sprintf("  Solver: %s\n", x$solver))
  cat(sprintf("  ATE estimate: %.4f\n", x$ate))

  if (is_cf) {
    cat(sprintf("  Fold ATEs: %s\n",
                paste(round(x$fold_ates, 4), collapse = ", ")))
  }

  # ESS (pooled across folds for cross-fitted)
  idx_t <- which(x$A == 1)
  idx_c <- which(x$A == 0)
  ess_t <- (sum(w[idx_t]))^2 / sum(w[idx_t]^2)
  ess_c <- (sum(w[idx_c]))^2 / sum(w[idx_c]^2)
  cat(sprintf("  ESS: treated = %.0f/%-d (%.0f%%)   control = %.0f/%-d (%.0f%%)\n",
              ess_t, n1, 100 * ess_t / n1,
              ess_c, n0, 100 * ess_c / n0))

  cat(strrep("-", 50), "\n")
  cat("Use summary() for covariate balance details.\n")
  invisible(x)
}


#' Summarize a forest_balance object
#'
#' Produces a detailed summary of the forest balance fit, including the ATE
#' estimate, covariate balance diagnostics (SMD, ESS, energy distance), and
#' kernel sparsity information.
#'
#' @param object A \code{forest_balance} object returned by
#'   \code{\link{forest_balance}}.
#' @param X.trans An optional matrix of nonlinear or transformed covariates
#'   (\eqn{n \times q}) on which to additionally assess balance (e.g.,
#'   interactions, squared terms).
#' @param threshold SMD threshold for flagging imbalanced covariates. Default is
#'   0.1.
#' @param energy.dist Logical; if \code{TRUE} (default), compute the weighted
#'   energy distance (skipped if \eqn{n > 5000}).
#' @param ... Ignored.
#'
#' @return Invisibly returns a list of class \code{"summary.forest_balance"}
#'   containing:
#' \describe{
#'   \item{ate}{The estimated ATE.}
#'   \item{balance_weighted}{A \code{forest_balance_diag} object for the
#'     weighted sample.}
#'   \item{balance_unweighted}{A \code{forest_balance_diag} object for the
#'     unweighted sample.}
#'   \item{kernel_sparsity}{Fraction of nonzero entries in the kernel matrix.}
#'   \item{n}{Total sample size.}
#'   \item{n1}{Number of treated.}
#'   \item{n0}{Number of control.}
#'   \item{num.trees}{Number of trees in the forest.}
#' }
#'
#' @examples
#' \donttest{
#' n <- 500; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
#' Y <- X[, 1] + rnorm(n)
#'
#' fit <- forest_balance(X, A, Y)
#' summary(fit)
#'
#' # With nonlinear balance assessment
#' X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
#' colnames(X.nl) <- c("X1^2", "X1*X2")
#' summary(fit, X.trans = X.nl)
#' }
#'
#' @export
summary.forest_balance <- function(object, X.trans = NULL, threshold = 0.1,
                                   energy.dist = TRUE, ...) {
  X <- object$X
  A <- object$A
  Y <- object$Y
  w <- object$weights
  n <- object$n

  # Unweighted balance (uniform weights)
  w_unif <- rep(1, n)
  bal_unwt <- compute_balance(X, A, w_unif, X.trans = X.trans,
                              energy.dist = energy.dist)

  # Weighted balance
  bal_wt <- compute_balance(X, A, w, X.trans = X.trans,
                            energy.dist = energy.dist)

  # Kernel sparsity (NULL when CG solver was used)
  K <- object$kernel
  if (is.null(K)) {
    kern_sparsity <- NA_real_
  } else if (inherits(K, "sparseMatrix")) {
    kern_sparsity <- length(K@x) / (as.numeric(n)^2)
  } else {
    kern_sparsity <- sum(K != 0) / (as.numeric(n)^2)
  }

  out <- list(
    ate                = object$ate,
    balance_weighted   = bal_wt,
    balance_unweighted = bal_unwt,
    kernel_sparsity    = kern_sparsity,
    n                  = n,
    n1                 = object$n1,
    n0                 = object$n0,
    num.trees          = object$forest[["_num_trees"]],
    threshold          = threshold
  )
  class(out) <- "summary.forest_balance"
  out
}


#' @rdname summary.forest_balance
#' @param x A \code{summary.forest_balance} object.
#' @export
print.summary.forest_balance <- function(x, ...) {
  threshold <- x$threshold
  bw <- x$balance_weighted
  bu <- x$balance_unweighted

  cat("Forest Kernel Energy Balancing\n")
  cat(strrep("=", 60), "\n")
  cat(sprintf("  n = %s  (n_treated = %d, n_control = %d)\n",
              format(x$n, big.mark = ","), x$n1, x$n0))
  cat(sprintf("  Trees: %d\n", x$num.trees))
  if (!is.na(x$kernel_sparsity)) {
    cat(sprintf("  Kernel density: %.1f%% nonzero\n",
                x$kernel_sparsity * 100))
  } else {
    cat("  Solver: CG (kernel not materialized)\n")
  }
  cat(sprintf("\n  ATE estimate: %.4f\n", x$ate))
  cat(strrep("=", 60), "\n\n")

  # --- Balance comparison table (covariates) ---
  cat("Covariate Balance (|SMD|)\n")
  cat(strrep("-", 60), "\n")

  smd_names <- names(bw$smd)
  smd_tbl <- data.frame(
    Covariate  = smd_names,
    Unweighted = sprintf("%.4f", bu$smd[smd_names]),
    Weighted   = sprintf("%.4f", bw$smd[smd_names]),
    stringsAsFactors = FALSE
  )
  # Flag imbalanced
  flag_uw <- ifelse(bu$smd[smd_names] > threshold, " *", "")
  flag_w  <- ifelse(bw$smd[smd_names] > threshold, " *", "")
  smd_tbl$Unweighted <- paste0(smd_tbl$Unweighted, flag_uw)
  smd_tbl$Weighted   <- paste0(smd_tbl$Weighted, flag_w)

  # Print table
  max_name_len <- max(nchar(smd_names), 10)
  fmt <- paste0("  %-", max_name_len, "s  %12s  %12s\n")
  cat(sprintf(fmt, "Covariate", "Unweighted", "Weighted"))
  cat(sprintf(fmt, strrep("-", max_name_len), strrep("-", 12), strrep("-", 12)))
  for (i in seq_len(nrow(smd_tbl))) {
    cat(sprintf(fmt, smd_tbl$Covariate[i], smd_tbl$Unweighted[i],
                smd_tbl$Weighted[i]))
  }
  cat(sprintf(fmt, strrep("-", max_name_len), strrep("-", 12), strrep("-", 12)))
  cat(sprintf(fmt, "Max |SMD|",
              sprintf("%.4f", bu$max_smd), sprintf("%.4f", bw$max_smd)))
  cat(sprintf("  (* indicates |SMD| > %.2f)\n", threshold))

  # --- Transformed covariates (if present) ---
  if (!is.null(bw$smd_trans)) {
    cat(sprintf("\nTransformed Covariate Balance (|SMD|)\n"))
    cat(strrep("-", 60), "\n")

    trans_names <- names(bw$smd_trans)
    max_tname_len <- max(nchar(trans_names), 10)
    fmt_t <- paste0("  %-", max_tname_len, "s  %12s  %12s\n")
    cat(sprintf(fmt_t, "Transform", "Unweighted", "Weighted"))
    cat(sprintf(fmt_t, strrep("-", max_tname_len), strrep("-", 12),
                strrep("-", 12)))
    for (nm in trans_names) {
      flag_uw_t <- ifelse(bu$smd_trans[nm] > threshold, " *", "")
      flag_w_t  <- ifelse(bw$smd_trans[nm] > threshold, " *", "")
      cat(sprintf(fmt_t, nm,
                  paste0(sprintf("%.4f", bu$smd_trans[nm]), flag_uw_t),
                  paste0(sprintf("%.4f", bw$smd_trans[nm]), flag_w_t)))
    }
    cat(sprintf(fmt_t, strrep("-", max_tname_len), strrep("-", 12),
                strrep("-", 12)))
    cat(sprintf(fmt_t, "Max |SMD|",
                sprintf("%.4f", bu$max_smd_trans),
                sprintf("%.4f", bw$max_smd_trans)))
  }

  # --- ESS ---
  cat("\nEffective Sample Size\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Treated: %.0f / %d  (%.0f%%)\n",
              bw$ess_treated * bw$n1, bw$n1, bw$ess_treated * 100))
  cat(sprintf("  Control: %.0f / %d  (%.0f%%)\n",
              bw$ess_control * bw$n0, bw$n0, bw$ess_control * 100))

  # --- Energy distance ---
  if (!is.na(bw$energy_dist)) {
    cat("\nEnergy Distance\n")
    cat(strrep("-", 60), "\n")
    cat(sprintf("  Unweighted: %.4f\n", bu$energy_dist))
    cat(sprintf("  Weighted:   %.4f\n", bw$energy_dist))
  }

  cat(strrep("=", 60), "\n")
  invisible(x)
}
