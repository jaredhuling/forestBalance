#' Estimate ATE using forest-based kernel energy balancing
#'
#' Fits a multivariate random forest that jointly models the relationship
#' between covariates, treatment, and outcome, computes a random forest
#' proximity kernel, and then uses kernel energy balancing to produce weights
#' for estimating the average treatment effect (ATE).
#'
#' @param X A numeric matrix or data frame of covariates (\eqn{n \times p}).
#' @param A A binary (0/1) vector of treatment assignments.
#' @param Y A numeric vector of outcomes.
#' @param num.trees Number of trees to grow in the forest. Default is 1000.
#' @param min.node.size Minimum number of observations per leaf node. Default
#'   is 10.
#' @param scale.outcomes If \code{TRUE} (default), the joint outcome matrix
#'   \code{cbind(A, Y)} is column-standardized before fitting the forest. This
#'   ensures that treatment and outcome contribute equally to the splits.
#' @param ... Additional arguments passed to
#'   \code{\link[grf]{multi_regression_forest}}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{ate}{The estimated average treatment effect (Hajek estimator).}
#'   \item{weights}{The balancing weights from \code{\link{kernel_balance}}.}
#'   \item{kernel}{The \eqn{n \times n} forest proximity kernel matrix.}
#'   \item{forest}{The trained \code{\link[grf]{multi_regression_forest}} object.}
#' }
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
#' @references
#' De, B. and Huling, J. (2024). Kernel Energy Balancing with Random Forest
#' Similarity. \emph{arXiv preprint arXiv:2512.18069}.
#'
#' @examples
#' \donttest{
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
#' Y <- X[, 1] + rnorm(n)  # true ATE = 0
#'
#' result <- forest_balance(X, A, Y)
#' result$ate
#' }
#'
#' @importFrom grf multi_regression_forest
#' @export
forest_balance <- function(X, A, Y,
                           num.trees = 1000,
                           min.node.size = 10,
                           scale.outcomes = TRUE,
                           ...) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (length(A) != n || length(Y) != n) {
    stop("X, A, and Y must have the same number of observations.")
  }
  if (!all(A %in% c(0, 1))) {
    stop("Treatment vector A must be binary (0/1).")
  }

  # Joint bivariate response: (treatment, outcome)
  response <- cbind(A, Y)
  if (scale.outcomes) {
    response <- scale(response)
  }

  # Fit multivariate random forest on covariates with joint response
  forest <- grf::multi_regression_forest(
    X, Y = response,
    num.trees = num.trees,
    min.node.size = min.node.size,
    ...
  )

  # Compute the random forest proximity kernel
  K <- forest_kernel(forest, newdata = X)

  # Compute kernel energy balancing weights
  bal <- kernel_balance(trt = A, kern = K)
  w <- bal$weights

  # Hajek (ratio) ATE estimator
  ate <- weighted.mean(Y[A == 1], w[A == 1]) -
         weighted.mean(Y[A == 0], w[A == 0])

  list(
    ate     = ate,
    weights = w,
    kernel  = K,
    forest  = forest
  )
}
