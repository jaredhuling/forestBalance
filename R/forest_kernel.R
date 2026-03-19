#' Compute random forest proximity kernel from a GRF forest
#'
#' Extracts the leaf node assignments for each observation across all trees
#' in a trained GRF forest, then computes the \eqn{n \times n} proximity kernel
#' matrix where entry \eqn{(i,j)} is the proportion of trees in which
#' observations \eqn{i} and \eqn{j} share a leaf node.
#'
#' @param forest A trained forest object from the \pkg{grf} package.
#' @param newdata A numeric matrix of observations. If \code{NULL} (default),
#'   uses the original training data.
#'
#' @return A symmetric numeric matrix of dimension \eqn{n \times n}.
#'
#' @details
#' This is a convenience function that calls \code{\link{get_leaf_node_matrix}}
#' followed by \code{\link{leaf_node_kernel}}. If you need both the leaf matrix
#' and the kernel, it is more efficient to call them separately.
#'
#' @examples
#' \donttest{
#' library(grf)
#' n <- 100
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- cbind(X[, 1] + rnorm(n), X[, 2] + rnorm(n))
#' forest <- multi_regression_forest(X, Y, num.trees = 50)
#'
#' K <- forest_kernel(forest)
#' }
#'
#' @export
forest_kernel <- function(forest, newdata = NULL) {
  leaf_mat <- get_leaf_node_matrix(forest, newdata)
  leaf_node_kernel(leaf_mat)
}
