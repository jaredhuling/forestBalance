#' Extract leaf node membership matrix from a GRF forest
#'
#' For each observation and each tree in a trained GRF forest, determines which
#' leaf node the observation falls into. The tree traversal is implemented in
#' C++ for speed, directly reading the forest's internal tree structure and
#' avoiding the overhead of \code{\link[grf]{get_tree}} and
#' \code{\link[grf]{get_leaf_node}}.
#'
#' @param forest A trained forest object from the \pkg{grf} package
#'   (e.g., from \code{\link[grf]{multi_regression_forest}},
#'   \code{\link[grf]{regression_forest}}, etc.).
#' @param newdata A numeric matrix of observations to predict leaf membership
#'   for. Must have the same number of columns as the training data. If
#'   \code{NULL} (default), uses the original training data stored in the forest.
#'
#' @return An integer matrix of dimension \code{nrow(newdata)} by
#'   \code{num.trees}, where entry \code{[i, b]} is the internal node index
#'   (1-based) of the leaf that observation \code{i} falls into in tree
#'   \code{b}.
#'
#' @details
#' This function reads the internal tree vectors stored in a GRF forest object
#' (\code{_child_nodes}, \code{_split_vars}, \code{_split_values},
#' \code{_root_nodes}, \code{_send_missing_left}) and traverses each tree in
#' compiled C++ code. This is dramatically faster than the per-observation
#' R-level loop in \code{grf::get_leaf_node}.
#'
#' @examples
#' \donttest{
#' library(grf)
#' n <- 200
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- cbind(X[, 1] + rnorm(n), X[, 2] + rnorm(n))
#' forest <- multi_regression_forest(X, Y, num.trees = 100)
#'
#' # Leaf membership for training data
#' leaf_mat <- get_leaf_node_matrix(forest)
#'
#' # Leaf membership for new data
#' X.test <- matrix(rnorm(50 * p), 50, p)
#' leaf_mat_test <- get_leaf_node_matrix(forest, X.test)
#' }
#'
#' @useDynLib forestBalance, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
get_leaf_node_matrix <- function(forest, newdata = NULL) {
  if (!inherits(forest, "grf")) {
    stop("'forest' must be a trained grf forest object.")
  }

  if (is.null(newdata)) {
    newdata <- forest[["X.orig"]]
  }
  newdata <- as.matrix(newdata)

  # Strip X.orig to avoid passing it into C++ (it's large and unneeded)
  forest_short <- forest[grep("^_", names(forest), value = TRUE)]
  forest_short[["_num_trees"]] <- forest[["_num_trees"]]

  get_leaf_nodes_cpp(forest_short, newdata)
}
