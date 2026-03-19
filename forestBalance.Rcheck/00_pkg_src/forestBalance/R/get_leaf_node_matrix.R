#' Extract leaf node membership matrix from a GRF forest
#'
#' For each observation and each tree in a trained GRF forest, determines which
#' leaf node the observation falls into. This function accesses the forest's
#' internal tree structure directly, avoiding the overhead of
#' \code{\link[grf]{get_tree}} and \code{\link[grf]{get_leaf_node}}, and
#' vectorizes the tree traversal across all observations simultaneously.
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
#' \code{_root_nodes}, \code{_send_missing_left}) and traverses each tree in a
#' vectorized fashion. All \eqn{n} observations are moved through the tree
#' simultaneously at each depth level, resulting in \eqn{O(d)} vectorized
#' operations per tree, where \eqn{d} is the tree depth (typically
#' \eqn{\sim \log n}). This is dramatically faster than the per-observation
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
#' @export
get_leaf_node_matrix <- function(forest, newdata = NULL) {
  if (!inherits(forest, "grf")) {
    stop("'forest' must be a trained grf forest object.")
  }

  if (is.null(newdata)) {
    newdata <- forest[["X.orig"]]
  }
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  ntrees <- forest[["_num_trees"]]

  result <- matrix(0L, nrow = n, ncol = ntrees)

  for (tr in seq_len(ntrees)) {
    # Access internal tree structure directly (0-based C++ indices)
    root   <- forest[["_root_nodes"]][[tr]]      # 0-based root node index
    left   <- forest[["_child_nodes"]][[tr]][[1]] # 0-based left child indices
    right  <- forest[["_child_nodes"]][[tr]][[2]] # 0-based right child indices
    svars  <- forest[["_split_vars"]][[tr]]        # 0-based split variable indices
    svals  <- forest[["_split_values"]][[tr]]      # split threshold values
    sml    <- forest[["_send_missing_left"]][[tr]] # send-missing-left flags

    # All observations start at the root (convert to 1-based R index)
    node <- rep(root + 1L, n)

    # Vectorized traversal: move all active observations one level per iteration.
    # Tree depth is O(log n), so this loop runs ~20 iterations at most.
    repeat {
      # Identify observations at internal (non-leaf) nodes.
      # In grf, a leaf node has left child == 0 AND right child == 0.
      at_internal <- (left[node] != 0L) | (right[node] != 0L)
      if (!any(at_internal)) break

      idx <- which(at_internal)

      # Split variable and value for each active observation's current node
      sv   <- svars[node[idx]] + 1L  # convert 0-based var index to 1-based column
      sval <- svals[node[idx]]
      smli <- sml[node[idx]]

      # Feature value for each active observation at its split variable
      x_val <- newdata[cbind(idx, sv)]

      # Decision logic matching grf's Tree.cpp::find_leaf_node:
      #   go left if: (value <= split_value) OR
      #               (value is NA and send_missing_left) OR
      #               (split_value is NA and value is NA)
      na_x <- is.na(x_val)
      na_s <- is.na(sval)
      go_left <- (!na_x & !na_s & x_val <= sval) |
                 (na_x & smli) |
                 (na_x & na_s)

      # Move to child node (convert 0-based child index to 1-based R index)
      node[idx] <- ifelse(go_left, left[node[idx]] + 1L, right[node[idx]] + 1L)
    }

    result[, tr] <- node
  }

  result
}
