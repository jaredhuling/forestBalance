#' Compute random forest proximity kernel from a leaf node matrix
#'
#' Given a matrix of leaf node assignments (observations x trees), computes the
#' \eqn{n \times n} kernel matrix where entry \eqn{(i,j)} is the proportion of
#' trees in which observations \eqn{i} and \eqn{j} fall into the same leaf node.
#'
#' @param leaf_matrix An integer matrix of dimension \eqn{n \times B}, where
#'   \code{leaf_matrix[i, b]} is the leaf node ID for observation \code{i} in
#'   tree \code{b}. Typically produced by \code{\link{get_leaf_node_matrix}}.
#' @param sparse Logical; if \code{TRUE} (default), returns a sparse
#'   \code{\linkS4class{dgCMatrix}} from the \pkg{Matrix} package. If
#'   \code{FALSE}, returns a dense base R matrix.
#'
#' @return A symmetric \eqn{n \times n} matrix (sparse or dense depending on
#'   \code{sparse}) where entry \eqn{(i,j)} equals
#'   \deqn{K(i,j) = \frac{1}{B} \sum_{b=1}^{B} \mathbf{1}(L_b(i) = L_b(j)),}
#'   i.e., the fraction of trees where \eqn{i} and \eqn{j} share a leaf.
#'
#' @details
#' For each tree \eqn{b}, the leaf assignment defines a sparse \eqn{n \times
#' L_b} indicator matrix \eqn{Z_b}. This function stacks all \eqn{B} indicator
#' matrices column-wise into a single sparse matrix
#' \eqn{Z = [Z_1 | Z_2 | \cdots | Z_B]} of dimension
#' \eqn{n \times \sum_b L_b}. The kernel is then obtained via a single sparse
#' cross-product: \eqn{K = Z Z^\top / B}. Leaf ID remapping is done in C++ for
#' speed.
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
#' leaf_mat <- get_leaf_node_matrix(forest)
#' K <- leaf_node_kernel(leaf_mat)               # sparse
#' K_dense <- leaf_node_kernel(leaf_mat, sparse = FALSE)  # dense
#' }
#'
#' @importFrom Matrix sparseMatrix tcrossprod
#' @export
leaf_node_kernel <- function(leaf_matrix, sparse = TRUE) {
  n <- nrow(leaf_matrix)
  B <- ncol(leaf_matrix)

  # Remap leaf IDs to consecutive integers with per-tree offsets (C++)
  remapped <- remap_leaves_cpp(leaf_matrix)
  total_cols <- max(remapped)

  # Build sparse indicator matrix Z: n rows, total_cols columns, n*B nonzeros
  row_idx <- rep(seq_len(n), times = B)
  col_idx <- as.integer(remapped)

  Z <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1,
    dims = c(n, total_cols)
  )

  # Single sparse tcrossprod: K = Z Z^T / B
  K <- Matrix::tcrossprod(Z) / B

  if (!sparse) {
    K <- as.matrix(K)
  }

  K
}
