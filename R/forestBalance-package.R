#' forestBalance: Forest Kernel Energy Balancing for Causal Inference
#'
#' Estimates average treatment effects (ATE) using kernel energy balancing with
#' random forest similarity kernels. A multivariate random forest jointly models
#' covariates, treatment, and outcome to build a proximity kernel that captures
#' confounding structure. Balancing weights are obtained via a closed-form
#' kernel energy distance solution.
#'
#' @section Main function:
#' \code{\link{forest_balance}} is the primary interface. It fits the forest,
#' constructs the kernel, computes balancing weights, and estimates the ATE.
#' By default it uses K-fold cross-fitting and an adaptive leaf size to minimize
#' overfitting bias.
#'
#' @section Key features:
#' \itemize{
#'   \item Adaptive \code{min.node.size} that scales with \eqn{n} and \eqn{p}
#'   \item K-fold cross-fitting to reduce kernel overfitting bias
#'   \item Rcpp-accelerated leaf node extraction
#'   \item Sparse kernel construction via single \code{tcrossprod}
#'   \item Conjugate gradient solver for large \eqn{n} (avoids forming the
#'     kernel matrix entirely)
#' }
#'
#' @section Lower-level interface:
#' For more control, the pipeline can be run step by step:
#' \code{\link{get_leaf_node_matrix}}, \code{\link{leaf_node_kernel}},
#' \code{\link{kernel_balance}}.
#'
#' @references
#' De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for causal
#' effect estimation for high dimensional data.
#' \emph{arXiv preprint arXiv:2512.18069}.
#'
#' @useDynLib forestBalance, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
