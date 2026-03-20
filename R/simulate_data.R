#' Simulate observational study data with confounding
#'
#' Generates data from a design with nonlinear confounding, where covariates
#' jointly influence both treatment assignment and the outcome through non-trivial
#' functions. The true average treatment effect is known, allowing evaluation of
#' estimator performance.
#'
#' @param n Sample size. Default is 500.
#' @param p Number of covariates. Must be at least 5 for \code{dgp = 1} and at
#'   least 8 for \code{dgp = 2}. Default is 10.
#' @param ate True average treatment effect. Default is 0.
#' @param rho Correlation parameter for the AR(1) covariance structure among
#'   covariates: \eqn{\Sigma_{jk} = \rho^{|j-k|}}. Default is \eqn{-0.25}.
#' @param sigma Noise standard deviation for the outcome. Default is 1.
#' @param dgp Integer selecting the data generating process. Default is 1.
#'   See Details.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{X}{The \eqn{n \times p} covariate matrix.}
#'   \item{A}{Binary (0/1) treatment assignment vector.}
#'   \item{Y}{Observed outcome vector.}
#'   \item{propensity}{True propensity scores \eqn{P(A=1 \mid X)}.}
#'   \item{ate}{The true ATE used in the simulation.}
#'   \item{n}{Sample size.}
#'   \item{p}{Number of covariates.}
#'   \item{dgp}{The DGP that was used.}
#' }
#'
#' @details
#' Both DGPs generate covariates from \eqn{X \sim N(0, \Sigma)} where
#' \eqn{\Sigma_{jk} = \rho^{|j-k|}}.
#'
#' \strong{DGP 1} (default): Confounding through \eqn{X_1} via a Beta density.
#'
#' \itemize{
#'   \item Propensity: \eqn{P(A=1|X) = 0.25(1 + B(X_1; 2, 4))} where
#'     \eqn{B} is the Beta(2,4) density.
#'   \item Outcome: \eqn{Y = 2(X_1-1) + 2 B(X_1;2,4) + X_2 + 2 B(X_5;2,4)
#'     + \tau A + \varepsilon}.
#' }
#'
#' \strong{DGP 2}: Rich outcome surface with moderate confounding. Designed to
#' illustrate the benefit of the augmented estimator. Confounding operates
#' through \eqn{X_1} and \eqn{X_2}, while the outcome depends on
#' \eqn{X_1, \ldots, X_8} with interactions and nonlinearities.
#'
#' \itemize{
#'   \item Propensity: \eqn{P(A=1|X) = \mathrm{logit}^{-1}(0.6 X_1 - 0.4 X_2
#'     + 0.2 X_1 X_2)}.
#'   \item Outcome: \eqn{Y = 2 X_1 + X_2^2 - 1.5 X_3 + \sin(2 X_4) +
#'     X_5 X_1 + 0.8 X_6 - \cos(X_7) + 0.5 X_8 + \tau A + \varepsilon}.
#' }
#'
#' @examples
#' dat1 <- simulate_data(n = 500, p = 10, dgp = 1)
#' dat2 <- simulate_data(n = 500, p = 20, dgp = 2)
#'
#' @importFrom stats dbeta plogis rbinom rnorm
#' @export
simulate_data <- function(n = 500, p = 10, ate = 0, rho = -0.25, sigma = 1,
                          dgp = 1) {
  if (!dgp %in% c(1, 2)) stop("dgp must be 1 or 2.")
  if (dgp == 1 && p < 5) stop("p must be at least 5 for dgp = 1.")
  if (dgp == 2 && p < 8) stop("p must be at least 8 for dgp = 2.")

  # AR(1) covariance
  Sig <- abs(rho)^abs(outer(1:p, 1:p, FUN = "-"))
  X <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  colnames(X) <- paste0("X", seq_len(p))

  if (dgp == 1) {
    # DGP 1: Beta density confounding through X1
    propensity <- 0.25 * (1 + dbeta(X[, 1], 2, 4))
    mu <- 2 * (X[, 1] - 1) + 2 * dbeta(X[, 1], 2, 4) +
          X[, 2] + 2 * dbeta(X[, 5], 2, 4)

  } else {
    # DGP 2: moderate confounding, rich outcome surface
    propensity <- plogis(0.6 * X[, 1] - 0.4 * X[, 2] + 0.2 * X[, 1] * X[, 2])
    mu <- 2 * X[, 1] + X[, 2]^2 - 1.5 * X[, 3] +
          sin(2 * X[, 4]) + X[, 5] * X[, 1] +
          0.8 * X[, 6] - cos(X[, 7]) + 0.5 * X[, 8]
  }

  A <- rbinom(n, 1, propensity)
  Y <- mu + ate * A + rnorm(n, sd = sigma)

  list(
    X          = X,
    A          = A,
    Y          = Y,
    propensity = propensity,
    ate        = ate,
    n          = n,
    p          = p,
    dgp        = dgp
  )
}
