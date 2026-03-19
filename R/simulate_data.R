#' Simulate observational study data with confounding
#'
#' Generates data from a design with nonlinear confounding, where covariates
#' jointly influence both treatment assignment and the outcome through non-trivial
#' functions. The true average treatment effect is known, allowing evaluation of
#' estimator performance.
#'
#' @param n Sample size. Default is 500.
#' @param p Number of covariates. Must be at least 5. Default is 10.
#' @param ate True average treatment effect. Default is 0.
#' @param rho Correlation parameter for the AR(1) covariance structure among
#'   covariates: \eqn{\Sigma_{jk} = \rho^{|j-k|}}. Default is \eqn{-0.25}.
#' @param sigma Noise standard deviation for the outcome. Default is 1.
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
#' }
#'
#' @details
#' The data generating process is:
#'
#' **Covariates:** \eqn{X \sim N(0, \Sigma)} where
#' \eqn{\Sigma_{jk} = \rho^{|j-k|}}.
#'
#' **Propensity score:** \eqn{P(A = 1 \mid X) = 0.25(1 + B(X_1; 2, 4))}, where
#' \eqn{B(\cdot; 2, 4)} is the Beta(2,4) density. This creates a nonlinear
#' relationship between \eqn{X_1} and treatment assignment.
#'
#' **Outcome model:**
#' \deqn{Y = 2(X_1 - 1) + 2\,B(X_1; 2, 4) + X_2 + 2\,B(X_5; 2, 4)
#'         + \tau \cdot A + \varepsilon,}
#' where \eqn{\varepsilon \sim N(0, \sigma^2)} and \eqn{\tau} is the ATE.
#' Confounding arises because \eqn{X_1} affects both the propensity score and
#' the outcome nonlinearly.
#'
#' @examples
#' dat <- simulate_data(n = 500, p = 10, ate = 0)
#' str(dat)
#'
#' # True ATE is 0, naive estimate is biased
#' mean(dat$Y[dat$A == 1]) - mean(dat$Y[dat$A == 0])
#'
#' @importFrom stats dbeta rbinom rnorm
#' @export
simulate_data <- function(n = 500, p = 10, ate = 0, rho = -0.25, sigma = 1) {
  if (p < 5) stop("p must be at least 5.")

  # AR(1) covariance
  Sig <- rho^abs(outer(1:p, 1:p, FUN = "-"))
  X <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = Sig)
  colnames(X) <- paste0("X", seq_len(p))

  # Propensity score: nonlinear in X1 via Beta(2,4) density
  propensity <- 0.25 * (1 + dbeta(X[, 1], 2, 4))

  # Treatment assignment
  A <- rbinom(n, 1, propensity)

  # Outcome with nonlinear confounding through X1 and X5
  mu <- 2 * (X[, 1] - 1) + 2 * dbeta(X[, 1], 2, 4) +
        X[, 2] + 2 * dbeta(X[, 5], 2, 4)
  Y <- mu + ate * A + rnorm(n, sd = sigma)

  list(
    X          = X,
    A          = A,
    Y          = Y,
    propensity = propensity,
    ate        = ate,
    n          = n,
    p          = p
  )
}
