#' GP Distribution Expected Information
#'
#' Calculates the expected information matrix for the GP distribution.
#'
#' @param scale,shape Numeric vectors. Respective values of the GP parameters
#'   scale parameter \eqn{\sigma} and shape parameter \eqn{\xi}. For
#'   `gpExpInfo`, `scale` and `shape` must have length 1.
#' @details `gpExpInfo` calculates, for a single pair of values
#'   \eqn{(\sigma, \xi) = } `(scale, shape)`, the expected information matrix for a
#'   single observation from a GP distribution with distribution function
#'   \deqn{F(x) = P(X \leq x) = 1 - \left[ 1 + \frac{\xi x}{\sigma} \right]_+^{-1/\xi}}
#'
#'   The other functions are vectorized and calculate the individual
#'   contributions to the expected information matrix. For example, `gev11e`
#'   calculates the expectation \eqn{i_{\mu\mu}} of the negated second
#'   derivative of the GEV log-density with respect to \eqn{\mu}, that is, each
#'   `1` indicates one derivative with respect to \eqn{\mu}. Similarly, `2`
#'   denotes one derivative with respect to \eqn{\sigma} and `3` one derivative
#'   with respect to \eqn{\xi}, so that, for example, `gev23e` calculates the
#'   expectation \eqn{i_{\sigma\xi}} of the negated GEV log-density after one
#'   taking one derivative with respect to \eqn{\sigma} and one derivative with
#'   respect to \eqn{\xi}. Note that \eqn{i_{\xi\xi}}, calculated using
#'   `gev33e`, depends only on \eqn{\xi}.
#'
#'   The expectation in `gev11e` can be calculated in a direct way for all
#'   \eqn{\xi > -0.5}. For the other components, direct calculation of the
#'   expectation is unstable when \eqn{\xi} is close to 0. Instead, we use
#'   a quadratic approximation over `(-eps, eps)`, from a Lagrangian
#'   interpolation of the values from the direct calculation for \eqn{\xi = }
#'   `-eps`, \eqn{0} and `eps`.
#' @returns `gevExpInfo` returns a 3 by 3 numeric matrix with row and column
#'   named `loc, scale, shape`. The other functions return a numeric vector of
#'   length equal to the maximum of the lengths of the arguments, excluding
#'   `eps`.
#' @examples
#' # Expected information matrices for ...
#' # ... scale = 1 and shape = 0.1
#' gpExpInfo(1, 0.1)
#' @export
gpExpInfo <- function(scale, shape) {
  if (shape <= -0.5) {
    stop("The GP expected information is undefined for shape <= -0.5")
  }
  val <- matrix(NA, 2, 2)
  val[1, 1] <- 1 / ((1 + 2 * shape) * scale ^ 2)
  val[2, 2] <- 2 / ((1 + shape) * (1 + 2 * shape))
  val[2, 1] <- val[1, 2] <- 1 / (scale * (1 + shape) * (1 + 2 * shape))
  dimnames(val) <- list(c("scale", "shape"), c("scale", "shape"))
  return(val)
}
