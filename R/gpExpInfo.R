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
#'   The matrix is of the form
#'   \deqn{\mathbf{I} = \frac{1}{\sigma (1 + 2 \xi)} \begin{bmatrix}
#'   \frac{1}{\sigma} & \frac{1}{1 + \xi} \\
#'   \frac{1}{1 + \xi} & \frac{2}{1 + \xi}
#'   \end{bmatrix}}
#'
#' @returns `gpExpInfo` returns a 2 by 2 numeric matrix with row and column
#'   named `scale, shape`.
#'
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
