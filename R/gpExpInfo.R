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
#'   \deqn{F(x) = P(X \leq x) = 1 - \left[ 1 + \displaystyle\frac{\xi x}{\sigma} \right]_+^{-1/\xi}}
#'   where \eqn{x_+ = \max(x, 0)}.
#'   The GP expected information is defined only for \eqn{\xi > -0.5}.
#'
#'   The matrix is of the form
#'   \deqn{\mathbf{I} = \displaystyle\frac{1}{1 + 2 \xi} \begin{bmatrix}
#'   \displaystyle\frac{1}{\sigma^2} & \displaystyle\frac{1}{\sigma (1 + \xi)} \\
#'   \displaystyle\frac{1}{\sigma (1 + \xi)} & \displaystyle\frac{2}{1 + \xi}
#'   \end{bmatrix}.}
#'
#'   The other functions are vectorized and calculate the individual
#'   contributions to the expected information matrix. For example, `gp11e`
#'   calculates the expectation \eqn{i_{\sigma\sigma}} of the negated second
#'   derivative of the GP log-density with respect to \eqn{\sigma}, that is,
#'   each `1` indicates one derivative with respect to \eqn{\sigma}. Similarly,
#'   `2` denotes one derivative with respect to \eqn{\xi}, so that, `gp12e`
#'   calculates the expectation \eqn{i_{\sigma\xi}} of the negated GP
#'   log-density after one taking one derivative with respect to \eqn{\sigma}
#'   and one derivative with respect to \eqn{\xi}.
#'
#' @returns `gpExpInfo` returns a 2 by 2 numeric matrix with row and column
#'   named `scale, shape`. The other functions return a numeric vector of
#'   length equal to the maximum of the lengths of the arguments.
#'
#' @examples
#' # Expected information matrices for ...
#' # ... scale = 1 and shape = 0.1
#' gpExpInfo(1, 0.1)
#'
#' # The individual components of this matrix
#' gp11e(1, 0.1)
#' gp12e(1, 0.1)
#' gp22e(1, 0.1)
#'
#' @name gpExpInfo
NULL
## NULL

#' @rdname gpExpInfo
#' @export
gp11e <- function(scale, shape) {
  m <- max(length(scale), length(shape))
  scale <- rep_len(scale, m)
  shape <- rep_len(shape, m)
  return(1 / ((1 + 2 * shape) * scale ^ 2))
}

#' @rdname gpExpInfo
#' @export
gp22e <- function(scale, shape) {
  m <- max(length(scale), length(shape))
  scale <- rep_len(scale, m)
  shape <- rep_len(shape, m)
  return(2 / ((1 + shape) * (1 + 2 * shape)))
}

#' @rdname gpExpInfo
#' @export
gp12e <- function(scale, shape) {
  m <- max(length(scale), length(shape))
  scale <- rep_len(scale, m)
  shape <- rep_len(shape, m)
  return(1 / (scale * (1 + shape) * (1 + 2 * shape)))
}

#' @rdname gpExpInfo
#' @export
gpExpInfo <- function(scale, shape) {
  if (any(shape <= -0.5)) {
    stop("The GP expected information is undefined for shape <= -0.5")
  }
  val <- matrix(NA, 2, 2)
  val[1, 1] <- gp11e(scale, shape)
  val[2, 2] <- gp22e(scale, shape)
  val[2, 1] <- val[1, 2] <- gp12e(scale, shape)
  dimnames(val) <- list(c("scale", "shape"), c("scale", "shape"))
  return(val)
}
