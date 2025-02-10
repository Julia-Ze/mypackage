#' Fit a flexible (?) threshold ...?
#'
#' Fits a flexible threshold to the provided data to extract exceedences
#' for further modelling using quantile regression. using boosting?
#' mention tau? using the function
#' [`mboost::mboost()`][`mboost::mboost`].
#'
#' @inheritParams mboost::mboost
#'
#' @param tau quantile level, level above which we would like the exceedences
#'
#' @param mstop an integer giving the number of initial boosting iterations.
#' If mstop = 0, the offset model is returned. (definition from boost_control)
#'
#' @return Returns a vector of fitted values which represent the threshold
#' for a given value of the covariate
#'
#' @details See [`mboost::mboost()`][`mboost::mboost`] for information about
#' the model and the fitting algorithm.
#'
#' @seealso what to put here?
#'
#' @examples
#'
#' set.seed(512024)
#' n <- 100
#' x <- stats::runif(n)
#' sigma <- 1 + 2 * x
#' xi <- 0.1
#' y <- nieve::rGPD2(n = 1, scale = sigma, shape = xi)
#' data <- data.frame(y = as.numeric(y), x = x)
#' plot(x, y)
#' mod_fitted <- setThreshold(y ~ bbs(x), data = data, mstop = 900)
#' points(x, mod_fitted, col = "blue", pch = 17)
#'
#' @import mboost
#' @export
setThreshold <- function(formula, data, tau = 0.9, mstop = 1000, ...) {

  mod <- mboost::mboost(formula = formula, data = data,
                       family = QuantReg(tau = tau),
#                      baselearner = "bbs",
                       control = boost_control(mstop = mstop, ...))
  return(fitted(mod))
}



