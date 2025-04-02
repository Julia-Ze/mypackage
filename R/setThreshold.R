#' Set an extreme value threshold
#'
#' Performs flexible quantile regression to set an extreme value threshold as
#' a user-supplied estimated quantile level of the response conditional on
#' the values of covariates.
#'
#' @inheritParams gamlss::gamlss
#'
#' @param tau A numeric scalar \eqn{\tau} between 0 and 1 passed as the
#'   argument `tau` to the [`Quantreg`][mboost::Family]. The threshold is set
#'   at an estimate of the \eqn{100\tau\%} conditional quantile of the response
#'   given the values of the covariates.
#' @param mstop The argument `mstop` to [`mboost::boost_control`].
#' @param replace A character vector. All instances of `pb()`, for example, in
#'   `formula` will be replaced by `bbs()` for use in setting the threshold
#'   using [`mboost::mboost`]. The default value of `replace`
#'   contains the P-spline fitting functions in [`gamlss::ps`].
#'
#'   The special cases of the cyclic spline functions `pbc(.)` and `cy(.)` are
#'   replaced by `bbs(., cyclic = TRUE)` and `bbs(., cyclic = TRUE)`
#'   respectively.
#' @param ... Additional arguments to be passed to [`mboost::boost_control`].
#' @details Uses the `QuantReg` family of the [`mboost::mboost`] function,
#'   which fits additive models using gradient boosting.
#' @return An object of class `mboost` returned by [`mboost::mboost`].
#'
#' To obtain a vector of the value of the threshold use `fitted(object)`, where
#' `object` is the name of the returned object.
#' @examples
#' ## A toy example
#'
#' Dat <- NULL
#' Dat$x <- rep(1:25, 20)
#' set.seed(1)
#' Dat$y <- SSlogis(Dat$x, 10, 12, 2) * rnorm(500, 1, 0.1)
#' # Non-cyclic spline function
#' fit <- setThreshold(y ~ pb(x), data = Dat)
#' plot(fit)
#' # Cyclic spline function
#' fit <- setThreshold(y ~ pbc(x), data = Dat)
#' plot(fit)
#'
#' ## The northern North Sea significant wave height data
#'
#' # NAO only
#' fit <- setThreshold(Hs ~ pb(NAO), data = waves)
#' plot(fit)
#' # The threshold superimposed on the data
#' threshold <- fitted(fit)
#' plot(waves$NAO, waves$Hs)
#' points(waves$NAO, threshold, col = "red", pch = 16)
#'
#' # Season only
#' fit <- setThreshold(Hs ~ pbc(season), data = waves)
#' threshold <- fitted(fit)
#' plot(waves$season, waves$Hs)
#' points(waves$season, threshold, col = "red", pch = 16)
#'
#' # Direction only
#' fit <- setThreshold(Hs ~ pbc(direction), data = waves)
#' threshold <- fitted(fit)
#' plot(waves$direction, waves$Hs)
#' points(waves$direction, threshold, col = "red", pch = 16)
#'
#' # Season and direction
#' fit <- setThreshold(Hs ~ pbc(season) + pbc(direction), data = waves)
#' par(mfrow = c(2, 1))
#' plot(fit)
#'
#' # NAO, season and direction
#' fit <- setThreshold(Hs ~ pb(NAO) + pbc(season) + pbc(direction), data = waves)
#' par(mfrow = c(2, 2))
#' plot(fit)
#'
#' @importFrom mboost bbs
#' @export
setThreshold <- function(formula, data, tau = 0.9, mstop = 1000,
                         replace = c("pb", "pbo", "pbp", "pbc", "cy", "pbm",
                                     "pbz", "ps", "pvc", "pvp"), ...) {
  # A function to replace pb() etc with bbs()
  bbsFormulaFun <- function(replace, formula) {
    # Replace the elements (function names) in replace
    # The cyclic formula functions pbc() and cy() need special treatment
    # because cy(x) and pbc(x) need to become bbs(x, cyclic = TRUE)
    if (replace == "pbc" || replace == "cy") {
      # Add the cyclic = TRUE argument(s)
      formula <- cy_and_pbc(formula, which = replace)
    } else {
      # Turn the formula into a character string
      formula <- deparse(formula)
    }
    # Replace all functions in replace with bbs()
    formula <- gsub(paste0(replace, "\\("), "bbs\\(", formula)
    # Convert back to a formula
    formula <- as.formula(formula)
    return(formula)
  }
  bbsFormula <- formula
  if (length(replace) > 0) {
    for (i in 1:length(replace)) {
      bbsFormula <- bbsFormulaFun(replace = replace[i], formula = bbsFormula)
    }
  }
  fit <- mboost::mboost(formula = bbsFormula, data = data,
                        family = mboost::QuantReg(tau = tau),
                        baselearner = "bbs",
                        control = mboost::boost_control(mstop = mstop, ...))
  return(fit)
}
