#' Methods for objects of class `"gamlssx"`
#'
#' Methods for objects of class `"gamlssx"` returned from [`fitThresholdGP`].
#'
#' @param x,object An object inheriting from class `"gamlssx"`, a result of a
#'   call to [`fitThresholdGP`].
#' @param ... Additional arguments.
#' @details These `print` and `summary` methods have been created to avoid
#'   printing to the console a huge statement of the call to the function
#'   `gamlss::gamlss`. Instead, the call to `fitThresholdGP` is printed.
#' @return See [`gamlss::print.gamlss`] and [`gamlss::summary.gamlss`].
#' @name gamlssx_methods
NULL
## NULL

# ================================ plot.gamlssx ============================= #

#' Plot method for objects of class `"gamlssx"`
#'
#' @rdname gamlssx_methods
#' @export
print.gamlssx <- function(x, ...) {
  # Create a copy of the input object to modify
  for_gamlss <- x
  for_gamlss$call <- x$GPcall
  # Remove the "gamlssx" part of the class to leave the "gamlss" part first
  class(for_gamlss) <- class(for_gamlss)[-1]
  print(for_gamlss, ...)
  return(x)
}

# =============================== summary.gamlssx =========================== #

#' Summary method for objects of class `"gamlssx"`
#'
#' @rdname gamlssx_methods
#' @export
summary.gamlssx <- function(object, ...) {
  # Create a copy of the input object to modify
  for_gamlss <- object
  for_gamlss$call <- object$GPcall
  # Remove the "gamlssx" part of the class to leave the "gamlss" part first
  class(for_gamlss) <- class(for_gamlss)[-1]
  summary(for_gamlss, ...)
  return(invisible())
}
