#' Generalized Pareto GAMLSS Model for threshold excesses
#'
#' Sets an extreme value threshold and then fits Generalized Additive Models
#' (GAMs) for Scale and Shape to threshold excesses with a GP response
#' distribution, using the function [`gamlss::gamlss()`][`gamlss::gamlss`].
#'
#' @param formula A formula object. See [`gamlss::gamlss`]. In this context,
#'   this could be named `sigma.formula` as it specifies the formula to be used
#'   for the GP scale parameter `sigma`. This is passed to [`gamlss::gamlss`]
#'   as the argument `formula`.
#' @param xi.formula A formula object. Specifies the formula to be used for
#'   the GP shape parameter `xi`. This is This is passed to [`gamlss::gamlss`]
#'   as the argument `sigma.formula` because `xi` is the second GP parameter.
#'   The default setting is that `xi.formula` is the same as `formula`, that
#'   is, the same formula is used for the GP scale `sigma` and shape `xi`.
#'   To force a single value of `xi` to estimated for all observations use
#'   `xi.formula = ~ 1`.
#' @inheritParams gamlss::gamlss
#' @param scoring A character scalar. If `scoring = "fisher"` then the weights
#'   used in the fitting algorithm are based on the expected Fisher
#'   information, that is, a Fisher's scoring algorithm is used.
#'   If `scoring = "quasi"` then these weights are based on the cross products
#'   of the first derivatives of the log-likelihood, leading to a quasi Newton
#'   scoring algorithm.
#' @param sigma.link,xi.link Character scalars to set the respective
#'   link functions for the scale (`sigma`) and shape (`xi`)
#'   parameters. The latter (`xi.link`) is passed to
#'   [`gamlss::gamlss()`][`gamlss::gamlss`] as `sigma.link` because `xi` is the
#'   second GP parameter.
#' @param stepLength A numeric vector of positive values. The initial
#'    values of the step lengths `sigma.step` and `nu.step` passed to
#'   [`gamlss::gamlss.control()`][`gamlss::gamlss.control`] in the first attempt
#'   to fit the model by calling [`gamlss::gamlss()`][`gamlss::gamlss`]. If
#'   `stepLength` has a length that is less than 2 then `stepLength` is
#'   recycled to have length 2.
#' @param stepAttempts A non-negative integer. If the first call to
#'   [`gamlss::gamlss()`][`gamlss::gamlss`] throws an error then we make
#'   `stepAttempts` further attempts to fit the model, each time dividing by 2
#'   the values of `sigma.step` and `nu.step` supplied to
#'   [`gamlss::gamlss.control()`][`gamlss::gamlss.control`]. If
#'   `stepAttempts < 1` then no further attempts are made.
#' @param stepReduce A number greater than 1. The factor by which the step
#'   lengths in `stepLength` are reduced for each extra attempt to fit the
#'   model. The default, `stepReduce = 2` means that the step lengths are
#'   halved for each extra attempt.
#' @param steps A logical scalar. Pass `steps = TRUE` to write to the
#'   console the current value of `stepLength` for each call to
#'   [`gamlss::gamlss()`].
#' @param ... Further arguments passed to [`gamlss::gamlss`], in particular
#'   `method`, which sets the fitting algorithm, with options `RS()`, `CG()`
#'   or `mixed()`. The default, `method = RS()` seems to work well, as does
#'   `method = mixed()`. In contrast, `method = CG()` often requires the step
#'   length to be reduced before convergence is achieved. [`fitGP`] attempts to
#'   do this automatically. See `stepAttempts`. Pass `trace = FALSE`
#'   (to [`gamlss::gamlss.control`]) to avoid writing to the console the global
#'   deviance after each outer iteration of the gamlss fitting algorithm.
#' @inheritParams setThreshold
#' @param quantile.formula A formula object to be passed to [`setThreshold`].
#'   If `quantile.formula` is supplied this this is the formula that is passed
#'   to [`mboost::mboost`] to determine the form of the quantile regression
#'   that sets the threshold. Otherwise, this formula is inferred from
#'   the `formula` argument to `fitThresholdGP`, as described in **Details**.
#'
#' @details The threshold is set using quantile regression. See
#'  [`setThreshold`], which uses the `QuantReg` family of the
#'  [`mboost::mboost`] function. If `quantile.formula` is missing then the
#'  form of the quantile regression model reflects the model formula in
#'  `formula` to be used for the GP GAM. See the description of the argument
#'  `replace`.
#'
#' The Generalised Pareto (GP) GAM is fitted using
#' [`gamlss::gamlss()`][`gamlss::gamlss`] based on the the [`GenPareto`] family.
#'
#' @return Returns a `gamlss` object. See the **Value** section of
#'   [`gamlss::gamlss()`][`gamlss::gamlss`]. The class of the returned object is
#'   `c("gamlssx", "gamlss", "gam", "glm", "lm")`.
#'
#'   The following components are added to the returned object list.
#'
#'   * `data`: the input data `data`.
#'   * `exc_data`: a subset on `data` containing only rows for which the
#'      response exceeds the threshold, that is, the threshold exceedances.
#'   * `u, exc_u`: numeric vectors of threshold values for all the observations
#'     (`u`) and only the threshold exceedances (`exc_u`).
#'   * `GPcall`: the call to `fitThresholdGP`.
#'   * `threshold`: the fitted threshold model object returned from [`mboost::mboost`].
#'   * `p_exc`: the probability of threshold exceedance, equal to `1 - tau`.
#'
#' @seealso [`GenPareto`], [`fitGP`], [`setThreshold`],
#'   [`gamlss.dist::gamlss.family()`][`gamlss.dist::gamlss.family`],
#'   [`gamlss::gamlss()`][`gamlss::gamlss`]
#' @examples
#' # Load gamlss
#' library(gamlss)
#'
#' ## Fit models to the waves data with 1 covariate
#'
#' # Seasonal covariate only
#' fit1 <- fitThresholdGP(Hs ~ pbc(season), data = waves)
#' summary(fit1)
#' # gamlss provides some diagnostic plots
#' plot(fit1)
#'
#' # Plot the data, the threshold and the GP fitted (mean) value
#' plot(waves$season, waves$Hs, xlab = "season", ylab = "Hs")
#' points(waves$season, fit1$u, col = "red", pch = 16)
#' newdata <- data.frame(season = waves$season)
#' sigma <- predict(fit1, newdata = newdata, what = "mu", type = "response")
#' xi <- predict(fit1, newdata = newdata, what = "sigma", type = "response")
#' fitted_median <- sigma * (2 ^ xi - 1) / xi
#' points(waves$season, fit1$u + fitted_median, col = "blue", pch = 16)
#' legend("top", legend = c("threshold u", "u + fitted GP median"), pch = 16,
#'        col = c("red", "blue"))
#'
#' # Plot the fitted values of the GP scale and shape parameters
#' # Recall that, for gamlss, sigma is mu and xi is sigma
#' plot(fit1$data$season, sigma, xlab = "season", ylab = "GP scale")
#' plot(fit1$data$season, xi, xlab = "season", ylab = "GP shape")
#'
#' # Directional covariate only
#' fit2 <- fitThresholdGP(Hs ~ pbc(direction), data = waves)
#' # We get a convergence warning, so start again from the estimates returned
#' fit2 <- fitThresholdGP(Hs ~ pbc(direction), data = waves,
#'                        mu.start = fitted(fit2, what = "mu"),
#'                        sigma.start = fitted(fit2, what = "sigma"))
#' summary(fit2)
#'
#' # Plot the data, the threshold and the GP fitted (mean) value
#' plot(waves$direction, waves$Hs, xlab = "direction", ylab = "Hs")
#' points(waves$direction, fit2$u, col = "red", pch = 16)
#' newdata <- data.frame(direction = waves$direction)
#' sigma <- predict(fit2, newdata = newdata, what = "mu", type = "response")
#' xi <- predict(fit2, newdata = newdata, what = "sigma", type = "response")
#' fitted_median <- sigma * (2 ^ xi - 1) / xi
#' points(waves$direction, fit2$u + fitted_median, col = "blue", pch = 16)
#' legend("top", legend = c("threshold u", "u + fitted GP median"), pch = 16,
#'        col = c("red", "blue"))
#'
#' # Plot the fitted values of the GP scale and shape parameters
#' # Recall that, for gamlss, sigma is mu and xi is sigma
#' plot(fit2$data$direction, sigma, xlab = "direction", ylab = "GP scale")
#' plot(fit2$data$direction, xi, xlab = "direction", ylab = "GP shape")
#'
#' ## Fit models to the waves data with 2 covariates
#'
#' # Main effects for Seasonal and directional covariates only
#'
#' # Note: when using this model we assume that the effect of season is the
#' # same for different values of direction and vice versa
#' fit3 <- fitThresholdGP(Hs ~ pbc(direction) + pbc(season), data = waves)
#' summary(fit3)
#'
#' directions <- 0:360
#'
#' # Estimated effect of direction in January
#' mid_january <- data.frame(direction = directions, season = 0.5 / 12)
#' xi_january <- predict(fit3, newdata = mid_january, what = "sigma",
#'                       type = "response")
#' plot(directions, xi_january, type = "l", xlab = "direction", ylab = "xi",
#'      main = "January", lwd = 2)
#'
#' # Estimated effect of direction in July
#' mid_july <- data.frame(direction = directions, season = 6.5 / 12)
#' xi_july <- predict(fit3, newdata = mid_july, what = "sigma",
#'                    type = "response")
#' plot(directions, xi_july, type = "l", xlab = "direction", ylab = "xi",
#'      main = "July", lwd = 2)
#'
#' # The rgl package is required for the following plot
#' # This is just for fun and in case it helps to see that direction has no
#' # influence on the estimated effect of season and vice versa
#' rgl::plot3d(waves$season, waves$direction, fit3$u)
#'
#' # How to calculate the threshold for a storm that arrives mid-January
#' # from direction 250
#' newdata <- data.frame(season = 0.5 / 12, direction = 250)
#' predict(fit3$threshold, newdata = newdata)
#'
#' # How to calculate the fitted values of the GP scale and shape parameters
#' # type = "response" accounts for the (log) link for the scale parameter
#' predict(fit3, what = "mu", newdata = newdata, type = "response")
#' predict(fit3, what = "sigma", newdata = newdata)
#'
#' # Load mgcv and gamlss.add because we need them for what follows
#' library(mgcv)
#' library(gamlss.add)
#'
#' # See ?ga and ?te
#' # More details in Section 9.5 of the gamlss book
#' # 2D effects of direction and season, constrained so that the effects of
#' # both are cyclical, that is, they are constrained on the boundary of a
#' # rectangle.
#' # We use quantile.formula to set the threshold using an equivalent approach
#' # We set tau = 0.5 because this is what Northrop, Jonathan and Randell (2016)
#' # did.
#' fit4 <- fitThresholdGP(Hs ~ ga(~te(direction, season, bs = "cc"),
#'                        knots=list(direction = c(0, 360), season = c(0, 1))),
#'                        data = waves, tau = 0.5,
#'                        quantile.formula =
#'                        Hs ~ bbs(direction, season, cyclic = TRUE))
#'
#' # I don't think that the type = "response" argument works here because
#' # some of the contours for "mu" (our "sigma") have negative values
#' vis.gam(getSmo(fit4, what = "mu"), type = "response", plot.type = "contour")
#' # I also don't believe the numbers of the plot for xi
#' vis.gam(getSmo(fit4, what = "sigma"), plot.type = "contour")
#'
#' # These numbers aren't correct either
#' term.plot(fit4)
#' term.plot(fit4, what = "sigma")
#'
#' # These numbers are correct I think
#' # Get fitted values for xi 'by hand'
#' p <- fitted(fit4, what = "mu")
#' summary(p)
#' p <- fitted(fit4, what = "sigma")
#' summary(p)
#' @export
fitThresholdGP <- function(formula, xi.formula = formula, data, tau = 0.75,
                           quantile.formula, mstop = 1000,
                           replace = c("pb", "pbo", "pbp", "pbc", "cy", "pbm",
                                       "pbz", "ps", "pvc", "pvp"),
                           scoring = c("fisher", "quasi"),
                           sigma.link = "log", xi.link = "identity",
                           stepLength = 1, stepAttempts = 2, stepReduce = 2,
                           steps = FALSE,  ...) {

  # Set a covariate-dependent threshold using mboost::mboost()
  # Set a threshold at the estimated 100tau% conditional quantile
  if (missing(quantile.formula)) {
    threshold <- setThreshold(formula = formula, data = data, tau = tau,
                              mstop = mstop, replace = replace)
  } else {
    threshold <- setThreshold(formula = quantile.formula, data = data,
                              tau = tau, mstop = mstop, replace = NULL)
  }
  # Extract just the fitted threshold values
  u <- fitted(threshold)

  # Get the name of the response variable from the formula, so that we can
  # extract the threshold exceedances and excesses from the data
  response_name <- all.vars(formula)[1]
  response_data <- data[, response_name]
  # Calculate the indicators of threshold exceedances and threshold excesses
  exceed_ind <- response_data > u
  excesses <- response_data[exceed_ind] - u[exceed_ind]

  # Save a copy of the input data
  all_data <- data

  # Retain in 'data' only the rows with threshold exceedances
  # Replace the response column with the threshold excesses
  data <- data[exceed_ind, ]
  data[, response_name] <- excesses

  # Check that one of the correct values of scoring has been supplied
  scoring <- match.arg(scoring)
  # Force stepLength to have length 2
  stepLength <- rep_len(stepLength, 2)
  # Set the scoring algorithm and links
  # For all the gamlss methods to work on the returned fitted model object, we
  # need the call to gamlss::gamlss to include explicitly the names of the
  # links for sigma and nu (xi here) as character scalars.
  # To achieve this, we create the internal function templateFit() and
  # modify the body of this function to include the names of the link
  # functions. This is rather clunky, but it works! Other attempts at passing
  # the link functions do not work completely. For example, vcov.gamlss(object)
  # does not work, because it performs calculations using the fitted model
  # object and needs to know the link functions.
  #
  # Fit using the supplied/default value of step length
  if (scoring == "fisher") {
    algor <- substitute(
      GPfisher(mu.link = sigma.link, sigma.link = xi.link)
    )
  } else {
    algor <- substitute(
      GPquasi(mu.link = sigma.link, sigma.link = xi.link)
    )
  }
  # Add the link functions to the call to gamlss() in templateFit()
  templateFit <- function(formula, sigma.formula, stepLength, data, ...) {
    dangerous <- NULL
    return(dangerous)
  }
  body(templateFit)[[2]] <- substitute(
    dangerous <- try(gamlss::gamlss(formula = formula,
                                    sigma.formula = sigma.formula,
                                    family = algor,
                                    mu.step = stepLength[1],
                                    sigma.step = stepLength[2], data = data,
                                    ...),
                     silent = TRUE)
  )
  if (steps) {
    cat("stepLength =", stepLength, "\n")
  }
  mod <- templateFit(formula = formula, sigma.formula = xi.formula,
                     stepLength = stepLength, data = data, ...)
  # If an error is thrown then try again stepAttempts times, each time reducing
  # the step length by a factor of a half
  isError <- inherits(mod, "try-error")
  while(isError & stepAttempts >= 1) {
    stepLength <- stepLength / stepReduce
    if (steps) {
      cat("stepLength =", stepLength, "\n")
    }
    # We need to update the value of stepLength in templateFit()
    body(templateFit)[[2]] <- substitute(
      dangerous <- try(gamlss::gamlss(formula = formula, family = algor,
                                      mu.step = stepLength[1],
                                      sigma.step = stepLength[2], data = data,
                                      ...),
                       silent = TRUE)
    )
    mod <- templateFit(formula = formula, stepLength = stepLength, data = data,
                       ...)
    isError <- inherits(mod, "try-error")
    stepAttempts <- stepAttempts - 1L
  }
  if (isError) {
    print(attr(mod, "condition"))
    stop("No convergence. An error was thrown from the last call to gamlss()")
  }
  # Add the data (raw and thresholded) and the threshold vector (for all
  # observations and just the threshold exceedances) to the returned object
  mod$data <- all_data
  mod$exc_data <- all_data[exceed_ind, ]
  mod$u <- u
  mod$exc_u <- u[exceed_ind]
  # Save the call to fitThresholdGP(), so that it can be printed by
  # print.gamlssx() and summary.gamlssx(), rather than the verbose call to
  # gamlss::gamlss()
  mod$GPcall <- match.call()
  mod$threshold <- threshold
  mod$p_exc <- 1- tau
  class(mod) <- c("gamlssx", class(mod))
  return(mod)
}
