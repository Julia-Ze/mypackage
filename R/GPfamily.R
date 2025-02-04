#' GP family distribution for fitting a GAMLSS
#'
#' The functions `GPfisher()` and `GPquasi()` each define the generalized
#' Pareto (GP) family distribution, a two parameter distribution, for
#' a [`gamlss.dist::gamlss.family()`][`gamlss.dist::gamlss.family`] object to
#' be used in GAMLSS fitting using the function
#' [`gamlss::gamlss()`][`gamlss::gamlss`]. The only difference
#' between `GPfisher()` and `GPquasi()` is the form of scoring method used to
#' define the weights used in the fitting algorithm. Fisher's scoring,
#' based on the expected Fisher information is used in `GPfisher()`, whereas
#' a quasi-Newton scoring, based on the cross products of the first derivatives
#' of the log-likelihood, is used in `GPquasi()`. The functions
#' `dGP`, `pGP`, `qGP` and `rGP` define the density, distribution function,
#' quantile function and random generation for the specific parameterization of
#' the generalized extreme value distribution given in **Details** below.
#'
#' @param mu.link Defines the `mu.link`, with `"log"` link as the default
#' for the `mu` parameter.
#' @param sigma.link Defines the `sigma.link`, with `"identity"` link as the default
#' for the `sigma` parameter.
#' @param x,q Vector of quantiles.
#' @param mu,sigma Vectors of scale and shape parameter values.
#' @param log,log.p Logical. If `TRUE`, probabilities `eqn{p}` are given as
#'   \eqn{\log(p)}.
#' @param lower.tail Logical. If `TRUE` (the default), probabilities are
#'   \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param p Vector of probabilities.
#' @param n Number of observations. If `length(n) > 1`, the length is taken to
#'   be the number required.
#'
#' @details For the purposes of using the gamlss package, the two parameters need to
#' be called \eqn{\mu} and \eqn{\sigma} for scale and shape respectively.
#' However, in mathematical notation the standard \eqn{\sigma} and \eqn{\xi} notation will be used.
#'
#' The distribution function of a GP distribution with parameters
#'  \code{scale} = \eqn{\sigma (> 0)} and \code{shape} = \eqn{\xi} (\eqn{= \nu}) is
#'   \deqn{F(x) = P(X \leq x) = 1 - \left( 1 + \frac{\xi x}{\sigma} \right)^{-1/\xi}.}
#'  If \eqn{\xi = 0} the Generalized Pareto distribution is equivalent to an
#'  Exponential distribution with parameter \eqn{1/\mu}.
#'
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x \geq 0}{x >= 0} for \eqn{\xi \geq 0};
#'  \eqn{0 \leq x \leq - \sigma / \xi}{0 <= x <= - \sigma / \xi} for \eqn{\xi < 0}{\xi < 0}.
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#'  and/or Chapter 4 of Coles (2001) for further information.
#'
#' For each observation in the data, the restriction that \eqn{\xi > -1/2} is
#' imposed, which is necessary for the usual asymptotic likelihood theory to be
#' applicable.
#'
#' @return `GPfisher()` and `GPquasi()` each return a
#'   [`gamlss.dist::gamlss.family()`][`gamlss.dist::gamlss.family`] object
#'   which can be used to fit a regression model with a GP response
#'   distribution using the
#'   [`gamlss::gamlss()`][`gamlss::gamlss`] function. `dGP()` gives the density,
#'   `pGP()` gives the distribution function, `qGP()` gives the quantile
#'   function, and `rGP()` generates random deviates.
#' @seealso [`fitGP`],
#'   [`gamlss.dist::gamlss.family()`][`gamlss.dist::gamlss.family`],
#'   [`gamlss::gamlss()`][`gamlss::gamlss`]
#' @references Coles, S. G. (2001) *An Introduction to Statistical
#'   Modeling of Extreme Values*, Springer-Verlag, London.
#'   Chapter 4: \doi{10.1007/978-1-4471-3675-0_4}
#' @section Examples:
#' See the examples in [`fitGP`].
#' @name GP
NULL
## NULL

#' @rdname GP
#' @export
GPfisher <- function(mu.link = "log", sigma.link = "identity") {

  mstats <- gamlss.dist::checklink("mu.link", "GP", substitute(mu.link),
                                   c("inverse", "log", "identity"))
  dstats <- gamlss.dist::checklink("sigma.link", "GP",substitute(sigma.link),
                                   c("inverse", "log", "identity"))

  structure(
    list(family = c("GP", "Generalized Pareto"),
         parameters = list(mu = TRUE, sigma = TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         dldm = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "scale"]
           return(dldm)
         },
         d2ldm2 = function(y, mu, sigma) {
           val <-  -gp11e(scale = mu, shape = sigma)
           return(val)
         },
         dldd = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "shape"]
           return(dldd)
         },
         d2ldd2 = function(y, mu, sigma) {
           val <-  -gp22e(scale = mu, shape = sigma)
           return(val)
         },
         d2ldmdd = function(y, mu, sigma) {
           val <-  -gp12e(scale = mu, shape = sigma)
           return(val)
         },
         G.dev.incr = function(y, mu, sigma,...) {
           val <- -2 * dGP(x = y, mu = mu, sigma = sigma, log = TRUE)
           return(val)
         },
         rqres = expression(rqres(pfun = "pGP", type = "Continuous",
                                  y = y, mu = mu, sigma = sigma)),
         mu.initial = expression(mu <- rep(mean(y), length(y))),
         sigma.initial = expression(sigma <- rep(0.1, length(y))),
         mu.valid = function(mu) all(mu > 0),
         sigma.valid = function(sigma) all(sigma > -0.5),
         y.valid = function(y) TRUE
    ),
    class = c("gamlss.family","family")
  )
}

#' @rdname GP
#' @export
GPquasi <- function(mu.link = "log", sigma.link = "identity") {

  mstats <- gamlss.dist::checklink("mu.link", "GP", substitute(mu.link),
                                   c("inverse", "log", "identity"))
  dstats <- gamlss.dist::checklink("sigma.link", "GP",substitute(sigma.link),
                                   c("inverse", "log", "identity"))

  structure(
    list(family = c("GP", "Generalized Pareto"),
         parameters = list(mu = TRUE, sigma = TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         dldm = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "scale"]
           return(dldm)
         },
         d2ldm2 = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "scale"]
           d2ldm2 <- -dldm * dldm
           return(d2ldm2)
         },
         dldd = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "shape"]
           return(dldd)
         },
         d2ldd2 = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "shape"]
           d2ldd2 <- -dldd * dldd
           return(d2ldd2)
         },
         d2ldmdd = function(y, mu, sigma) {
           dl <- nieve::dGPD2(x = y, scale = mu, shape = sigma,
                              log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "scale"]
           dldd <- attr(dl, "gradient")[, "shape"]
           d2ldmdd <- -dldm * dldd
           return(d2ldmdd)
         },
         G.dev.incr = function(y, mu, sigma,...) {
           val <- -2 * dGP(x = y, mu = mu, sigma = sigma, log = TRUE)
           return(val)
         },
         rqres = expression(rqres(pfun = "pGP", type = "Continuous",
                                  y = y, mu = mu, sigma = sigma)),
         mu.initial = expression(mu <- rep(mean(y), length(y))),
         sigma.initial = expression(sigma <- rep(0.1, length(y))),
         mu.valid = function(mu) all(mu > 0),
         sigma.valid = function(sigma) all(sigma > -0.5),
         y.valid = function(y) TRUE
    ),
    class = c("gamlss.family","family")
  )
}

#' @rdname GP
#' @export
dGP <- function(x, mu = 1, sigma = 0, log = FALSE) {
  return(nieve::dGPD2(x = x, scale = mu, shape = sigma, log = log))
}

#' @rdname GP
#' @export
pGP <- function(q, mu = 1, sigma = 0, lower.tail = TRUE, log.p = FALSE) {
  return(nieve::pGPD2(q = q, scale = mu, shape = sigma))
}

#' @rdname GP
#' @export
qGP <- function(p, mu = 1, sigma = 0, lower.tail = TRUE, log.p = FALSE) {
  return(nieve::qGPD2(p = p, scale = mu, shape = sigma))
}

#' @rdname GP
#' @export
rGP <- function(n, mu = 1, sigma = 0) {
  return(nieve::rGPD2(n = n, scale = mu, shape = sigma))
}
