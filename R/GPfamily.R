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
#' @param sigma.link Defines the `sigma.link`, with `"log"` link as the default
#' for the `sigma` parameter.
#' @param nu.link Defines the `nu.link`, with `"identity"` link as the default
#' for the `nu` parameter.
#' @param x,q Vector of quantiles.
#' @param sigma,nu Vectors of scale and shape parameter values.
#' @param log,log.p Logical. If `TRUE`, probabilities `eqn{p}` are given as
#'   \eqn{\log(p)}.
#' @param lower.tail Logical. If `TRUE` (the default), probabilities are
#'   \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param p Vector of probabilities.
#' @param n Number of observations. If `length(n) > 1`, the length is taken to
#'   be the number required.
#'
#' @details The distribution function of a GP distribution with parameters
#'  \code{scale} = \eqn{\sigma (> 0)} and \code{shape} = \eqn{\xi} (\eqn{= \nu}) is
#'   \deqn{F(x) = P(X \leq x) = 1 - \left( 1 + \frac{\xi x}{\sigma} \right)^{-1/\xi}.}
#'  If \eqn{\xi = 0} the distribution function is defined as the limit as
#'  \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x \leq \mu - \sigma / \xi}{x <= \mu - \sigma / \xi} for \eqn{\xi < 0};
#'  \eqn{x \geq \mu - \sigma / \xi}{x >= \mu - \sigma / \xi} for \eqn{\xi > 0};
#'  and \eqn{x} is unbounded for \eqn{\xi = 0}.
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
#'  and/or Chapter 3 of Coles (2001) for further information.
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
#'   Chapter 3: \doi{10.1007/978-1-4471-3675-0_3}
#' @section Examples:
#' See the examples in [`fitGP`].
#' @name GP
NULL
## NULL

#' @rdname GP
#' @export
GPfisher <- function(sigma.link = "log",
                      nu.link = "identity") {

  dstats <- gamlss.dist::checklink("sigma.link", "GP", substitute(sigma.link),
                                   c("inverse", "log", "identity"))
  vstats <- gamlss.dist::checklink("nu.link", "GP",substitute(nu.link),
                                   c("inverse", "log", "identity"))

  structure(
    list(family = c("GP", "Generalized Pareto"),
         parameters = list(sigma = TRUE, nu = TRUE),
         nopar = 2,
         type = "Continuous",
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,
         dldd = function(y, sigma, nu) {
           dl <- nieve::dGPD2(x = y, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "scale"]
           return(dldd)
         },
         d2ldd2 = function(y, sigma, nu) {
           dl <- nieve::dGPD2(x = y, scale = sigma, shape = nu,
                              log = TRUE, hessian = TRUE)
           dldd2 <- attr(dl, "hessian")[, "scale", "scale"]
           return(dldd2)
         },
         dldv = function(y, sigma, nu) {
           dl <- nieve::dGPD2(x = y, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldv <- attr(dl, "gradient")[, "shape"]
           return(dldv)
         },
         d2ldv2 = function(y, mu, sigma, nu) {
           dl <- nieve::dGPD2(x = y, scale = sigma, shape = nu,
                              log = TRUE, hessian = TRUE)
           dldd2 <- attr(dl, "hessian")[, "scale", "shape"]
           return(dldv2)
         },
         d2ldddv = function(y, mu, sigma, nu) {
           dl <- nieve::dGPD2(x = y, scale = sigma, shape = nu,
                              log = TRUE, hessian = TRUE)
           dldd2 <- attr(dl, "hessian")[, "shape", "shape"]
           return(dldddv)
         },
         G.dev.incr = function(y, mu, sigma, nu,...) {
           val <- -2 * dGEV(x = y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
           return(val)
         },
         rqres = expression(rqres(pfun = "pGEV", type = "Continuous",
                                  y = y, mu = mu, sigma = sigma, nu = nu)),
         # sqrt(6) / pi is approximately 0.78
         # 0.57722 * sqrt(6) / pi is approximately 0.45
         # The gamlss.dist::RGE() code had a typo in mu.initial + 0.45 should be - 0.45
         # The next (commented out) line starts from a crude stationary fit
         #         mu.initial = expression(mu <- rep(mean(y) - 0.45 * sd(y), length(y))),
         mu.initial = expression(mu <- y - 0.45 * sd(y)),
         sigma.initial = expression(sigma <- rep(0.78 * sd(y), length(y))),
         nu.initial = expression(nu <- rep(0.1, length(y))),
         mu.valid = function(mu) TRUE,
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > -0.5),
         y.valid = function(y) TRUE
    ),
    class = c("gamlss.family","family")
  )
}

#' @rdname GP
#' @export
GEVquasi <- function(mu.link = "identity", sigma.link = "log",
                     nu.link = "identity") {

  mstats <- gamlss.dist::checklink("mu.link", "GEV", substitute(mu.link),
                                   c("1/mu^2", "log", "identity"))
  dstats <- gamlss.dist::checklink("sigma.link", "GEV", substitute(sigma.link),
                                   c("inverse", "log", "identity"))
  vstats <- gamlss.dist::checklink("nu.link", "GEV",substitute(nu.link),
                                   c("inverse", "log", "identity"))

  structure(
    list(family = c("GEV", "Generalized Extreme Value"),
         parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
         nopar = 3,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,
         dldm = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "loc"]
           return(dldm)
         },
         d2ldm2 = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "loc"]
           dldm2 <- -dldm * dldm
           return(dldm2)
         },
         dldd = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "scale"]
           return(dldd)
         },
         d2ldd2 = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "scale"]
           dldd2 <- -dldd * dldd
           return(dldd2)
         },
         dldv = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldv <- attr(dl, "gradient")[, "shape"]
           return(dldv)
         },
         d2ldv2 = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldv <- attr(dl, "gradient")[, "shape"]
           dldv2 <- -dldv * dldv
           return(dldv2)
         },
         d2ldmdd = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "loc"]
           dldd <- attr(dl, "gradient")[, "scale"]
           dldmdd <- -dldm * dldd
           return(dldmdd)
         },
         d2ldmdv = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldm <- attr(dl, "gradient")[, "loc"]
           dldv <- attr(dl, "gradient")[, "shape"]
           dldmdv <- -dldm * dldv
           return(dldmdv)
         },
         d2ldddv = function(y, mu, sigma, nu) {
           dl <- nieve::dGEV(x = y, loc = mu, scale = sigma, shape = nu,
                             log = TRUE, deriv = TRUE)
           dldd <- attr(dl, "gradient")[, "scale"]
           dldv <- attr(dl, "gradient")[, "shape"]
           dldddv <- -dldd * dldv
           return(dldddv)
         },
         G.dev.incr = function(y, mu, sigma, nu,...) {
           val <- -2 * dGEV(x = y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
           return(val)
         },
         rqres = expression(rqres(pfun = "pGEV", type = "Continuous",
                                  y = y, mu = mu, sigma = sigma, nu = nu)),
         mu.initial = expression(mu <- y + 0.45 * sd(y)),
         sigma.initial = expression(sigma <- rep(0.78 * sd(y), length(y))),
         nu.initial = expression(nu <- rep(0.1, length(y))),
         mu.valid = function(mu) TRUE,
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > -0.5),
         y.valid = function(y) TRUE
    ),
    class = c("gamlss.family","family")
  )
}

#' @rdname GP
#' @export
dGP <- function(x, sigma = 1, nu = 0, log = FALSE) {
  return(nieve::dGPD2(x = x, scale = sigma, shape = nu,
                     log = log))
}

#' @rdname GP
#' @export
pGP <- function(q, sigma = 1, nu = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  return(nieve::pGPD2(q = q, scale = sigma, shape = nu))
}

#' @rdname GP
#' @export
qGP <- function(p, sigma = 1, nu = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  return(nieve::qGPD2(p = p, scale = sigma, shape = nu))
}

#' @rdname GP
#' @export
rGP <- function(n, sigma = 1, nu = 0) {
  return(nieve::rGPD2(n = n, scale = sigma, shape = nu))
}
