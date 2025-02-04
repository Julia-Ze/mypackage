setThreshold <- function(formula, data, tau, mstop = 1000) {

  mod <- mboost:mboost(y ~ bbs(x), data, family = QuantReg(tau = 0.9),
                       baselearner = "bbs")
  return(fitted(mod))

}

# Q1: formula is an argument, so the user doesn't need to know to call y ~ bb(x) [the function bb],
# because all of this is done by us behind the scenes so we implement it either way?
#
# Q2: how to connect data to y ~ bb(x), how do we extract the response and covariate, what is the form
# of the object data that we provide?
#
# Q3: import quantreg beforehand, the whole package?
