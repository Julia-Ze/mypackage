#' Storm Peak Significant Wave Heights in the Northern North Sea
#'
#' Data used for extreme value regression modelling in Chapter 5 of Dey and
#' Yan (2016). The raw data are hindcast time series for significant wave
#' height \eqn{H_S}, dominant wave direction \eqn{d} and seasonal degree
#' day \eqn{s} for three-hour sea states, at an undisclosed location in the
#' Northern North Sea, for the period September 1957 to December 2011.
#'
#' @format This data frame contains the following numeric variables:
#'
#' * `year`: year in which the storm peak occurred.
#' * `month`: month in which the storm peak occurred.
#' * `day`: day in which the storm peak occurred.
#' * `Hs`: storm peak significant wave height, in metres.
#' * `direction`: dominant wave direction, a measure of the direction from
#'   which the storm is coming, in degrees from North, increasing clockwise.
#' * `season`: day of the year, for a standardised year of 360 days.
#' * `NAO`: the NOAA's daily North Atlantic Oscillation index for the day on
#'    which the storm peak occurred.
#'
#' @references Dey, D.K. and Yan, J. (Eds.). (2016). Extreme Value Modeling and
#'   Risk Analysis: Methods and Applications (1st ed.). Chapman and Hall/CRC.
#'   \doi{10.1201/b19721}
#' @source [Computer Code Supplement](http://merlot.stat.uconn.edu/~jyan/docs/evrisk_r/evrisk.html)
#'   to Chapter 5 of Dey and Yan (2016). The daily NAO data were sourced from the
#'   [NOAA](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml)
#'   (National Oceanic and Atmospheric Administration).
#' @examples
#' # Summaries of data
#' head(waves)
#' tail(waves)
#' summary(waves)
#'
#' # Plots of Hs vs potential explanatory variables
#' plot(waves$season, waves$Hs)
#' plot(waves$NAO, waves$Hs)
#' plot(waves$direction, waves$Hs)
"waves"
