#' Simulated data set generated with fngendata
#'
#' Simulated data, used to demonstrate the functions of GCEstim.
#'
#' @format A \code{data.frame} containing:
#' \describe{
#'   \item{X001}{A N(0,1) independent variable.}
#'   \item{X002}{A N(0,1) independent variable.}
#'   \item{X003}{A N(0,1) independent variable.}
#'   \item{X004}{A N(0,1) independent variable.}
#'   \item{X005}{A N(0,1) independent variable.}
#'   \item{X006}{A N(0,1) independent variable.}
#'   \item{y}{A Dependent variable: y = 2.5 - 8 * X004 + 19 * X005 - 13 * X006 + error;
#'   the error follows a normal distribution with mean equal to zero and variance
#'   such that the signal to noise ratio is equal to 1.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(dataincRidGME.test)
#'
#' plot(dataincRidGME.test)

"dataincRidGME.test"
