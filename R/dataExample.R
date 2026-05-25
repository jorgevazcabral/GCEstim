#' Simulated data set generated with fngendata
#'
#' Simulated data, used to demonstrate the functions of GCEstim.
#'
#' @format A \code{data.frame} containing:
#' \describe{
#'   \item{X001}{A N(0,1) independent variable.}
#'   \item{X002}{A N(0,1) independent variable.}
#'   \item{y}{A Dependent variable: y = 1 - 6 * X001 + 9 * X002 +
#'   error; the error follows a normal distribution with mean equal to zero and
#'   variance such that the signal to noise ratio is equal to 1; N = 8.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(dataExample)
#'
#' plot(dataExample)

"dataExample"
