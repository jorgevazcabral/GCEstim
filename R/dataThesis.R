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
#'   \item{y}{A Dependent variable: y = -4 - 16 * X002 - 12 * X003 - 5 * X004 +
#'   error; the error follows a normal distribution with mean equal to zero and
#'   variance such that the signal to noise ratio is equal to 5; N = 75.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(dataThesis)
#'
#' plot(dataThesis)

"dataThesis"

#' Coefficients used in simulated data set generated with fngendata
#'
#' Coefficients used in simulated data, used to demonstrate the functions of GCEstim.
#'
#' @format A vector containing the coefficients used to generate dataThesis
#'
#' @keywords datasets
#'
#' @examples
#'
#' coef.dataThesis

"coef.dataThesis"
