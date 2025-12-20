#' Simulated data set generated with fngendata
#'
#' Simulated data, used to demonstrate the functions of GCEstim. The seed used
#' is the different from the one used to generate \code{dataGCE} but the
#' remaining parameters are the same.
#'
#' @format A \code{data.frame} containing:
#' \describe{
#'   \item{X001}{A N(0,1) independent variable.}
#'   \item{X002}{A N(0,1) independent variable.}
#'   \item{X003}{A N(0,1) independent variable.}
#'   \item{X004}{A N(0,1) independent variable.}
#'   \item{X005}{A N(0,1) independent variable.}
#'   \item{y}{A Dependent variable: y = 1 + 3 * X003 + 6 * X004 + 9 * X005 +
#'   error; the error follows a normal distribution with mean equal to zero and
#'   variance such that the signal to noise ratio is equal to 5.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(dataGCE.test)
#'
#' plot(dataGCE.test)

"dataGCE.test"
