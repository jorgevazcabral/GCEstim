#' Change the step from \code{\link{lmgce}} object
#'
#' Changes the number of GCE reestimations of a \code{\link{lmgce}} object
#'
#' @param object fitted \code{\link{lmgce}} object.
#' @param twosteps.n An integer that defines the number of GCE reestimations to
#' be used.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @return An \code{\link{lmgce}} object with the specified number of GCE
#' reestimations
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         twosteps.n = 10,
#'         seed = 230676)
#'
#' res_gce_package_change_step <- changestep(res_gce_package, 5)
#'
#' summary(res_gce_package)
#'
#' summary(res_gce_package_change_step)
#' }
#'
#' @export

changestep <- function(object, twosteps.n, verbose = 0)
{
  if (is.null(object$results$twosteps) | object$twosteps.n != twosteps.n) {

    tochange <- c(
      "coefficients",
      "residuals",
      "fitted.values",
      "nep",
      "nepk",
      "vcov",
      "error.measure",
      "p",
      "w",
      "lambda",
      "convergence",
      "p0",
      "w0",
      "error.measure.cv.mean",
      "nep.cv.mean",
      "error.measure.cv.sd",
      "nep.cv.sd")

    suppressWarnings(
      aux.object <-
        update(
          object,
          support.signal = object$support.matrix,
          twosteps.n = twosteps.n,
          verbose = verbose
        )
    )
    object[tochange] <- aux.object[tochange]
    object$results$twosteps <-  aux.object$results$twosteps
    object$results$bootstrap <- aux.object$results$bootstrap
  }

  return(object)

}
