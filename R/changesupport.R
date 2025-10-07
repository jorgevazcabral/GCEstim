#' Change the support from \code{\link{lmgce}} object
#'
#' Changes the support spaces of a \code{\link{lmgce}} object
#'
#' @param object fitted \code{\link{lmgce}} object.
#' @param support One of c("min", "1se", "elbow") or a chosen support from
#' \code{object$support.ok}.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @return An \code{\link{lmgce}} object with the specified support spaces
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' res_gce_package_change <- changesupport(res_gce_package, "min")
#'
#' summary(res_gce_package)
#'
#' summary(res_gce_package_change)
#'
#' }
#'
#' @export

changesupport <- function(object, support, verbose = 0)
{
  if (is.null(object$error.which)) {
    stop(cat('It is not possible to change the support since the lmgce object was
             computed for a specific support'),
         call. = FALSE)
  }

  newsupport <- NULL

  if (!(support %in% c("min", "1se", "elbow"))) {

    if (!(support %in% object$support.ok)) {
      stop(cat('argument `support` must be one of c("min", "1se", "elbow") or ',
               object$support.ok),
           call. = FALSE)
    } else
    {
      object$error.which <- "manual"
      object$support.signal.manual <- as.numeric(support)
      newsupport <- as.character(support)

    }
  } else {
  if (support != object$error.which) {

    object$error.which <- support

    object$support.signal.manual <- NULL

    support <- paste0("support.signal.", support)

    newsupport <- as.character(object[support])
  }
  }

  if (!(is.null(newsupport))) {
    tochange <- c(
      "coefficients",
      "residuals",
      "fitted.values",
      "nep",
      "nepk",
      "vcov",
      "error.measure",
      "support.matrix",
      "p",
      "w",
      "lambda",
      "convergence",
      "p0",
      "w0",
      "support.stdUL",
      "error.measure.cv.mean",
      "nep.cv.mean",
      "error.measure.cv.sd",
      "nep.cv.sd")

    object$results$nocv[tochange[1:15]] <-
      object$results$nocv$support.results[[newsupport]][tochange[1:15]]

      suppressWarnings(
        aux.object <-
          update(
            object,
            support.method = "standardized",
            support.signal = as.matrix(
              object$results$nocv$support.results[[newsupport]]$support.matrix),
            verbose = verbose
            )
      )
        object[tochange[c(1:12,16:19)]] <- aux.object[tochange[c(1:12,16:19)]]
        object$support.stdUL <- as.numeric(newsupport)
        object$results$bootstrap <- aux.object$results$bootstrap
  }

return(object)

}
