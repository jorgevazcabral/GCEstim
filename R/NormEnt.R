#' Normalized Entropy
#'
#' Returns the normalized entropy of the model or the normalized entropy of the
#' predictors.
#'
#' @param object fitted \code{\link{lmgce}} or \code{\link{tsbootgce}} object.
#' @param model Boolean value. If \code{model = TRUE}, the model's normalized
#' entropy is returned. If \code{model = FALSE} the normalized entropy of each
#' estimate is returned. The default is \code{model = TRUE}.
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing,
#' all parameters are considered.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @return the value of the normalized entropy of the model or parameters.
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#'
#' NormEnt(res_gce_package)
#'
#' @export

NormEnt <- function(object, model = TRUE, parm)
{
  if (!inherits(object, "lmgce") & !inherits(object, "tsbootgce"))
    stop("use only with \"lmgce\" or \"tsbootgce\" objects")

  if (isTRUE(model)) {
    object$nep
  } else {
    pnames <- names(object$nepk)
    if (missing(parm))
      parm <- pnames
    else if (is.numeric(parm))
      parm <- pnames[parm]
    if (!all(parm %in% pnames)) {
      stop("Invalid choice of parameter")
    }
    object$nepk[parm]
  }
}
