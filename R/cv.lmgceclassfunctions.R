# Generic functions for the \code{\link{cv.lmgce}} class

#' Print \code{\link{cv.lmgce}} object
#'
#' Print \code{\link{cv.lmgce}} object
#'
#' @param x fitted \code{\link{cv.lmgce}} object.
#' @param digits  significant digits in printout.
#' @param ... additional print arguments.
#'
#' @return A small summary of a \code{cv.lmgce} object is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res.cv.lmgce <-
#'   cv.lmgce(y ~ .,
#'            data = dataGCE)
#'
#' res.cv.lmgce
#' }
#'
#' @method print cv.lmgce
#' @export

print.cv.lmgce <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")

  cat("\nBest combination:\n",
      paste(x$support.signal.points.best,
            "points for the signal support;",
            x$support.noise.points.best,
            "points for the noise support;",
            "a weight of",
            x$weight.best,
            "for the noise support."
            ),
      "\n\n",
      sep = "")

  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2L,
                  quote = FALSE)
  }
  else
    cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' Extract \code{\link{cv.lmgce}} Coefficients
#'
#' Extract coefficients from a \code{\link{cv.lmgce}} object
#'
#' @param object Fitted \code{\link{cv.lmgce}} model object.
#' @param ... Additional arguments (not used).
#'
#' @return Returns the coefficients from a \code{\link{cv.lmgce}} object. The
#' coefficients are obtained from the \code{\link{lmgce}} object with best
#' performance. These coefficients are stored in \code{object$best$coefficients}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res.cv.lmgce <-
#'   cv.lmgce(y ~ .,
#'            data = dataGCE)
#'
#' coef(res.cv.lmgce)
#' }
#'
#' @method coef cv.lmgce
#' @importFrom stats coef
#' @export

coef.cv.lmgce <- function(object, ...)
{
  object$best$coefficients
}

#' Extract \code{\link{cv.lmgce}} Coefficients
#'
#' Extract coefficients from a \code{\link{cv.lmgce}} object
#'
#' @param object Fitted \code{\link{cv.lmgce}} model object.
#' @param ... Additional arguments (not used).
#'
#' @return Returns the coefficients from a \code{\link{cv.lmgce}} object. The
#' coefficients are obtained from the \code{\link{lmgce}} object with best
#' performance. These coefficients are stored in \code{object$best$coefficients}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname coefficients.cv.lmgce
#' @examples
#' \donttest{
#' res.cv.lmgce <-
#'   cv.lmgce(y ~ .,
#'            data = dataGCE)
#'
#' coefficients(res.cv.lmgce)
#' }
#'
#' @method coefficients cv.lmgce
#' @importFrom stats coefficients
#' @export

coefficients.cv.lmgce <- coef.cv.lmgce

#' Plot Diagnostics for a \code{\link{cv.lmgce}} Object
#'
#' One plot (selectable by \code{which}) is currently available to
#' evaluate a \code{\link{cv.lmgce}} object. The plot depicts the error change
#' with the combination of different arguments of \code{\link{cv.lmgce}}.
#'
#' @param x Fitted \code{cv.lmgce} model object.
#' @param which A subset of the numbers 1:1.
#' @param ncol Number of columns of the plot (see
#' \code{\link[ggplot2]{facet_wrap}}).
#' @param scales One of c("free", "fixed") (see
#' \code{\link[ggplot2]{facet_wrap}}).
#' @param ... additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @seealso \code{\link[GCEstim]{cv.lmgce}}
#'
#' @examples
#' \donttest{
#' res.cv.lmgce <-
#'   cv.lmgce(y ~ .,
#'            data = dataGCE)
#'
#' plot(res.cv.lmgce)
#' }
#'
#' @method plot cv.lmgce
#' @importFrom rlang .data
#' @export

plot.cv.lmgce <-
  function (x,
            which = 1,
            ncol = 1,
            scales = "free",
            ...) {
    if (!inherits(x, "cv.lmgce"))
      stop("use only with \"cv.lmgce\" objects")

    show <- rep(FALSE, 1)
    show[which] <- TRUE

    if (show[1L]) {
    ggplot2::ggplot(x$results,
                    ggplot2::aes(x = .data$support.signal.points,
                                 y = .data$error.measure.cv.mean,
                                 group = as.factor(.data$support.noise.points),
                                 colour = as.factor(.data$support.noise.points))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~weight, ncol = ncol, scales = scales) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlab("Number of points of the signal support") +
      ggplot2::ylab(paste0("CV-", x$best$error)) +
      ggplot2::guides(
        colour =
          ggplot2::guide_legend(title = "Number of points of the noise support"))

    }
  }
