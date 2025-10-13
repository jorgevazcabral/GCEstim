# Generic functions for the \code{\link{ridgetrace}} class

#' Print a \code{\link{ridgetrace}} object
#'
#' Concise summary of a \code{\link{ridgetrace}} object
#'
#' @param x fitted \code{\link{ridgetrace}} object.
#' @param digits  significant digits in printout.
#' @param ... additional print arguments.
#'
#' @return A small summary of a \code{\link{ridgetrace}} object is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res.ridgetrace <-
#'   ridgetrace(
#'     formula = y ~ X001 + X002 + X003 + X004 + X005,
#'     data = dataGCE)
#'
#' res.ridgetrace
#'
#' @method print ridgetrace
#' @export

print.ridgetrace <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
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

#' Extract \code{\link{ridgetrace}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{ridgetrace}} object
#'
#' @param object Fitted \code{\link{ridgetrace}} model object.
#' @param which One of \code{c("min.error", "max.abs")}.  If
#' \code{which = "min.error"}, the default, the coefficients that produced the
#' lowest error cross-validation error (\code{cv = TRUE}),or in sample error are
#'  returned (\code{cv = FALSE}). If \code{which = "max.abs"} then the maximum
#'  absolute coefficients are returned.
#' @param ... Additional arguments (not used).
#'
#' @return Returns the coefficients from a \code{ridgetrace} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res.ridgetrace <-
#'   ridgetrace(
#'     formula = y ~ X001 + X002 + X003 + X004 + X005,
#'     data = dataGCE)
#'
#' coef(res.ridgetrace)
#'
#' @method coef ridgetrace
#' @importFrom stats coef
#' @export

coef.ridgetrace <- function(object, which = "min.error",...)
{
  if (which == "max.abs") {
    object$max.abs.coef} else if (which == "min.error") {
      if (is.null(object$error.lambda.cv)) {
        object$coef.lambda[, which.min(object$error.lambda)]
      } else {
        object$coef.lambda[, which.min(apply(object$error.lambda.cv, 2, mean))[[1]]]
      }
    }
}

#' Extract \code{\link{ridgetrace}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{ridgetrace}} object
#'
#' @param object Fitted \code{\link{ridgetrace}} model object.
#' @param which One of \code{c("min.error", "max.abs")}.  If
#' \code{which = "min.error"}, the default, the coefficients that produced the
#' lowest error cross-validation error (\code{cv = TRUE}),or in sample error are
#'  returned (\code{cv = FALSE}). If \code{which = "max.abs"} then the maximum
#'  absolute coefficients are returned.
#' @param ... Additional arguments.
#'
#' @return Returns the coefficients from a \code{ridgetrace} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname coefficients.ridgetrace
#' @examples
#' res.ridgetrace <-
#'   ridgetrace(
#'     formula = y ~ X001 + X002 + X003 + X004 + X005,
#'     data = dataGCE)
#'
#' coefficients(res.ridgetrace)
#'
#' @method coefficients ridgetrace
#' @importFrom stats coefficients
#' @export

coefficients.ridgetrace <- coef.ridgetrace

#' Plot Diagnostics for a \code{\link{ridgetrace}} Object
#'
#' @param x Fitted \code{\link{ridgetrace}} model object.
#' @param coef A vector of true coefficients if available.
#' @param ... additional arguments.
#'
#' @return Supports are returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' res.ridgetrace <-
#'   ridgetrace(
#'     formula = y ~ X001 + X002 + X003 + X004 + X005,
#'     data = dataGCE)
#'
#' plot(res.ridgetrace)
#'
#' @method plot ridgetrace
#' @importFrom rlang .data
#' @export

plot.ridgetrace <- function(x, coef = NULL, ...){

  if (!inherits(x, "ridgetrace"))
    stop("use only with \"ridgetrace\" objects")

  col.coef.all <- viridis::turbo(length(x$max.abs.coef))

  p1.data <-
    data.frame(variable = row.names(x$coef.lambda),
               lambda = rep(log(x$lambda, base = x$lambda[2] / x$lambda[1]),
                            each = nrow(x$coef.lambda)),
               coefficient = as.vector(x$coef.lambda))

  p1 <-
    ggplot2::ggplot(p1.data,
      ggplot2::aes(x = .data$lambda,
                   y = .data$coefficient,
                   group = .data$variable,
                   colour = .data$variable)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = col.coef.all) +
    ggplot2::theme_minimal() +
    ggplot2::ylab("Coefficients") +
    ggplot2::xlab("log(lambda)") +
    ggplot2::ggtitle("RIDGE TRACE")

  if (!is.null(coef)) {
    p1 <- p1 +
      ggplot2::geom_hline(yintercept = coef,
                          linetype = "dashed",
                          colour = col.coef.all)
  }

  p1
}
