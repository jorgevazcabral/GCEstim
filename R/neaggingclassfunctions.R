# Generic functions for the \code{\link{neagging}} class

#' Extract \code{\link{neagging}} Coefficients
#'
#' Extract coefficients from a \code{\link{neagging}} object
#'
#' @param object Fitted \code{\link{neagging}} model object.
#' @param which Number of aggregated models. The coefficients returned are by
#' default the ones that produced the lowest in sample error.
#' @param ... Additional arguments.
#'
#' @return Returns the coefficients from a \code{\link{neagging}} object
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
#' res_neagging <- neagging(res_gce_package)
#' coef(res_neagging)
#' coef(res_neagging, which = ncol(res_neagging$matrix))
#' }
#'
#' @method coef neagging
#' @importFrom stats coef
#' @export

coef.neagging <- function(object, which = which.min(object$error)[[1]], ...)
{
  object$matrix[, which]
}

#' Extract \code{\link{neagging}} Coefficients
#'
#' Extract coefficients from a \code{\link{neagging}} object
#'
#' @param object Fitted \code{\link{neagging}} model object.
#' @param which Number of aggregated models. The coefficients returned are by
#' default the ones that produced the lowest in sample error.
#' @param ... Additional arguments.
#'
#' @return Returns the coefficients from a \code{\link{neagging}} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname coefficients.neagging
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' res_neagging <- neagging(res_gce_package)
#' coefficients(res_neagging)
#' coefficients(res_neagging, which = ncol(res_neagging$matrix))
#' }
#'
#' @method coefficients neagging
#' @importFrom stats coefficients
#' @export

coefficients.neagging <- coef.neagging

#' Plot Diagnostics for a \code{\link{neagging}} Object
#'
#' Two plots (selectable by \code{which}) are currently available to
#' evaluate a \code{\link{neagging}} object: plots of the estimates and in
#' sample error against the number of bootstrap samples aggregated.
#'
#' @param x Fitted \code{neagging} model object.
#' @param which Numbers 1 or 2.
#' @param ... additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @seealso \code{\link[GCEstim]{lmgce}}, \code{\link[GCEstim]{tsbootgce}} and
#' \code{\link[GCEstim]{neagging}}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' res_neagging <- neagging(res_gce_package)
#'
#' plot(res_neagging)
#' }
#'
#' @method plot neagging
#' @importFrom rlang .data
#' @export

plot.neagging <-
  function (x,
            which = 1,
            ...) {
    object <- x
    if (!inherits(object, "neagging"))
      stop("use only with \"neagging\" objects")

    vect.coef <- NULL

    for (i in 1:ncol(object$matrix)) {
      vect.coef <- c(vect.coef, object$matrix[, i])
    }

    dt.coef <-
      as.data.frame(
        cbind(rep(1:ncol(object$matrix), each = nrow(object$matrix)),
              names(vect.coef),
              vect.coef))
    names(dt.coef) <- c("B","predictor","estimate")
    row.names(dt.coef) <- NULL
    dt.coef$estimate <- as.numeric(dt.coef$estimate)
    dt.coef$B <- as.numeric(dt.coef$B)

    col.coef.all <- viridis::turbo(length(object$coefficients))
    if ("(Intercept)" %in% names(object$coefficients))
      col.coef <- col.coef.all[-1] else
        col.coef <- col.coef.all

    if (which == 1) {
      dt.error <- data.frame(B = 1:length(object$error),
                             error = object$error)

      ggplot2::ggplot(data = dt.error,
                      ggplot2::aes(x = .data$B,
                                   y = .data$error)) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = which.min(object$error)[[1]],
                            linetype = "dashed") +
        ggplot2::xlab("Bootstrapp sample number") +
        ggplot2::ylab(paste0("in sample ",
                             names(object$error.object))) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.y = ggplot2::element_text(size = 12,
                                              colour = "black"),
          axis.text.x.bottom = ggplot2::element_text(size = 12,
                                                     colour = "black"),
          axis.title.x = ggplot2::element_text(size = 12,
                                               colour = "black"),
          legend.title = ggplot2::element_blank()) +
        ggplot2::geom_hline(yintercept = object$error.object,
                            linetype = "dotted")
    } else  if (which == 2) {
      ggplot2::ggplot(data = dt.coef,
                      ggplot2::aes(x = .data$B,
                                   y = .data$estimate,
                                   group = .data$predictor,
                                   colour = .data$predictor)) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = which.min(object$error)[[1]],
                            linetype = "dashed") +
        ggplot2::xlab("Bootstrapp sample number") +
        ggplot2::ylab("Estimates") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.y = ggplot2::element_text(size = 12,
                                              colour = "black"),
          axis.text.x.bottom = ggplot2::element_text(size = 12,
                                                     colour = "black"),
          axis.title.x = ggplot2::element_text(size = 12,
                                               colour = "black"),
          legend.title = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(values = col.coef.all) +
        ggplot2::geom_hline(yintercept = object$coefficients.object,
                            color = col.coef.all,
                            linetype = "dotted")
    }
  }
