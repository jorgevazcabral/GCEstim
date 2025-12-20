#' Normalized Entropy Aggregation for Inhomogeneous Large-Scale Data - Neagging
#'
#' Computes the estimates for the Normalized Entropy Aggregation
#'
#' @param object Fitted \code{\link{lmgce}} or \code{\link{tsbootgce}} model
#' object.
#' @param boot.B To use with a \code{\link{lmgce}} object. A single positive
#' integer greater or equal to 10 for the number of bootstrap replicates for the
#' computation of the Normalized Entropy Aggregation estimate(s), to be used
#' when \code{object} was created with \code{boot.B = 0}. The default is
#' \code{boot.B = 100} when the \code{object} has no previous sampling
#' information and \code{boot.B = object$boot.B} otherwise, which corresponds
#' to the \code{boot.B} given to \code{lmgce} when the \code{object} was
#' created.
#' @param boot.method To use with a \code{\link{lmgce}} object. Method used for
#' bootstrapping. One of \code{c("residuals", "cases", "wild")} which
#' corresponds to resampling on residuals, on individual cases or on residuals
#' multiplied by a N(0,1) variable, respectively. The default is
#'  \code{boot.method = object$boot.method}.
#' @param error Loss function (error) to be used for the selection
#' of the support spaces. One of
#' c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE"). The default is
#' \code{boot.method = object$error}.
#'
#' @return
#'  An object of \code{\link[base]{class}} \code{neagging} is a list containing
#'  at least the following components:
#'
#' \item{matrix}{a matrix where each column contains sequentially the aggregated
#'  estimates.}
#' \item{error}{a named vector with the in sample error for each aggregated set
#' of estimates.}
#' \item{error.object}{the in sample error of the \code{object}.}
#' \item{coefficients}{the aggregated coefficients that produced the lowest in
#' sample error.}
#' \item{coefficients.object}{the coefficients of the \code{object}.}
#'
#' @seealso
#' The generic functions \code{\link{plot.neagging}} and
#' \code{\link{coef.neagging}}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @references
#' da Conceição Costa, M. and Macedo, P. (2019).
#' \emph{Normalized Entropy Aggregation for Inhomogeneous Large-Scale Data.}
#' In O. Valenzuela, F. Rojas, H. Pomares, & I. Rojas (Eds.), Theory and
#' Applications of Time Series Analysis (pp. 19–29).
#' Springer International Publishing.
#' \doi{10.1007/978-3-030-26036-1_2}
#'
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' neagging(res_gce_package, boot.method = "cases")
#'
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' neagging(res.tsbootgce)
#' }
#'
#' @export

neagging <- function(object,
                     boot.B = ifelse(object$boot.B == 0,
                                     100,
                                     object$boot.B),
                     boot.method = object$boot.method,
                     error = object$error)
{
  if (!inherits(object, "lmgce") & !inherits(object, "tsbootgce"))
    stop("use only with \"lmgce\" or \"tsbootgce\" objects")

  res <- list(matrix = NA,
              error = NA,
              error.object = NA,
              coefficients = NA,
              coefficients.object = NA)

  if (inherits(object, "lmgce")) {
  if (is.null(object$y) ||
      !("x" %in% names(object)) ||
      is.null(object$results$bootstrap$coefficients) ||
      boot.B != object$boot.B ||
      boot.method != object$boot.method) {

    udobj <- update(
      object,
      y = TRUE,
      x = TRUE,
      support.signal = object$support.matrix,
      boot.B = boot.B,
      boot.method = boot.method,
      verbose = 0
    )
    object$x <- udobj$x
    object$y <- udobj$y
    object$results$bootstrap <-
      udobj$results$bootstrap
    }

    res$coefficients.object <- coef(object)
    res$error.object <- object$error.measure
    names(res$error.object) <- object$error

  res$matrix <-
    data.frame(
      matrix(NA,
             ncol = ncol(object$results$bootstrap$nepk),
             nrow = nrow(object$results$bootstrap$nepk)))

  row.names(res$matrix) <- row.names(object$results$bootstrap$nepk)

  res$matrix[, 1] <- object$results$bootstrap$coefficients[, 1]
  res$error[1] <- accmeasure(object$x %*% res$matrix[, 1],
                             object$y,
                             error)
  for (i in 2:ncol(object$results$bootstrap$nepk)) {
    nepkweights <-
        t(apply(
          (1 - object$results$bootstrap$nepk[, 1:i]),
          1,
          function(x) {x/sum(x)}))

    res$matrix[, i] <-
      diag(
        as.matrix(
          object$results$bootstrap$coefficients[, 1:i]) %*% t(nepkweights))

    res$error[i] <- accmeasure(object$x %*% res$matrix[, i],
                               object$y,
                               error)
  }

  names(res$error) <- 1:boot.B
  res$coefficients <- res$matrix[, i]
  names(res$coefficients) <- row.names(res$matrix)
  res$matrix <- as.matrix(res$matrix)

  } else {

    res$coefficients.object <- coef(object)
    res$error.object <- object$error.measure
    names(res$error.object) <- object$error

    res$matrix <-
      data.frame(
        matrix(NA,
               ncol = ncol(object$results$bootstrap$nepk.matrix),
               nrow = nrow(object$results$bootstrap$nepk.matrix)))

    row.names(res$matrix) <- row.names(object$results$bootstrap$nepk.matrix)

    res$matrix[, 1] <- object$results$bootstrap$coef.matrix[, 1]
    res$error[1] <- accmeasure(object$lmgce$x %*% res$matrix[, 1],
                                        object$lmgce$y,
                                        error)

    for (i in 2:ncol(object$results$bootstrap$nepk.matrix)) {
      nepkweights <-
        t(apply(
          (1 - object$results$bootstrap$nepk.matrix[, 1:i]),
          1,
          function(x) {x/sum(x)}))

      res$matrix[, i] <-
        diag(
          as.matrix(
            object$results$bootstrap$coef.matrix[, 1:i]) %*% t(nepkweights))

      res$error[i] <- accmeasure(object$lmgce$x %*% res$matrix[, i],
                                          object$lmgce$y,
                                          error)
    }
    names(res$error) <- 1:ncol(object$results$bootstrap$nepk.matrix)
    res$coefficients <- res$matrix[, i]
    names(res$coefficients) <- row.names(res$matrix)
    res$matrix <- as.matrix(res$matrix)
  }
  class(res) <- "neagging"

  return(res)

}
