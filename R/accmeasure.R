#' Accuracy measures
#'
#' Function that allows to calculate different types of errors for point predictions:
#' \enumerate{
#' \item MAE - Mean Absolute Error,
#' \item MAD - Mean Absolute Deviation,
#' \item MSE - Mean Squared Error,
#' \item RMSE - Root Mean Squared Error,
#' \item MAPE - Mean Absolute Percentage Error,
#' \item sMAPE - symmetric Mean Absolute Percentage Error,
#' \item MASE - Mean Absolute Scaled Error (Hyndman & Koehler, 2006)
#' }
#'
#' @param y_pred fitted values.
#' @param y_true observed values.
#' @param which one of c("RMSE", "MAPE", "sMAPE", "MAE", "MAD", "MASE")
#'
#'
#' @return The value of the chosen error is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @references
#' Hyndman, R. J., & Koehler, A. B. (2006)
#' \emph{Another look at measures of forecast accuracy.}
#' International Journal of Forecasting, 22(4), 679–688.
#' \doi{10.1016/j.ijforecast.2006.03.001}\cr
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
#' accmeasure(fitted(res_gce_package), dataGCE$y, which = "MSE")
#'
#' @export

accmeasure <- function (y_pred,
                        y_true,
                        which = c("RMSE", "MSE", "MAPE", "sMAPE", "MAE", "MAD", "MASE"))
{
  which <- match.arg(which)
  switch(
    which,
    RMSE = {
      accm <- sqrt(mean((y_true - y_pred)^2))
    },
    MSE = {
      accm <- mean((y_true - y_pred)^2)
    },
    MAPE = {
      accm <- 100 / length(y_true) * sum(abs((y_true - y_pred) / y_true))
    },
    sMAPE = {
      accm <- 100 / length(y_true) * sum(abs(y_true - y_pred) / (abs(y_true) + abs(y_pred)))
    },
    MAE = {
      accm <- 1 / length(y_true) * sum(abs(y_true - y_pred))
    },
    MAD = {
      accm <- 1 / length(y_true) * sum(abs(y_true - mean(y_true)))
    },
    MASE = {
      accm <- (1 / length(y_true) * sum(abs(y_true - y_pred))) / (1 / length(y_true) * sum(abs(y_true - mean(y_true))))
    },
    stop(
      'which must be one of c("RMSE", "MAPE", "sMAPE", "MAE", "MAD", "MASE")'
    )
  )
  return(accm)
}
