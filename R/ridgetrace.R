# flexiRidGCEw estimation

#' Function to obtain the ridge trace and choose the support limits given a
#' formula
#'
#' Function to obtain the ridge trace and choose the support limits given a
#' formula
#'
#' @param formula An object of class \code{\link[stats]{formula}} (or one that
#'  can be coerced to that class): a symbolic description of the model to be
#'  fitted.
#' @param data A data frame (or object coercible by
#' \code{\link[base]{as.data.frame}} to a data frame) containing the variables
#'  in the model.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s. The default is set by the \code{na.action} setting of
#' \code{\link[base]{options}}, and is \code{\link[stats]{na.fail}} if that is
#' unset. The ‘factory-fresh’ default is \code{\link[stats]{na.omit}}. Another
#' possible value is \code{NULL}, no action. Value
#'  \code{\link[stats]{na.exclude}} can be useful.
#' @param offset this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting. This should be \code{NULL}
#' or a numeric vector or matrix of extents matching those of the response. One
#' or more \code{\link[stats]{offset}} terms can be included in the formula
#' instead or as well, and if more than one are specified their sum is used.
#' See \code{\link[stats]{model.offset}}.
#' @param contrasts An optional list. See the \code{contrasts.arg} of
#' \code{\link[stats]{model.matrix.default}}.
#' @param lambda Ridge parameter. The default is \code{lambda = NULL} and a
#' lambda logarithmic sequence will be computed based on
#'  \code{lambda.n}, \code{lambda.min} and \code{lambda.max}. Supplying a lambda
#'  sequence overrides this.
#' @param lambda.min Minimum value for the \code{lambda} sequence. The default
#' id \code{lambda.min = 10^-3}. To be used when \code{lambda = NULL}.
#' @param lambda.max Maximum value for the \code{lambda} sequence. The default
#' id \code{lambda.max = 10^3}. To be used when \code{lambda = NULL}.
#' @param lambda.n The number of lambda values. The default is
#' \code{lambda.n = 100}. To be used when \code{lambda = NULL}.
#' @param standardize Boolean value. If \code{TRUE}, the default, then: i)
#' centering is done by subtracting the column means of x and y from their
#' corresponding columns; ii) scaling is done by dividing the (centered) columns
#'  of x and y by their standard deviations.
#' @param penalize.intercept Boolean value. If \code{TRUE}, the default, the
#' intercept will be penalized when \code{standardize = FALSE}.
#' @param cv Boolean value. If \code{TRUE} the error, \code{errormeasure},
#' will be computed using cross-validation. If \code{FALSE} the error will be
#' computed in sample. The default is \code{cv = TRUE}.
#' @param cv.nfolds number of folds used for cross-validation when
#' \code{cv = TRUE}. The default is \code{cv.nfolds = 5}.
#' @param errormeasure Loss function (error) to be used for the selection
#' of the support spaces. One of
#' c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE"). The default is
#' \code{errormeasure = "RMSE"}.
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = 230676}.
#'
#' @return
#'  An object of \code{\link[base]{class}} \code{ridgetrace} is a list
#'  containing at least the following components:
#'
#' \item{lambda}{the lambda sequence used}
#' \item{min.coef}{a named vector of coefficients (minimum
#' coefficients)}
#' \item{max.coef}{a named vector of coefficients (maximum
#' coefficients)}
#' \item{max.abs.coef}{a named vector of coefficients (maximum absolute
#' coefficients)}
#' \item{max.abs.residual}{the maximum absolute residual}
#' \item{coef.lambda}{a data.frame with the coefficients for each lambda tested}
#' \item{coef.lambda.cv}{a list of length \code{cv.nfolds} containing
#' data.frames with the coefficients for each lambda tested in each
#' cross-validation fold}
#' \item{error.lambda}{a vector with the in sample error}
#' \item{error.lambda.cv}{a data.frame with cross-validation errors}
#' \item{call}{the matched call}
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
#' @export

ridgetrace <- function(formula,
                       data,
                       subset,
                       na.action,
                       offset,
                       contrasts = NULL,
                       lambda = NULL,
                       lambda.min = 10^-3,
                       lambda.max = 10^3,
                       lambda.n = 100,
                       standardize = TRUE,
                       penalize.intercept = TRUE,
                       errormeasure = c("RMSE","MSE", "MAE",
                                        "MAPE", "sMAPE", "MASE"),
                       cv = TRUE,
                       cv.nfolds = 5,
                       seed = 230676)
{
  errormeasure <- match.arg(errormeasure)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")


  if (is.empty.model(mt))
    stop("Empty model specified.")

  y <- model.response(mf, "numeric")

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(
        gettextf(
          "number of offsets is %d, should equal %d (number of observations)",
          length(offset),
          NROW(y)
        ),
        domain = NA
      )
  }

  X <- model.matrix(mt, mf, contrasts)

  if (is.null(n <- nrow(X)))
    stop("'x' must be a matrix")

  if (n == 0L)
    stop("0 (non-NA) cases")

  if (ncol(X) == 0L)
    stop("NULL model")

  if (!is.null(offset))
    y <- y - offset

  if (NROW(y) != n)
    stop("incompatible dimensions")

  res <-
    ridgetrace.Xy(
      X = X,
      y = y,
      lambda = lambda,
      lambda.min = lambda.min,
      lambda.max = lambda.max,
      lambda.n = lambda.n,
      standardize = standardize,
      penalize.intercept = penalize.intercept,
      errormeasure = errormeasure,
      cv = cv,
      cv.nfolds = cv.nfolds,
      seed = seed
    )

  res$call <- cl
  return(res)
}

#' Function to obtain the ridge trace and choose the support limits given X
#' and y
#'
#' Function to obtain the ridge trace and choose the support limits given X
#' and y
#'
#' @param X A model matriz.
#' @param y A vector containing the response.
#' @param lambda Ridge parameter. The default is \code{lambda = NULL} and a
#' lambda logarithmic sequence will be computed based on \code{lambda.n},
#' \code{lambda.min} and \code{lambda.max}. Supplying a lambda sequence
#' overrides this.
#' @param lambda.min Minimum value for the \code{lambda} sequence. The default
#' id \code{lambda.min = 10^-3}.
#' @param lambda.max Maximum value for the \code{lambda} sequence. The default
#' id \code{lambda.max = 10^3}.
#' @param lambda.n The number of lambda values. The default is
#' \code{lambda.n = 100}.
#' @param standardize Boolean value. If \code{TRUE}, the default, then: i)
#' centering is done by subtracting the column means of x and y from their
#' corresponding columns; ii) scaling is done by dividing the (centered) columns
#'  of x and y by their standard deviations.
#' @param penalize.intercept Boolean value. If \code{TRUE}, the default, the
#' intercept will be penalized when \code{standardize = FALSE}.
#' @param cv Boolean value. If \code{TRUE} the error, \code{errormeasure},
#' will be computed using cross-validation. If \code{FALSE} the error will be
#' computed in sample. The default is \code{cv = TRUE}.
#' @param cv.nfolds number of folds used for cross-validation when
#' \code{cv = TRUE}. The default is \code{cv.nfolds = 5}.
#' @param errormeasure Loss function (error) to be used for the selection
#' of the support spaces. One of
#' c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE"). The default is
#' \code{errormeasure = "RMSE"}.
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = 230676}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ridgetrace.Xy <- function(X,
                          y,
                          lambda = NULL,
                          lambda.min = 10^-3,
                          lambda.max = 10^3,
                          lambda.n = 100,
                          standardize = TRUE,
                          penalize.intercept = TRUE,
                          errormeasure = c("RMSE","MSE", "MAE",
                                           "MAPE", "sMAPE", "MASE"),
                          cv = TRUE,
                          cv.nfolds = 5,
                          seed = 230676){

  if (is.null(lambda)) {
    lambda <-
      lambda.min * (lambda.max / lambda.min) ^ (
        (0:(lambda.n - 1)) / (lambda.n - 1))
  } else {
    lambda.n <- length(lambda)
  }

  k <- ncol(X)
  n <- nrow(X)

  X.intercept <- "(Intercept)" %in% colnames(X)

  coef_lambda <- matrix(0, ncol = lambda.n, nrow = k)
  resid_lambda <- fitted_lambda <- matrix(0, ncol = lambda.n, nrow = n)
  error.lambda.cv <- NULL
  coef.cv <- NULL #vector("list", cv.nfolds)

  if (standardize) {
    y_tilde <- scale(y)
    if (X.intercept) {
      X_tilde <- scale(X[, colnames(X) != "(Intercept)"])
      penalty <- diag(rep(1, k - 1))
      coef_lambda_tilde <- matrix(0, ncol = lambda.n, nrow = k - 1)
    } else {
      X_tilde <- scale(X)
      coef_lambda_tilde <- matrix(0, ncol = lambda.n, nrow = k)
    }
  } else {
    if (X.intercept) {
      penalty <- diag(c(as.numeric(penalize.intercept), rep(1, k - 1)))
    } else {
      penalty <- diag(rep(1, k))
      }
  }

  if (isTRUE(cv)){
    coef.cv <- vector("list", cv.nfolds)
    names(coef.cv) <- paste0("fold", 1:cv.nfolds)
    for (i in 1:cv.nfolds) {
      coef.cv[[i]] <- coef_lambda
      colnames(coef.cv[[i]]) <- paste0("lambda", round(lambda, 8))
      rownames(coef.cv[[i]]) <- colnames(X)
    }
    if (!is.null(seed))
      set.seed(seed)
    auxfolds = cut(seq(1, nrow(X)),
                   breaks = cv.nfolds,
                   labels = FALSE)
    change_order <- sample(nrow(X))

    error.lambda.cv <- data.frame(matrix(NA, ncol = lambda.n, nrow = cv.nfolds))
    }

  for (i in 1:lambda.n) {
    if (standardize) {
      coef_lambda_tilde[, i] <-
        solve(
          t(X_tilde) %*% X_tilde + lambda[i] * penalty) %*% (
            t(X_tilde) %*% y_tilde)
      if (X.intercept) {
        coef_lambda[, i] <-
          scalebackcoef(X_tilde,
                        y_tilde,
                        c(0, coef_lambda_tilde[, i]),
                        intercept = TRUE)
      } else {
        coef_lambda[, i] <-
          scalebackcoef(X_tilde,
                        y_tilde,
                        coef_lambda_tilde[, i],
                        intercept = FALSE)
        }
    } else {
      coef_lambda[, i] <-
        solve(t(X) %*% X + lambda[i] * penalty) %*% (t(X) %*% y)
    }

    fitted_lambda[, i] <- X %*% coef_lambda[, i]
    resid_lambda[, i] <- y - fitted_lambda[, i]

    if (isTRUE(cv)){
      for (cv.n in 1:cv.nfolds) {
        y.cv = y[change_order][auxfolds != cv.n]
        X.cv = X[change_order, ][auxfolds != cv.n,]

        if (standardize) {
          y_tilde.cv = y_tilde[change_order][auxfolds != cv.n, drop = FALSE]
          attr(y_tilde.cv, "scaled:center") <-
            attr(y_tilde, "scaled:center")
          attr(y_tilde.cv, "scaled:scale")  <-
            attr(y_tilde, "scaled:scale")

          X_tilde.cv = X_tilde[change_order, ][auxfolds != cv.n, ,drop = FALSE]
          attr(X_tilde.cv, "scaled:center") <-
            attr(X_tilde, "scaled:center")
          attr(X_tilde.cv, "scaled:scale")  <-
            attr(X_tilde, "scaled:scale")

          coef_lambda_cv <-
            solve(t(X_tilde.cv) %*% X_tilde.cv + lambda[i] * penalty) %*% (
              t(X_tilde.cv) %*% y_tilde.cv)
          if (X.intercept) {
            coef_lambda_cv <-
              scalebackcoef(X_tilde.cv,
                            y_tilde.cv,
                            c(0, coef_lambda_cv),
                            intercept = TRUE)
          } else {
            coef_lambda_cv <-
              scalebackcoef(X_tilde.cv,
                            y_tilde.cv,
                            coef_lambda_cv,
                            intercept = FALSE)
          }
        } else {
          coef_lambda_cv <-
            solve(t(X.cv) %*% X.cv + lambda[i] * penalty) %*% (t(X.cv) %*% y.cv)
          }

        coef.cv[[cv.n]][, i] <- coef_lambda_cv
        error.lambda.cv[cv.n, i] <-
            accmeasure(y[change_order][auxfolds == cv.n],
                       X[change_order, ][auxfolds == cv.n,] %*% coef_lambda_cv,
                       errormeasure)
      }
      colnames(error.lambda.cv) <- paste0("lambda", round(lambda, 8))
      rownames(error.lambda.cv) <- paste0("fold", 1:cv.nfolds)
    }
  }

  min.coef <- apply(coef_lambda, 1, function(x){min(x)})
  max.coef <- apply(coef_lambda, 1, function(x){max(x)})
  max.abs.coef <- apply(coef_lambda, 1, function(x){max(abs(x))})
  max.abs.residual <- max(abs(resid_lambda))

  error.lambda <- apply(fitted_lambda, 2, accmeasure, y, errormeasure)

  names(max.coef) <- colnames(X)
  names(min.coef) <- colnames(X)
  names(max.abs.coef) <- colnames(X)
  rownames(coef_lambda) <- colnames(X)
  colnames(coef_lambda) <- paste0("lambda", round(lambda, 8))

  res <- list(lambda = lambda,
              max.coef = max.coef,
              min.coef = min.coef,
              max.abs.coef = max.abs.coef,
              max.abs.residual = max.abs.residual,
              coef.lambda = coef_lambda,
              coef.lambda.cv = coef.cv,
              error.lambda = error.lambda,
              error.lambda.cv = error.lambda.cv
              )

  class(res) <- "ridgetrace"

  return(res)
  }
