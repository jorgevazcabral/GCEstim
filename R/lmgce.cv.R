#' Generalized Cross entropy estimation
#'
#' Internal function used to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @inheritParams lmgce.assign.noci
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
lmgce.cv <- function(y,
                     X,
                     offset,
                     y.test = NULL,
                     X.test = NULL,
                     offset.test = NULL,
                     cv.nfolds = 5,
                     cv.repeats = 1,
                     errormeasure = "RMSE",
                     min.coef = NULL,
                     max.coef = NULL,
                     max.abs.residual = NULL,
                     max.abs.coef = NULL,
                     support.signal = NULL,
                     support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                     support.noise = NULL,
                     support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                     weight = 0.5,
                     method = "dual.lbfgsb3c",
                     caseGLM = "D",
                     seed = NULL,
                     verbose = 0) {
  if (verbose >= 2)
    cat(0, "% ", sep = "")

  res <- list(
    cvresults = list(),
    nep.cv.mean = NULL,
    nep.cv.sd = NULL,
    error.measure.cv.mean = NULL,
    error.measure.cv.sd = NULL
  )

  res <-
    lmgce.fit(
      y,
      X,
      offset,
      y.test,
      X.test,
      offset.test,
      errormeasure,
      min.coef,
      max.coef,
      max.abs.residual,
      support.signal,
      support.signal.points,
      support.noise,
      support.noise.points,
      weight,
      method,
      caseGLM)

  if (verbose >= 2) {
    cat(round(1 / (cv.repeats * cv.nfolds + 1) * 100, 0), "% ", sep = "")
  }

  for (r in 1:cv.repeats) {
    res$cvresults[[r]] <- list()

    if (!is.null(seed))
      set.seed(seed)

    auxfolds = cut(seq(1, nrow(X)),
                   breaks = cv.nfolds,
                   labels = FALSE)
    change_order <- sample(nrow(X))

    for (cv in 1:cv.nfolds) {
      res$cvresults[[r]][[cv]] <- list()

      res$cvresults[[r]][[cv]] <-
        lmgce.fit(
          y = y[change_order][auxfolds != cv],
          X = X[change_order, ][auxfolds != cv, ],
          offset[change_order][auxfolds != cv],
          y.test = y[change_order][auxfolds == cv],
          X.test = X[change_order, ][auxfolds == cv, ],
          offset.test[change_order][auxfolds == cv],
          errormeasure,
          min.coef,
          max.coef,
          max.abs.residual,
          support.signal,
          support.signal.points,
          support.noise,
          support.noise.points,
          weight,
          method,
          caseGLM
        )

      res$error.measure.cv.mean <-
        c(res$error.measure.cv.mean,
          res$cvresults[[r]][[cv]]$error.measure[[1]])

      res$nep.cv.mean <-
        c(res$nep.cv.mean,
          res$cvresults[[r]][[cv]]$nep[[1]])

      if (verbose >= 2)
        cat(
          round(
            ((r - 1) * cv.nfolds + cv + 1) / (cv.repeats * cv.nfolds + 1) * 100,
            0),
          "% ",
          sep = "")
    }

    names(res$cvresults[[r]]) <- paste0("fold", 1:cv.nfolds)
  }

  names(res$cvresults) <- paste0("repeats", 1:cv.repeats)

  res$cvresults$error.measure.cv <- res$error.measure.cv.mean
  res$cvresults$error.measure.cv.mean <- mean(res$error.measure.cv.mean)
  res$cvresults$error.measure.cv.sd <- sd(res$error.measure.cv.mean)
  res$error.measure.cv.sd <- res$cvresults$error.measure.cv.sd
  res$error.measure.cv.mean <- res$cvresults$error.measure.cv.mean

  res$cvresults$nep.cv <- res$nep.cv.mean
  res$cvresults$nep.cv.mean <- mean(res$nep.cv.mean)
  res$cvresults$nep.cv.sd <- sd(res$nep.cv.mean)
  res$nep.cv.sd <- res$cvresults$nep.cv.sd
  res$nep.cv.mean <- res$cvresults$nep.cv.mean

  return(res)

}
