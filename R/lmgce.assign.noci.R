#' Generalized Cross entropy estimation
#'
#' Internal function to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @inheritParams lmgce.assign.ci
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
lmgce.assign.noci <- function(y,
                              X,
                              offset,
                              cv = TRUE,
                              cv.nfolds = 5,
                              cv.repeats = 1,
                              errormeasure = "RMSE",
                              errormeasure.which = "min",
                              max.abs.coef = NULL,
                              support.signal = NULL,
                              support.signal.vector = NULL,
                              support.signal.vector.min = 0.5,
                              support.signal.vector.max = 50,
                              support.signal.vector.n = 20,
                              support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                              support.noise = NULL,
                              support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                              weight = 0.5,
                              method = "dual.lbfgsb3c",
                              caseGLM = "D",
                              seed = NULL,
                              verbose = 0) {
  X.test <- NULL
  y.test <- NULL
  offset.test <- NULL

  sweep <- is.null(support.signal)

  if (sweep) {
    if (is.null(support.signal.vector)) {
      support.signal.vector <-
        round(
          exp(seq(log(support.signal.vector.min),
                  log(support.signal.vector.max),
                  length.out = support.signal.vector.n)),
          4)
    } else {
      support.signal.vector <- sort(unique(support.signal.vector), decreasing = FALSE)
    }
    support.signal <- support.signal.vector[1]
  }

  if (cv) {
    if (sweep) {
      res <-
        lmgce.sscv(
          y,
          X,
          offset,
          y.test,
          X.test,
          offset.test,
          cv.nfolds,
          cv.repeats,
          errormeasure,
          errormeasure.which,
          max.abs.coef,
          support.signal,
          support.signal.vector,
          support.signal.points,
          support.noise,
          support.noise.points,
          weight,
          method,
          caseGLM,
          seed,
          verbose)

      names(res$results)[names(res$results) == "cvresults"] <- "cv"

    } else {
      res <- list(results = list(nocv = list(support.results = list()),
                                 cv = NULL,
                                 bootstrap = NULL),
                  support.signal.manual = NULL,
                  support.signal.min = NULL,
                  support.signal.1se = NULL,
                  support.signal.elbow = NULL,
                  support = ifelse(length(support.signal) == 1, support.signal, "interval"),
                  support.ok = ifelse(length(support.signal) == 1, support.signal, "interval"))
      res <-
        c(
          res,
          lmgce.cv(
            y,
            X,
            offset,
            y.test,
            X.test,
            offset.test,
            cv.nfolds,
            cv.repeats,
            errormeasure,
            max.abs.coef,
            support.signal,
            support.signal.points,
            support.noise,
            support.noise.points,
            weight,
            method,
            caseGLM,
            seed,
            verbose
          )
        )
      res$results$cv <- res$cvresults
      res <- res[names(res) != "cvresults"]
      res$results$nocv$support.results[[1]] <- res[names(res$results$cv[[1]][[1]])]
      names(res$results$nocv$support.results) <- res$support
    }
  }
  else {
    if (sweep) {
      res <- list(results = list(nocv = NULL,
                                 cv = NULL,
                                 bootstrap = NULL),
                  nep.cv.mean = NULL,
                  nep.cv.sd = NULL,
                  error.measure.cv.mean = NULL,
                  error.measure.cv.sd = NULL)
      res <-
        c(
          res,
          lmgce.ss(
            y,
            X,
            offset,
            y.test,
            X.test,
            offset.test,
            errormeasure,
            errormeasure.which,
            max.abs.coef,
            support.signal,
            support.signal.vector,
            support.signal.points,
            support.noise,
            support.noise.points,
            weight,
            method,
            caseGLM,
            verbose
          )
        )
      res$results$nocv$support.results <- res$support.results
      res <- res[names(res) != "support.results"]
    } else {
      res <- list(results = list(nocv = list(support.results = list()),
                                 cv = NULL,
                                 bootstrap = NULL),
                  support.signal.manual = NULL,
                  support.signal.min = NULL,
                  support.signal.1se = NULL,
                  support.signal.elbow = NULL,
                  nep.cv.mean = NULL,
                  nep.cv.sd = NULL,
                  error.measure.cv.mean = NULL,
                  error.measure.cv.sd = NULL,
                  support = ifelse(length(support.signal) == 1, support.signal, "interval"),
                  support.ok = ifelse(length(support.signal) == 1, support.signal, "interval"))

      res$results$nocv$support.results[[1]] <-
        lmgce.fit(
          y,
          X,
          offset,
          y.test,
          X.test,
          offset.test,
          errormeasure,
          max.abs.coef,
          support.signal,
          support.signal.points,
          support.noise,
          support.noise.points,
          weight,
          method,
          caseGLM
        )
      res <-
        c(res, res$results$nocv$support.results[[1]])
      names(res$results$nocv$support.results) <- res$support
    }
  }

  return(res)
}
