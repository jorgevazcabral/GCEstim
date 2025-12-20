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
lmgce.sscv <- function(y,
                       X,
                       offset,
                       y.test = NULL,
                       X.test = NULL,
                       offset.test = NULL,
                       cv.nfolds = 5,
                       cv.repeats = 1,
                       errormeasure = "RMSE",
                       errormeasure.which = "min",
                       min.coef = NULL,
                       max.coef = NULL,
                       max.abs.residual = NULL,
                       support.signal = NULL,
                       support.signal.vector = NULL,
                       support.signal.points =
                         c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                       support.noise = NULL,
                       support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                       weight = 0.5,
                       method = "dual.lbfgsb3c",
                       caseGLM = "D",
                       seed = NULL,
                       verbose = 0) {
  if (verbose >= 2)
    cat("[0%]\n", sep = "")

  names.supports <- support.signal.vector
  nsupports <- length(names.supports)

  aux.res <-
    lmgce.ss(
      y,
      X,
      offset,
      y.test,
      X.test,
      offset.test,
      errormeasure,
      errormeasure.which,
      min.coef,
      max.coef,
      max.abs.residual,
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

  res <- list(
    results = list(nocv = aux.res,
                   cvresults = list(),
                   bootstrap = NULL),
    support = names.supports,
    support.ok = NULL
  )

  mean.measure.ss <- as.data.frame(matrix(
    data = 0,
    nrow = cv.repeats * cv.nfolds,
    ncol = nsupports,
    dimnames = list(c(
      apply(expand.grid(
        paste0("fold", 1:cv.nfolds),
        paste0("repeat", 1:cv.repeats)
      ),
      1,
      paste,
      collapse = "_")
    ),
    names.supports)
  ))

  mean.nep.ss <- mean.measure.ss

  if (verbose >= 2)
    cat("\n[",
        round(1 / (cv.repeats * cv.nfolds + 1) * 100, 0),
        "%]\n",
        sep = "")

  for (r in 1:cv.repeats) {
    res$results$cvresults[[r]] <- list()

    if (!is.null(seed))
      set.seed(seed)

    auxfolds = cut(seq(1, nrow(X)),
                   breaks = cv.nfolds,
                   labels = FALSE)
    change_order <- sample(nrow(X))

    for (cv in 1:cv.nfolds) {
      res$results$cvresults[[r]][[cv]] <- list()

      res_k <-
        lmgce.ss(
          y = y[change_order][auxfolds != cv],
          X = X[change_order, ][auxfolds != cv,],
          offset[change_order][auxfolds != cv],
          y.test = y[change_order][auxfolds == cv],
          X.test = X[change_order, ][auxfolds == cv,],
          offset.test[change_order][auxfolds == cv],
          errormeasure,
          errormeasure.which,
          min.coef,
          max.coef,
          max.abs.residual,
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

      #aux.mean.measure.ss <- NULL
      #aux.mean.nep.ss <- NULL

      for (amm.ss in 1:nsupports) {
        mean.measure.ss[(r - 1) * cv.nfolds + cv, amm.ss] <-
          res_k$support.results[[amm.ss]]$error.measure[[1]]
        mean.nep.ss[(r - 1) * cv.nfolds + cv, amm.ss] <-
          res_k$support.results[[amm.ss]]$nep[[1]]
      }

      res$results$cvresults[[r]][[cv]] <- res_k

      if (verbose >= 2)
        cat("\n[",
          round(((r - 1) * cv.nfolds + cv + 1) / (cv.repeats * cv.nfolds + 1) *
                           100, 0), "%]\n", sep = ""
        )
    }

    names(res$results$cvresults[[r]]) <- paste0("fold", 1:cv.nfolds)

  }

  names(res$results$cvresults) <- paste0("repeats", 1:cv.repeats)

  res$results$cvresults$error.measure.cv <- mean.measure.ss

  res$results$cvresults$error.measure.cv.mean <-
    colMeans(mean.measure.ss, na.rm = TRUE)

  names(res$results$cvresults$error.measure.cv.mean) <- names.supports

  res$results$cvresults$error.measure.cv.sd <-
    apply(mean.measure.ss, 2, sd, na.rm = TRUE)

  names(res$results$cvresults$error.measure.cv.sd) <- names.supports

  res$results$cvresults$nep.cv <- mean.nep.ss

  res$results$cvresults$nep.cv.mean <-
    colMeans(mean.nep.ss, na.rm = TRUE)

  names(res$results$cvresults$nep.cv.mean) <- names.supports

  res$results$cvresults$nep.cv.sd <-
    apply(mean.nep.ss, 2, sd, na.rm = TRUE)

  names(res$results$cvresults$nep.cv.sd) <- names.supports

  res$support.ok <-
    as.numeric(
     intersect(
      names(
        which(
          !sapply(
            lapply(res$results$cvresults$error.measure.cv, is.na),
            any))),
    names(
      which(
        !sapply(
          lapply(res$results$nocv$support.results, is.na),
          isTRUE)))))

  if (length(res$support.ok) < 2)
    stop('Estimation successfully completed for less than 2 supports.
         \nPlease choose different supports.',
         call. = FALSE)

  res$support.signal.manual <- NULL
  res$support.signal.min <-
    as.numeric(names(which.min(res$results$cvresults$error.measure.cv.mean)))

  aux.support.signal.1se <-
    which.min(res$results$cvresults$error.measure.cv.mean)
  res$support.signal.1se <-
    (res$results$cvresults$error.measure.cv.mean[[aux.support.signal.1se]] +
       res$results$cvresults$error.measure.cv.sd[[
         aux.support.signal.1se]] / cv.nfolds) -
    res$results$cvresults$error.measure.cv.mean
  res$support.signal.1se[res$support.signal.1se > 0] <- -Inf
  res$support.signal.1se[
    (which.min(
      res$results$cvresults$error.measure.cv.mean) + 1):length(
        res$results$cvresults$error.measure.cv.mean)] <- -Inf

  res$support.signal.1se <- as.numeric(names(which.max(res$support.signal.1se)))

  res$support.signal.elbow <-
    pathviewr::find_curve_elbow(
      na.omit(
        data.frame(support.signal.vector,
                   res$results$cvresults$error.measure.cv.mean)),
      export_type = "")[[1]]

  if (errormeasure.which == "min")
    chosen_support <- res$support.signal.min
  else if (errormeasure.which == "1se")
    chosen_support <- res$support.signal.1se
  else
    chosen_support <- res$support.signal.elbow

  chosen_support <- which(support.signal.vector == chosen_support)

  res <- c(res,
           res$results$nocv$support.results[[chosen_support]],
           res$results$cvresults$nep.cv.mean[chosen_support],
           res$results$cvresults$nep.cv.sd[chosen_support],
           res$results$cvresults$error.measure.cv.mean[chosen_support],
           res$results$cvresults$error.measure.cv.sd[chosen_support]
           )

  names(res)[length(res) - 3] <- "nep.cv.mean"
  names(res)[length(res) - 2] <- "nep.cv.sd"
  names(res)[length(res) - 1] <- "error.measure.cv.mean"
  names(res)[length(res)] <- "error.measure.cv.sd"

  return(res)

}
