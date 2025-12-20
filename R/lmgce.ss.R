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
lmgce.ss <- function(y,
                     X,
                     offset,
                     y.test = NULL,
                     X.test = NULL,
                     offset.test = NULL,
                     errormeasure = "RMSE",
                     errormeasure.which = "min",
                     min.coef = NULL,
                     max.coef = NULL,
                     max.abs.residual = NULL,
                     support.signal = NULL,
                     support.signal.vector = NULL,
                     support.signal.points = c(1/5, 1/5, 1/5, 1/5, 1/5),
                     support.noise = NULL,
                     support.noise.points = c(1/3, 1/3, 1/3),
                     weight = 0.5,
                     method = "dual.lbfgsb3c",
                     caseGLM = "D",
                     verbose = 0) {
  if (verbose >= 3) {
    cat("0% ", sep = "")
  }

  nsupports <- length(support.signal.vector)

  res <- list(
    support.results = list(list()),
    support = support.signal.vector,
    support.ok = NULL
  )

  for (su in 1:nsupports) {

    res$support.results[[su]] <-
      tryCatch({
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
          support.signal = support.signal.vector[su],
          support.signal.points,
          support.noise,
          support.noise.points,
          weight,
          method,
          caseGLM)
      },
      error = function(e) {
        return(list(error.measure = NA,
                    nep = NA))
      })

    if (verbose >= 3) {
      cat(round((su) / nsupports * 100, 0), "% ", sep = "")
    }
  }

  names(res$support.results) <- support.signal.vector

  error.measure.ss <- NULL
  nep.ss <- NULL

  for (ms in 1:nsupports) {
    error.measure.ss <- c(error.measure.ss,
                          res$support.results[[ms]]$error.measure)
    nep.ss <- c(nep.ss,
                res$support.results[[ms]]$nep)
  }

  names(error.measure.ss) <- support.signal.vector
  names(nep.ss) <- support.signal.vector

  res$support.results$error.measure.ss <- error.measure.ss
  res$support.results$nep.ss <- nep.ss

  res$support.ok <-
    as.numeric(
      names(
        which(!is.na(error.measure.ss))))

  if (length(res$support.ok) < 2)
    stop('Estimation successfully completed for less than 2 supports.\nPlease
         choose different supports.',
         call. = FALSE)

  coef.matrix.ss <- matrix(NA,
                           nrow = ncol(X) * 2L,
                           ncol = length(res$support.ok))
  rownames(coef.matrix.ss) <- rep(colnames(X), 2L)
  colnames(coef.matrix.ss) <- res$support.ok

  for (i in 1:length(res$support.ok)) {
    coef.matrix.ss[,i] <-
      c(res$support.results[[
        which(
          names(res$support.results) %in% res$support.ok[i])]]$coefficients,
        res$support.results[[
          which(
            names(res$support.results) %in% res$support.ok[i])]]$nepk)
  }

  res$support.results$coef.matrix.ss <- coef.matrix.ss

  res$support.signal.manual <- NULL
  res$support.signal.min <- as.numeric(names(which.min(error.measure.ss)))
  res$support.signal.1se <- NULL
  res$support.signal.elbow <-
    pathviewr::find_curve_elbow(
      na.omit(data.frame(support.signal.vector,
                         error.measure.ss)),
      export_type = "")[[1]]

  if (errormeasure.which == "min")
    chosen_support <- res$support.signal.min
  else if (errormeasure.which == "elbow")
    chosen_support <- res$support.signal.elbow
  else if (errormeasure.which == "1se")
    chosen_support <- NULL

  if (!is.null(chosen_support)) {
  chosen_support <- which(names(res$support.results) == chosen_support)

  res <- c(res, res$support.results[[chosen_support]])
  }

  return(res)

}
