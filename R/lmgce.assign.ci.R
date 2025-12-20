#' Generalized Cross entropy estimation
#'
#' Internal function to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @param y A vector of observations of length n.
#' @param X A design matrix of dimension n * k.
#' @param object Fitted \code{lmgce} model object.
#' @inheritParams lmgce
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
lmgce.assign.ci <- function(y,
                            X,
                            offset,
                            errormeasure = "RMSE",
                            support.signal = NULL,
                            support.signal.points =
                              c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                            support.noise = NULL,
                            support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                            weight = 0.5,
                            method = "dual.lbfgsb3c",
                            caseGLM = "D",
                            seed = NULL,
                            object,
                            boot.B = 10,
                            boot.method = "residuals",
                            verbose = 0)
{
  res <- object

  res$results$bootstrap$coefficients <-
    data.frame(matrix(
      NA,
      ncol = boot.B,
      nrow = length(res$coefficients)
    ))

  row.names(res$results$bootstrap$coefficients) <-
    names(res$coefficients)

  res$results$bootstrap$nep <- NA

  res$results$bootstrap$nepk <- res$results$bootstrap$coefficients

  b = 0L
  attempts_b <- 0L
  complete_b <- NULL

  if (verbose >= 1L) cat("0% ")

  while (b < boot.B & attempts_b < boot.B * 3L) {
    attempts_b <- attempts_b + 1L
    b <- b + 1L
    if (is.null(complete_b)) {b <- 1L} else
    {if (!((b - 1L) %in% complete_b)) {
      b <- b - 1L}}

    tryCatch({
      #ci_seed <- ifelse(is.null(seed), 0L, seed) + attempts_b
      #set.seed(ci_seed)
      ci_ind <- sample(nrow(X), nrow(X), replace = TRUE)
      if (boot.method == "cases") {
        y.aux <- y[ci_ind]
        X.aux <- X[ci_ind,]
        offset.aux <- offset[ci_ind]
      } else if (boot.method == "residuals") {
        y.aux <- fitted(res) + resid(res)[ci_ind]
        X.aux <- X
        offset.aux <- offset
      } else if (boot.method == "wild") {
        y.aux <- fitted(res) + resid(res) * rnorm(length(resid(res)))
        X.aux <- X
        offset.aux <- offset
      }

      res_ci <-
        lmgce.assign.noci(
          y.aux,
          X.aux,
          offset.aux,
          FALSE,
          5,
          1,
          errormeasure,
          "min",
          NULL,
          NULL,
          NULL,
          support.signal,
          NULL,
          0.5,
          50,
          20,
          support.signal.points,
          support.noise,
          support.noise.points,
          weight,
          method,
          caseGLM,
          seed = seed,#ci_seed,
          verbose
        )

      res$results$bootstrap$coefficients[, b] <- res_ci$coefficients
      res$results$bootstrap$nep[b] <- res_ci$nep
      res$results$bootstrap$nepk[, b] <- res_ci$nepk

      complete_b <- c(complete_b, b)
    },
    error = function(e) {
      b <<- b - 1L
    })

    if (verbose >= 1L && (b %% 10 == 0)) cat(round(b / boot.B * 100, 0),
                                             "% ",
                                             sep = "")
  }

  if (b < boot.B & attempts_b == boot.B * 3L) {
    res$results$bootstrap$coefficients <- NULL
    res$results$bootstrap$nep <- NULL
    res$results$bootstrap$nepk <- NULL
    warning("\nNot possible to determine bootstrapp CI.
            Number of attemps exceeded.")
  }

  res$results$bootstrap$attempts <- attempts_b

  return(res)

}
