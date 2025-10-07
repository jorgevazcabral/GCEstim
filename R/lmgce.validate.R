#' Validation of \code{\link{lmgce}} arguments
#'
#' This function validates the arguments for \code{\link{lmgce}}
#'
#' @inheritParams lmgce
#'
#' @return
#' The process stops if arguments are not valid.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

lmgce.validate <- function(formula,
                           data,
                           subset,
                           na.action,
                           offset,
                           contrasts = NULL,
                           model = TRUE,
                           x = FALSE,
                           y = FALSE,
                           cv = TRUE,
                           cv.nfolds = 5,
                           cv.repeats = 1,
                           errormeasure = "RMSE",
                           errormeasure.which = "min",
                           support.method = c("standardized", "ridge"),
                           support.method.ridge.penalize.intercept = TRUE,
                           support.signal = NULL,
                           support.signal.vector = NULL,
                           support.signal.vector.min = 0.5,
                           support.signal.vector.max = 50,
                           support.signal.vector.n = 20,
                           support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                           support.noise = NULL,
                           support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                           weight = 0.5,
                           twosteps.n = 1,
                           method = "dual.lbfgsb3c",
                           caseGLM = "D",
                           boot.B = 0,
                           boot.method = "residuals",
                           seed = NULL,
                           OLS = TRUE,
                           verbose = 0) {

  if (missing(formula)) stop("argument `formula` is required.",call. = FALSE)
  if (missing(data)) stop("argument `data` is required.",call. = FALSE)

  if (support.method == "ridge" & (length(support.signal) > 1))
    stop("argument `support.signal` must be NULL or an integer when `support.method='ridge'`.",
         call. = FALSE)

  if (is.null(method) ||
      !(method %in% c("primal.solnl", "primal.solnp", "dual",
                      "dual.BFGS", "dual.CG", "dual.L-BFGS-B",
                      "dual.Rcgmin", "dual.bobyqa", "dual.newuoa",
                      "dual.nlminb", "dual.nlm",
                      "dual.lbfgs",
                      "dual.lbfgsb3c",
                      "dual.optimParallel"))) {
    stop(
      'argument `method` must be one of c("primal.solnl", "primal.solnp", "dual",
      "dual.BFGS", "dual.CG", "dual.L-BFGS-B", "dual.Rcgmin", "dual.bobyqa",
      "dual.newuoa", "dual.nlminb", "dual.nlm",
      "dual.lbfgs",
      "dual.lbfgsb3c",
      "dual.optimParallel")',
      call. = FALSE
    )
  }

  if (method == "dual" & length(unique(support.signal.points)) != 1) {
    stop('`method` = "dual" supports only GME.',
         call. = FALSE)
  }

  if (method == "dual" & weight != 0.5) {
    stop('`method` = "dual" supports only equal weights for the signal and noise.',
         call. = FALSE)
  }

  if (is.null(support.signal)) {
    if(any(support.signal.vector <= 0)) {
      stop('argument `support.signal.vector` must be an positive number or a vector of positive numbers',
           call. = FALSE)
      }
    }

  if (is.null(errormeasure) ||
      (!errormeasure %in% c("RMSE","MSE", "MAE",
                                    "MAPE", "sMAPE", "MASE")))
    stop('argument `errormeasure` must be one of
         c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE")',
         call. = FALSE)

  if (is.null(errormeasure.which) ||
      (!errormeasure.which %in% c("min", "1se", "elbow")))
    stop('argument `errormeasure.which` must be one of
         c("min", "1se", "elbow")',
         call. = FALSE)
  if (errormeasure.which == "1se" & !cv)
    stop('argument `errormeasure.which` must be one of
         c("min", "elbow") when cv = FALSE',
         call. = FALSE)

  logical.param <- c(
    model = model,
    x = x,
    y = y,
    cv = cv,
    OLS = OLS,
    support.method.ridge.penalize.intercept = support.method.ridge.penalize.intercept
  )

  if (!is.logical(logical.param) |
      length(logical.param) != 6) {
    stop(
      "The following arguments must be TRUE or FALSE:
          `model`, `x`, `y`, `cv`, `OLS`, `support.method.ridge.penalize.intercept`",
      call. = FALSE
    )
  }

  if (is.null(verbose) || !(verbose %in% c(0, 1, 2, 3))) {
    stop("argument `verbose` must be one of 0, 1, 2 or 3.",
         call. = FALSE)
  }

  if (is.null(cv.nfolds) || cv.nfolds < 3 || cv.nfolds %% 1 != 0)
    stop("argument `cv.nfolds` must be an integer equal or greater than 3.",
         call. = FALSE)

  if (is.null(cv.repeats) || cv.repeats < 1 || cv.repeats %% 1 != 0)
    stop("argument `cv.repeats` must be a positive integer.", call. = FALSE)

  if (!is.numeric(support.signal.points)) {
    stop("argument `support.signal.points` must be a numeric vector.", call. = FALSE)
  } else if (length(support.signal.points) != 1) {
    if (is.vector(support.signal.points) && !all.equal(sum(support.signal.points), 1)) {
    stop("the sum of the elements of argument `support.signal.points` must be equal to 1.",
         call. = FALSE)
    } else {
      if (is.matrix(support.signal.points) &&
          !all(apply(matrix(apply(support.signal.points, 1, sum), ncol = 1), 1, all.equal, 1))
      ) {
        stop("the sum of the elements by row of argument `support.signal.points` must be
             equal to 1.",
             call. = FALSE)
    }
    }
  }

  if (!is.numeric(support.noise.points)) {
    stop("argument `support.noise.points` must be a numeric vector.", call. = FALSE)
  } else if (length(support.noise.points) != 1) {
    if (is.vector(support.noise.points) && !all.equal(sum(support.noise.points), 1)) {
      stop("the sum of the elements of argument `support.noise.points` must be equal to 1.",
           call. = FALSE)
    } else {
      if (is.matrix(support.noise.points) &&
          !all(apply(matrix(apply(support.noise.points, 1, sum), ncol = 1), 1, all.equal, 1))
          ) {
        stop("the sum of the elements by row of argument `support.noise.points` must be
             equal to 1.",
             call. = FALSE)
      }
    }
  }

  if (is.null(caseGLM) ||
      (!caseGLM %in% c("D", "M", "NM")))
    stop('argument `caseGLM` must be one of
         c("D", "M", "NM")',
         call. = FALSE)

  if (!is.numeric(boot.B) || boot.B < 0 || (boot.B > 0 & boot.B < 10))
    stop("argument `boot.B` must be 0 or a positive integer greater or equal to 10.",
         call. = FALSE)

  if (is.null(boot.method) ||
      (!boot.method %in% c("residuals", "cases", "wild")))
    stop('argument `boot.method` must be one of c("residuals", "cases", "wild")',
         call. = FALSE)

}
