#' Cross-validation for \code{\link{lmgce}}
#'
#' Performs k-fold cross-validation for some of the \code{\link{lmgce}}
#' parameters.
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
#' @param model Boolean value. if \code{TRUE}, the model frame used is returned.
#' The default is \code{model = TRUE}.
#' @param x Boolean value. if \code{TRUE}, the model matrix used is returned.
#' The default is \code{x = FALSE}.
#' @param y Boolean value. if \code{TRUE}, the response used is returned.
#' The default is \code{y = FALSE}.
#' @param cv Boolean value. If \code{TRUE} the error, \code{errormeasure},
#' will be computed using cross-validation. If \code{FALSE} the error will be
#' computed in sample. The default is \code{cv = TRUE}.
#' @param cv.nfolds number of folds used for cross-validation when
#' \code{cv = TRUE}. The default is \code{cv.nfolds = 5} and the smallest value
#' allowable is \code{cv.nfolds = 3}.
#' @param errormeasure Loss function (error) to be used for the selection
#' of the support spaces. One of c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE").
#' The default is \code{errormeasure = "RMSE"}.
#' @param errormeasure.which Which value of \code{errormeasure}
#' to be used for selecting a support space upper limit from \code{support.signal.vector}.
#' One of \code{c("min", "1se", "elbow")} where \code{"min"} corresponds to the
#' support spaces that produced the lowest error, \code{"1se"} corresponds to
#' the support spaces such that error is within 1 standard error of the CV error
#' for \code{"min"} and \code{"elbow"} corresponds to the elbow point of the error
#' curve (the point that maximizes the distance between each observation, i.e,
#' the pair composed by the upper limit of the support space and the error, and
#' the line between the first and last observations, i.e., the lowest and the
#' highest upper limits of the support space respectively. See
#' \code{\link[pathviewr]{find_curve_elbow}}). The default is
#' \code{errormeasure.which = "1se"}.
#' @param support.method One of c("standardized", "ridge"). If
#' \code{support.method = "standardized}, the default, standardized coefficients
#' are used to define the signal support spaces. If
#' \code{support.method = "ridge} the signal support spaces are define by the
#' ridge trace.
#' @param support.method.ridge.lambda Ridge parameter. The default is
#' \code{support.method.ridge.lambda = NULL} and a lambda base
#' \code{support.method.ridge.base} logarithmic sequence will be computed based
#' on \code{support.method.ridge.lambda.n}, \code{support.method.ridge.lambda.min}
#'  and \code{support.method.ridge.lambda.max}. Supplying a lambda sequence
#' overrides this. To be used when \code{support.method = "ridge"}.
#' @param support.method.ridge.lambda.base Value for the base of logarithmic
#' sequence of ridge parameters. The default is
#' \code{support.method.ridge.lambda.base = 10}. To be used when
#' \code{support.method = "ridge"} and \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.lambda.min Minimum value for the
#' \code{support.method.ridge.lambda} sequence. The default is
#' \code{support.method.ridge.lambda.min = 10^-3}. To be used when
#' \code{support.method = "ridge"} and \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.lambda.max Maximum value for the
#' \code{support.method.ridge.lambda} sequence. The default is
#' \code{support.method.ridge.lambda.max = 10^3}. To be used when
#' \code{support.method = "ridge"} and \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.lambda.n The number of ridge parameters values.
#' The default is \code{support.method.ridge.lambda.n = 100}. To be used when
#'  \code{support.method = "ridge"} and \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.standardize Boolean value. If \code{TRUE}, the
#' default, then: i) centering is done by subtracting the column means of x and
#' y from their corresponding columns; ii) scaling is done by dividing the
#' (centered) columns of x and y by their standard deviations. To be used when
#' \code{support.method = "ridge"}.
#' @param support.method.ridge.penalize.intercept Boolean value. if \code{TRUE},
#' the default, the intercept will be penalized. To be used when
#' \code{support.method = "ridge"} and
#' \code{support.method.ridge.standardize = FALSE}.
#' @param support.method.ridge.symm Boolean value. If \code{TRUE}, the default,
#' signal supports will be symmetrical and the upper limit will be the maximum
#'  absolute values of the estimated ridge coefficients for
#'  \code{support.method.ridge.lambda}. If \code{FALSE}, the lower and upper
#'  limits will be, respectively, the minimum and maximum values of the
#'  estimated ridge coefficients.
#' @param support.method.ridge.maxresid Boolean value. if \code{TRUE}, the
#' default, noise supports will symmetrical and the upper limit will be the
#' maximum absolute value of the residuals of ridge estimation for
#'  \code{support.method.ridge.lambda}. If \code{FALSE} limits are computed
#'  using the empirical three-sigma rule (Pukelsheim (1994)).
#' @param support.signal \code{NULL} or fixed positive upper limit (L) for the
#' support spaces (-L,L) on standardized data (when
#' \code{support.method = "standardized"}); \code{NULL} or fixed positive factor
#'  to be multiplied by the maximum absolute value of the ridge trace for each
#'  coefficient (when \code{support.method = "ridge"}); a pair (LL,UL) or a
#'  matrix ((k+1) x 2) for the support spaces on original data. The default is
#'  \code{support.signal = NULL}.
#' @param support.signal.vector NULL or a vector of positive values when
#' \code{support.signal = NULL}. If \code{support.signal.vector = NULL},
#' the default, a vector
#' \code{c(support.signal.vector.min,...,support.signal.vector.max)} of dimension
#'  \code{support.signal.vector.n} and logarithmically equally spaced will be
#' generated. Each value represents the upper limits for the standardized support
#'  spaces, when \code{support.method = "standardized"} or the factor to be
#'  multiplied by the maximum absolute value of the ridge trace for each
#'  coefficient, when \code{support.method = "ridge"}.
#' @param support.signal.vector.min A positive value for the lowest limit of the
#' \code{support.signal.vector} when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.min = 0.3}.
#' @param support.signal.vector.max A positive value for the highest limit of the
#' \code{support.signal.vector} when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.max = 20}.
#' @param support.signal.vector.n A positive integer for the number of support
#' spaces to be used when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.n = 20}.
#' @param support.signal.points A vector of positive integers defining the number
#'  of points for the signal support to be tested .The default is
#' \code{support.signal.points = c(3, 5, 7, 9)}.
#' @param support.noise An interval, preferably centered around zero, given in the form
#' \code{c(LL,UL)}. If \code{support.noise = NULL}, the default, then a vector
#' \code{c(-L,L)} is computed using the empirical three-sigma rule
#' Pukelsheim (1994).
#' @param support.noise.points A vector of positive integers defining the number
#'  of points for the noise support to be tested .The default is
#' \code{support.noise.points = c(3, 5, 7, 9)}.
#' @param weight a vector of values between zero and one representing the
#' prediction-precision loss trade-off. The default is
#' \code{weight = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)}.
#' @param twosteps.n Number of GCE reestimations using a previously estimated
#' vector of signal probabilities.
#' @param method  Use \code{"primal.solnl"} (GCE using Sequential Quadratic
#' Programming (SQP) method; see \code{\link[NlcOptim]{solnl}}) or
#' \code{"primal.solnp"} (GCE using the augmented Lagrange multiplier method
#' with an SQP interior algorithm; see \code{\link[Rsolnp]{solnp}}) for primal
#' form of the optimization problem and \code{"dual"} (GME), \code{"dual.CG"}
#' (GCE using a conjugate gradients method; see \code{\link[stats]{optim}}),
#' \code{"dual.BFGS"} (GCE using Broyden-Fletcher-Goldfarb-Shanno quasi-Newton
#' method; see \code{\link[stats]{optim}}), \code{"dual.L-BFGS-B"} (GCE using a
#' box-constrained optimization with limited-memory modification of the BFGS
#' quasi-Newton method; see \code{\link[stats]{optim}}), \code{dual.Rcgmin}
#' (GCE using an update of the conjugate gradient algorithm; see
#' \code{\link[optimx]{optimx}}),
#' \code{dual.bobyqa} (GCE using a derivative-free optimization by quadratic
#' approximation; see \code{\link[optimx]{optimx}} and
#' \code{\link[minqa]{bobyqa}}), \code{dual.newuoa} (GCE using a
#' derivative-free optimization by quadratic approximation; see
#' \code{\link[optimx]{optimx}} and \code{\link[minqa]{newuoa}}),
#' \code{dual.nlminb} (GCE; see \code{\link[optimx]{optimx}} and
#' \code{\link[stats]{nlminb}}), \code{dual.nlm} (GCE; see
#' \code{\link[optimx]{optimx}} and \code{\link[stats]{nlm}}),
#' \code{dual.lbfgs} (GCE using the Limited-memory
#' Broyden-Fletcher-Goldfarb-Shanno; see \code{\link[lbfgs]{lbfgs}}),
#' \code{dual.lbfgsb3c} (GCE using L-BFSC-B implemented in Fortran code and with
#' an Rcpp interface; see \code{\link[lbfgsb3c]{lbfgsb3c}}) or
#' \code{dual.optimParallel} (GCE using parallel version of the L-BFGS-B; see
#' \code{\link[optimParallel]{optimParallel}}) for dual form. The
#' default is \code{method = "dual.BFGS"}.
#' @param caseGLM special cases of the generic general linear model. One of
#' \code{c("D", "M", "NM")}, where "D" stands for data, "M" for moment and
#'  "NM" for normed-moment The default is
#' \code{caseGLM = "D"}.
#' @param boot.B A single positive integer greater or equal to 10 for the number
#' of bootstrap replicates to be used for the computation of the bootstrap
#' confidence interval(s). Zero value will generate no replicate. The default
#' is \code{boot.B = 0}.
#' @param boot.method Method to be use for bootstrapping. One of
#' \code{c("residuals", "cases", "wild")} which corresponds to resampling on
#' residuals, on individual cases or on residuals multiplied by a N(0,1) variable,
#' respectively. The default is \code{boot.method = "residuals"}.
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = 230676}.
#' @param OLS Boolean value. if \code{TRUE}, the default, OLS estimation is
#' performed.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @details
#'
#' The \code{cv.lmgce} function fits several linear regression models via
#'  generalized cross according to the defined arguments. In particular,
#'  \code{support.signal.points}, \code{support.noise.points} and
#'  \code{weight} can be defined as vectors.
#'
#' @return
#' \code{cv.lmgce} returns an object of \code{\link[base]{class}} \code{cv.lmgce}.
#'
#'  An object of \code{\link[base]{class}} \code{cv.lmgce} is a list containing at
#'  least the following components:
#'
#' \item{results}{a \eqn{C \times 8} \code{data.frame}, where C is the number of
#' combinations of the arguments \code{support.signal.points},
#' \code{support.noise.points} and \code{weight}. Contains information about the
#' arguments, error, convergence of the optimization method and time of
#' computation.}
#' \item{best}{a \code{\link{lmgce}} object obtained with the combination of
#' arguments that produced the lowest cross-validation error.}
#' \item{support.signal.points}{a vector of the \code{support.signal.points}
#' tested.}
#' \item{support.signal.points.best}{the value of \code{support.signal.points}
#' that produced the lowest cross-validation error.}
#' \item{support.noise.points}{a vector of the \code{support.noise.points}
#' tested.}
#' \item{support.noise.points.best}{the value of \code{support.noise.points}
#' that produced the lowest cross-validation error.}
#' \item{weight}{a vector of the \code{weight} tested.}
#' \item{weight.best}{the value of \code{weight} that produced the lowest
#' cross-validation error.}
#' @param coef A vector of the true coefficients, when available.
#'
#' @seealso
#' See the generic functions \code{\link{plot.cv.lmgce}},
#' \code{\link{print.cv.lmgce}} and \code{\link{coef.cv.lmgce}}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @references
#' Golan, A., Judge, G. G. and Miller, D. (1996)
#' \emph{Maximum entropy econometrics : robust estimation with limited data.}
#' Wiley.\cr
#' Golan, A. (2008).
#' \emph{Information and Entropy Econometrics — A Review and Synthesis.}
#' Foundations and Trends® in Econometrics, 2(1–2), 1–145.
#' \doi{10.1561/0800000004}\cr
#' Golan, A. (2017)
#' \emph{Foundations of Info-Metrics: Modeling, Inference, and Imperfect Information (Vol. 1).}
#' Oxford University Press.
#' \doi{10.1093/oso/9780199349524.001.0001}\cr
#' Pukelsheim, F. (1994)
#' \emph{The Three Sigma Rule.}
#' The American Statistician, 48(2), 88–91.
#' \doi{10.2307/2684253}
#'
#' @examples
#' \donttest{
#' res.cv.lmgce <-
#'   cv.lmgce(y ~ .,
#'            data = dataGCE)
#'
#' res.cv.lmgce
#' }
#'
#' @export

cv.lmgce <- function(formula,
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
                     errormeasure = c("RMSE", "MSE", "MAE", "MAPE", "sMAPE", "MASE"),
                     errormeasure.which =
                       {
                         if (isTRUE(cv))
                           c("1se", "min", "elbow")
                         else
                           c("min", "elbow")
                       },
                     support.method = c("standardized", "ridge"),
                     support.method.ridge.lambda = NULL,
                     support.method.ridge.base = 10,
                     support.method.ridge.lambda.min = 10^-3,
                     support.method.ridge.lambda.max = 10^3,
                     support.method.ridge.lambda.n = 100,
                     support.method.ridge.standardize = TRUE,
                     support.method.ridge.penalize.intercept = TRUE,
                     support.method.ridge.symm = TRUE,
                     support.method.ridge.maxresid = TRUE,
                     support.signal = NULL,
                     support.signal.vector = NULL,
                     support.signal.vector.min = 0.3,
                     support.signal.vector.max = 20,
                     support.signal.vector.n = 20,
                     support.signal.points = c(3, 5, 7, 9),
                     support.noise = NULL,
                     support.noise.points = c(3, 5, 7, 9),
                     weight = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                     twosteps.n = 1,
                     method = c(
                       "dual.lbfgsb3c",
                       "dual.BFGS",
                       "dual",
                       "primal.solnl",
                       "primal.solnp",
                       "dual.CG",
                       "dual.L-BFGS-B",
                       "dual.Rcgmin",
                       "dual.bobyqa",
                       "dual.newuoa",
                       "dual.nlminb",
                       "dual.nlm",
                       "dual.lbfgs",
                       "dual.optimParallel"
                     ),
                     caseGLM = c("D", "M", "NM"),
                     boot.B = 0,
                     boot.method = c("residuals", "cases", "wild"),
                     seed = 230676,
                     OLS = TRUE,
                     verbose = 0,
                     coef = NULL) {

  support.method <- match.arg(support.method)
  errormeasure <- match.arg(errormeasure)
  errormeasure.which <- match.arg(errormeasure.which)
  method <- match.arg(method)
  caseGLM <- match.arg(caseGLM)
  boot.method <- match.arg(boot.method)

  cv.repeats <- 1

  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  if (is.empty.model(mt))
    stop("Empty model specified.")

  y <- model.response(mf, "numeric")

  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }

  X <- model.matrix(mt, mf, contrasts)

  if (is.null(n <- nrow(X)))
    stop("'x' must be a matrix")

  if(n == 0L)
    stop("0 (non-NA) cases")

  if (ncol(X) == 0L)
    stop("NULL model")

  if (support.method == "standardized") {
  if (is.null(support.signal) || length((support.signal) == 1)){
    if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
      if (any(as.numeric(apply(cbind(X[, -1], y), 2, sd)) == 0)) {
        stop("standardization not possible because some predictors are constants.",
             call. = FALSE)
      }
    } else {
      if (any(as.numeric(apply(cbind(X, y), 2, sd)) == 0)) {
        stop("standardization not possible because some predictors are constants.",
             call. = FALSE)
      }
    }
  }
  }

  if(!is.null(offset))
    y <- y - offset

  if (NROW(y) != n)
    stop("incompatible dimensions")

  max.abs.coef <- NULL
  if (support.method == "ridge") {
    max.abs.coef <-
      ridgetrace.Xy(
        X,
        y,
        lambda = support.method.ridge.lambda,
        lambda.base = support.method.ridge.lambda.base,
        lambda.min = support.method.ridge.lambda.min,
        lambda.max = support.method.ridge.lambda.max,
        lambda.n = support.method.ridge.lambda.n,
        standardize = support.method.ridge.standardize,
        penalize.intercept = support.method.ridge.penalize.intercept,
        cv = FALSE)$max.abs.coef
  }

  res <- list(best = list(),
              results = NULL,
              support.signal.points = support.signal.points,
              support.signal.points.best = NULL,
              support.noise.points = support.noise.points,
              support.noise.points.best = NULL,
              weight = weight,
              weight.best = NULL)

  parameter_grid <- expand.grid(
    support.signal.points = support.signal.points,
    support.noise.points = support.noise.points,
    weight = weight
  )

  parameter_grid$error.measure <- NA
  parameter_grid$error.measure.cv.mean <- NA
  parameter_grid$error.measure.cv.mean.beta <- NA
  parameter_grid$convergence <- NA
  parameter_grid$time <- NA

  for (i in 1:nrow(parameter_grid)) {
    if (verbose >= 1L) {
      cat("\n", round(i/nrow(parameter_grid) * 100, 3), "%", sep = "")
      cat("\nEstimation\n", sep = "")
    }

      start.time <- Sys.time()

      res.aux <-
        lmgce.assign.noci(
          y,
          X,
          offset,
          cv,
          cv.nfolds,
          cv.repeats,
          errormeasure,
          errormeasure.which,
          max.abs.coef,
          support.signal,
          support.signal.vector,
          support.signal.vector.min,
          support.signal.vector.max,
          support.signal.vector.n,
          support.signal.points = parameter_grid[i, 1],
          support.noise,
          support.noise.points = parameter_grid[i, 2],
          weight = parameter_grid[i, 3],
          method,
          caseGLM,
          seed,
          verbose
        )

    res.aux$results$twosteps <- list()

    if (twosteps.n > 0) {
      if (verbose >= 1L)
        cat("\n\nTwo steps\n", sep = "")

      for (ts in 1:twosteps.n) {
        res.aux$results$twosteps[[ts]] <-
          lmgce.assign.noci(
            y,
            X,
            offset,
            cv,
            cv.nfolds,
            cv.repeats,
            errormeasure,
            errormeasure.which,
            max.abs.coef,
            {if (ts == 1)
              res.aux$support.matrix
              else
                res.aux$results$twosteps[[ts - 1]]$support.matrix},
            support.signal.vector,
            support.signal.vector.min,
            support.signal.vector.max,
            support.signal.vector.n,
            {if (ts == 1)
              res.aux$p
              else
                res.aux$results$twosteps[[ts - 1]]$p},
            support.noise,
            support.noise.points = parameter_grid[i, 2],
            weight = parameter_grid[i, 3],
            method,
            caseGLM,
            seed,
            {if (verbose == 0) verbose else 1}
          )
      }

      names(res.aux$results$twosteps) <-
        sprintf(paste0("ts_%0",floor(log10(twosteps.n)) + 2,"d"), 1:twosteps.n)


      tochange <- c(
        "coefficients",
        "residuals",
        "fitted.values",
        "nep",
        "nepk",
        "vcov",
        "error.measure",
        "p",
        "w",
        "lambda",
        "convergence",
        "p0",
        "w0",
        "error.measure.cv.mean",
        "nep.cv.mean",
        "error.measure.cv.sd",
        "nep.cv.sd")
      res.aux[tochange] <- res.aux$results$twosteps[[twosteps.n]][tochange]

    }

    parameter_grid$error.measure[i] <- res.aux$error.measure
    parameter_grid$error.measure.cv.mean[i] <- res.aux$error.measure.cv.mean
    if (!is.null(coef)) {
    parameter_grid$error.measure.cv.mean.beta[i] <-
      mean(sapply(res.aux$results$cv$repeats1,
                  function(x){accmeasure(x$coefficients, coef, errormeasure)}))}
    parameter_grid$convergence[i] <- res.aux$convergence
    parameter_grid$time[i] <-
      difftime(Sys.time(),
               start.time,
               units = "secs")

  if (i == 1) {
    res$best <- res.aux
    res$best$weight <- parameter_grid[i, 3]
    res$weight.best <- res$best$weight
    res$support.signal.points.best <- parameter_grid[i, 1]
    res$support.noise.points.best <- parameter_grid[i, 2]
  } else
    if (res.aux$error.measure.cv.mean < res$best$error.measure.cv.mean) {
      res$best <- res.aux
      res$best$weight <- parameter_grid[i, 3]
      res$weight.best <- res$best$weight
      res$support.signal.points.best <- parameter_grid[i, 1]
      res$support.noise.points.best <- parameter_grid[i, 2]
    }
  }

  res$results <- parameter_grid

  res$best$results$bootstrap <-
    list(
      coefficients = NULL,
      nep = NULL,
      nepk = NULL,
      attempts = NULL)

  boot.B <- as.integer(boot.B)

  if (boot.B != 0L) {
    if (verbose >= 1L)
      cat("\n\nConfidence intervals\n", sep = "")
    res$best <-
      lmgce.assign.ci(
        y = y,
        X = X,
        offset = offset,
        errormeasure = errormeasure,
        support.signal = res$best$support.matrix,
        support.signal.points = res$best$p0,
        support.noise = support.noise,
        support.noise.points = res$best$w0,
        weight = res$best$weight,
        method = method,
        caseGLM = caseGLM,
        seed = seed,
        res$best,
        boot.B = boot.B,
        boot.method = boot.method,
        verbose = verbose
      )
  }

  # OLS ##

  res$best$results$OLS <- list(res = NULL,
                          matrix.coef = NULL,
                          error = NULL)

  if (isTRUE(OLS)) {
    res.OLS <-
      lm(formula,
         data[row.names(data) %in% row.names(X), ]
         # , subset,
         # na.action,
         # offset,
         # contrasts
      ) ### CAREFUL ##

    if (isTRUE(cv)) {

      coef.OLS <- matrix(NA,
                         ncol = cv.repeats * cv.nfolds,
                         nrow = ncol(X))
      error.OLS <- NA

      for (r in 1:cv.repeats) {

        if (!is.null(seed)) set.seed(seed)

        auxfolds = cut(seq(1, nrow(X)),
                       breaks = cv.nfolds,
                       labels = FALSE)
        change_order <- sample(nrow(X))

        for (cv.n in 1:cv.nfolds) {

          res.OLS.cv <-
            lm(as.formula(formula),
               data[row.names(data) %in% row.names(X), ][change_order, ][auxfolds != cv.n, ]
               # , subset,
               # na.action,
               # offset,
               # contrasts
            ) ### CAREFUL ##

          coef.OLS[, (r - 1) * cv.nfolds + cv.n] <- coef(res.OLS.cv)
          error.OLS[(r - 1) * cv.nfolds + cv.n] <-
            accmeasure(predict(res.OLS.cv,
                               data[row.names(data) %in% row.names(X), ][change_order, ][auxfolds == cv.n, ]),
                       y[change_order][auxfolds == cv.n],
                       errormeasure)
        }
      }
    } else {
      coef.OLS <- coef(res.OLS)
      error.OLS <- accmeasure(fitted(res.OLS), y, errormeasure)
    }

    res$best$results$OLS <- list(res = res.OLS,
                            matrix.coef = coef.OLS,
                            error = error.OLS)
  }

  res$best$support.method <- support.method
  res$best$twosteps.n <- twosteps.n
  res$best$boot.B <- boot.B
  res$best$boot.method <- boot.method
  res$best$caseGLM <- caseGLM
  res$best$error <- errormeasure
  res$best$error.which <- {if (is.null(support.signal)) errormeasure.which else NULL}
  res$best$df.residual <- nrow(X) - ncol(X)
  res$best$offset <- offset
  res$best$na.action <- attr(mf, "na.action")
  res$best$contrasts <- attr(X, "contrasts")
  res$best$xlevels <- .getXlevels(mt, mf)
  res$best$call <- cl
  res$best$terms <- mt
  if (model)
    res$best$model <- mf
  if (ret.x)
    res$best$x <- X
  if (ret.y)
    res$best$y <- y

  res$best <- res$best[order(names(res$best))]

  class(res$best) <- "lmgce"
  class(res) <- "cv.lmgce"

  return(res)
}
