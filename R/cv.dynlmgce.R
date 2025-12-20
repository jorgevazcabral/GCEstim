#' Cross-validation for \code{\link{dynlmgce}}
#'
#' Performs k-fold cross-validation for some of the \code{\link{dynlmgce}}
#' parameters.
#'
#' @param formula a "formula" describing the linear model to be fit. For details
#' see \code{\link[stats]{lm}} and \code{\link[dynlm]{dynlm}}.
#' @param data A \code{\link[base]{data.frame}} (or object coercible by
#' \code{\link[base]{as.data.frame}} to a data frame) or time series object
#' (e.g., \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}), containing the
#' variables in the model.
#' @param start The time of the first observation. Either a single number
#'  or a vector of two numbers (the second of which is an integer), which
#'  specify a natural time unit and a (1-based) number of samples into the time
#'  unit (see \code{\link[stats]{ts}}).
#' @param end The time of the last observation, specified in the same way as
#' \code{start} (see \code{\link[stats]{ts}}).
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
#' of the support spaces. One of
#' c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE"). The default is
#' \code{errormeasure = "RMSE"}.
#' @param errormeasure.which Which value of \code{errormeasure} to be used for
#' selecting a support space upper limit from \code{support.signal.vector}.
#' One of \code{c("min", "1se", "elbow")} where \code{"min"} corresponds to the
#' support spaces that produced the lowest error, \code{"1se"} corresponds to
#' the support spaces such that error is within 1 standard error of the CV error
#' for \code{"min"} and \code{"elbow"} corresponds to the elbow point of the
#' error curve (the point that maximizes the distance between each observation,
#' i.e, the pair composed by the upper limit of the support space and the error,
#' and the line between the first and last observations, i.e., the lowest and
#' the highest upper limits of the support space respectively. See
#' \code{\link[pathviewr]{find_curve_elbow}}). The default is
#' \code{errormeasure.which = "1se"}.
#' @param support.method One of c("standardized", "ridge"). If
#' \code{support.method = "standardized}, the default, standardized coefficients
#' are used to define the signal support spaces. If
#' \code{support.method = "ridge} the signal support spaces are define by the
#' ridge trace.
#' @param support.method.ridge.lambda Ridge parameter. The default is
#' \code{support.method.ridge.lambda = NULL} and a lambda logarithmic sequence
#' will be computed based on \code{support.method.ridge.lambda.n},
#' \code{support.method.ridge.lambda.min} and
#' \code{support.method.ridge.lambda.max}. Supplying a lambda sequence
#' overrides this. To be used when \code{support.method = "ridge"}.
#' @param support.method.ridge.lambda.min Minimum value for the
#' \code{support.method.ridge.lambda} sequence. The default is
#' \code{support.method.ridge.lambda.min = 10^-3}. To be used when
#' \code{support.method = "ridge"} and
#' \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.lambda.max Maximum value for the
#' \code{support.method.ridge.lambda} sequence. The default is
#' \code{support.method.ridge.lambda.max = 10^3}. To be used when
#' \code{support.method = "ridge"} and
#' \code{support.method.ridge.lambda = NULL}.
#' @param support.method.ridge.lambda.n The number of ridge parameters values.
#' The default is \code{support.method.ridge.lambda.n = 100}. To be used when
#'  \code{support.method = "ridge"} and
#'  \code{support.method.ridge.lambda = NULL}.
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
#' \code{c(support.signal.vector.min,...,support.signal.vector.max)} of
#' dimension \code{support.signal.vector.n} and logarithmically equally spaced
#' will be generated. Each value represents the upper limits for the
#' standardized support spaces, when \code{support.method = "standardized"} or
#' the factor to be multiplied by the maximum absolute value of the ridge trace
#' for each coefficient, when \code{support.method = "ridge"}.
#' @param support.signal.vector.min A positive value for the lowest limit of the
#' \code{support.signal.vector} when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.min = 0.3}.
#' @param support.signal.vector.max A positive value for the highest limit of
#' the \code{support.signal.vector} when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.max = 20}.
#' @param support.signal.vector.n A positive integer for the number of support
#' spaces to be used when \code{support.signal = NULL} and
#' \code{support.signal.vector = NULL}. The default is
#' \code{support.signal.vector.n = 20}.
#' @param support.signal.points A vector of positive integers defining the
#' number of points for the signal support to be tested .The default is
#' \code{support.signal.points = c(3, 5, 7, 9)}.
#' @param support.noise An interval, preferably centered around zero, given in
#' the form \code{c(LL,UL)}. If \code{support.noise = NULL}, the default, then a
#'  vector \code{c(-L,L)} is computed using the empirical three-sigma rule
#' Pukelsheim (1994).
#' @param support.noise.points A vector of positive integers defining the number
#'  of points for the noise support to be tested. The default is
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
#' residuals, on individual cases or on residuals multiplied by a N(0,1)
#' variable, respectively. The default is \code{boot.method = "residuals"}.
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
#' The \code{cv.dynlmgce} function fits several dynamic linear regression
#' models via generalized cross according to the defined arguments.
#' In particular, \code{support.signal.points}, \code{support.noise.points} and
#'  \code{weight} can be defined as vectors.
#'
#' @return
#' \code{cv.dynlmgce} returns an object of \code{\link[base]{class}}
#' \code{cv.lmgce} containing at least the following components:
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
#' \emph{Foundations of Info-Metrics: Modeling, Inference, and Imperfect
#' Information (Vol. 1).}
#' Oxford University Press.
#' \doi{10.1093/oso/9780199349524.001.0001}\cr
#' Pukelsheim, F. (1994)
#' \emph{The Three Sigma Rule.}
#' The American Statistician, 48(2), 88–91.
#' \doi{10.2307/2684253}
#'
#' @examples
#' \donttest{
#' res.cv.dynlmgce <-
#'   cv.dynlmgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' res.cv.dynlmgce
#' }
#'
#' @importFrom stats density ts
#' @importFrom utils head tail
#' @importFrom zoo merge.zoo
#' @export

cv.dynlmgce <- function(formula,
                     data,
                     subset,
                     na.action,
                     offset,
                     contrasts = NULL,
                     start = NULL,
                     end = NULL,
                     cv = TRUE,
                     cv.nfolds = 5,
                     errormeasure = c("RMSE", "MSE", "MAE",
                                      "MAPE", "sMAPE", "MASE"),
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
                       "dual.BFGS",
                       "dual.lbfgsb3c",
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

  cv.repeats <- 1

  if(!is.ts(data))
    stop("argument `data` must be a \"ts\" object")

  #suppressMessages(require(zoo, quietly = TRUE))

  x = TRUE
  y = TRUE
  model = TRUE
  support.method <- match.arg(support.method)
  errormeasure <- match.arg(errormeasure)
  errormeasure.which <- match.arg(errormeasure.which)
  method <- match.arg(method)
  caseGLM <- match.arg(caseGLM)
  boot.method <- match.arg(boot.method)

  Zenv <- new.env(parent = environment(formula))
  assign("dynformula",
         function(x) structure(x, class = unique(c("dynformula",
                                                   oldClass(x)))), envir = Zenv)
  assign("L", function(x, k = 1) {
    if (length(k) > 1) {
      res <- lapply(k, function(i) lag(x, k = -i))
      res <- if (inherits(x, "ts"))
        do.call("ts.intersect", res)
      else do.call("merge", c(res, list(all = FALSE)))
      colnames(res) <- k
    }
    else {
      res <- lag(x, k = -k)
    }
    return(res)
  }, envir = Zenv)
  assign("d", function(x, lag = 1) diff(x, lag = lag), envir = Zenv)
  assign("season", function(x, ref = NULL) {
    freq <- frequency(x)
    stopifnot(freq > 1 && identical(all.equal(freq, round(freq)),
                                    TRUE))
    freq <- ofreq <- round(freq)
    freq <- if (freq == 12)
      month.abb
    else if (freq == 4)
      paste("Q", 1:4, sep = "")
    else 1:freq
    res <- factor(zoo::coredata(cycle(x)), labels = freq)
    if (!is.null(ref))
      res <- relevel(res, ref = ref)
    res <- zoo::zoo(res, zoo::index(x), ofreq)
    return(res)
  }, envir = Zenv)
  assign("trend", function(x, scale = TRUE) {
    freq <- ofreq <- if (inherits(x, "ts"))
      frequency(x)
    else attr(x, "frequency")
    if (is.null(freq) | !scale)
      freq <- 1
    stopifnot(freq >= 1 && identical(all.equal(freq, round(freq)),
                                     TRUE))
    freq <- round(freq)
    res <- zoo::zoo(seq_along(zoo::index(x))/freq,
                    zoo::index(x),
                    frequency = ofreq)
    return(res)
  }, envir = Zenv)
  assign("harmon", function(x, order = 1) {
    freq <- frequency(x)
    stopifnot(freq > 1 && identical(all.equal(freq, round(freq)),
                                    TRUE))
    freq <- round(freq)
    order <- round(order)
    stopifnot(order <= freq/2)
    res <- outer(2 * pi * zoo::index(x), 1:order)
    res <- cbind(apply(res, 2, cos), apply(res, 2, sin))
    colnames(res) <- if (order == 1) {
      c("cos", "sin")
    }
    else {
      c(paste("cos", 1:order, sep = ""), paste("sin",
                                               1:order, sep = ""))
    }
    if ((2 * order) == freq)
      res <- res[, -(2 * order)]
    return(res)
  }, envir = Zenv)
  assign("model.frame.dynformula",
         function(formula,
                  data = NULL,
                  subset = NULL,
                  na.action = na.omit,
                  drop.unused.levels = FALSE,
                  xlev = NULL, ...) {
    if (is.null(data)) {
      data <- as.list(parent.frame())
      data <- data[!sapply(data, inherits, "function")]
    }
    if (!is.list(data))
      data <- as.list(data)
    args <- as.list(attr(terms(formula, data = data), "variables"))[-1]
    args$retclass <- "list"
    args$all <- FALSE
    formula <- terms(formula, data = data)
    attr(formula, "predvars") <- as.call(append(merge.zoo,
                                                args))
    attr(formula, "predvars")[[1]] <- as.name("merge.zoo")
    NextMethod("model.frame", formula = formula, ...)
  }, envir = Zenv)
  if (missing(data))
    data <- Zenv
  orig.class <- if (is.data.frame(data) || is.environment(data))
    class(eval(attr(terms(formula, data = data), "variables")[[2]],
               data, Zenv))
  else class(data)

  ret.x <- x
  ret.y <- y

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]],
                                            as.name("|"))) {
    twostage <- TRUE
    ff <- formula
    mf$formula[[3]][1] <- call("+")
    ff1 <- . ~ .
    ff2 <- ~.
    ff1[[2]] <- ff[[2]]
    ff1[[3]] <- ff[[3]][[2]]
    ff2[[3]] <- ff[[3]][[3]]
    ff2[[2]] <- NULL
  }
  else {
    twostage <- FALSE
  }
  mf[[1]] <- as.name("model.frame")
  mf[[2]] <- as.call(list(as.name("dynformula"), mf[[2]]))
  mf[[2]] <- eval(mf[[2]], envir = Zenv)
  environment(mf[[2]]) <- Zenv
  mf <- eval(mf, envir = Zenv)
  mfna <- attr(mf, "na.action")
  if (length(zoo::index(mf[, 1])) > nrow(mf)) {
    for (i in 1:NCOL(mf)) {
      attr(mf[, i], "index") <-
        attr(mf[, i], "index")[-as.vector(mfna)]
      }
  }
  is.zoofactor <- function(x) !is.null(attr(x, "oclass")) &&
    attr(x, "oclass") == "factor"
  for (i in 1:NCOL(mf)) if (is.zoofactor(mf[, i]))
    mf[, i] <- zoo::coredata(mf[, i])
  mf1 <- mf[, 1]

  if (is.null(start))
    start <- start(data)

  if (length(start) > 1)
    start <- start[1] + (start[2] - 1)/frequency(mf1)

  start <- min(which(zoo::index(mf1) >= start))

  if (is.null(end))
    end <- end(data)

  if (length(end) > 1)
    end <- end[1] + (end[2] - 1)/frequency(mf1)

  end <- max(which(zoo::index(mf1) <= end))

  if (end < start) {
    warning("empty model frame specified")
    mf <- mf[0, ]
    mf1 <- mf1[0, ]
  }
  else {
    mf <- mf[start:end, , drop = FALSE]
    mf1 <- mf1[start:end]
    if (!is.null(mfna))
      attr(mf, "na.action") <-
      structure(mfna[as.vector(mfna) >= start & as.vector(mfna) <= end],
                class = class(mfna))
  }
  if ("ts" %in% orig.class && zoo::is.regular(mf1, strict = TRUE)) {
    for (i in 1:ncol(mf)) if (!is.factor(mf[, i])) {
      mf[, i] <- astszoo(mf[, i])
    }
  }
  if (all(orig.class == "numeric")) {
    for (i in 1:ncol(mf)) if (!is.factor(mf[, i]))
      mf[, i] <- as.vector(mf[, i])
  }
  rownames(mf) <- zoo::index2char(zoo::index(mf1), frequency(mf1))

  mt <- attr(mf, "terms")
  attr(mt, "predvars") <- NULL
  attr(mt, "dataClasses") <- NULL
  y <- model.response(mf, "numeric")

  offset <- model.offset(mf)
  if (!is.null(offset)) {
    offset <- as.numeric(offset)
    if (length(offset) != NROW(y))
      stop("Number of offsets is ", length(offset), ", should equal ",
           NROW(y), " (number of observations)")
  }
  if (is.empty.model(mt)) {
    x <- NULL
    res <- list(coefficients = numeric(0), residuals = y,
                fitted.values = 0 * y)
    if (!is.null(offset))
      res$fitted.values <- offset
  }
  else {
    ##### Improve ############
    if (twostage) {
      stop("two stage algorithm not implemented")
    }
    else {
      x <- model.matrix(mt, mf, contrasts)
    }

    res <- list()

    ## start lmgce

    if ("(Intercept)" %in% colnames(x)) {
      aux.x <- x[, - 1]
    } else {
      aux.x <- x
    }

    data.lmgce <-
      data.frame(y = y, aux.x)

    if ("(Intercept)" %in% colnames(x)) {
      colnames(data.lmgce)[-1] <- colnames(x)[-1]
    } else {
      colnames(data.lmgce)[-1] <- colnames(x)
    }

    res <-
      cv.lmgce(
      formula = y ~ .,
      data = data.lmgce,
      model = model,
      x = ret.x,
      y = ret.y,
      cv = cv,
      cv.nfolds = cv.nfolds,
      errormeasure = errormeasure,
      errormeasure.which = errormeasure.which,
      support.method = support.method,
      support.method.ridge.lambda = support.method.ridge.lambda,
      support.method.ridge.lambda.min = support.method.ridge.lambda.min,
      support.method.ridge.lambda.max = support.method.ridge.lambda.max,
      support.method.ridge.lambda.n = support.method.ridge.lambda.n,
      support.method.ridge.standardize = support.method.ridge.standardize,
      support.method.ridge.penalize.intercept =
        support.method.ridge.penalize.intercept,
      support.method.ridge.symm = support.method.ridge.symm,
      support.method.ridge.maxresid = support.method.ridge.maxresid,
      support.signal = support.signal,
      support.signal.vector = support.signal.vector,
      support.signal.vector.min = support.signal.vector.min,
      support.signal.vector.max = support.signal.vector.max,
      support.signal.vector.n = support.signal.vector.n,
      support.signal.points = support.signal.points,
      support.noise = support.noise,
      support.noise.points = support.noise.points,
      weight = weight,
      twosteps.n = twosteps.n,
      method = method,
      caseGLM = caseGLM,
      boot.B = boot.B,
      boot.method = boot.method,
      seed = seed,
      OLS = OLS,
      verbose = verbose
    )
    res$index <- zoo::index(mf1)
    res$frequency <- frequency(mf1)
    res$call <- cl
  }
  return(res)
}
