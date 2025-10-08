#' Generalized Cross entropy estimation
#'
#' This generic function fits a linear regression model via generalized cross
#' entropy. Initial support spaces can be provided or computed.
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
#' \code{support.method.ridge.lambda = NULL} and a lambda logarithmic sequence
#' will be computed based on \code{support.method.ridge.lambda.n},
#' \code{support.method.ridge.lambda.min} and
#' \code{support.method.ridge.lambda.max}. Supplying a lambda sequence overrides
#'  this. To be used when \code{support.method = "ridge"}.
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
#' @param support.signal.points A positive integer, a vector or a matrix. Prior
#' weights for the signal. If not a positive integer then the sum of weights by
#' row must be equal to 1. The default is
#' \code{support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5)}.
#' @param support.noise An interval, preferably centered around zero, given in
#' the form \code{c(LL,UL)}. If \code{support.noise = NULL}, the default, then a
#'  vector \code{c(-L,L)} is computed using the empirical three-sigma rule
#' (Pukelsheim (1994)).
#' @param support.noise.points A positive integer, a vector or a matrix. Prior
#' weights for the noise. If not a positive integer then the sum of weights by
#' row must be equal to 1. The default is
#' \code{support.noise.points = c(1 / 3, 1 / 3, 1 / 3)}.
#' @param weight a value between zero and one representing the
#' prediction-precision loss trade-off. If \code{weight = 0.5}, the default,
#' equal weight is placed on the signal and noise entropies. A higher than 0.5
#' value places more weight on the noise entropy whereas a lower than 0.5 value
#' places more weight on the signal entropy.
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
#'  "NM" for normed-moment The default is \code{caseGLM = "D"}.
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
#' The \code{lmgce} function fits a linear regression model via generalized cross
#' entropy. Models for \code{lmgce} are specified symbolically. A typical model has the
#' form response ~ terms where response is the (numeric) response vector and
#' terms is a series of terms which specifies a linear predictor for response.
#'  \code{lmgce} calls the lower level functions \code{lmgce.validate},
#'  \code{lmgce.assign.ci}, \code{lmgce.assign.noci}, \code{lmgce.sscv},
#'  \code{lmgce.ss}, \code{lmgce.cv} and \code{lmgce.fit}.
#'
#' @return
#' \code{lmgce} returns an object of \code{\link[base]{class}} \code{lmgce}.
#'  The function \code{\link{summary.lmgce}} is used to obtain and print a
#'  summary of the results. The generic accessory functions
#'  \code{\link{coef.lmgce}}, \code{\link{fitted.values.lmgce}},
#'  \code{\link{residuals.lmgce}} and \code{\link{df.residual.lmgce}}, extract
#'  various useful features of the value returned by \code{object} of class
#'  \code{lmgce}.
#'
#'  An object of \code{\link[base]{class}} \code{lmgce} is a list containing at
#'  least the following components:
#'
#' \item{coefficients}{a named vector of coefficients.}
#' \item{residuals}{the residuals, that is response minus fitted values.}
#' \item{fitted.values}{the fitted mean values.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{call}{the matched call.}
#' \item{terms}{the \code{\link[stats]{terms}} object used.}
#' \item{contrast}{(only where relevant) the contrasts used.}
#' \item{xlevels}{(only where relevant) a record of the levels of the factors
#' used in fitting.}
#' \item{offset}{the offset used (missing if none were used).}
#' \item{y}{if requested (the default), the response used.}
#' \item{x}{if requested (the default), the model matrix used.}
#' \item{model}{if requested (the default), the model frame used.}
#' \item{na.action}{(where relevant) information returned by
#' \code{\link[stats]{model.frame}} on the special handling of \code{NA}s.}
#' \item{boot.B}{number of bootstrap replicates used.}
#' \item{boot.method}{method used for bootstrapping.}
#' \item{caseGLM}{case of the generic general linear model used.}
#' \item{convergence}{an integer code. 0 indicates successful
#' optimization completion. Other numbers indicate different errors. See
#' \code{\link[stats]{optim}}, \code{\link[optimx]{optimx}},
#' \code{\link[NlcOptim]{solnl}}, \code{\link[Rsolnp]{solnp}},
#' \code{\link[lbfgs]{lbfgs}}) and \code{\link[lbfgsb3c]{lbfgsb3c}}).}
#' \item{error}{loss function (error) used for the selection of the support
#' spaces.}
#' \item{error.measure}{in sample error for the selected support space.}
#' \item{error.measure.cv.mean}{cross-validation mean error for the selected
#' support space.}
#' \item{error.measure.cv.sd}{standard deviation of the cross-validation error
#' for the selected support space.}
#' \item{error.which}{which criterion/standardized/factor support was used}
#' \item{support.signal.1se}{upper limit of the standardized support space or
#' factor that produced the error within one standard error from the minimum
#' error.}
#' \item{support.signal.elbow}{upper limit of the standardized support space or
#' factor that produced the error correspondent to the elbow of the error curve.}
#' \item{support.signal.min}{upper limit of the standardized support space or
#' factor that produced the minimum error.}
#' \item{p0}{vector of prior weights used for the signal.}
#' \item{p}{estimated probabilities associated with the signal.}
#' \item{w0}{vector of prior weights used for the noise.}
#' \item{w}{estimated probabilities associated with the noise.}
#' \item{lambda}{estimated Lagrange multipliers.}
#' \item{nep}{normalized entropy of the signal of the model.}
#' \item{nep.cv.mean}{cross-validation normalized entropy of the signal of
#' the model.}
#' \item{nep.cv.sd}{standard deviation of the cross-validation normalized
#' entropy of the signal of the model.}
#' \item{nepk}{normalized entropy of the signal of each coefficient.}
#' \item{results}{results from the different support spaces with or without
#' cross-validation, and from bootstrap replicates, namely number of attempts
#' (if the number of attempts is greater than three times the
#' number of bootstrap replicates the bootstrapping process stops), coefficients
#' and normalized entropies (nep - model, and nepk - coefficients), when
#' applicable; results from OLS estimation if \code{OLS = TRUE}; results from
#' GCE reestimation if \code{twosteps.n} is greater than 0.}
#' \item{support}{vector of given positive upper limits for the
#' support spaces on standardized data or factors, when
#' \code{support.signal = NULL} or \code{support.signal = L}, or
#' \code{"interval"} otherwise.}
#' \item{support.matrix}{matrix with the support spaces used for estimation on
#' original data.}
#' \item{support.method}{method chosen for the support's limits}
#' \item{support.ok}{vector of successful positive upper limits for the
#' support spaces on standardized data (\code{support.method = "standardized"})
#' or factors (\code{support.method = "ridge"}), when \code{support.signal = NULL}
#'  or \code{support.signal = L}, or \code{"interval"} otherwise.}
#' \item{support.stdUL}{when applicable, the upper limit of the standardized
#' support chosen, when \code{support.method = "standardized"} or the factor used
#'  when \code{support.method = "ridge"}.}
#' \item{vcov}{variance-covariance matrix of the coefficients.}
#'
#' @seealso
#' \code{\link{summary.lmgce}} for more detailed summaries.
#' The generic functions \code{\link{plot.lmgce}}, \code{\link{print.lmgce}},
#'  \code{\link{coef.lmgce}} and \code{\link{confint.lmgce}}.
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
#' \doi{10.2307/2684253}\cr
#' Macedo, P., Cabral, J., Afreixo, V., Macedo, F., Angelelli, M. (2025)
#' \emph{RidGME estimation and inference in ill-conditioned models.}
#' In: Gervasi O, Murgante B, Garau C, et al., eds. Computational Science and
#' Its Applications – ICCSA 2025 Workshops. Springer Nature Switzerland; 2025:300-313.
#' \doi{10.1007/978-3-031-97589-9_21}
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
#' res_gce_package
#'
#' @importFrom stats .checkMFClasses .getXlevels as.formula cycle delete.response frequency
#' is.empty.model is.ts lag lm median model.frame model.offset model.response na.omit na.pass
#' napredict naprint naresid optim pchisq pnorm printCoefmat qnorm quantile rbinom relevel
#' rnorm rt sd symnum terms time update var
#' @export

lmgce <- function(formula,
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
                  errormeasure = c("RMSE","MSE", "MAE", "MAPE", "sMAPE", "MASE"),
                  errormeasure.which =
                    {if (isTRUE(cv)) c("1se", "min", "elbow") else c("min", "elbow")},
                  support.method = c("standardized", "ridge"),
                  support.method.ridge.lambda = NULL,
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
                  support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
                  support.noise = NULL,
                  support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
                  weight = 0.5,
                  twosteps.n = 1,
                  method = c("dual.BFGS",
                             "dual.lbfgsb3c",
                             "dual", "primal.solnl", "primal.solnp",
                             "dual.CG", "dual.L-BFGS-B",
                             "dual.Rcgmin", "dual.bobyqa", "dual.newuoa",
                             "dual.nlminb", "dual.nlm",
                             "dual.lbfgs",
                             "dual.optimParallel"),
                  caseGLM = c("D", "M", "NM"),
                  boot.B = 0,
                  boot.method = c("residuals", "cases", "wild"),
                  seed = 230676,
                  OLS = TRUE,
                  verbose = 0) {

  support.method <- match.arg(support.method)
  errormeasure <- match.arg(errormeasure)
  errormeasure.which <- match.arg(errormeasure.which)
  method <- match.arg(method)
  caseGLM <- match.arg(caseGLM)
  boot.method <- match.arg(boot.method)

  cv.repeats <- 1

  lmgce.validate(
    formula,
    data,
    subset,
    na.action,
    offset,
    contrasts,
    model,
    x,
    y,
    cv,
    cv.nfolds,
    cv.repeats,
    errormeasure,
    errormeasure.which,
    support.method,
    support.method.ridge.penalize.intercept,
    support.signal,
    support.signal.vector,
    support.signal.vector.min,
    support.signal.vector.max,
    support.signal.vector.n,
    support.signal.points,
    support.noise,
    support.noise.points,
    weight,
    twosteps.n,
    method,
    caseGLM,
    boot.B,
    boot.method,
    seed,
    OLS,
    verbose
  )

  if (verbose >= 1L)
    cat("Estimation\n", sep = "")

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
        lambda.min = support.method.ridge.lambda.min,
        lambda.max = support.method.ridge.lambda.max,
        lambda.n = support.method.ridge.lambda.n,
        standardize = support.method.ridge.standardize,
        penalize.intercept = support.method.ridge.penalize.intercept,
        errormeasure = errormeasure,
        cv = FALSE)$max.abs.coef
  }
  res <-
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
      support.signal.points,
      support.noise,
      support.noise.points,
      weight,
      method,
      caseGLM,
      seed,
      verbose
    )

  res$results$twosteps <- list()

  if (twosteps.n > 0) {
  if (verbose >= 1L)
    cat("\n\nTwo steps\n", sep = "")

  for (ts in 1:twosteps.n) {
    res$results$twosteps[[ts]] <-
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
          res$support.matrix
          else
            res$results$twosteps[[ts - 1]]$support.matrix},
        support.signal.vector,
        support.signal.vector.min,
        support.signal.vector.max,
        support.signal.vector.n,
        {if (ts == 1)
          res$p
          else
            res$results$twosteps[[ts - 1]]$p},
        support.noise,
        support.noise.points,
        weight,
        method,
        caseGLM,
        seed,
        {if (verbose == 0) verbose else 1}
      )
  }

  names(res$results$twosteps) <-
      sprintf(paste0("ts_%",floor(log10(twosteps.n)) + 2,"d"), 1:twosteps.n)


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
  res[tochange] <- res$results$twosteps[[twosteps.n]][tochange]

  }

  res$results$bootstrap <-
    list(
      coefficients = NULL,
      nep = NULL,
      nepk = NULL,
      attempts = NULL)

  boot.B <- as.integer(boot.B)

  if (boot.B != 0L) {
    if (verbose >= 1L)
      cat("\n\nConfidence intervals\n", sep = "")
    res <-
      lmgce.assign.ci(
        y = y,
        X = X,
        offset = offset,
        errormeasure = errormeasure,
        support.signal = res$support.matrix,
        support.signal.points = res$p0,
        support.noise = support.noise,
        support.noise.points = support.noise.points,
        weight = weight,
        method = method,
        caseGLM = caseGLM,
        seed = seed,
        res,
        boot.B = boot.B,
        boot.method = boot.method,
        verbose = verbose
      )
  }

  # OLS ##

  res$results$OLS <- list(res = NULL,
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

  res$results$OLS <- list(res = res.OLS,
                          matrix.coef = coef.OLS,
                          error = error.OLS)
}

  #### ###

  res$support.method <- support.method
  res$twosteps.n <- twosteps.n
  res$boot.B <- boot.B
  res$boot.method <- boot.method
  res$weight <- weight
  res$caseGLM <- caseGLM
  res$error <- errormeasure
  res$error.which <- {if (is.null(support.signal)) errormeasure.which else NULL}
  res$df.residual <- nrow(X) - ncol(X)
  res$offset <- offset
  res$na.action <- attr(mf, "na.action")
  res$contrasts <- attr(X, "contrasts")
  res$xlevels <- .getXlevels(mt, mf)
  res$call <- cl
  res$terms <- mt
  if (model)
    res$model <- mf
  if (ret.x)
    res$x <- X
  if (ret.y)
    res$y <- y

  res <- res[order(names(res))]

  if (is.null(support.signal) && max(res$support.ok) == res$support.signal.min) {
    #cat("\n")
    warning("\n\nThe minimum error was found for the highest upper limit of the support. Confirm if higher values should be tested.")
  }

  if (is.null(support.signal) && cv && min(res$support.ok) == res$support.signal.1se) {
    #cat("\n")
    warning("\n\nThe 1se error was found for the lowest upper limit of the support. Confirm if lower values should be tested.")
    }

  class(res) <- "lmgce"

  return(res)

}
