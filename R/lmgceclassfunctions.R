# Generic functions for the \code{\link{lmgce}} class

#' Print a \code{\link{lmgce}} object
#'
#' Concise summary of a \code{\link{lmgce}} object
#'
#' @param x fitted \code{\link{lmgce}} object.
#' @param digits  significant digits in printout.
#' @param ... additional print arguments.
#'
#' @return A small summary of a \code{\link{lmgce}} object is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' res_gce_package
#'
#' @method print lmgce
#' @export

print.lmgce <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2L,
                  quote = FALSE)
  }
  else
    cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' Extract design matrix from \code{\link{lmgce}} object
#'
#' Returns the design matrix used to fit \code{\link{lmgce}} object.
#'
#' @param object fitted \code{\link{lmgce}} object.
#' @param ... additional arguments.
#'
#' @return A numeric matrix with one row for each observation and one column for
#'  each parameter in the model.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         x = TRUE,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' model.matrix(res_gce_package)
#' }
#'
#' @method model.matrix lmgce
#' @importFrom stats model.matrix
#' @export

model.matrix.lmgce <- function(object, ...)
{
  if (n_match <- match("x", names(object), 0L))
    object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    NextMethod("model.matrix",
               data = data,
               contrasts.arg = object$contrasts)
  }
}

#' Extract Model Formula from \code{\link{lmgce}} object
#'
#' Returns the model used to fit \code{\link{lmgce}} object.
#'
#' @param x fitted \code{\link{lmgce}} object.
#' @param ... additional arguments.
#'
#' @return An object of class \code{formula} representing the model formula.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#'
#' formula(res_gce_package)
#' }
#'
#' @method formula lmgce
#' @importFrom stats formula
#' @export

formula.lmgce <- function(x, ...)
{
  formula(x$terms)
}

#' Extract \code{\link{lmgce}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{lmgce}} object
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... Additional arguments (not used).
#'
#' @return Returns the coefficients from a \code{\link{lmgce}} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' coef(res_gce_package)
#'
#' @method coef lmgce
#' @importFrom stats coef
#' @export

coef.lmgce <- function(object, ...)
{
  object$coefficients
}

#' Extract \code{\link{lmgce}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{lmgce}} object
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... Additional arguments (not used).
#'
#' @return Returns the coefficients from a \code{\link{lmgce}} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname coefficients.lmgce
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         seed = 230676)
#' }
#' coefficients(res_gce_package)
#'
#'
#' @method coefficients lmgce
#' @importFrom stats coefficients
#' @export

coefficients.lmgce <- coef.lmgce

#' Extract \code{\link{lmgce}} Model Residuals
#'
#' \code{residuals} is a function which extracts model residuals from
#' \code{\link{lmgce}} objects.
#' The abbreviated form \code{resid} is an alias for \code{residuals}.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return Returns the residuals from a \code{\link{lmgce}} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' residuals(res_gce_package)
#'
#' @method residuals lmgce
#' @importFrom stats residuals
#' @export

residuals.lmgce <- function(object, ...)
{
  res <- object$residuals
  res <- naresid(object$na.action, res)
  res
}

#' Extract \code{\link{lmgce}} Model Residuals
#'
#'\code{resid} is a function which extracts model residuals from
#'\code{\link{lmgce}} objects.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return Returns the residuals from a \code{\link{lmgce}} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname resid.lmgce
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' resid(res_gce_package)
#'
#' @method resid lmgce
#' @importFrom stats resid
#' @export

resid.lmgce <- residuals.lmgce

#' Extract \code{\link{lmgce}} Model's Variance-Covariance Matrix
#'
#' Returns the variance-covariance matrix of the main parameters of a
#' \code{\link{lmgce}} object
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return A matrix of the estimated covariances between the parameter estimates
#'  in the linear predictor of the \code{\link{lmgce}} model.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' vcov(res_gce_package)
#'
#' @method vcov lmgce
#' @importFrom stats vcov
#' @export

vcov.lmgce <- function(object, ...)
{
  object$vcov
}

#' Variable Names of \code{\link{lmgce}} Fitted Models
#'
#' Simple utility returning variable names.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return  A character vector containing the names of the variables in the
#'   \code{\link{lmgce}} model object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' variable.names(res_gce_package)
#'
#' @method variable.names lmgce
#' @importFrom stats variable.names
#' @export

variable.names.lmgce <- function(object, ...) {
  names(object$coefficients)
}

#' Extract the Number of Observations from a \code{\link{lmgce}} model fit
#'
#' Extract the number of ‘observations’ from a \code{\link{lmgce}} model fit.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return An integer scalar representing the number of observations (rows) used
#'  in fitting the \code{\link{lmgce}} model object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' nobs(res_gce_package)
#'
#' @method nobs lmgce
#' @importFrom stats nobs
#' @export

nobs.lmgce <- function(object, ...) {
  NROW(object$residuals)
}

#' Case Names of \code{\link{lmgce}} Fitted Models
#'
#' Simple utility returning case names.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... Additional arguments (not used).
#'
#' @return A character vector containing the names or labels of the cases
#' (observations) in the \code{\link{lmgce}} model object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' case.names(res_gce_package)
#'
#' @method case.names lmgce
#' @importFrom stats case.names
#' @export

case.names.lmgce <- function(object, ...) {
  names(residuals.lmgce(object))
}

#' Calculate \code{\link{lmgce}} Fitted Values
#'
#' The fitted values for the linear model represented by a \code{\link{lmgce}}
#' object are extracted.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return Returns a vector with the fitted values for the linear model
#' represented by a \code{\link{lmgce}} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' fitted(res_gce_package)
#'
#' @method fitted lmgce
#' @importFrom stats fitted
#' @export

fitted.lmgce <- function(object, ...)
{
  xx <- object$fitted.values
  napredict(object$na.action, xx)
}
#' Calculate \code{\link{lmgce}} Fitted Values
#'
#' The fitted values for the linear model represented by a \code{\link{lmgce}}
#' object are extracted.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return Returns a vector with the fitted values for the linear model
#' represented by a \code{\link{lmgce}} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname fitted.values.lmgce
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' fitted.values(res_gce_package)
#'
#' @method fitted.values lmgce
#' @importFrom stats fitted.values
#' @export

fitted.values.lmgce <- fitted.lmgce

#' Residual Degrees-of-Freedom
#'
#' Returns the residual degrees-of-freedom extracted from a fitted model
#' \code{\link{lmgce}} object.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param ... additional arguments.
#'
#' @return The value of the residual degrees-of-freedom extracted from a
#' \code{\link{lmgce}} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' df.residual(res_gce_package)
#'
#' @method df.residual lmgce
#' @importFrom stats df.residual
#' @export

df.residual.lmgce <- function(object, ...) {
  object$df.residual
}

#' Confidence Intervals for \code{\link{lmgce}} Model Parameters and
#' Normalized Entropy
#'
#' Computes confidence intervals for one or more parameters or Normalized Entropy
#' in a \code{\link{lmgce}} fitted model.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing,
#' all parameters are considered.
#' @param level the confidence level required. The default is
#'  \code{level = 0.95}.
#' @param which One of \code{c("estimates", "NormEnt")}. The default is
#' \code{which = "estimates"}.
#' @param method method used to compute the interval. One of
#' \code{c("z","percentile", "basic")}. The default is \code{method = "z"} and
#' is only valid for the parameters.
#' @param boot.B A single positive integer greater or equal to 10 for the number
#' of bootstrap replicates for the computation of the bootstrap confidence
#' interval(s), to be used when \code{method = c("percentile", "basic")} and
#' when \code{object} was created with \code{boot.B = 0}. The default is
#' \code{boot.B = 100} when the \code{object} has no previous sampling information
#' and \code{boot.B = object$boot.B} otherwise, which corresponds to
#' the \code{boot.B} given to \code{lmgce} when the \code{object} was created.
#' @param boot.method Method used for bootstrapping. One of
#' \code{c("residuals", "cases", "wild")} which corresponds to resampling on
#' residuals, on individual cases or on residuals multiplied by a N(0,1) variable,
#' respectively. The default is \code{boot.method = object$boot.method}.
#' @param ... additional arguments.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter. These will be labelled as (1-level)/2 and
#' 1 - (1-level)/2 in percentage (by default 2.5 percent and 97.5 percent).
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' confint(res_gce_package, method = "percentile")
#'
#' confint(res_gce_package, which = "NormEnt", level = 0.99)
#'
#' confint(res_gce_package, parm = c("X005"), level = 0.99)
#'
#' @method confint lmgce
#' @importFrom stats confint
#' @export

confint.lmgce <- function(object,
                          parm,
                          level = 0.95,
                          which = c("estimates", "NormEnt"),
                          method = {if (which == "estimates") {
                            c("z","percentile", "basic")} else {c("percentile", "basic")}},
                          boot.B = ifelse(object$boot.B == 0,
                                          100,
                                          object$boot.B),
                          boot.method = object$boot.method,
                          ...)
{
  which <- match.arg(which)
  method <- match.arg(method)

  if (which == "estimates")
    cf <- coef(object) else
      cf <- NormEnt(object, model = FALSE)
  pnames <- names(cf)
  if (missing(parm)) parm <- pnames
  else if (is.numeric(parm)) parm <- pnames[parm]
  if (!all(parm %in% pnames)) {
    stop("Invalid choice of parameter")
  }
  if (!is.numeric(level)) {
    stop("Non numeric level!")
  } else {
    if (level > 1 | level < 0) {
      stop("level must be greater than 0 smaller than 1")
    }
  }
  a <- (1 - level)/2
  a <- c(a, 1 - a)

  if (method == "z") {

      aux.conf.int <-
        matrix(
          cbind(
            cf[parm] - sqrt(diag(vcov(object))[parm])*qnorm(a[2]),
            cf[parm] + sqrt(diag(vcov(object))[parm])*qnorm(a[2])
          ),
          ncol = 2
        )
      colnames(aux.conf.int) <- paste0(a*100,"%")
      rownames(aux.conf.int) <- parm
    return(aux.conf.int)
  } else {
    if (which == "estimates")
      mcf <- object$results$bootstrap$coefficients
      else
        mcf <- object$results$bootstrap$nepk
    if (is.null(mcf) ||
        boot.B != object$boot.B ||
        boot.method != object$boot.method) {
      object$results$bootstrap <-
        update(
          object,
          support.signal = object$support.matrix,
          boot.B = boot.B,
          boot.method = boot.method,
          verbose = 0
        )$results$bootstrap

      if (which == "estimates")
        mcf <- object$results$bootstrap$coefficients
      else
        mcf <- object$results$bootstrap$nepk
      }
      if (method == "basic") {
        aux.conf.int <-
          2 * cf[parm] -
          t(apply(mcf[parm, ], 1, quantile, rev(a)))
        colnames(aux.conf.int) <- paste0(a*100,"%")

        return(aux.conf.int)
        } else if (method == "percentile") {
          return(t(apply(mcf[parm, ], 1, quantile, a)))
        }
      }
}

#' Summarise a linear regression model via generalized cross entropy fit
#'
#' summary method for class \code{\link{lmgce}}. Function used to produce
#' summary information from a fitted linear regression model via generalized
#' cross entropy as represented by \code{object} of class \code{\link{lmgce}}.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param call Boolean value. if \code{TRUE}, the call used is returned.
#' The default is \code{model = TRUE}.
#' @param correlation Boolean value. if \code{TRUE}, the correlation
#' matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor Boolean value. if \code{TRUE}, print the correlations in
#' a symbolic form (see \code{\link[stats]{symnum}}) rather than as numbers.
#' @param ci.level the confidence level (0,1) required to compute the confidence
#' interval. The default is \code{ci.level = NULL} which results in the omission
#'  of the confidence interval.
#' @param ci.method method used to compute a confidence interval. One of
#' c("z","percentile", "basic"). The default is \code{ci.method = "z"}.
#' @param boot.B A single positive integer greater or equal to 10 for the number
#' of bootstrap replicates for the computation of the bootstrap confidence
#' interval(s), to be used when \code{method = c("percentile", "basic")} and
#' when \code{object} was created with \code{boot.B = 0}. The default is
#' \code{boot.B = 100} when the \code{object} has no previous sampling information
#' and \code{boot.B = object$boot.B} otherwise, which corresponds to
#' the \code{boot.B} given to \code{lmgce} when the \code{object} was created.
#' @param boot.method Method used for bootstrapping. One of
#' \code{c("residuals", "cases", "wild")} which corresponds to resampling on
#' residuals, on individual cases or on residuals multiplied by a N(0,1) variable,
#' respectively. The default is \code{boot.method = object$boot.method}.
#' @param ... additional arguments.
#'
#' @return The function \code{summary.lmgce} computes and returns a list of
#' summary statistics of the fitted \code{\link{lmgce}} linear model given in
#' \code{object}, using the components (list elements) "call" and "terms" from
#' its argument, plus
#'
#' \item{residuals}{the residuals, that is response minus fitted values.}
#' \item{coefficients}{a \eqn{p \times 4} matrix, where \eqn{p} is the number of
#' non-aliased coefficients, with columns for the estimated coefficient, its
#' standard error, z-statistic and corresponding (two-sided)
#' p-value. Aliased coefficients are omitted.}
#' \item{support}{a \eqn{p \times 3} matrix with columns for the normalized
#' entropy (NormEnt), and lower (LL) and upper (UL) limits for each of the
#' \eqn{K+1} support spaces.}
#' \item{aliased}{named logical vector showing if the original coefficients are
#' aliased.}
#' \item{sigma}{the square root of the estimated variance of the random error.}
#' \item{df}{degrees of freedom, a 3-vector \eqn{(p, n - p)} the first being the
#' number of non-aliased coefficients, the last being the \eqn{p} minus the
#' number of included individuals \eqn{n}.}
#' \item{r.squared}{\eqn{R^2}, the ‘fraction of variance explained by the model’}
#' \item{adj.r.squared}{the above \eqn{R^2} statistic ‘adjusted’, penalizing for
#' higher \eqn{p}.}
#' \item{cov.unscaled}{a \eqn{p \times p} matrix of covariances of the
#' \eqn{\hat \beta}}
#' \item{support.stdUL}{when applicable, the upper limit of the standardized
#' support chosen, when \code{support.method = "standardized"} or the factor used
#'  when \code{support.method = "ridge"}.}
#' \item{support.method}{method chosen for the support's limits}
#' \item{nep}{the normalized entropy of the model.}
#' \item{nep.cv.mean}{the cross-validation normalized entropy of the model.}
#' \item{nep.cv.sd}{the standard deviation of the cross-validation normalized
#' entropy of the model.}
#' \item{error}{the error measure chosen}
#' \item{error.which}{which criterion/standardized/factor support was used}
#' \item{error.measure}{the value of the error measure}
#' \item{error.measure.cv.mean}{the cross-validation value of the error measure}
#' \item{error.measure.cv.sd}{the standard deviation of the cross-validation
#' value of the error measure}
#' \item{correlation}{the correlation matrix corresponding to the above
#' cov.unscaled, if \code{correlation = TRUE} is specified.}
#' \item{symbolic.cor}{(only if \code{correlation = TRUE}) The value of the
#' argument \code{symbolic.cor}.}
#' \item{na.action}{from object, if present there.}
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' sm_res_gce_package <- summary(res_gce_package)
#'
#' str(sm_res_gce_package)
#'
#' sm_res_gce_package$coefficients
#'
#' @method summary lmgce
#' @export

summary.lmgce <- function(object,
                          call = TRUE,
                          correlation = FALSE,
                          symbolic.cor = FALSE,
                          ci.level = NULL,
                          ci.method = c("z", "percentile", "basic"),
                          boot.B = ifelse(object$boot.B == 0,
                                          100,
                                          object$boot.B),
                          boot.method = object$boot.method,
                          ...) {
  if (is.null(object$terms))
    stop("invalid 'lmgce' object:  no 'terms' component")
  if (!inherits(object, "lmgce"))
    warning("calling summary.lmgce(<fake-lmgce-object>) ...")

  ci.method <- match.arg(ci.method)

  if (is.null(object$coefficients)) {
    ans <- object[c("call", "terms")]
    ans$aliased <- is.na(coef(object))
    class(ans) <- "summary.lmgce"
    ans
  } else {
  rdf <- object$df.residual
  r <- as.numeric(object$residuals)
  n <- length(r)
  p <- length(object$coefficients)
  f <- object$fitted.values
  if (!is.null(object$offset)) {
    f <- f - object$offset
  }

  mss <- if (attr(object$terms, "intercept"))
    sum((f - mean(f))^2) else sum(f^2)
  rss <- sum(r^2)

  resvar <- rss / rdf

  if (!is.null(object$offset)) {
    f <- f - object$offset
  }

  est <- coef(object)
  R <- vcov(object)
  std <- sqrt(diag(R))[names(est)] #sqrt(diag(R) * resvar)

  zval <- est/std

  p1 <- 1L:p
  ans <- object[c("call", "terms")]
  if (!isTRUE(call))
    ans$call <- NULL
  ans$residuals <- r

  ans$coefficients <-
    cbind(Estimate = est,
          `Std. Deviation` = std,
          `z value` = zval,
          `Pr(>|t|)` = 2 * pnorm(abs(zval), lower.tail = FALSE))

  if (!is.null(ci.level)) {
    ans$conf.int.est <-
      confint(object,
              level = ci.level,
              method = ci.method,
              boot.B = boot.B,
              boot.method = boot.method)
    if (ci.method != "z") {
    ans$conf.int.NormEnt <-
      confint(object,
              which = "NormEnt",
              level = ci.level,
              method = ci.method,
              boot.B = boot.B,
              boot.method = boot.method)}
    else {
      ans$conf.int.NormEnt <- NULL
      }
    } else {
    ans$conf.int.est <- NULL}

  ans$support <-
    cbind(`NormEnt` = object$nepk,
          `SupportLL` = object$support.matrix[,1],
          `SupportUL` = object$support.matrix[,2])

  ans$aliased <- is.na(est)
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf)
  if (p != attr(object$terms, "intercept")) {
    df.int <- if (attr(object$terms, "intercept")) 1L else 0L
    ans$r.squared <- mss / (mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
  } else ans$r.squared <- ans$adj.r.squared <- 0

  ans$cov.unscaled <- R

  if (!is.null(object$na.action))
    ans$na.action <- object$na.action

  if (correlation) {
    Invstd <- solve(diag(na.omit(std)))
    ans$correlation <- Invstd %*% R %*% Invstd
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }

  ans$support.method <- object$support.method
  ans$support.stdUL <- object$support.stdUL
  ans$nep <- object$nep
  ans$nep.cv.mean <- object$nep.cv.mean
  ans$nep.cv.sd <- object$nep.cv.sd
  ans$error <- object$error
  ans$error.which <- object$error.which
  ans$error.measure <- object$error.measure
  ans$error.measure.cv.mean <- object$error.measure.cv.mean
  ans$error.measure.cv.sd <- object$error.measure.cv.sd

  class(ans) <- "summary.lmgce"
  ans
  }
}

#' Print Summary of \code{\link{lmgce}} Model Fits
#'
#' \code{print.summary} method for class \code{\link{lmgce}}.
#'
#' @param x an object of class \code{\link{summary.lmgce}}, usually, a result of
#'  a call to \code{\link{summary.lmgce}}.
#' @param digits The number of significant digits to use when printing.
#' @param symbolic.cor Boolean value. if \code{TRUE}, print the correlations in
#' a symbolic form (see \code{\link[stats]{symnum}}) rather than as numbers.
#' @param signif.stars Boolean value. if \code{TRUE}, ‘significance stars’ are
#' printed for each coefficient.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function \code{print.summary.lmgce} prints the information in a
#' \code{\link{summary.lmgce}} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' summary(res_gce_package)
#'
#' summary(res_gce_package, ci.level = 0.90, ci.method = "basic")
#'
#' @method print summary.lmgce
#' @export

print.summary.lmgce <-
  function(x,
           digits = max(3L, getOption("digits") - 3L),
           symbolic.cor = x$symbolic.cor,
           signif.stars = getOption("show.signif.stars"),
           ...) {
    if (!is.null(x$call)) {
      cat("\nCall:\n",
          paste(deparse(x$call), sep = "\n", collapse = "\n"),
          "\n\n",
          sep = "")
    }
    if (length(x$aliased) != 0L) {
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat("Residuals:\n", sep = "")

    if (rdf > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2L)
        structure(
          apply(t(resid), 1L, quantile), dimnames = list(nam,
                                                         dimnames(resid)[[2L]]))
      else {
        zz <- zapsmall(quantile(resid), digits + 1L)
        structure(zz, names = nam)
      }
      print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
      print(resid, digits = digits, ...)
    }
    else {
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
      cat("\n")
    }
    }
    if (length(x$aliased) == 0L) {
      cat("\nNo Coefficients\n")
    }
    else {
      cat("\nCoefficients:\n")
      coefs <- x$coefficients
      if (!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(
          NA,
          length(aliased),
          4,
          dimnames = list(cn,colnames(coefs)))
        coefs[!aliased,] <- x$coefficients
      }

      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)

      if (!is.null(x$conf.int.est)) {
        cat("\n")
        print(round(x$conf.int.est,
                    max(1L, getOption("digits") - 1L)))
      }
    }
    if (length(x$aliased) != 0L) {
      cat("\nNormalized Entropy:")
      cat("\n")
      print(round(x$support,
                  max(1L, getOption("digits") - 1L)))
      if (!is.null(x$conf.int.NormEnt)) {
        cat("\n")
        print(round(x$conf.int.NormEnt,
                    max(1L, getOption("digits") - 1L)))
      }
      cat("\nResidual standard error:",
          format(signif(x$sigma, digits)),
          "on",
          rdf,
          "degrees of freedom")
      cat("\n")
      if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep = "")
      cat(
        {if (x$support.method == "standardized") {
          "Chosen Upper Limit for Standardized Supports: "} else {
            "Chosen factor for the Upper Limit of the Supports: "
          }},
        x$support.stdUL,
        ", Chosen Error: ",
        x$error.which,
        sep = ""
      )
      cat("\n")
      cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
      cat(", Adjusted R-squared:",
          formatC(x$adj.r.squared, digits = digits))
      cat("\n")
      cat("NormEnt:", formatC(x$nep, digits = digits))
      cat(
        ", CV-NormEnt: ",
        formatC(x$nep.cv.mean, digits = digits),
        " (",
        formatC(x$nep.cv.sd, digits = digits),
        ")",
        sep = ""
      )
      cat("\n")
      cat(x$error, ": ", formatC(x$error.measure, digits = digits), sep = "")
      cat(
        ", CV-",
        x$error,
        ": ",
        formatC(x$error.measure.cv.mean, digits = digits),
        " (",
        formatC(x$error.measure.cv.sd, digits = digits), ")", sep = ""
      )
      correl <- x$correlation
      if (!is.null(correl)) {
        cat("\n")
        p <- NCOL(correl)
        if (p > 1L) {
          cat("\nCorrelation of Coefficients:\n")
          if (is.logical(symbolic.cor) && symbolic.cor) {
            print(symnum(correl, abbr.colnames = NULL))
          }
          else {
            correl <- format(round(correl, 2), nsmall = 2, digits = digits)
            correl[!lower.tri(correl)] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
          }
        }
      }
    }
    cat("\n")

  }

#' Plot Diagnostics for a \code{\link{lmgce}} Object
#'
#' Seven plots (selectable by \code{which}) are currently available to
#' evaluate a \code{\link{lmgce}} object: a plot of the Estimates and confidence
#' intervals; four plots of supports against Prediction Error, Estimates,
#' Normalized Entropy and Precision Error; two plots of GCE reestimation against
#' Prediction and Precision Errors. Note that plots regarding Precision Error are
#' only produced if the argument \code{coef} is not \code{NULL}.
#'
#' @param x Fitted \code{\link{lmgce}} model object.
#' @param type One of \code{c("ggplot2", "plotly")}. "ggplot2" is used
#' by default.
#' @param which A subset of the numbers 1:7.
#' @param ci.level the confidence level (0,1) required to compute the confidence
#' interval.
#' @param ci.method the method used to compute the confidence interval. One of
#' c("z","percentile", "basic"). The default is \code{method = "z"}.
#' @param boot.B A single positive integer greater or equal to 10 for the number
#' of bootstrap replicates for the computation of the bootstrap confidence
#' interval(s), to be used when \code{method = c("percentile", "basic")} and
#' when \code{object} was created with \code{boot.B = 0}. The default is
#' \code{boot.B = 100} when the \code{object} has no previous sampling information
#' and \code{boot.B = object$boot.B} otherwise, which corresponds to
#' the \code{boot.B} given to \code{lmgce} when the \code{object} was created.
#' @param boot.method Method used for bootstrapping. One of
#' \code{c("residuals", "cases", "wild")} which corresponds to resampling on
#' residuals, on individual cases or on residuals multiplied by a N(0,1) variable,
#' respectively. The default is \code{boot.method = object$boot.method}.
#' @param coef A vector of true coefficients to be used when \code{which = c(5,7)}.
#' @param OLS Boolean value. if \code{TRUE}, the default, plots the OLS results.
#' @param NormEnt Boolean value. if \code{TRUE}, the default, plots the
#' normalized entropy.
#' @param caption Captions to appear above the plots;
#' \code{\link[base]{character}} vector or \code{\link[base]{list}}
#' of valid graphics annotations, see \code{\link[grDevices]{as.graphicsAnnot}},
#' of length 7, the j-th entry corresponding to which[j]. Can be set to "" or
#' \code{NA} to suppress all captions.
#' @param ... additional arguments.
#'
#' @return A named list of \code{ggplot} or \code{plotly} objects, each
#' representing a separate plot.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @seealso \code{\link{lmgce}}
#'
#' @examples
#' \donttest{
#' res_gce_package <-
#'   lmgce(y ~ .,
#'         data = dataGCE,
#'         boot.B = 50,
#'         seed = 230676)
#' }
#' plot(res_gce_package)
#'
#' @method plot lmgce
#' @importFrom rlang .data
#' @importFrom grDevices as.graphicsAnnot
#' @importFrom graphics abline arrows axis legend lines matplot mtext par points title
#' @export

plot.lmgce <-
  function (x,
            type = c("ggplot2", "plotly"),
            which = 1:7,
            ci.level = 0.95,
            ci.method = c("z", "percentile", "basic"),
            boot.B = ifelse(x$boot.B == 0,
                            100,
                            x$boot.B),
            boot.method = x$boot.method,
            coef = NULL,
            OLS = TRUE,
            NormEnt = TRUE,
            caption = list(
              paste0("Estimates (", ci.method[1], " ", ci.level * 100, "% CI)"),
              "Prediction Error vs supports",
              "Estimates vs supports",
              "Normalized Entropy vs supports",
              "Precision Error vs supports",
              "Prediction Error vs GCE reestimation",
              "Precision Error vs GCE reestimation"
            ),
            ...) {
    object <- x
    if (!inherits(object, "lmgce"))
      stop("use only with \"lmgce\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 7))
      stop("'which' must be in 1:7")

    type <- match.arg(type)
    ci.method  <- match.arg(ci.method)
    show <- rep(FALSE, 7)
    show[which] <- TRUE

    plots <- list(p1 = NULL,
                  p2 = NULL,
                  p3 = NULL,
                  p4 = NULL,
                  p5 = NULL,
                  p6 = NULL,
                  p7 = NULL)

    getCaption <- function(k) {
      if (length(caption) < k)
        ""
    else
      as.graphicsAnnot(caption[[k]])
    }

    manual <- object$support.signal.manual
    stdULmin <- object$support.signal.min
    stdUL1se <- object$support.signal.1se
    elbow <- object$support.signal.elbow

    error.OLS <- NULL

    col.coef.all <- viridis::turbo(length(object$coefficients))
    if ("(Intercept)" %in% names(object$coefficients))
      col.coef <- col.coef.all[-1] else
      col.coef <- col.coef.all

    # plot 1 ####

    if (show[1L]) {
      coefs <- coef(object)
      var.names = variable.names(object)
      est.perc <- confint(object,
                          level = ci.level,
                          method = ci.method,
                          boot.B = boot.B,
                          boot.method = boot.method)
      ncoefs <- length(coefs)
      idx <- seq(1, ncoefs)
      k <- 1 / ncoefs

      data.plot.1 <-
          data.frame(
            predictor = length(coefs):1,
            labels = variable.names(object),
            estimate = coefs,
            LL = est.perc[, 1],
            UL = est.perc[, 2])

        plots$p1 <-
        ggplot2::ggplot(data.plot.1,
                        ggplot2::aes(y = .data$predictor,
                                     x = .data$estimate)) +
          ggplot2::geom_point(size = 1.5, colour = col.coef.all) +
          ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$LL,
                                               xmax = .data$UL),
                                 orientation = "y",
                                 colour = col.coef.all, width = 0) +
          ggplot2::geom_vline(
            xintercept = 0,
            color = "black",
            linetype = "dashed",
            linewidth = 0.7,
            alpha = 0.75
          ) +
          ggplot2::scale_y_continuous(
            name = "",
            breaks = length(coefs):1,
            labels = data.plot.1$labels,
            trans = "reverse"
          ) +
          ggplot2::xlab("Estimates") +
          ggplot2::ylab("") +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
          ggplot2::ggtitle(getCaption(1))

        if (length(object$results$OLS$error) != 0 &&
            isTRUE(OLS))
          plots$p1 <-
          plots$p1 +
          ggplot2::geom_point(
            data = data.frame(x = coef(object$results$OLS$res), y = length(coefs):1),
            ggplot2::aes(x = .data$x, y = .data$y),
            shape = 1,
            size = 1.5,
            colour = col.coef.all
          )

        if (!is.null(coef))
          plots$p1 <-
          plots$p1 +
          ggplot2::geom_point(
            data = data.frame(x = coef, y = length(coefs):1),
            ggplot2::aes(x = .data$x, y = .data$y),
            shape = 4,
            size = 1.5,
            colour = col.coef.all
          )

          if (type == "plotly")
         plots$p1 <- plotly::ggplotly(plots$p1)
    }

    # plot 2 ####

    if (show[2L] &
      (length(object$support.ok) > 1)) {
    x <- round(as.numeric(object$support.ok), 8)

    if (!is.null(object$results$cv)) {
      ycv <- object$results$cv$error.measure.cv.mean[as.character(object$support.ok)]
      sderror <- object$results$cv$error.measure.cv.sd[as.character(object$support.ok)]

      nepcv <- object$results$cv$nep.cv.mean[as.character(object$support.ok)]
      sdnep <- object$results$cv$nep.cv.sd[as.character(object$support.ok)]

      if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
        error.OLS <- mean(object$results$OLS$error)

        data.plot.2 <-
          data.frame(
            support = x,
            error = ycv,
            error.LL = ycv - sderror,
            error.UL = ycv + sderror,
            nep = nepcv,
            nep.LL = nepcv - sdnep,
            nep.UL = nepcv + sdnep
          )

        plots$p2 <-
          ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                    x = .data$support)) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$error.LL, ymax = .data$error.UL),
            width = 0,
            colour = "darkgrey"
          ) +
          ggplot2::geom_line(colour = "red", linetype = "dashed") +
          ggplot2::geom_point(size = 1.5, colour = "red") +
          ggplot2::geom_point(
            data = data.frame(
              x = object$support.stdUL,
              y = object$error.measure.cv.mean
            ),
            ggplot2::aes(x = .data$x, y = .data$y),
            size = 1.5,
            colour = "red4"
          )

        if (isTRUE(NormEnt)) {
          ylim.left <- c(min(data.plot.2$error.LL), max(data.plot.2$error.UL))
          ylim.right <- c(min(data.plot.2$nep.LL), max(data.plot.2$nep.UL))

          b <- diff(ylim.left) / diff(ylim.right)
          a <- ylim.left[1] - b * ylim.right[1]

          plots$p2 <-
            plots$p2 +
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = a + .data$nep.LL * b,
                           ymax = a + .data$nep.UL * b),
              width = 0,
              colour = "black"
            ) +
            ggplot2::geom_point(ggplot2::aes(y = a + .data$nep * b,
                                             x = .data$support),
                                size = 1.5,
                                colour = "green") +
            ggplot2::geom_line(
              ggplot2::aes(y = a + .data$nep * b,
                           x = .data$support),
              colour = "green",
              linetype = "dashed"
            ) +
            ggplot2::geom_point(
              data = data.frame(
                x = object$support.stdUL,
                y = a + object$nep.cv.mean * b
              ),
              ggplot2::aes(x = .data$x, y = .data$y),
              size = 1,
              colour = "green4"
            ) +
            ggplot2::scale_y_continuous(
              paste0("CV-", object$error),
              sec.axis =  ggplot2::sec_axis( ~ (. - a) / b, name = "CV Normalized Entropy")
            )
        }

        plots$p2 <-
          plots$p2 +
          ggplot2::xlab("Upper limit of the support spaces") +
          ggplot2::ylab(paste0("CV-", object$error)) +
          ggplot2::geom_vline(
            xintercept = c(stdULmin, stdUL1se, elbow, manual),
            linetype = c(3, {
              if (is.null(stdUL1se))
                NULL
              else
                2
            }, 3, {
              if (is.null(manual))
                NULL
              else
                3
            }),
            colour = c(1, {
              if (is.null(stdUL1se))
                NULL
              else
                1
            }, 2, {
              if (is.null(manual))
                NULL
              else
                3
            })
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black")
          ) +
          ggplot2::ggtitle(getCaption(2))

        if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
          plots$p2 <- plots$p2 +
            ggplot2::geom_hline(yintercept = error.OLS,
                                linetype = "dotted")}
        if (type == "plotly")
          plots$p2 <- plotly::ggplotly(plots$p2)

    } else {
      y <- object$results$nocv$support.results$error.measure.ss[as.character(object$support.ok)]
      nep <- object$results$nocv$support.results$nep.ss[as.character(object$support.ok)]

      if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
        error.OLS <- object$results$OLS$error

        data.plot.2 <-
          data.frame(support = x,
                     error = y,
                     nep = nep)

        plots$p2 <-
          ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                    x = .data$support)) +
          ggplot2::geom_point(size = 1.5, colour = "red") +
          ggplot2::geom_line(colour = "red", linetype = "dashed") +
          ggplot2::geom_point(
            data = data.frame(x = object$support.stdUL, y = object$error.measure),
            ggplot2::aes(x = .data$x, y = .data$y),
            size = 1.5,
            colour = "red4"
          )

        if (isTRUE(NormEnt)) {
          ylim.left <- c(min(data.plot.2$error), max(data.plot.2$error))
          ylim.right <- c(min(data.plot.2$nep), max(data.plot.2$nep))

          b <- diff(ylim.left) / diff(ylim.right)
          a <- ylim.left[1] - b * ylim.right[1]

          plots$p2 <-
            plots$p2 +
            ggplot2::geom_point(ggplot2::aes(y = a + .data$nep * b,
                                             x = .data$support),
                                size = 1.5,
                                colour = "green") +
            ggplot2::geom_point(
              data = data.frame(x = object$support.stdUL, y = a + object$nep * b),
              ggplot2::aes(x = .data$x, y = .data$y),
              size = 1.5,
              colour = "green4"
            ) +
            ggplot2::geom_line(
              ggplot2::aes(y = a + .data$nep * b, x = .data$support),
              colour = "green",
              linetype = "dashed"
            ) +
            ggplot2::scale_y_continuous(
              paste0("CV-", object$error),
              sec.axis =  ggplot2::sec_axis( ~ (. - a) / b, name = "Normalized Entropy")
            )
        }

        plots$p2 <-
          plots$p2 +
          ggplot2::xlab("Upper limit of the support spaces") +
          ggplot2::ylab(paste0("CV-", object$error))+
          ggplot2::geom_vline(
            xintercept = c(stdULmin, stdUL1se, elbow, manual),
            linetype = c(3, {
              if (is.null(stdUL1se))
                NULL
              else
                2
            }, 3, {
              if (is.null(manual))
                NULL
              else
                3
            }),
            colour = c(1, {
              if (is.null(stdUL1se))
                NULL
              else
                1
            }, 2, {
              if (is.null(manual))
                NULL
              else
                3
            })
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
          ggplot2::ggtitle(getCaption(2))

        if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
          plots$p2 <- plots$p2 +
            ggplot2::geom_hline(yintercept = error.OLS,
                                linetype = "dotted")}

        if (type == "plotly")
          plots$p2 <- plotly::ggplotly(plots$p2)

    }
  }

    # plot 3 and 4 ####

  if ((show[3L] || show[4L]) &
      (length(object$support.ok) > 1)) {
    coef.matrix.wide <- object$results$nocv$support.results$coef.matrix.ss
    coef.matrix.long <- data.frame(predictor = NA,
                                   support = NA,
                                   estimate = NA,
                                   nepk = NA,
                                   nepk_nep = NA)

    for (i in 1:ncol(coef.matrix.wide)) {
      coef.matrix.long <-
        rbind(coef.matrix.long,
              cbind(predictor = row.names(coef.matrix.wide)[1:(nrow(coef.matrix.wide) / 2)],
                    support = colnames(coef.matrix.wide)[i],
                    estimate =coef.matrix.wide[1:(nrow(coef.matrix.wide) / 2), i],
                    nepk = coef.matrix.wide[(nrow(coef.matrix.wide) / 2 + 1):nrow(coef.matrix.wide), i],
                    nepk_nep = coef.matrix.wide[(nrow(coef.matrix.wide) / 2 + 1):nrow(coef.matrix.wide), i] / object$results$nocv$support.results$nep.ss[i])
              )
    }
    coef.matrix.long <- coef.matrix.long[-1, ]
    coef.matrix.long$support <- as.numeric(coef.matrix.long$support)
    coef.matrix.long$estimate <- as.numeric(coef.matrix.long$estimate)
    coef.matrix.long$nepk <- as.numeric(coef.matrix.long$nepk)
    coef.matrix.long$nepk_nep <- as.numeric(coef.matrix.long$nepk_nep)

    if (show[3L]) {

        plots$p3 <-
          ggplot2::ggplot(data = coef.matrix.long,
                          ggplot2::aes(x = .data$support,
                                       y = .data$estimate)) +
          ggplot2::geom_line(ggplot2::aes(group = .data$predictor,
                                          colour = .data$predictor)) +
          ggplot2::geom_point(data = data.frame(x = rep(object$support.stdUL,
                                                        length(coef(object))),
                                                y = coef(object)),
                              ggplot2::aes(x = .data$x, y = .data$y),
                              size = 1.5,
                              colour = col.coef.all) +
          ggplot2::geom_vline(
            xintercept = c(stdULmin, stdUL1se, elbow, manual),
            linetype = c(3, {
              if (is.null(stdUL1se))
                NULL
              else
                2
            }, 3, {
              if (is.null(manual))
                NULL
              else
                3
            }),
            colour = c(1, {
              if (is.null(stdUL1se))
                NULL
              else
                1
            }, 2, {
              if (is.null(manual))
                NULL
              else
                3
            })
          ) +
          ggplot2::xlab("Upper limit of the support spaces") +
          ggplot2::ylab("Estimates") +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black"),
            legend.title = ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = col.coef.all) +
          ggplot2::ggtitle(getCaption(3))

        if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
          plots$p3 <-
            plots$p3 +
            ggplot2::geom_hline(yintercept = coef(object$results$OLS$res),
                                linetype = "dotted",
                                colour = col.coef.all)

        if (!is.null(coef))
          plots$p3 <-
            plots$p3 +
            ggplot2::geom_hline(yintercept = coef,
                                linetype = "dotdash",
                                colour = col.coef.all)

        if (type == "plotly") {
          plots$p3 <- plotly::ggplotly(plots$p3)
        }
      }

    if (show[4L]) {

      plots$p4 <-
        ggplot2::ggplot(data = na.omit(coef.matrix.long),
                        ggplot2::aes(x = .data$support,
                                     y = .data$nepk)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$predictor,
                                        colour = .data$predictor)) +
        ggplot2::geom_point(data = data.frame(x = rep(object$support.stdUL,
                                                      length(coef(object))),
                                              y = NormEnt(object, model = FALSE)),
                            ggplot2::aes(x = .data$x, y = .data$y),
                            size = 1.5,
                            colour = col.coef.all) +
        ggplot2::geom_vline(
          xintercept = c(stdULmin, stdUL1se, elbow, manual),
          linetype = c(3, {
            if (is.null(stdUL1se))
              NULL
            else
              2
          }, 3, {
            if (is.null(manual))
              NULL
            else
              3
          }),
          colour = c(1, {
            if (is.null(stdUL1se))
              NULL
            else
              1
          }, 2, {
            if (is.null(manual))
              NULL
            else
              3
          })
        ) +
        ggplot2::xlab("Upper limit of the support spaces") +
        ggplot2::ylab("Normalized Entropy") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
          axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
          axis.title.x = ggplot2::element_text(size = 12, colour = "black"),
          legend.title = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(values = col.coef.all) +
        ggplot2::ggtitle(getCaption(4))

      if (type == "plotly") {
        plots$p4 <- plotly::ggplotly(plots$p4)
      }
    }
  }

    # plot 5 ####

    if ((show[5L] | show[7]) &
        !is.null(coef) &
         (length(object$results$twosteps) != 0)) {
      if (!is.null(object$results$cv)) {
      beta.error.cv.step <-
        matrix(NA,
               ncol = length(object$results$twosteps) + 1,
               nrow = length(object$results$cv$repeats1))

      for (nfolds in 1:length(object$results$cv$repeats1)) {
        beta.matrix <-
          matrix(NA,
                 ncol = length(object$results$twosteps) + 1,
                 nrow = length(coef(object)))

        if (is.null(object$results$cv$repeats1[[nfolds]]$support.results$coef.matrix.ss)) {
          beta.matrix[, 1] <- object$results$cv$repeats1$fold1$coefficients} else {
            beta.matrix[, 1] <-
              object$results$cv$repeats1[[nfolds]]$support.results$coef.matrix.ss[1:length(coef(object)),
                                                                                  colnames(object$results$cv$repeats1[[nfolds]]$support.results$coef.matrix.ss) == object$support.stdUL]
          }

        for (i_mat in 2:(length(object$results$twosteps) + 1)) {
          beta.matrix[, i_mat] <-
            object$results$twosteps[[i_mat - 1]]$results$cv$repeats1[[nfolds]]$coefficients
        }
        beta.error.cv.step[nfolds, ] <- apply(beta.matrix, 2, accmeasure, coef)
      }
      } else {
        beta.matrix.step <-
          matrix(NA,
                 ncol = length(object$results$twosteps) + 1,
                 nrow = length(coef(object)))

        if (length(object$results$nocv$support.results) == 1) {
          beta.matrix.step[, 1] <- object$results$nocv$support.results[[1]]$coefficients
        } else {
          beta.matrix.step[, 1] <- object$results$nocv$support.results[[as.character(object$support.stdUL)]]$coefficients
        }
        for (i_mat in 2:(length(object$results$twosteps) + 1)) {
          beta.matrix.step[, i_mat] <-
            object$results$twosteps[[i_mat - 1]]$coefficients
        }
    }}

    if (show[5L] && (!is.null(coef)) &&
        (length(object$support.ok) > 1)) {
      x <- round(as.numeric(object$support.ok), 8)

      if (!is.null(object$results$cv)) {
      beta.error.cv <-
        matrix(NA,
               ncol = length(object$results$cv$repeats1$fold1$support.results) - 3,
               nrow = length(object$results$cv$repeats1))

      for (i in 1:length(object$results$cv$repeats1)) {
        beta.error.cv[i,] <-
          apply(
            object$results$cv$repeats1[[i]]$support.results$coef.matrix.ss[1:length(coef(object)),],
            2,
            accmeasure,
            coef
          )
      }

      if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
        error.OLS <- mean(apply(object$results$OLS$matrix.coef,
                                2,
                                accmeasure,
                                coef))

      ycv <- apply(beta.error.cv,2,mean)
      sderror <- apply(beta.error.cv,2,sd)

      nepcv <- object$results$cv$nep.cv.mean[as.character(object$support.ok)]
      sdnep <- object$results$cv$nep.cv.sd[as.character(object$support.ok)]

        data.plot.5 <-
          data.frame(
            support = x,
            error = ycv,
            error.LL = ycv - sderror,
            error.UL = ycv + sderror,
            nep = nepcv,
            nep.LL = nepcv - sdnep,
            nep.UL = nepcv + sdnep
          )

        plots$p5 <-
          ggplot2::ggplot(data.plot.5, ggplot2::aes(y = .data$error,
                                                    x = .data$support)) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$error.LL,
                                              ymax = .data$error.UL),
                                 width = 0,
                                 colour = "darkgrey") +
          ggplot2::geom_point(size = 1.5, colour = "red") +
          ggplot2::geom_line(colour = "red",
                             linetype = "dashed") +
          ggplot2::geom_point(data = data.frame(x = object$support.stdUL,
                                                y = {if (length(object$results$twosteps) == 0) {
                                                  ycv[object$support.ok == object$support.stdUL]
                                                } else {
                                                  mean(beta.error.cv.step[ ,ncol(beta.error.cv.step)])
                                                }}),
                              ggplot2::aes(x = .data$x, y = .data$y),
                              size = 1.5,
                              colour = "red4")

        if (isTRUE(NormEnt)) {
          ylim.left <- c(min(data.plot.5$error.LL), max(data.plot.5$error.UL))
          ylim.right <- c(min(data.plot.5$nep.LL), max(data.plot.5$nep.UL))

          b <- diff(ylim.left) / diff(ylim.right)
          a <- ylim.left[1] - b * ylim.right[1]

          plots$p5 <-
            plots$p5 +
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = a + .data$nep.LL * b,
                           ymax = a + .data$nep.UL * b),
              width = 0,
              colour = "black"
            ) +
            ggplot2::geom_point(ggplot2::aes(y = a + .data$nep * b,
                                             x = .data$support),
                                size = 1.5,
                                colour = "green") +
            ggplot2::geom_line(
              ggplot2::aes(y = a + .data$nep * b,
                           x = .data$support),
              colour = "green",
              linetype = "dashed"
            ) +
            ggplot2::geom_point(
              data = data.frame(
                x = object$support.stdUL,
                y = a + object$nep.cv.mean * b
              ),
              ggplot2::aes(x = .data$x, y = .data$y),
              size = 1.5,
              colour = "green4"
            ) +
            ggplot2::scale_y_continuous(
              paste0("CV-", object$error),
              sec.axis =  ggplot2::sec_axis( ~ (. - a) /
                                               b, name = "CV Normalized Entropy")
            )
        }

        plots$p5 <-
          plots$p5 +
          ggplot2::xlab("Upper limit of the support spaces") +
          ggplot2::ylab(paste0("CV-", object$error)) +
          ggplot2::geom_vline(
            xintercept = c(stdULmin, stdUL1se, elbow, manual),
            linetype = c(3, {
              if (is.null(stdUL1se))
                NULL
              else
                2
            }, 3, {
              if (is.null(manual))
                NULL
              else
                3
            }),
            colour = c(1, {
              if (is.null(stdUL1se))
                NULL
              else
                1
            }, 2, {
              if (is.null(manual))
                NULL
              else
                3
            })
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
          ggplot2::ggtitle(getCaption(5))

        if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
          plots$p5 <- plots$p5 +
            ggplot2::geom_hline(yintercept = error.OLS,
                                linetype = "dotted")}

        if (type == "plotly")
          plots$p5 <- plotly::ggplotly(plots$p5)

    } else {

      y <-
        apply(
          sapply(object$results$nocv$support.results[1:(length(object$results$nocv$support.results) - 3)],
                 function(x) {x$coefficients}
          ),
          2,
          accmeasure,
          coef
        )

      nep <-
          object$results$nocv$support.results$nep.ss[as.character(object$support.ok)]

      if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
        error.OLS <- accmeasure(object$results$OLS$matrix.coef, coef)

        data.plot.5 <-
          data.frame(support = x,
                     error = y,
                     nep = nep)

        plots$p5 <-
          ggplot2::ggplot(data.plot.5, ggplot2::aes(y = .data$error,
                                                    x = .data$support)) +
          ggplot2::geom_point(size = 1, colour = "red") +
          ggplot2::geom_line(colour = "red",
                             linetype = "dashed") +
          ggplot2::geom_point(data = data.frame(x = object$support.stdUL,
                                                y = {if (length(object$results$twosteps) == 0) {
                                                  y[object$support.ok == object$support.stdUL]
                                                } else {
                                                  accmeasure(beta.matrix.step[ ,ncol(beta.matrix.step)],
                                                             coef)
                                                }}),
                              ggplot2::aes(x = .data$x, y = .data$y),
                              size = 1.5,
                              colour = "red4")

        if (isTRUE(NormEnt)) {
          ylim.left <- c(min(data.plot.5$error), max(data.plot.5$error))
          ylim.right <- c(min(data.plot.5$nep), max(data.plot.5$nep))

          b <- diff(ylim.left) / diff(ylim.right)
          a <- ylim.left[1] - b * ylim.right[1]

          plots$p5 <-
            plots$p5 +
            ggplot2::geom_point(ggplot2::aes(y = a + .data$nep * b,
                                             x = .data$support),
                                size = 1.5,
                                colour = "green") +
            ggplot2::geom_point(
              data = data.frame(x = object$support.stdUL, y = a + object$nep * b),
              ggplot2::aes(x = .data$x, y = .data$y),
              size = 1.5,
              colour = "green4"
            ) +
            ggplot2::geom_line(
              ggplot2::aes(y = a + .data$nep * b,
                           x = .data$support),
              colour = "green",
              linetype = "dashed"
            ) +
            ggplot2::scale_y_continuous(
              paste0("CV-", object$error),
              sec.axis =  ggplot2::sec_axis( ~ (. - a) /
                                               b, name = "Normalized Entropy")
            )
        }
        plots$p5 <-
          plots$p5 +
          ggplot2::xlab("Upper limit of the support spaces") +
          ggplot2::ylab(paste0("CV-", object$error)) +
          ggplot2::geom_vline(
            xintercept = c(stdULmin, stdUL1se, elbow, manual),
            linetype = c(3, {
              if (is.null(stdUL1se))
                NULL
              else
                2
            }, 3, {
              if (is.null(manual))
                NULL
              else
                3
            }),
            colour = c(1, {
              if (is.null(stdUL1se))
                NULL
              else
                1
            }, 2, {
              if (is.null(manual))
                NULL
              else
                3
            })
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"),
            axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
            axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
            axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
          ggplot2::ggtitle(getCaption(5))

        if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
          plots$p5 <- plots$p5 +
            ggplot2::geom_hline(yintercept = error.OLS,
                                linetype = "dotted")}

        if (type == "plotly")
          plots$p5 <- plotly::ggplotly(plots$p5)

    }
    }

    # plot 6 and 7 ####
    if ((show[6L] | show[7L]) & (length(object$results$twosteps) != 0)) {
      x <- 0:length(object$results$twosteps)

      if (!is.null(object$results$cv)) {
        if (show[6L]) {

          if (length(object$results$cv$error.measure.cv.mean) == 1) {
            ycv <- c(object$results$cv$error.measure.cv.mean,
                     sapply(object$results$twosteps,
                            function(x) {x$error.measure.cv.mean}))
            sderror <- c(object$results$cv$error.measure.cv.sd,
                         sapply(object$results$twosteps,
                                function(x) {x$error.measure.cv.sd}))
          } else {
          ycv <- c(object$results$cv$error.measure.cv.mean[as.character(object$support.stdUL)],
                   sapply(object$results$twosteps,
                          function(x) {x$error.measure.cv.mean}))
          sderror <- c(object$results$cv$error.measure.cv.sd[as.character(object$support.stdUL)],
                       sapply(object$results$twosteps,
                              function(x) {x$error.measure.cv.sd}))
          }

          if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
            error.OLS <- mean(object$results$OLS$error)

          data.plot.2 <-
            data.frame(
              support = x,
              error = ycv,
              error.LL = ycv - sderror,
              error.UL = ycv + sderror)

          plots$p6 <-
            ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                      x = .data$support)) +
            ggplot2::geom_line(colour = "orange", linetype = "dashed") +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$error.LL,
                                                ymax = .data$error.UL),
                                   width = 0,
                                   colour = "darkgrey") +
            ggplot2::geom_point(size = 1.5,
                                colour = c("red",
                                           {if (length(object$results$twosteps) > 1) {
                                             rep("orange",
                                                 length(object$results$twosteps) - 1)}},
                                           "red4")) +
            ggplot2::xlab("number of reestimations") +
            ggplot2::ylab(paste0("CV-", object$error)) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              panel.background = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.line = ggplot2::element_line(colour = "black"),
              axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
              axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
              axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
            ggplot2::ggtitle(getCaption(6)) +
            ggplot2::scale_x_continuous(breaks = x, labels = x)

          if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
            plots$p6 <- plots$p6 +
              ggplot2::geom_hline(yintercept = error.OLS,
                                  linetype = "dotted")}

          if (type == "plotly")
            plots$p6 <- plotly::ggplotly(plots$p6)

        }

        if (show[7L] && (!is.null(coef))) {

          ycv <- apply(beta.error.cv.step, 2, mean)
          sderror <- apply(beta.error.cv.step, 2, sd)

          if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
            error.OLS <-  mean(apply(object$results$OLS$matrix.coef,
                                     2,
                                     accmeasure,
                                     coef))

            data.plot.2 <-
              data.frame(
                support = x,
                error = ycv,
                error.LL = ycv - sderror,
                error.UL = ycv + sderror)

            plots$p7 <-
              ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                        x = .data$support)) +
              ggplot2::geom_line(colour = "orange", linetype = "dashed") +
              ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$error.LL,
                                                  ymax = .data$error.UL),
                                     width = 0,
                                     colour = "darkgrey") +
              ggplot2::geom_point(size = 1.5,
                                  colour = c("red",
                                             {if (length(object$results$twosteps) > 1) {
                                               rep("orange",
                                                   length(object$results$twosteps) - 1)}},
                                             "red4")) +
              ggplot2::xlab("number of reestimations") +
              ggplot2::ylab(paste0("CV-", object$error)) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black"),
                axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
                axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
                axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
              ggplot2::ggtitle(getCaption(7)) +
              ggplot2::scale_x_continuous(breaks = x, labels = x)

            if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
              plots$p7 <- plots$p7 +
                ggplot2::geom_hline(yintercept = error.OLS,
                                    linetype = "dotted")}

            if (type == "plotly")
              plots$p7 <- plotly::ggplotly(plots$p7)

          }
        } else {

        if (show[6L]) {
          if (length(object$results$nocv$support.results) == 1) {
            y <- c(object$results$nocv$support.results[[1]]$error.measure,
                   sapply(object$results$twosteps,
                          function(x) {x$error.measure}))
            } else {
            y <- c(object$results$nocv$support.results[[as.character(object$support.stdUL)]]$error.measure,
                   sapply(object$results$twosteps,
                          function(x) {x$error.measure}))
            }

          if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
            error.OLS <- object$results$OLS$error

          data.plot.2 <-
            data.frame(support = x,
                       error = y)

          plots$p6 <-
            ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                      x = .data$support)) +
            ggplot2::geom_line(colour = "orange", linetype = "dashed") +
            ggplot2::geom_point(size = 1.5,
                                colour = c("red",
                                           {if (length(object$results$twosteps) > 1) {
                                             rep("orange",
                                                 length(object$results$twosteps) - 1)}},
                                           "red4")) +
            ggplot2::xlab("number of reestimations") +
            ggplot2::ylab(object$error) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              panel.background = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.line = ggplot2::element_line(colour = "black"),
              axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
              axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
              axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
            ggplot2::ggtitle(getCaption(6)) +
            ggplot2::scale_x_continuous(breaks = x, labels = x)

          if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
            plots$p6 <- plots$p6 +
              ggplot2::geom_hline(yintercept = error.OLS,
                                  linetype = "dotted")}

          if (type == "plotly")
            plots$p6 <- plotly::ggplotly(plots$p6)
        }

          if (show[7L] && (!is.null(coef))) {

            y <- apply(beta.matrix.step, 2, accmeasure, coef)

            if (length(object$results$OLS$error) != 0 && isTRUE(OLS))
              error.OLS <-  accmeasure(object$results$OLS$matrix.coef, coef)

            data.plot.2 <-
              data.frame(support = x,
                         error = y)

            plots$p7 <-
              ggplot2::ggplot(data.plot.2, ggplot2::aes(y = .data$error,
                                                        x = .data$support)) +
              ggplot2::geom_line(colour = "orange", linetype = "dashed") +
              ggplot2::geom_point(size = 1.5,
                                  colour = c("red",
                                             {if (length(object$results$twosteps) > 1) {
                                               rep("orange",
                                                   length(object$results$twosteps) - 1)}},
                                             "red4")) +
              ggplot2::xlab("number of reestimations") +
              ggplot2::ylab(object$error) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black"),
                axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
                axis.text.x.bottom = ggplot2::element_text(size = 12, colour = "black"),
                axis.title.x = ggplot2::element_text(size = 12, colour = "black")) +
              ggplot2::ggtitle(getCaption(7)) +
              ggplot2::scale_x_continuous(breaks = x, labels = x)

            if (length(object$results$OLS$error) != 0 && isTRUE(OLS)) {
              plots$p7 <- plots$p7 +
                ggplot2::geom_hline(yintercept = error.OLS,
                                    linetype = "dotted")}

            if (type == "plotly")
              plots$p7 <- plotly::ggplotly(plots$p7)

        }
        }
    }

    plots[!sapply(plots, is.null)]

}

#' Predict method for \code{\link{lmgce}} Linear Model Fits
#'
#' Predicted values based on a fitted model \code{\link{lmgce}} object.
#'
#' @param object Fitted \code{\link{lmgce}} model object.
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used.
#' @param interval One of \code{c("none", "confidence")}. Type of interval
#' calculation. Can be abbreviated.
#' @param type One of \code{c("response", "terms")}. Type of prediction
#' (response or model term). Can be abbreviated.
#' @param level Tolerance/confidence level (0,1).
#' @param terms if \code{type = "terms"}, which terms (default is all terms),
#' a \code{\link[base]{character}} vector.
#' @param na.action function determining what should be done with missing values
#' in \code{newdata}. The default is to predict \code{NA}.
#' @param ... additional arguments.
#'
#' @return \code{predict.lmgce} produces a vector of predictions.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' predict(res_gce_package, dataGCE.test)
#'
#' @method predict lmgce
#' @importFrom stats predict
#' @export

predict.lmgce <-
  function(object,
           newdata,
           interval = c("none", "confidence"),
           type = c("response", "terms"),
           level = 0.95,
           terms = NULL,
           na.action = na.pass,
           ...)
  {
    tt <- terms(object)
    if (!inherits(object, "lmgce"))
      warning("calling predict.lmgce(<fake-lmgce-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
      mm <- X <- model.matrix(object)
      mmDone <- TRUE
      offset <- object$offset
    }
    else {
      Terms <- delete.response(tt)
      m <- model.frame(Terms,
                       newdata,
                       na.action = na.action,
                       xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
      offset <- rep(0, nrow(X))
      if (!is.null(off.num <- attr(tt, "offset")))
        for (i in off.num)
          offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
      if (!is.null(object$call$offset))
        offset <- offset + eval(object$call$offset, newdata)
      mmDone <- FALSE
    }

    n <- length(object$residuals)
    beta <- object$coefficients
    p <- length(object$coefficients)
    p1 <- seq_len(p)
    piv <- p1

    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])

    if (!is.null(offset))
      predictor <- predictor + offset

    interval <- match.arg(interval)
    predictor.ci <- NULL
    if (interval == "confidence") {
      a <- (1 - level) / 2
      a <- c(a, 1 - a)

      predictor.ci <-
        data.frame(matrix(
          NA,
          nrow = nrow(X),
          ncol = ncol(object$results$bootstrap$coefficients)
        ))

      for (b in 1:ncol(predictor.ci)) {
        predictor.ci[, b] <-
          drop(X[, piv, drop = FALSE] %*% object$results$bootstrap$coefficients[, b])
        #X %*% object$results$bootstrap$coefficients[, b]
        if (!is.null(offset))
          predictor.ci[, b] <- predictor.ci[, b] + offset
      }
    }

    type <- match.arg(type)

    if (type == "terms") {
      if (!mmDone) {
        mm <- model.matrix(object)
        mmDone <- TRUE
      }
      aa <- attr(mm, "assign")
      ll <- attr(tt, "term.labels")
      hasintercept <- attr(tt, "intercept") > 0L
      if (hasintercept)
        ll <- c("(Intercept)", ll)
      aaa <- factor(aa, labels = ll)
      asgn <- split(order(aa), aaa)
      if (hasintercept) {
        asgn$"(Intercept)" <- NULL
        avx <- colMeans(mm)
        termsconst <- sum(avx[piv] * beta[piv])
      }
      nterms <- length(asgn)
      if (nterms > 0) {
        predictor <- matrix(ncol = nterms, nrow = NROW(X))
        dimnames(predictor) <- list(rownames(X), names(asgn))

        if (hasintercept)
          X <- sweep(X, 2L, avx, check.margin = FALSE)
        unpiv <- rep.int(0L, NCOL(X))
        unpiv[piv] <- p1

        for (i in seq.int(1L, nterms, length.out = nterms)) {
          iipiv <- asgn[[i]]
          ii <- unpiv[iipiv]
          iipiv[ii == 0L] <- 0L
          predictor[, i] <-
            if (any(iipiv > 0L))
              X[, iipiv, drop = FALSE] %*% beta[iipiv]
          else
            0
        }
        if (!is.null(terms)) {
          predictor <- predictor[, terms, drop = FALSE]
        }
      } else {
        predictor <- ip <- matrix(0, n, 0L)
      }
      attr(predictor, 'constant') <- if (hasintercept)
        termsconst
      else
        0
    }

    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
      predictor <- napredict(na.act, predictor)
    }

    if (interval != "none") {
      predictor <-
        cbind(predictor, t(apply(predictor.ci, 1, quantile, a, na.rm = TRUE)))
      colnames(predictor) <- c("fit", "lwr", "upr")
      predictor
    } else {
      predictor
    }
  }
