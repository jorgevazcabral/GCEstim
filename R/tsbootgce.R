#' Time series bootstrap Cross entropy estimation
#'
#' This generic function fits a linear regression model using bootstrapped time
#' series via generalized cross entropy.
#'
#' @param formula a "formula" describing the linear model to be fit. For details
#' see \code{\link[stats]{lm}} and \code{\link[dynlm]{dynlm}}.
#' @param data A \code{\link[base]{data.frame}} (or object coercible by
#' \code{\link[base]{as.data.frame}} to a data frame) or time series object
#' (e.g., \code{\link[stats]{ts}} or \code{\link[zoo]{zoo}}), containing the
#' variables in the model.
#' @param trim The trimming proportion (see \code{\link[meboot]{meboot}}).
#' The default is \code{trim = 0.05}.
#' @param reps The number of replicates to generate (see
#' \code{\link[meboot]{meboot}}). The default is \code{reps = 1000}.
#' @param start The time of the first observation. Either a single number
#'  or a vector of two numbers (the second of which is an integer), which
#'  specify a natural time unit and a (1-based) number of samples into the time
#'  unit (see \code{\link[stats]{ts}}).
#' @param end The time of the last observation, specified in the same way as
#' \code{start} (see \code{\link[stats]{ts}}).
#' @param coef.method Method used to estimate the coefficients. One of
#' \code{c("mode", "median")}. for \code{"mode"} see \code{\link[hdrcde]{hdr}}
#' @inheritParams lmgce
#'
#' @details
#'
#' The \code{tsbootgce} function fits several linear regression models via
#' generalized cross entropy in replicas of time series obtained using
#' \code{\link[meboot]{meboot}}. Models for \code{\link{tsbootgce}} are specified
#' symbolically (see \code{\link[stats]{lm}} and \code{\link[dynlm]{dynlm}}).
#'
#' @return
#' \code{tsbootgce} returns an object of \code{\link[base]{class}} \code{tsbootgce}.
#' The generic accessory functions \code{\link{coef.tsbootgce}},
#'  \code{\link{confint.tsbootgce}} and \code{\link{plot.tsbootgce}} extract
#'  various useful features of the value returned by \code{object} of class
#'  \code{tsbootgce}.
#'
#'  An object of \code{\link[base]{class}} \code{tsbootgce} is a list containing at
#'  least the following components:
#'
#' \item{call}{the matched call.}
#' \item{coefficients}{a named data frame of coefficients determined by
#'  \code{coef.method}.}
#' \item{data.ts}{\code{ts} object.}
#' \item{error}{loss function (error) used for the selection of the support
#' spaces.}
#' \item{error.measure}{in sample error for the selected support space.}
#' \item{fitted.values}{the fitted mean values.}
#' \item{frequency}{see \code{link[zoo]{zoo}}.}
#' \item{index}{see \code{link[zoo]{zoo}}.}
#' \item{lmgce}{\code{lmgce} object.}
#' \item{meboot}{\code{meboot} replicates.}
#' \item{model}{the model frame used.}
#' \item{nep}{normalized entropy of the signal of the model.}
#' \item{nepk}{normalized entropy of the signal of each coefficient.}
#' \item{residuals}{the residuals, that is response minus fitted values.}
#' \item{results}{a list containing the bootstrap results: "coef.matrix", a named
#'  data frame of all the coefficients; "nepk.matrix", a named data frame of all
#'  the normalized entropy values of each parameter; "nep.vector", a vector of
#'  all the normalized entropy values of the model.}
#' \item{seed}{the seed used.}
#' \item{terms}{the \code{\link[stats]{terms}} object used.}
#' \item{x}{if requested (the default), the model matrix used.}
#' \item{xlevels}{(only where relevant) a record of the levels of the factors
#' used in fitting.}
#' \item{y}{if requested (the default), the response used.}
#'
#' @seealso
#' The generic functions \code{\link{plot.tsbootgce}}, \code{\link{print.tsbootgce}},
#'  and \code{\link{coef.tsbootgce}}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @references
#' Golan, A., Judge, G. G. and Miller, D. (1996)
#' \emph{Maximum entropy econometrics : robust estimation with limited data.}
#' Wiley.\cr
#'
#' Golan, A. (2008)
#' \emph{Information and Entropy Econometrics — A Review and Synthesis.}
#' Foundations and Trends® in Econometrics, 2(1–2), 1–145.
#' \doi{10.1561/0800000004}\cr
#'
#' Golan, A. (2017)
#' \emph{Foundations of Info-Metrics: Modeling, Inference, and Imperfect Information (Vol. 1).}
#' Oxford University Press.
#' \doi{10.1093/oso/9780199349524.001.0001}\cr
#'
#' Hyndman, R.J. (1996)
#' \emph{Computing and graphing highest density regions.}
#'  American Statistician, 50, 120-126.
#'  \doi{10.2307/2684423}\cr
#'
#' Pukelsheim, F. (1994)
#' \emph{The Three Sigma Rule.}
#' The American Statistician, 48(2), 88–91.
#' \doi{10.2307/2684253}\cr
#'
#' Vinod, H. D., & Lopez-de-Lacalle, J. (2009).
#' \emph{Maximum Entropy Bootstrap for Time Series: The meboot R Package.}
#'  Journal of Statistical Software, 29(5), 1–19.
#' \doi{10.18637/jss.v029.i05}
#'
#' @examples
#' \donttest{
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' res.tsbootgce
#' }
#'
#' @importFrom stats density ts
#' @importFrom utils head tail
#' @importFrom zoo merge.zoo
#' @export

tsbootgce <- function(formula,
                      data,
                      subset,
                      na.action,
                      offset,
                      contrasts = NULL,
                      trim = 0.05,
                      reps = 1000,
                      start = NULL,
                      end = NULL,
                      coef.method = c("mode", "median"),
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
                      support.method.ridge.lambda.min = 10^-3,
                      support.method.ridge.lambda.max = 10^3,
                      support.method.ridge.lambda.n = 100,
                      support.method.ridge.standardize = TRUE,
                      support.method.ridge.penalize.intercept = TRUE,
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

  if(!is.ts(data))
    stop("argument `data` must be a \"ts\" object")

  #suppressMessages(require(zoo, quietly = TRUE))

  x = TRUE
  y = TRUE
  model = TRUE
  coef.method <- match.arg(coef.method)
  support.method <- match.arg(support.method)
  errormeasure <- match.arg(errormeasure)
  errormeasure.which <- match.arg(errormeasure.which)
  method <- match.arg(method)
  caseGLM <- match.arg(caseGLM)
  boot.method <- match.arg(boot.method)

  Zenv <- new.env(parent = environment(formula))
  assign("dynformula", function(x) structure(x, class = unique(c("dynformula",
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
    res <- zoo::zoo(seq_along(zoo::index(x))/freq, zoo::index(x), frequency = ofreq)
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
  assign("model.frame.dynformula", function(formula, data = NULL,
                                            subset = NULL, na.action = na.omit, drop.unused.levels = FALSE,
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
    for (i in 1:NCOL(mf)) attr(mf[, i], "index") <- attr(mf[,
                                                            i], "index")[-as.vector(mfna)]
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
      attr(mf, "na.action") <- structure(mfna[as.vector(mfna) >=
                                                start & as.vector(mfna) <= end], class = class(mfna))
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

  res$lmgce <-
    lmgce(
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
      support.method.ridge.penalize.intercept = support.method.ridge.penalize.intercept,
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

  ## end lmgce

  mf_ts <- mf

  me_mf_ts <- list()

  if (!is.null(seed))
    set.seed(seed)

  for (i in 1:ncol(mf_ts)) {
    me_mf_ts[[i]] <-
      data.frame(meboot::meboot(x = mf_ts[, i], reps = reps, trim = trim)$ensemble)
  }

  coef.df <- data.frame(matrix(
    NA,
    nrow = ncol(mf) - 1 + attr(mt, "intercept"),
    ncol = reps
  ))
  colnames(coef.df) <- sprintf(paste0("rep_%0", floor(log10(reps)) + 1, "d"), 1:reps)
  rownames(coef.df) <- c(ifelse(attr(mt, "intercept") == 1, "(Intercept)"), attr(mt, "term.labels"))

  nepk.df <- coef.OLS.df <- coef.df

  coefnepk.list <- list(coef = vector(mode = "list",
                                      length = cv.nfolds),
                        nepk = vector(mode = "list",
                                      length = cv.nfolds))

  coef.OLS.list <- list(coef = vector(mode = "list",
                                      length = cv.nfolds))

  for (i in 1:cv.nfolds){
    coefnepk.list$coef[[i]] <- coef.df
    coefnepk.list$nepk[[i]] <- nepk.df
    coef.OLS.list$coef[[i]] <- coef.df
  }

  nep <- NULL

  for (i in 1:reps) {
    data_me_mf_ts <- data.frame(sapply(me_mf_ts, function(x){x[, i]}))

    res.aux.lmgce <-
      suppressWarnings(
        lmgce(
          formula = as.formula(paste0("X1 ~",
                                      ifelse(attr(mt,
                                                  "intercept") == 1,
                                             " 1 + ",
                                             " -1 + "),
                                      ".")),
          data = data_me_mf_ts,
          contrasts = contrasts,
          model = TRUE,
          x = TRUE,
          y = TRUE,
          cv = cv,
          cv.nfolds = cv.nfolds,
          errormeasure = errormeasure,
          errormeasure.which =
            {if (errormeasure.which == "1se") {
                "min"
              } else {errormeasure.which}
            },
          support.method = "standardized",#support.method,
          support.method.ridge.penalize.intercept = support.method.ridge.penalize.intercept,
          support.signal =
            #{if (is.null(support.signal)) {
            #  res$lmgce$support.matrix
            #  } else {
            #    support.signal}},
            res$lmgce$support.matrix,
          support.signal.vector = support.signal.vector,
          support.signal.vector.min = support.signal.vector.min,
          support.signal.vector.max = support.signal.vector.max,
          support.signal.vector.n = support.signal.vector.n,
          support.signal.points = support.signal.points,
          support.noise = support.noise,
          support.noise.points = support.noise.points,
          weight = weight,
          twosteps.n =
            {if (is.null(support.signal)) {
              0
            } else {
              twosteps.n}},
          method = method,
          caseGLM = caseGLM,
          boot.B = 0,
          boot.method = boot.method,
          seed = seed,
          OLS = TRUE,
          verbose = verbose
        ))

    coef.df[, i] <- coef(res.aux.lmgce)
    nepk.df[, i] <- NormEnt(res.aux.lmgce, model = FALSE)
    nep[i] <- NormEnt(res.aux.lmgce, model = TRUE)

    coef.OLS.df[, i] <- coef(res.aux.lmgce$results$OLS$res)

    if (isTRUE(cv)) {
    for (j in 1:cv.nfolds){
      coefnepk.list$coef[[j]][, i] <-
        res.aux.lmgce$results$cv$repeats1[[j]]$coefficients
      coefnepk.list$nepk[[j]][, i] <-
        res.aux.lmgce$results$cv$repeats1[[j]]$nepk
      coef.OLS.list$coef[[j]][, i] <- res.aux.lmgce$results$OLS$matrix.coef[, j]
    }
    }
  }

  if (coef.method == "mode") {
    if (!is.null(seed))
      set.seed(seed)
    coef.matrix.cv <-
      sapply(coefnepk.list$coef,
             function(x){apply(x,
                               1,
                               function(x) {
                                 hdrcde::hdr(x = as.numeric(x), prob = 95)[[2]]} )})
    coef.OLS.matrix.cv <-
      sapply(coef.OLS.list$coef,
             function(x){apply(x,
                               1,
                               function(x) {
                                 hdrcde::hdr(x = as.numeric(x), prob = 95)[[2]]} )})
  } else {
    coef.matrix.cv <-
      sapply(coefnepk.list$coef,
             function(x){apply(x,
                               1,
                               median)})
    coef.OLS.matrix.cv <-
      sapply(coef.OLS.list$coef,
             function(x){apply(x,
                               1,
                               median)})
  }

  if (!is.null(seed))
    set.seed(seed)

  auxfolds = cut(seq(1, nrow(x)),
                 breaks = cv.nfolds,
                 labels = FALSE)
  change_order <- sample(nrow(x))

  error.cv <- error.OLS.cv <- NULL

  for (cv.n in 1:cv.nfolds) {
    y.cv = y[change_order][auxfolds != cv.n]
    X.cv = x[change_order, ][auxfolds != cv.n,]
    error.cv[cv.n] <- accmeasure(y.cv,
                                 X.cv %*% coef.matrix.cv[, cv.n],
                                 errormeasure)
    error.OLS.cv[cv.n] <- accmeasure(y.cv,
                                     X.cv %*% coef.OLS.matrix.cv[, cv.n],
                                     errormeasure)
  }
  res$na.action <- attr(mf, "na.action")
  res$offset <- offset
  res$contrasts <- attr(x, "contrasts")
  res$xlevels <- .getXlevels(mt, mf)
  res$call <- cl
  res$terms <- mt
  if (model)
    res$model <- mf
  if (ret.x)
    res$x <- x
  if (ret.y)
    res$y <- y
  res$index <- zoo::index(mf1)
  res$frequency <- frequency(mf1)
  # res$twostage <- twostage
  # if (twostage) {
  #   res$formula <- ff
  #res$residuals <- y1 - as.vector(x1 %*% res$coefficients)
  #res$fitted.values <- y1 - res$residuals
  # }
  res$results$bootstrap$coef.matrix <- coef.df
  res$results$bootstrap$nepk.matrix <- nepk.df
  res$results$bootstrap$nep.vector <- nep
  res$results$bootstrap$coef.matrix.OLS <- coef.OLS.df
  if (coef.method == "mode") {
    if (!is.null(seed)) set.seed(seed)
    res$coefficients <- apply(res$results$bootstrap$coef.matrix, 1, function(x) {
      hdrcde::hdr(x = as.numeric(x), prob = 95)[[2]]
    })
    res$coefficients.OLS <- apply(res$results$bootstrap$coef.matrix.OLS, 1, function(x) {
      hdrcde::hdr(x = as.numeric(x), prob = 95)[[2]]
    })
  } else {
    res$coefficients <- apply(res$results$bootstrap$coef.matrix, 1, median)
    res$coefficients.OLS <- apply(res$results$bootstrap$coef.matrix.OLS, 1, median)
  }
  res$coefficients.cv <- coef.matrix.cv
  res$res.aux.lmgce <- res.aux.lmgce
  res$fitted.values <- model.matrix(mt, mf, contrasts) %*% res$coefficients
  res$residuals <- model.response(mf, "numeric") - res$fitted.values
  res$fitted.values.OLS <- model.matrix(mt, mf, contrasts) %*% res$coefficients.OLS
  res$residuals.OLS <- model.response(mf, "numeric") - res$fitted.values.OLS
  res$nep <- median(res$results$bootstrap$nep.vector)
  res$nepk <- apply(res$results$bootstrap$nepk.matrix, 1, median)
  res$data.ts <- mf_ts
  names(me_mf_ts) <- colnames(mf)
  res$meboot <- me_mf_ts
  res$error <- errormeasure
  res$error.measure <-
    accmeasure(model.response(mf, "numeric"),
                        res$fitted.values,
                        errormeasure)
  res$error.measure.OLS <-
    accmeasure(model.response(mf, "numeric"),
                        res$fitted.values.OLS,
                        errormeasure)
  res$error.measure.cv.mean <- mean(error.cv)
  res$error.measure.cv.sd <- sd(error.cv)
  res$error.measure.OLS.cv.mean <- mean(error.OLS.cv)
  res$error.measure.OLS.cv.sd <- sd(error.OLS.cv)
  res$seed <- seed
  }
  class(res) <- "tsbootgce"
  return(res)
}

#' Auxiliary function
#'
#' Returns time series (from zoo)
#'
#' @param x .
#' @param ... Additional arguments.
#'
#' @return time series
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

astszoo <- function (x, ...)
{
  if (zoo::is.regular(x)) {
    attr(x, "frequency") <- frequency(x)
    return(astszooreg(x, ...))
  }
  else {
    warning(paste(sQuote("x"), "does not have an underlying regularity"))
    return(ts(zoo::coredata(x)))
  }
}

#' Auxiliary function
#'
#' Returns time series (from zoo)
#'
#' @param x .
#' @param ... Additional arguments.
#'
#' @return time series
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

astszooreg <- function (x, ...)
{
  freq <- frequency(x)
  deltat <- 1/freq
  tt <- as.numeric(time(x))
  round. <- if (isTRUE(all.equal(c(deltat, tt), round(c(deltat,
                                                        tt))))) {
    function(x) floor(x + 0.5)
  }
  else {
    function(x) deltat * floor(x/deltat + 0.5)
  }
  tt <- round.(tt)
  tt2 <- round.(seq(head(tt, 1), tail(tt, 1), deltat))
  fill <- list(...)$fill
  if (is.null(fill))
    fill <- NA
  xx <- merge(zoo::zoo(zoo::coredata(x), tt), zoo::zoo(, tt2), fill = fill)
  ts(zoo::coredata(xx), start = tt[1], frequency = freq)
}

#' Auxiliary function
#'
#' Invert Box Cox transformation (from hdrcde)
#'
#' @param x a numeric vector or time series
#' @param lambda transformation parameter
#'
#' @return a numeric vector of the same length as x.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
InvertBoxCox <- function (x, lambda)
{
  if (is.list(x))
    x <- x[[1]]
  if (lambda == 0)
    exp(x)
  else (x * lambda + 1)^(1/lambda)
}
