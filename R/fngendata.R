#' Data generating function
#'
#' Generates data
#'
#' @param n Number of individuals.
#' @param bin.k Number of binary variables not used for generating y.
#' @param bin.prob A vector of probabilities with length equal to \code{bin.k}.
#' @param cont.k Number of continuous variables not used for generating y.
#' @param y.gen.bin.k Number of binary variables used for generating y.
#' @param y.gen.bin.beta A vector of coefficients with length equal to
#' \code{bin.k}used to generate y.
#' @param y.gen.bin.prob A vector of probabilities with length equal
#' to \code{y.gen.bin.k}.
#' @param y.gen.cont.beta A vector of coefficients with length equal to
#' \code{cont.k} used to generate y.
#' @param y.gen.cont.mod.k Experimental
#' @param y.gen.cont.mod.beta Experimental
#' @param y.gen.bin.mod.prob Experimental
#' @param y.gen.cont.sp.k Experimental
#' @param y.gen.cont.sp.groups Experimental
#' @param y.gen.cont.sp.rho Experimental
#' @param y.gen.cont.sp.dif Experimental
#' @param intercept.beta Value for the constant used to generate y.
#' @param Xgenerator.method Method used to generate X data ( \code{"simstudy"}
#' or \code{"svd"}).
#' @param corMatrix A positive number for alphad
#' (see \code{\link[clusterGeneration]{rcorrmatrix}}), NULL or a correlation
#'  matrix to be used when \code{Xgenerator} is \code{"simstudy"}.
#' @param rho Correlation coefficient, \code{-1 <= rho <= 1}. Use when
#' \code{Xgenerator} is \code{"simstudy"} and \code{corMatrix} is NULL.
#' @param corstr correlation structure (\code{"ind"}, \code{"cs"} or
#' \code{"ar1"}) (see \code{\link[simstudy]{genCorData}}) to be used when
#' \code{Xgenerator} is \code{"simstudy"} and \code{corMatrix} is NULL.
#' @param condnumber A value for the condition number of the X matrix to be used
#'  when \code{Xgenerator} is \code{"svd"}.
#' @param mu The mean of the variables. To be used when all variables have the
#' same mean.
#' @param muvect A vector of means. To be used when variables have different
#' means. The length of \code{muvect} must be \code{k}.
#' @param sd Standard deviation of the variables. To be used when all variables
#' have the same standard deviation.
#' @param sdvect A vector of standard deviations. To be used when variables have
#' different standard deviations. The length of \code{sdvect} must be \code{k}.
#' @param error.dist Distribution of the error. \code{"normal"} for normal
#' distribution or \code{"t"} for t-student distribution.
#' @param error.dist.mean Mean value used when \code{error.dist} is
#' \code{"normal"}.
#' @param error.dist.sd Standard deviation value used when \code{error.dist} is
#' \code{"normal"}.
#' @param error.dist.snr Signal to noise ratio. If not \code{NULL}, the value of
#' \code{error.dist.sd} will be ignored and it will be determined accordingly.
#' @param error.dist.df Degrees of freedom used when \code{error.dist} is
#' \code{"t"}.
#' @param dataframe Logical. If \code{TRUE}, the default, returns a
#' \code{data.frame} else returns a \code{list}.
#' @param seed A seed for reproducibility.
#'
#' @return A \code{data.frame} or a \code{list} composed of a matrix of
#' independent variables values (X), a vector of the dependent variable values
#' (y), a vector of coefficient values (coefficients), a vector of non-zero
#' coefficients (y.coefficients), and a vector of the error values (epsilon).
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#'
#' dataGCEstim <- fngendata(
#'   n = 100, cont.k = 2,
#'   y.gen.cont.beta = c(3, 6, 9),
#'   intercept.beta = 1,
#'   Xgenerator.method = "svd", condnumber = 50,
#'   mu = 0, sd = 1,
#'   error.dist = "normal", error.dist.mean = 0, error.dist.snr = 5,
#'   dataframe = TRUE, seed = 230676)
#'
#' summary(dataGCEstim)
#'
#' @export
#' @importFrom simstudy genCorData
#' @importFrom clusterGeneration rcorrmatrix

fngendata <- function(n,
                      bin.k = 0,
                      bin.prob = NULL,
                      cont.k = 5,
                      y.gen.bin.k = 0,
                      y.gen.bin.beta = NULL,
                      y.gen.bin.prob = NULL,
                      y.gen.cont.beta = c(2, 4, 6, 8, 10),
                      y.gen.cont.mod.k = 0,
                      y.gen.cont.mod.beta = matrix(c(-2, 2),
                                                   1,
                                                   2,
                                                   byrow = TRUE),
                      y.gen.bin.mod.prob = c(0.5),
                      y.gen.cont.sp.k = 0,
                      y.gen.cont.sp.groups = 2,
                      y.gen.cont.sp.rho = 0.2,
                      y.gen.cont.sp.dif = 1,
                      intercept.beta = 0,
                      Xgenerator.method = "simstudy",
                      corMatrix = 100,
                      rho = NULL,
                      corstr = NULL,
                      condnumber = 1,
                      mu = 0,
                      muvect = NULL,
                      sd = 1,
                      sdvect = NULL,
                      error.dist = "normal",
                      error.dist.mean = 0,
                      error.dist.sd = 1,
                      error.dist.snr = NULL,
                      error.dist.df = 2,
                      dataframe = TRUE,
                      seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- intercept.beta

  y.gen.cont.k <- length(y.gen.cont.beta)
  if ((cont.k + y.gen.cont.k) != 0) {
    betas.cont <- NULL

    if (cont.k != 0) {
      betas.cont <- c(betas.cont, rep(0, cont.k))
    }
    if (y.gen.cont.k != 0) {
      betas.cont <- c(betas.cont, y.gen.cont.beta)
    }

    if (is.null(muvect)) {muvect <- rep(mu, cont.k + y.gen.cont.k)}

    if (is.null(sdvect)) {sdvect <- rep(sd, cont.k + y.gen.cont.k)}

    if (Xgenerator.method == "svd") {
      if (cont.k + y.gen.cont.k > 1) {
        rm <- matrix(0, n, cont.k + y.gen.cont.k, byrow = FALSE)

        for (i in 1:(cont.k + y.gen.cont.k)) {
          # if (!is.null(seed)) {
          #   set.seed(seed + i)
          # }
          rm[, i] <- rnorm(n, mean = muvect[i], sd = sdvect[i])
        }

        res_svd <- svd(rm, nu = n, nv = (cont.k + y.gen.cont.k))

        d <- matrix(0, n, cont.k + y.gen.cont.k)

        diag(d) <- c(2 / (1 + condnumber),
                     rep(1, min(n, cont.k + y.gen.cont.k) - 2),
                     2 * condnumber / (1 + condnumber))

        X.cont <- res_svd$u %*% d %*% t(res_svd$v)
      } else {
        X.cont <- rnorm(n, mean = mu, sd = sd)
      }

    } else if (Xgenerator.method == "simstudy") {

      if (!is.null(corMatrix)) {
        if (length(corMatrix) == 1) {
          # if (!is.null(seed)) {
          #   set.seed(seed)
          # }
          corr_aux <-
            clusterGeneration::rcorrmatrix(cont.k + y.gen.cont.k,
                                           alphad = as.numeric(corMatrix))
        } else {
          corr_aux <- corMatrix
        }
      } else {
        corr_aux <- NULL
      }

      # if (!is.null(seed)) {
      #   set.seed(seed)
      # }
      X.cont <- as.matrix(
        simstudy::genCorData(
          n = n,
          mu = muvect,
          sigma = sdvect,
          corMatrix = corr_aux,
          rho = rho,
          corstr = corstr
        )[, -1]
      )

      attr(X.cont, "dimnames") <- NULL
    }

    if (cont.k + y.gen.cont.k > 1) {
      y <- y + X.cont %*% betas.cont
    } else {
      y <- y + X.cont * betas.cont
    }

  } else {
    X.cont <- NULL
    betas.cont <- NULL
  }

  if ((bin.k + y.gen.bin.k) != 0) {
    betas.bin <- NULL

    if (bin.k != 0) {
      betas.bin <- c(betas.bin, rep(0, bin.k))
    }
    if (y.gen.bin.k != 0) {
      betas.bin <- c(betas.bin, y.gen.bin.beta)
    }

    X.bin <- matrix(0, n, bin.k + y.gen.bin.k)

    for (i in 1:(bin.k + y.gen.bin.k)) {
      # if (!is.null(seed)) {
      #   set.seed(seed + i)
      # }
      X.bin[, i] <-
        rbinom(n, 1, ifelse(i <= bin.k,
                            bin.prob[i],
                            y.gen.bin.prob[i - bin.k]))
    }

    y <- y + X.bin %*% betas.bin

  } else {
    X.bin <- NULL
    betas.bin <- NULL
  }

  if (y.gen.cont.mod.k != 0) {
    betas.cont.mod <- rep(0, y.gen.cont.mod.k)
    betas.bin.mod <- rep(0, y.gen.cont.mod.k)
    X.cont.mod <- matrix(0, n, y.gen.cont.mod.k)
    X.bin.mod <- matrix(0, n, y.gen.cont.mod.k)
    for (i in 1:y.gen.cont.mod.k) {
      # if (!is.null(seed)) {
      #   set.seed(seed + 1000 * i)
      # }
      X.cont.mod[, i] <- rnorm(n, 0, 1)
      X.bin.mod[, i] <- rbinom(n, 1, y.gen.bin.mod.prob[i])
    }

    y <-
      y + X.cont.mod * ifelse(X.bin.mod[, 1] == 0,
                              y.gen.cont.mod.beta[1, 1],
                              y.gen.cont.mod.beta[1, 2])

  } else {
    X.cont.mod <- NULL
    X.bin.mod <- NULL
    betas.cont.mod <- NULL
    betas.bin.mod <- NULL
  }

  if (y.gen.cont.sp.k != 0) {
    X.cont.sp <-
      bayestestR::simulate_simpson(
        n = n / y.gen.cont.sp.groups,
        r = y.gen.cont.sp.rho,
        groups = y.gen.cont.sp.groups,
        difference = y.gen.cont.sp.dif,
        group_prefix = ""
      )

    coef.sp <- coef(lm(V2 ~ ., X.cont.sp))

    X.bin.sp <- as.numeric(X.cont.sp[, 3]) - 1
    X.cont.sp <- X.cont.sp[, 1]

    intercept.beta.sp <- coef.sp[[1]]
    betas.cont.sp <- coef.sp[[2]]
    betas.bin.sp <- coef.sp[[3]]
    y <-
      y + intercept.beta.sp + betas.cont.sp * X.cont.sp +
      betas.bin.sp * X.bin.sp

  } else {
    X.cont.sp <- NULL
    X.bin.sp <- NULL
    intercept.beta.sp <- 0
    betas.cont.sp <- NULL
    betas.bin.sp <- NULL
  }

  X <-
    cbind(X.cont, X.bin, X.bin.mod, X.cont.mod, X.bin.sp, X.cont.sp)
  colnames(X) <- paste0("X", sprintf('%0.3d', 1:ncol(X)), sep = "")

  if (error.dist == "normal") {
    if (is.null(error.dist.snr)) {
      epsilon <- rnorm(n, error.dist.mean, error.dist.sd)
    } else {
      epsilon <- rnorm(n, error.dist.mean, sqrt(var(y) / error.dist.snr))
    }
  } else {
    epsilon <- rt(n, error.dist.df)
  }

  y <- y + epsilon

  if (isTRUE(dataframe)) {
    return(data.frame(X, y))
  } else {
    return(list(
      X = X,
      y = y,
      coefficients = c(
        intercept.beta + intercept.beta.sp,
        betas.cont,
        betas.bin,
        betas.bin.mod,
        betas.cont.mod,
        betas.bin.sp,
        betas.cont.sp
      ),
      y.coefficients = c(
        intercept.beta + intercept.beta.sp,
        if (y.gen.cont.k == 0) {
          NULL
        } else {
          y.gen.cont.beta
        },
        betas.bin,
        betas.bin.mod,
        betas.bin.sp,
        betas.cont.sp
      ),
      epsilon = epsilon
    ))
  }
}
