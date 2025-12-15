# Generic functions for the \code{\link{tsbootgce}} class

#' Print \code{\link{tsbootgce}} object
#'
#' Print \code{\link{tsbootgce}} object
#'
#' @param x fitted \code{lmgce} object.
#' @param digits significant digits in printout.
#' @param ... additional print arguments.
#'
#' @return A small summary of a \code{\link{tsbootgce}} object is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
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
#' @method print tsbootgce
#' @export

print.tsbootgce <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")
  if (length(coef(x))) {

    res.print <-
      as.matrix(rbind(coef(x),
                      coef(x$lmgce),
                      coef(x$lmgce$results$OLS$res)),
                nrow = 2)

    dimnames(res.print) <- list(c("tsbootgce", "lmgce", "lm"),
                                row.names(x$results$bootstrap$coef.matrix))

    cat("Coefficients:\n")
    print.default(format(res.print, digits = digits),
                  print.gap = 2L,
                  quote = FALSE)
  }
  else
    cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' Extract \code{\link{tsbootgce}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{tsbootgce}} object
#'
#' @param object Fitted \code{\link{tsbootgce}} model object.
#' @param which  The default is \code{which = NULL} and returns the coefficients
#' defined in the argument \code{coef.method} from the \code{tsbootgce} object.
#' Can be set as "mode" or "median" and the mode and median coefficients will be
#' computed, respectively (see \code{\link[hdrcde]{hdr}}).
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = object$seed}.
#' @param ... Additional arguments.
#'
#' @return Returns the coefficients from a \code{tsbootgce} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' coef(res.tsbootgce)
#' }
#'
#' @method coef tsbootgce
#' @importFrom stats coef
#' @export

coef.tsbootgce <- function(object, which = NULL, seed = object$seed, ...)
{
  if (is.null(which))
    object$coefficients
  else
    if (which == "mode") {
      if (!is.null(seed))
        set.seed(seed)
      apply(object$results$bootstrap$coef.matrix, 1, function(x) {
        hdrcde::hdr(x = as.numeric(x), prob = 95)[[2]]
      })
    } else {
      apply(object$results$bootstrap$coef.matrix, 1, median)
    }
}

#' Extract \code{\link{tsbootgce}} Model Coefficients
#'
#' Extract coefficients from a \code{\link{tsbootgce}} object
#'
#' @param object Fitted \code{\link{tsbootgce}} model object.
#' @param which  The default is \code{which = NULL} and returns the coefficients
#' defined in the argument \code{coef.method} from the \code{tsbootgce} object.
#' Can be set as "mode" or "median" and the mode and median coefficients will be
#' computed, respectively (see \code{\link[hdrcde]{hdr}}).
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = object$seed}.
#' @param ... Additional arguments.
#'
#' @return Returns the coefficients from a \code{tsbootgce} object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @rdname coefficients.tsbootgce
#' @examples
#' \donttest{
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' coefficients(res.tsbootgce)
#' }
#'
#' @method coefficients tsbootgce
#' @importFrom stats coefficients
#' @export

coefficients.tsbootgce <- coef.tsbootgce

#' Confidence Intervals for \code{\link{tsbootgce}} Model Parameters and Normalized
#'  Entropy
#'
#' Computes confidence intervals for one or more parameters or Normalized Entropy
#' in a \code{\link{tsbootgce}} fitted model.
#'
#' @param object Fitted \code{\link{tsbootgce}} model object.
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing,
#' all parameters are considered.
#' @param level the confidence level required. The default is
#'  \code{level = 0.95}.
#' @param which One of \code{c("estimates", "NormEnt")}. The default is
#' \code{which = "estimates"}.
#' @param method method used to compute the interval. One of
#' c("hdr", "percentile", "basic"). The default is \code{method = "hdr"}
#' (see \code{\link[hdrcde]{hdr}}).
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = object$seed}.
#' @param ... additional arguments.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter. Generally, these will be labelled as (1-level)/2 and
#' 1 - (1-level)/2 in percentage (by default 2.5 percent and 97.5 percent).
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' \donttest{
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' confint(res.tsbootgce, method = "percentile")
#'
#' confint(res.tsbootgce, which = "NormEnt", level = 0.99)
#'
#' confint(res.tsbootgce, parm = c("L(GDP, 1)"), level = 0.99)
#' }
#'
#' @method confint tsbootgce
#' @importFrom stats confint
#' @export

confint.tsbootgce <- function(object,
                              parm,
                              level = 0.95,
                              which = c("estimates", "NormEnt"),
                              method = c("hdr", "percentile", "basic"),
                              seed = object$seed,
                              ...)
{
  if (!inherits(object, "tsbootgce"))
    stop("use only with \"tsbootgce\" objects")

  #which <- match.arg(which)
  which <- which[1]
  #method <- match.arg(method)
  method <- method[1]

  if (which == "estimates")
    cf <- coef(object)
  else
    cf <- NormEnt(object, model = FALSE)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
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

  if (which == "estimates")
    mcf <- object$results$bootstrap$coef.matrix
  else
    mcf <- object$results$bootstrap$nepk.matrix

  if (method == "hdr") {
    if (!is.null(seed))
      set.seed(seed)

    res.conf <-
      apply(mcf, 1, function(x) {
        hdrcde::hdr(x = as.numeric(x), prob = level * 100)[[1]]})

    if (is.list(res.conf)) {
    confint.matrix <-
      data.frame(matrix(NA,
                        ncol = max(sapply(res.conf, function(x){length(x)})),
                        nrow = length(res.conf)))
    for (i in 1:length(res.conf)) {
      for (j in 1:length(res.conf[[i]])) {
        confint.matrix[i, j] <- res.conf[[i]][j]
      }
    }

    rownames(confint.matrix) <- names(res.conf)
    colnames(confint.matrix) <- paste(rep(paste0(a * 100, "%"),
                                          ncol(confint.matrix)/2),
                                      if (ncol(confint.matrix) > 2) {
                                      rep(c(1:(ncol(confint.matrix)/2)),
                                          each = 2)})
    return(as.matrix(confint.matrix))
    } else {
      res.conf <- t(res.conf)
      colnames(res.conf) <- paste0(a * 100, "%")
      return(res.conf)
    }
  }
  else
    if (method == "basic") {
      aux.conf.int <-
        2 * cf[parm] -
        t(apply(mcf[parm, ], 1, quantile, rev(a)))
      colnames(aux.conf.int) <- paste0(a * 100, "%")
      return(aux.conf.int)
    } else if (method == "percentile") {
      return(t(apply(mcf[parm, ], 1, quantile, a)))
    }
}


#' Plot Diagnostics for a \code{\link{tsbootgce}} object
#'
#' Three plots (selectable by \code{which}) are currently available to
#' evaluate a \code{\link{tsbootgce}} object.
#'
#' @param x Fitted \code{tsbootgce} object.
#' @param which Integers from 1 to 3. The default is \code{which = c(1,2)}.
#' @param group Boolean value. If \code{group = TRUE}, the default, plots are
#' grouped in one image.
#' @param group.ncol Number of columns (see \code{\link[ggpubr]{ggarrange}}).
#'  The default is \code{group.ncol = NULL}.
#' @param group.nrow Number of rows. (see \code{\link[ggpubr]{ggarrange}}).
#'  The default is \code{group.nrow = NULL}.
#' @param ci.levels the confidence levels (maximum of 4) required to compute the
#'  confidence interval. The default is \code{ci.levels = c(0.90, 0.95, 0.99)}.
#' @param ci.method One of \code{c("hdr", "basic", "percentile")}. The default
#' is \code{ci.method = "hdr"} (see \code{\link[hdrcde]{hdr}}).
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL} for randomness. The default is \code{seed = object$seed}.
#' @param lambda Box-Cox transformation parameter. Value between 0 and 1. The
#' default is \code{lambda = 1} (see \code{\link[hdrcde]{hdr}}).
#' @param col Vector of colors for regions. The default is \code{col = NULL}.
#' @param plot.lines Boolean. The default is \code{plot.lines = TRUE}.
#' @param legend.position The default is \code{legend.position = "bottom"}.
#' @param ... additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @seealso \code{\link[GCEstim]{tsbootgce}}
#'
#' @examples
#' \donttest{
#' res.tsbootgce <-
#'   tsbootgce(
#'     formula = CO2 ~ 1 + L(GDP, 1) + L(EPC, 1) + L(EU, 1),
#'     data = moz_ts)
#'
#' plot(res.tsbootgce, which = 2, group = TRUE)
#' }
#'
#' @method plot tsbootgce
#' @importFrom rlang .data
#' @export

plot.tsbootgce <-
  function (x,
            which = c(1, 2),
            group = TRUE,
            group.ncol = NULL,
            group.nrow = NULL,
            ci.levels = c(0.90, 0.95, 0.99),
            ci.method = c("hdr", "basic", "percentile"),
            seed = object$seed,
            lambda = 1,
            col = NULL,
            plot.lines = TRUE,
            legend.position = "bottom",
            ...) {
    object <- x
    if (!inherits(object, "tsbootgce"))
      stop("use only with \"tsbootgce\" objects")
    if (length(ci.levels) > 4)
      stop("argument `ci.levels` must have length not greater than 4.", call. = FALSE)

    ci.method <- match.arg(ci.method)

    if (!is.null(seed))
      set.seed(seed)

    show <- rep(FALSE, 3)
    show[which] <- TRUE

    plots <- list(p1 = NULL,
                  p2 = NULL,
                  p3 = NULL)

    if (is.null(col)) {
      col <- c("#DF536B", "#61D04F", "#2297E6", "#28E2E5")
    }

    ci.levels <- sort(ci.levels)
    prob <- ci.levels * 100

    if (show[1L]) {

      for (i in 1:nrow(object$results$bootstrap$coef.matrix)) {

        p1 <- NULL
        LL <- UL <- Median <- NULL
        for (k in 1:length(ci.levels)) {
          LL <- c(LL, as.numeric(apply(object$meboot[[i]], 1, function(x){quantile(x, probs = (1 - ci.levels[k])/2)})))
          UL <- c(UL, as.numeric(apply(object$meboot[[i]], 1, function(x){quantile(x, probs = (1 + ci.levels[k])/2)})))
          Median <- c(Median, as.numeric(apply(object$meboot[[i]], 1, median)))
        }

        plot_data <-
          data.frame(Time = rep(as.numeric(time(object$data.ts[, 1])), length(ci.levels)),
                     y = rep(as.matrix(object$data.ts[, i]), length(ci.levels)),
                     LL,
                     UL,
                     Median,
                     ".width" = rep(ci.levels, each = length(time(object$data.ts[, 1]))))

        p1 <-
          ggplot2::ggplot(data = plot_data,
                          ggplot2::aes(x = .data$Time,
                                       y = .data$y,
                                       ymin = .data$LL,
                                       ymax = .data$UL)) +
          ggdist::geom_lineribbon() +
          ggplot2::scale_fill_manual(values = col[1:length(prob)],
                                     labels = paste0(prob,"%")) +
          ggplot2::geom_line(ggplot2::aes(y = Median),
                             linetype = "dashed",
                             colour = "yellow") +
          ggplot2::theme_minimal() +
          ggplot2::ylab(names(object$meboot)[i]) +
          ggplot2::theme(legend.title = ggplot2::element_blank())

        plots$p1[[i]] <- p1
      }
    }

    if (show[2L]) {

      for (i in 1:nrow(object$results$bootstrap$coef.matrix)) {
        x <- as.numeric(object$results$bootstrap$coef.matrix[i, ])
        h = hdrcde::hdrbw(hdrcde::BoxCox(x, lambda), mean(prob))
        if (lambda == 1) {
          den <- density(x, bw = h, n = 1001)
        } else {
          den.y <- hdrcde::BoxCox(x, lambda)
          den <- density(den.y, bw = h, n = 1001)
          den.j <- den$x > 0.1 - 1/lambda
          den$y <- den$y[j]
          den$x <- den$x[j]
          xgrid <- InvertBoxCox(den$x, lambda)
          den$y <- c(0, den$y * xgrid^(lambda - 1))
          den$x <- c(0, xgrid)
        }

        maxden <- max(den$y)
        stepy <- maxden * 0.02

        ylim <- c((1 - length(prob)) * stepy, maxden)
        ylim[1] <- min(ylim[1], -length(prob) * stepy)

        rangex <- range(x, na.rm = TRUE)
        minx <- rangex[1]
        maxx <- rangex[2]
        rangex <- maxx - minx

        if (ci.method == "hdr") {
          hd <- hdrcde::hdr(x = x, prob = prob, den = den, h = h)
        } else {
          hd <- list(hdr = matrix(NA,
                                  nrow = length(prob),
                                  ncol = 2,
                                  dimnames = list(paste0(rev(prob), "%"))),
                     mode = NA,
                     falpha =
                       {falpha <- rep(NA, length(prob))
                       names(falpha) <- paste0(rev(100-prob), "%")
                       falpha
                       },
                     falphamax =
                       {falphamax <- rep(NA, length(prob))
                       names(falphamax) <- paste0(rev(100-prob), "%")
                       falphamax
                       }
                     )

          for (ci.l in 1:length(prob)) {
            hd$hdr[ci.l,] <-
              confint(object,
                      level = rev(prob/100)[ci.l],
                      which = "estimates",
                      method = ci.method,
                      parm = i)
            hd$falpha[ci.l] <- den$y[which.min(abs(den$x - hd$hdr[ci.l, 1]))]
            hd$falphamax[ci.l] <- den$y[which.min(abs(den$x - hd$hdr[ci.l, 2]))]
          }
          hd$mode <- median(x)
          }

        nregions <- nrow(hd$hdr)
        leng <- dim(hd$hdr)[[2]]

        data.plot <- data.frame(x = den$x,
                                y = den$y)

        p2 <-
          ggplot2::ggplot(data = data.plot,
                          ggplot2::aes(x = .data$x, y = .data$y)) +
          ggplot2::geom_line() +
          ggplot2::xlab(paste0("N = ",
                               den$n,
                               "; Bandwidth = ",
                               round(den$bw, 5),
                               {if(ci.method == "hdr") ", Mode = " else", Median = "},
                               round(hd$mode, 5))) +
          ggplot2::ylab("Density") +
          ggplot2::ggtitle(
            latex2exp::TeX(
              paste0(
                "HDR of ",
                "$b_",
                i-ifelse("(Intercept)" %in% rownames(object$results$bootstrap$coef.matrix),
                         1,
                         0),
                "$",
                " estimates - ",
                rownames(object$results$bootstrap$coef.matrix)[i]))) +
          ggplot2::ylim(ylim) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.line = ggplot2::element_line(color='black'),
                         plot.background = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.border = ggplot2::element_blank())

        if (plot.lines) {
          for (r in 1:nregions) {
            if (ci.method == "hdr") {
              # p2 <- p2 +
              #   ggplot2::geom_hline(yintercept = hd$falpha[r],
              #                       col = col[r],
              #                       linetype = "dashed")
            for (j in 1:length(hd$hdr[r, ])) {
              data.segment <- data.frame(x = hd$hdr[r, j],
                                         xend = hd$hdr[r, j],
                                         y = stepy * (r - nregions),
                                         yend = hd$falpha[r])
              p2 <- p2 +
                ggplot2::geom_segment(data = data.segment,
                                      ggplot2::aes(x = .data$x,
                                                   xend = .data$xend,
                                                   y = .data$y,
                                                   yend = .data$yend),
                                      colour = col[r],
                                      linetype = "dashed")

            }
            } else {
              for (j in 1:length(hd$hdr[r, ])) {
                data.segment <- data.frame(x = hd$hdr[r, j],
                                           xend = hd$hdr[r, j],
                                           y = stepy * (r - nregions),
                                           yend = hd[[2 + j]][r])
                p2 <- p2 +
                  ggplot2::geom_segment(data = data.segment,
                                        ggplot2::aes(x = .data$x,
                                                     xend = .data$xend,
                                                     y = .data$y,
                                                     yend = .data$yend),
                                        colour = col[r],
                                        linetype = "dashed")

              }
              }
          }
        }

        for (r in 1:nregions) {

          hdr <- hd$hdr[r, ]
          nint <- length(hdr[!is.na(hdr)])/2
          if (nint != 0) {
            for (j in 1:nint) {
              l <- j * 2 - 1
              data.interval <- data.frame(x = hdr[l],
                                          xend = hdr[l + 1],
                                          y = (r - nregions - 0.5) * stepy,
                                          yend = (r - nregions - 0.5) * stepy)

              p2 <- p2 +
                ggplot2::geom_segment(data = data.interval,
                                      ggplot2::aes(x = .data$x,
                                                   xend = .data$xend,
                                                   y = .data$y,
                                                   yend = .data$yend),
                                      colour = col[r],
                                      linewidth = 2)
            }
          }
        }

        legend_df <- data.frame(
          x = rep(minx, length(prob)),
          y = rep(0, length(prob)),
          Level = factor(prob, levels = prob))

        p2 <- p2 +
          ggplot2::geom_vline(xintercept = hd$mode,
                              linetype = "dotted") +
          ggplot2::geom_point(data = legend_df,
                              ggplot2::aes(x = .data$x,
                                           y = .data$y,
                                           color = .data$Level),
                              size = 3, alpha = 0) +
          ggplot2::scale_color_manual(name = "",
                                      labels = paste0(prob,"%"),
                                      values = col[1:length(prob)]) +
          ggplot2::guides(
            color = ggplot2::guide_legend(
              override.aes = list(shape = 15,
                                  size = 5,
                                  alpha = 1)))

        if ((0 >= min(den$x)) && (0 <= max(den$x))) {

          p2 <- p2 +
            ggplot2::geom_segment(data = data.frame(x = 0,
                                                    xend = 0,
                                                    y = ylim[1],
                                                    yend = 0),
                                  ggplot2::aes(x = .data$x,
                                               xend = .data$xend,
                                               y = .data$y,
                                               yend = .data$yend),
                                  colour = "black",
                                  linewidth = 0.5)
        }

        plots$p2[[i]] <- p2
      }
    }

    if (show[3L]) {

      for (i in 1:nrow(object$results$bootstrap$nepk.matrix)) {
        x <- as.numeric(object$results$bootstrap$nepk.matrix[i, ])
        h = hdrcde::hdrbw(hdrcde::BoxCox(x, lambda), mean(prob))
        if (lambda == 1) {
          den <- density(x, bw = h, n = 1001)
        } else {
          den.y <- hdrcde::BoxCox(x, lambda)
          den <- density(den.y, bw = h, n = 1001)
          den.j <- den$x > 0.1 - 1/lambda
          den$y <- den$y[j]
          den$x <- den$x[j]
          xgrid <- InvertBoxCox(den$x, lambda)
          den$y <- c(0, den$y * xgrid^(lambda - 1))
          den$x <- c(0, xgrid)
        }
        maxden <- max(den$y)
        stepy <- maxden * 0.02

        ylim <- c((1 - length(prob)) * stepy, maxden)
        ylim[1] <- min(ylim[1], -length(prob) * stepy)

        rangex <- range(x, na.rm = TRUE)
        minx <- rangex[1]
        maxx <- rangex[2]
        rangex <- maxx - minx

        if (ci.method == "hdr") {
          hd <- hdrcde::hdr(x = x, prob = prob, den = den, h = h)
        } else {
          hd <- list(hdr = matrix(NA,
                                  nrow = length(prob),
                                  ncol = 2,
                                  dimnames = list(paste0(rev(prob), "%"))),
                     mode = NA,
                     falpha = {falpha <- rep(NA, length(prob))
                     names(falpha) <- paste0(rev(100-prob), "%")
                     falpha})

          for (ci.l in 1:length(prob)) {
            hd$hdr[ci.l,] <-
              confint(object,
                      level = rev(prob/100)[ci.l],
                      which = "NormEnt",
                      method = ci.method,
                      parm = i)
          }
          hd$mode <- median(x)
        }

        nregions <- nrow(hd$hdr)
        leng <- dim(hd$hdr)[[2]]

        data.plot <- data.frame(x = den$x,
                                y = den$y)

        p3 <-
          ggplot2::ggplot(data = data.plot,
                          ggplot2::aes(x = .data$x, y = .data$y)) +
          ggplot2::geom_line() +
          ggplot2::xlab(paste0("N = ",
                               den$n,
                               "; Bandwidth = ",
                               round(den$bw, 5),
                               {if(ci.method == "hdr") ", Mode = " else", Median = "},
                               round(hd$mode, 5))) +
          ggplot2::ylab("Density") +
          ggplot2::ggtitle(
            latex2exp::TeX(
              paste0(
                "HDR of ",
                "$b_",
                i-ifelse("(Intercept)" %in% rownames(object$results$bootstrap$coef.matrix),
                         1,
                         0),
                "$",
                " estimates - ",
                rownames(object$results$bootstrap$coef.matrix)[i]))) +
          ggplot2::ylim(ylim) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.line = ggplot2::element_line(color='black'),
                         plot.background = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.border = ggplot2::element_blank())

        if (plot.lines) {
          for (r in 1:nregions) {
            if (ci.method == "hdr") {
              # p3 <- p3 +
              #   ggplot2::geom_hline(yintercept = hd$falpha[r],
              #                       col = col[r],
              #                       linetype = "dashed")
              for (j in 1:length(hd$hdr[r, ])) {
                data.segment <- data.frame(x = hd$hdr[r, j],
                                           xend = hd$hdr[r, j],
                                           y = stepy * (r - nregions),
                                           yend = hd$falpha[r])
                p3 <- p3 +
                  ggplot2::geom_segment(data = data.segment,
                                        ggplot2::aes(x = .data$x,
                                                     xend = .data$xend,
                                                     y = .data$y,
                                                     yend = .data$yend),
                                        colour = col[r],
                                        linetype = "dashed")

              }
            } else {
              for (j in 1:length(hd$hdr[r, ])) {
                data.segment <- data.frame(x = hd$hdr[r, j],
                                           xend = hd$hdr[r, j],
                                           y = stepy * (r - nregions),
                                           yend = hd[[2 + j]][r])
                p3 <- p3 +
                  ggplot2::geom_segment(data = data.segment,
                                        ggplot2::aes(x = .data$x,
                                                     xend = .data$xend,
                                                     y = .data$y,
                                                     yend = .data$yend),
                                        colour = col[r],
                                        linetype = "dashed")

              }
            }
          }
        }

        for (r in 1:nregions) {

          hdr <- hd$hdr[r, ]
          nint <- length(hdr[!is.na(hdr)])/2
          if (nint != 0) {
            for (j in 1:nint) {
              l <- j * 2 - 1
              data.interval <- data.frame(x = hdr[l],
                                          xend = hdr[l + 1],
                                          y = (r - nregions - 0.5) * stepy,
                                          yend = (r - nregions - 0.5) * stepy)

              p3 <- p3 +
                ggplot2::geom_segment(data = data.interval,
                                      ggplot2::aes(x = .data$x,
                                                   xend = .data$xend,
                                                   y = .data$y,
                                                   yend = .data$yend),
                                      colour = col[r],
                                      linewidth = 2)
            }
          }
        }

        legend_df <- data.frame(
          x = rep(minx, length(prob)),
          y = rep(0, length(prob)),
          Level = factor(prob, levels = prob))

        p3 <- p3 +
          ggplot2::geom_vline(xintercept = hd$mode,
                              linetype = "dotted") +
          ggplot2::geom_point(data = legend_df,
                              ggplot2::aes(x = .data$x,
                                           y = .data$y,
                                           color = .data$Level),
                              size = 3, alpha = 0) +
          ggplot2::scale_color_manual(name = "",
                                      labels = paste0(prob,"%"),
                                      values = col[1:length(prob)]) +
          ggplot2::guides(
            color = ggplot2::guide_legend(
              override.aes = list(shape = 15,
                                  size = 5,
                                  alpha = 1)))

        if ((1 >= min(den$x)) && (1 <= max(den$x))) {

          p3 <- p3 +
            ggplot2::geom_segment(data = data.frame(x = 1,
                                                    xend = 1,
                                                    y = ylim[1],
                                                    yend = 0),
                                  ggplot2::aes(x = .data$x,
                                               xend = .data$xend,
                                               y = .data$y,
                                               yend = .data$yend),
                                  colour = "black",
                                  linewidth = 0.5)
        }

        plots$p3[[i]] <- p3
      }
    }

      if (group) {
        plotsg <- list(p1 = NULL,
                       p2 = NULL,
                       p3 = NULL)

        if (show[1L]) {
          plotsg$p1 <-
            ggpubr::ggarrange(plotlist = plots$p1,
                              common.legend = TRUE,
                              legend = legend.position,
                              ncol = group.ncol,
                              nrow = group.nrow)
        }
        if (show[2L]) {
          plotsg$p2 <-
            ggpubr::ggarrange(plotlist = plots$p2,
                              common.legend = TRUE,
                              legend = legend.position,
                              ncol = group.ncol,
                              nrow = group.nrow)
        }
        if (show[3L]) {
          plotsg$p3 <-
            ggpubr::ggarrange(plotlist = plots$p3,
                              common.legend = TRUE,
                              legend = legend.position,
                              ncol = group.ncol,
                              nrow = group.nrow)
        }

        plotsg[!sapply(plotsg, is.null)]
      } else {
        plots[!sapply(plots, is.null)]
      }
  }
