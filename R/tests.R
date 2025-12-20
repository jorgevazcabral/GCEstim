#' Entropy Ratio test
#'
#' The Entropy Ratio test - which corresponds to the likelihood ratio, or
#' empirical ratio, test - measures the entropy discrepancy between the
#' constrained and the unconstrained models.
#'
#' @param object fitted \code{\link{lmgce}} object.
#'
#' @return A matrix with the X-squared statistics, degrees of freedom and
#' p-value for each parameter.
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
#' ER.test(res_gce_package)
#'
#' @export

ER.test <- function(object)
{

  p <- object$p

  # pnames <- names(object$nepk) #row.names(p)
  # if (missing(parm)) parm <- pnames
  # else if (is.numeric(parm)) parm <- pnames[parm]
  # if (!all(parm %in% pnames)) {
  #   stop("Invalid choice of parameter")
  # }


  digits = max(3L, getOption("digits") - 2L)
  dig.tst = max(1L, min(5L, digits - 1L))
  eps.Pvalue = .Machine$double.eps

  k <- nrow(p)

  ER <- rep(0, k)
  for (k_aux in 1:k) {
    ER[k_aux] <- 2 * (log(ncol(p)) + sum(p[k_aux, ] * log(p[k_aux, ])))
  }

  res <- data.frame(
    cbind(ER,
          1,
          format.pval(pchisq(ER, 1, lower.tail = FALSE),
                      digits = dig.tst,
                      eps = eps.Pvalue),
          format(symnum(pchisq(ER,
                               1,
                               lower.tail = FALSE),
                        corr = FALSE, na = FALSE,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", " "))))
  )

  row.names(res) <- row.names(p)
  colnames(res) <- c("X-squared", "df","p-value")

  return(res)
}
