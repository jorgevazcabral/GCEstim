# Functions for optimization

#' Objective function for primal optimization formulation using solnp
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param X .
#' @param n .
#' @param k .
#' @param m .
#' @param j .
#' @param p0 .
#' @param w0 .
#' @param s1 .
#' @param S .
#' @param weight .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ObjFunGCE.primal.solnp <- function(x0, X, n, k, m, j, p0, w0,
                                   s1, S, weight, ...) {
  p <- x0[1:(k*m)]
  w <- x0[(k*m + 1):length(x0)]

  return(
   (1 - weight) * sum(p * log(p / as.numeric(p0))) +
    weight * sum(w * log(w / as.numeric(w0)))
  )
}

#' Constraint function for primal optimization formulation using solnp
#'
#' Returns y value and probabilities
#'
#' @inheritParams ObjFunGCE.primal.solnp
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ConstFunGCE.primal.solnp <- function(x0, X, n, k, m, j,
                                     p0, w0, s1, S, weight) {
  p <- matrix(x0[1:(k * m)], nrow = k, ncol = m)
  w <- matrix(x0[(k * m + 1):length(x0)], nrow = n, ncol = j)

  signal_term <- X %*% rowSums(p * s1)
  error_term <- rowSums(w * S)

  y_constraint <- signal_term + error_term

  p_add <- rowSums(p)
  w_add <- rowSums(w)

  c(y_constraint, p_add, w_add)
}

#' Objective function for primal optimization formulation using solnl
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param m .
#' @param k .
#' @param p0 .
#' @param n .
#' @param w0 .
#' @param weight .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ObjFunGCE.primal.solnl <- function(x0, m, k, p0, n, w0, weight) {
  p <- x0[1:(k * m)]
  w <- x0[(k * m + 1):ncol(x0)]
  term_p <- sum(p * (log(p + 1e-12) - log(p0)))
  term_w <- sum(w * (log(w + 1e-12) - log(w0)))

  return((1 - weight) * term_p + weight * term_w)
}

#' LSE
#'
#' Returns LSE
#'
#' @param x .
#'
#' @return LSE
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

logsumexp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

#' Row LSE
#'
#' Returns row LSE
#'
#' @param x .
#'
#' @return row LSE
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

row_logsumexp <- function(mat) {
  apply(mat, 1, logsumexp)
}


#' Objective function for dual optimization formulation
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param y .
#' @param X .
#' @param s1 .
#' @param s2 .
#' @param p0 .
#' @param w0 .
#' @param n .
#' @param k .
#' @param weight .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ObjFunGCE.dual.optim <- function(x0, y, X, s1, s2, p0,
                                 w0, n, k, weight,...)
{
  temp_scalar <- t(X) %*% x0 / (1 - weight)
  temp_scalar_mat <- temp_scalar %*% matrix(1, ncol = ncol(s1))

  exponent_Omega <- s1 * temp_scalar_mat
  Omega_log <- row_logsumexp(log(p0) + exponent_Omega)

  exponent_Psi <- x0 %*% t(s2) / weight
  Psi_log <- row_logsumexp(log(w0) + exponent_Psi)

  -sum(x0 * y) + (1 - weight) * sum(Omega_log) + weight * sum(Psi_log)
}


#' Gradient function for dual optimization formulation
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param y .
#' @param X .
#' @param s1 .
#' @param s2 .
#' @param p0 .
#' @param w0 .
#' @param n .
#' @param k .
#' @param m.optim .
#' @param j .
#' @param weight .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

GradFunGCE.dual.optim <- function(x0, y, X, s1, s2, p0, w0,
                                  n, k, m.optim, j, weight) {

  temp_scalar <- t(X) %*% x0 / (1 - weight)
  temp_scalar_mat <- temp_scalar %*% matrix(1, ncol=ncol(s1))

  exponent_Omega <- s1 * temp_scalar_mat
  Omega_lse <- row_logsumexp(log(p0) + exponent_Omega)

  p <- exp(log(p0) + exponent_Omega - Omega_lse)

  beta <- matrix(rowSums(s1 * p), ncol = 1)

  exponent_Psi <- x0 %*% t(s2) / weight

  w <- exp(log(w0) + exponent_Psi - row_logsumexp(log(w0) +
        exponent_Psi) %*% t(rep(1, length(s2))))

  epsilon <- w %*% matrix(s2, ncol = 1)

  -y + X %*% beta + epsilon
}
