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

ObjFunGCE.primal.solnp <- function(x0, X, n, k, m, j, p0, w0, s1, S, weight) {

  p <- as.matrix(as.numeric(abs(matrix(x0[1:(k * m)], k, m))), ncol = 1)
  w <- as.matrix(as.numeric(abs(matrix(x0[(k * m + 1):length(x0)], n, j))), ncol = 1)

  ## CHECK #####
  # p <- round(p, 8)
  # w <- round(w, 8)
  # p[p == 0] <- 10^-8
  # w[w == 0] <- 10^-8
  #### ###

  return(as.numeric(
    (1 - weight) * (t(p) %*% log(p) - t(p) %*% log(p0)) +  weight * (t(w) %*% log(w) - t(w) %*% log(w0))
  ))
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

ConstFunGCE.primal.solnp <- function(x0, X, n, k, m, j, p0, w0, s1, S, weight) {
  aux.x0.X <- matrix(x0[1:(k * m)], k, m)
  aux.x0.error <- matrix(x0[(k * m + 1):length(x0)], n, j)
  aux.y = as.matrix(X) %*% (rowSums(s1 * aux.x0.X)) +
    rowSums(S * aux.x0.error)
  aux.p = c(rowSums(aux.x0.X), rowSums(aux.x0.error))
  return(c(aux.y, aux.p))
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
  p <- as.matrix(x0[1:(k * m)])
  w <- as.matrix(x0[(k * m + 1):ncol(x0)])

  ## CHECK #####
  # p <- round(p, 8)
  # w <- round(w, 8)
  # p[p == 0] <- 10^-8
  # w[w == 0] <- 10^-8
  #### ###

  return((1 - weight) * (t(p) %*% log(p) - t(p) %*% log(p0)) + weight * (t(w) %*% log(w) - t(w) %*% log(w0)))
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

ObjFunGCE.dual.optim <- function(x0, y, X, s1, s2, p0, w0, n, k, weight,...)
{
  # (7.3.9) page 112 Golan 1996
  Omega <- rep(0, k)
  for (k_aux in 1:k) {
    Omega[k_aux] <- sum(p0[k_aux, ] * exp(s1[k_aux, ] * sum(x0 * X[, k_aux]) * (1 / (1 - weight))))
  }

  # (7.3.10) page 112 Golan 1996 for gamma = 0.5.
  Psi <- rep(0, n)
  for (n_aux in 1:n) {
    Psi[n_aux] <- sum(w0[n_aux, ] * exp(x0[n_aux] * s2 * (1 / weight)))
  }

  # (7.3.11) page 112 Golan 1996
  return(-sum(x0 * y) + (1 - weight) * sum(log(Omega)) + weight * sum(log(Psi)))
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

GradFunGCE.dual.optim <- function(x0, y, X, s1, s2, p0, w0, n, k, m.optim, j, weight) {

  m <- m.optim
  p <- matrix(0, k, m)

  Omega <- rep(0, k)

  for (k_aux in 1:k) {
    temp <- sum(x0 * X[, k_aux])
    for (m_aux in 1:m) {
      p[k_aux, m_aux] <- p0[k_aux, m_aux] * exp(s1[k_aux, m_aux] * temp * (1 / (1 - weight)))
    }
    Omega[k_aux] <- sum(p[k_aux, ])
    p[k_aux, ] <- p[k_aux, ] / Omega[k_aux]
  }

  ## CHECK #####
  p <- round(p, 8)
  p[p == 0] <- 10^-8
  ### ###

  beta <- matrix(apply(s1 * p, 1, sum), ncol = 1)

  w <- matrix(0, n, j)
  Psi <- rep(0, n)

  for (n_aux in 1:n) {
    for (j_aux in 1:j) {
      w[n_aux, j_aux] <- w0[n_aux, j_aux] * exp(s2[j_aux] * x0[n_aux] * (1 / weight))
    }
    Psi[n_aux]  <- sum(w[n_aux, ])
    w[n_aux, ] <- w[n_aux, ] / Psi[n_aux]
  }

  ## CHECK #####
  w <- round(w, 8)
  w[w == 0] <- 10^-8
  ### ###

  epsilon <- w %*% matrix(s2, ncol = 1)

  # (7.3.17) page 113 Golan 1996 for gamma = 0.5.
  return(-y + X %*% beta + epsilon)
}
