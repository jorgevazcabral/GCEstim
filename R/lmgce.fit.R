#' Generalized Cross entropy estimation
#'
#' Internal function used to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @inheritParams lmgce.assign.noci
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
lmgce.fit <-
  function(y,
           X,
           offset,
           y.test = NULL,
           X.test = NULL,
           offset.test = NULL,
           errormeasure = "RMSE",
           min.coef = NULL,
           max.coef = NULL,
           max.abs.residual = NULL,
           support.signal = NULL,
           support.signal.points = c(1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5),
           support.noise = NULL,
           support.noise.points = c(1 / 3, 1 / 3, 1 / 3),
           weight = 0.5,
           method = "dual.lbfgsb3c",
           caseGLM = "D")
{
  X_model <- X
  y_model <- y

  var_beta <- NULL
  lambda_hat <- NULL

  if (caseGLM %in% c("M", "NM")) {
    y_model <- t(X_model) %*% y_model
    X_model <- t(X_model) %*% X_model
    if (caseGLM == "NM") {
      y_model <- y_model/nrow(X)
      X_model <- X_model/nrow(X)
    }
  }

  y_model <- as.matrix(y_model)
  n <- nrow(X_model)
  k <- ncol(X_model)

  if (length(support.signal.points) == 1) {
    m <- support.signal.points
    p0 <- matrix(1/support.signal.points, nrow = k, ncol = m, byrow = TRUE)
  } else if (is.vector(support.signal.points)) {
    m <- length(support.signal.points)
    p0 <- matrix(support.signal.points, nrow = k, ncol = m, byrow = TRUE)
  } else {
    m <- ncol(support.signal.points)
    p0 <- support.signal.points
  }

  if (length(support.noise.points) == 1) {
    j <- support.noise.points
    w0 <- matrix(1/support.noise.points, nrow = n, ncol = j, byrow = TRUE)
  } else if (is.vector(support.noise.points)) {
    j <- length(support.noise.points)
    w0 <- matrix(support.noise.points, nrow = n, ncol = j, byrow = TRUE)
  } else {
    j <- ncol(support.noise.points)
    w0 <- support.noise.points
  }

  if (length(support.signal) == 1 & is.null(max.coef)) {
    X_scaled <- X_model
    y_scaled <- y_model
    k_scaled <- k
    if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
      X_scaled <- X_scaled[, -1]
      k_scaled <- k_scaled - 1
      if (caseGLM %in% c("M", "NM")) {
        X_scaled <- X_scaled[-1, ]
        y_scaled <- y_scaled[-1, ]
      }
    }
    X_scaled <- scale(X_scaled)
    y_scaled <- scale(y_scaled)
  }

  if (length(support.signal) == 1) {
    if (is.null(max.coef)) {
    zlim <- matrix(c(-support.signal, support.signal),
                   k_scaled,
                   2,
                   byrow = TRUE)

    zlim <-
      as.matrix(data.frame(
        "LL" =
          {
            if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
              scalebackcoef(X_scaled,
                            y_scaled,
                            c(0, zlim[, 1]),
                            intercept = TRUE)
            } else {
              scalebackcoef(X_scaled,
                            y_scaled,
                            zlim[, 1],
                            intercept = FALSE)
            }
          },
        "UL" =
          {
            if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
              scalebackcoef(X_scaled,
                            y_scaled,
                            c(0, zlim[, 2]),
                            intercept = TRUE)
            } else {
              scalebackcoef(X_scaled,
                            y_scaled,
                            zlim[, 2],
                            intercept = FALSE)
            }
          }
      ))

    if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
      zlim[1 , ] <-
        c(-max(abs(zlim[1, ])), max(abs(zlim[1, ])))
    }
  } else {
    if (is.null(min.coef)) {
      zlim <- support.signal * as.matrix(cbind(-max.coef, max.coef))
    } else {
      zlim <-
        as.matrix(cbind(
          (max.coef + min.coef) / 2 + support.signal * (max.coef - min.coef) * (-0.5),
          (max.coef + min.coef) / 2 + support.signal * (max.coef - min.coef) * (0.5)))
      }
  }
    } else if (length(support.signal) == 2) {
    zlim <- matrix(sort(support.signal), k, 2, byrow = TRUE)
  } else {
    zlim <- t(apply(support.signal, 1, sort))
  }

  colnames(zlim) <- c("SupportLL","SupportUL")
  row.names(zlim) <- colnames(X)

  if (is.null(support.noise)) {
    if (is.null(max.abs.residual)) {
      vlim <- c(-3 * sd(y_model), 3 * sd(y_model))
    } else {
      vlim <- c(-max.abs.residual, max.abs.residual)
    }
      } else {
    vlim <- support.noise
  }

  s1 <- matrix(0, k, m)
  Z <- matrix(0, k, k * m)
  for (i in 1:k) {
    s1[i,] <- seq(zlim[i, 1],
                  zlim[i, 2],
                  by = (zlim[i, 2] - zlim[i, 1]) / (m - 1))
    Z[i, ((i - 1) * m + 1):(i * m)] <- s1[i, 1:m]
  }

  s2 <- seq(vlim[1], vlim[2], by = (vlim[2] - vlim[1]) / (j - 1))
  S <- matrix(rep(s2, n), ncol = length(s2), byrow = TRUE)
  V <- matrix(0, n, n * j)
  for (i in 1:n) {
    V[i, ((i - 1) * j + 1):(i * j)] <- s2
  }

  #method.maxfeval = 1e+04
  method.maxiter = 5000
  method.tol = 1e-06

  if (method == "dual.optimParallel") {
    cl <- parallel::makeCluster(2)
    parallel::setDefaultCluster(cl = cl)

    res.opt <-
    optimParallel::optimParallel(
      par = rep(1e-8, n),
      fn = ObjFunGCE.dual.optim,
      gr = GradFunGCE.dual.optim,
      y = y_model,
      X = X_model,
      s1 = s1,
      s2 = s2,
      p0 = p0,
      w0 = w0,
      n = n,
      k = k,
      m.optim = m,
      j = j,
      weight = weight,
      hessian = TRUE
    )

    parallel::setDefaultCluster(cl = NULL)
    parallel::stopCluster(cl)
  }
  else if (method %in% c("dual.lbfgsb3c")) {

    res.opt <-
      lbfgsb3c::lbfgsb3c(
        par = rep(1e-8, n),
        fn = ObjFunGCE.dual.optim,
        gr = GradFunGCE.dual.optim,
        y = y_model,
        X = X_model,
        s1 = s1,
        s2 = s2,
        p0 = p0,
        w0 = w0,
        n = n,
        k = k,
        m.optim = m,
        j = j,
        weight = weight
      )

  } else

  if (method %in% c("dual.lbfgs")) {

    res.opt <- suppressWarnings(
      lbfgs::lbfgs(
        ObjFunGCE.dual.optim,
        GradFunGCE.dual.optim,
        vars = rep(1e-8, n),
        y = y_model,
        X = X_model,
        s1 = s1,
        s2 = s2,
        p0 = p0,
        w0 = w0,
        n = n,
        k = k,
        m.optim = m,
        j = j,
        weight = weight,
        invisible = 1
      )
    )

  } else

  if (method %in% c("dual.Rcgmin", "dual.bobyqa", "dual.newuoa",
                    "dual.nlminb", "dual.nlm")) {

    res.opt <- suppressWarnings(
      optimx::optimx(
        par = rep(1e-8, n),
        ObjFunGCE.dual.optim,
        GradFunGCE.dual.optim,
        y = y_model,
        X = X_model,
        s1 = s1,
        s2 = s2,
        p0 = p0,
        w0 = w0,
        n = n,
        k = k,
        m.optim = m,
        j = j,
        weight = weight,
        method = gsub("dual.", "", method),
        #control = list(maxit = 200, trace=0),
        hessian = TRUE
      )
    )

    res.opt <-
      list(par = as.numeric(res.opt[1:n]),
           convergence = res.opt$convcode)
  } else

  if (method %in% c("dual.BFGS", "dual.CG", "dual.L-BFGS-B")) {

    res.opt <- optim(
      par = rep(1e-8, n),
      ObjFunGCE.dual.optim,
      GradFunGCE.dual.optim,
      y = y_model,
      X = X_model,
      s1 = s1,
      s2 = s2,
      p0 = p0,
      w0 = w0,
      n = n,
      k = k,
      m.optim = m,
      j = j,
      weight = weight,
      method = gsub("dual.", "", method),
      #control = list(maxit = 200, trace=0),
      hessian = TRUE)

  } else

  if (method == "primal.solnp") {
    res.opt = Rsolnp::solnp(
      pars = c(rep(1 / m, k * m), rep(1 / j, n * j)),
      fun = ObjFunGCE.primal.solnp,
      eqfun = ConstFunGCE.primal.solnp,
      eqB = c(y_model, rep(1, n + k)),
      LB = rep(1e-5, k * m + n * j),
      UB = rep(1, k * m + n * j),
      control =
        list(
        #tol = method.tol,
        #inner.iter = method.maxiter,
        trace = 0
      ),
      X = X_model,
      n = n,
      k = k,
      m = m,
      j = j,
      p0 = as.numeric(p0),
      w0 = as.numeric(w0),
      s1 = s1,
      S = S,
      weight = weight
    )

    aux.convergence <- res.opt$convergence

    p <- matrix(res.opt$pars[1:(k*m)],
                nrow = k,
                ncol = m)
    w <- matrix(res.opt$pars[(k*m + 1):length(res.opt$pars)],
                nrow = n,
                ncol = j)

    beta_hat <- Z %*% as.numeric(t(p))

    #e_hat <- V %*% as.numeric(t(w))

  } else if (method == "primal.solnl") {

    Aeq <-
      rbind(cbind(as.matrix(X_model) %*% Z, V),
            cbind(kronecker(diag(k),
                            matrix(1, 1, m)),
                  matrix(0, k, n * j)),
            cbind(matrix(0, n, k * m),
                  kronecker(diag(n),
                            matrix(1, 1, j))))

    res.opt <- pracma::fmincon(
      x0 = t(c(rep(1 / m, k * m), rep(1 / j, n * j))),
      fn = ObjFunGCE.primal.solnl,
      Aeq = Aeq,
      beq = c(y_model, rep(1, n + k)),
      lb = rep(1e-5, k * m + n * j),
      ub = rep(1, k * m + n * j),
      m = m,
      k = k,
      p0 = as.numeric(t(p0)),
      n = n,
      w0 = as.numeric(t(w0)),
      weight = weight
      #,tol = method.tol,
      #maxfeval = method.maxfeval,
      #maxiter = method.maxiter
    )

    aux.convergence <- res.opt$convergence

    p <- res.opt$par[1:(k * m)]
    beta_hat <- Z %*% p
    p <- matrix(p, ncol = m, nrow = k, byrow = TRUE)

    w <- res.opt$par[(k * m + 1):length(res.opt$par)]
    #e_hat <- V %*% w
    w <- matrix(w, ncol = j, nrow = n, byrow = TRUE)

   } else if (method == "dual") {
    dimZ <- ncol(Z)
    t <- 1
    u <- 1
    lambda <- increm <- matrix(0, n, 1)
    iter <- 0
    while (u > method.tol & method.maxiter > iter) {
      iter <- iter + 1
      lambda <- lambda + increm
      newz <- exp(-t(X_model) %*% lambda %*% matrix(1, 1, dimZ) * Z)
      p_m2 <- newz / (newz %*% matrix(1, dimZ, dimZ))
      newv <- exp(-lambda %*% matrix(1, 1, j) * S)
      w9 <- newv / (newv %*% matrix(1, j, j))
      g9 <-
        y_model - as.matrix(X_model) %*% ((Z * p_m2) %*% matrix(1, dimZ, 1)) -
        ((S * w9) %*% matrix(1, j, 1))
      inv_z <-
        diag((apply((p_m2 * (
          Z ^ 2
        )), 1, sum) -
          apply((p_m2 * Z), 1, sum) ^ 2) ^ (-1))
      inv_v <-
        diag((apply((w9 * (
          S ^ 2
        )), 1, sum) -
          apply((w9 * S), 1, sum) ^ 2) ^ (-1))
      temp <- inv_v %*% as.matrix(X_model)
      inv_H = -inv_v +
        ((temp %*% solve(inv_z + t(as.matrix(
          X_model
        )) %*% temp)) %*% t(temp))
      increm <- inv_H %*% g9
      t0 <- t
      t <- t(g9) %*% increm
      u <- abs(t - t0)
    }

    aux.convergence <- ifelse(iter > method.maxiter , 1, 0)

    beta_hat <- apply(p_m2 * Z, 1, sum)

    p <- matrix(0, 1, k * m)

    for (i in 1:k) {
      pos <- (i - 1) * m + 1
      p[1,pos:(pos + m - 1)] <-
        p_m2[i, pos:(pos + m - 1)] + (sum(p_m2[i,-c(pos:(pos + m - 1))])/m)
    }
    p <- t(p)[,1]

    p <- matrix(p, ncol = m, nrow = k, byrow = TRUE)
    }

  if (method %in% c("dual.optimParallel",
                    "dual.BFGS", "dual.CG", "dual.L-BFGS-B",
                    "dual.Rcgmin", "dual.bobyqa", "dual.newuoa",
                    "dual.nlminb", "dual.nlm",
                    "dual.lbfgs", "dual.lbfgsb3c")) {

    aux.convergence <- res.opt$convergence
    lambda_hat <- res.opt$par

    p <- matrix(0, k, m)
    Omega <- rep(0, k)

    for (k_aux in 1:k) {
      temp <- sum(lambda_hat * X_model[, k_aux])

      for (m_aux in 1:m) {
        p[k_aux, m_aux] <- p0[k_aux, m_aux] * exp(s1[k_aux, m_aux] * temp * (1 / (1 - weight)))
      }

      Omega[k_aux] <- sum(p[k_aux, ])
      p[k_aux, ] <- p[k_aux, ] / Omega[k_aux]
    }

    p <- round(p, 8)
    p[p == 0] <- 10^-8

    beta_hat <- matrix(apply(s1 * p, 1, sum), ncol = 1)

    Psi <- rep(0, n)

    w <- matrix(0, n, j)
    for (n_aux in 1:n) {
      for (j_aux in 1:j) {
        w[n_aux, j_aux] <-
          w0[n_aux, j_aux] * exp(s2[j_aux] * lambda_hat[n_aux] * (1 / weight))
      }
      Psi[n_aux]  <- sum(w[n_aux, ])
      w[n_aux, ] <- w[n_aux, ] / Psi[n_aux]
    }

    w <- round(w, 8)
    w[w == 0] <- 10^-8

    sigma2_zeta <- rep(0, n)
    for (n_aux in 1:n) {
      sigma2_zeta[n_aux] <-
        sum((s2 * s2) * w[n_aux, ]) - (sum(s2 * w[n_aux, ]))^2
    }
    var_beta <-
      ((sum(lambda_hat * lambda_hat) / n) / ((sum(1 / sigma2_zeta) / n)^2)) * solve(t(X_model) %*% X_model)
  }

  if (method %in% c("primal.solnp")){

    lambda_hat <- res.opt$lagrange[1:n]
    p <- round(p, 8)
    p[p == 0] <- 10^-8
    w <- round(w, 8)
    w[w == 0] <- 10^-8

    sigma2_zeta <- rep(0, n)
    for (n_aux in 1:n) {
      sigma2_zeta[n_aux] <-
        sum((s2 * s2) * w[n_aux, ]) - (sum(s2 * w[n_aux, ]))^2
    }
    var_beta <-
      ((sum(lambda_hat * lambda_hat) / n) / ((sum(1 / sigma2_zeta) / n)^2)) * solve(t(X_model) %*% X_model)

  }

  nep <- sum(p * log(p)) / (- k * log(m))

  nepk <- matrix(0, k, 1)

  for (k_aux in 1:k) {

    nepk[k_aux, 1] <-
        sum(p[k_aux, ] * log(p[k_aux, ])) / (- log(m))
    }

  if (is.null(y.test) || is.null(X.test)) {
    y.fitted <- as.matrix(X) %*% beta_hat
    y.values <- as.matrix(y)
  } else {
    y.fitted <- as.matrix(X.test) %*% beta_hat
    y.values <- as.matrix(y.test)
  }

  aux_resid <- as.numeric(y.values - y.fitted)
  names(aux_resid) <- attr(y.values,"dimnames")[[1]]

  y.fitted <- as.numeric(y.fitted)
  names(y.fitted) <- attr(y.values,"dimnames")[[1]]

  if (is.null(y.test) ||
      is.null(X.test)) {
    if (!is.null(offset)) {
      y.fitted <- y.fitted + offset
    }
  } else {
    if (!is.null(offset.test)) {
      y.fitted <- y.fitted + offset.test
    }
  }

  beta_hat <- as.numeric(beta_hat)
  names(beta_hat) <- colnames(X)

  nepk <- as.numeric(nepk)
  names(nepk) <- colnames(X)

  row.names(p) <- names(nepk)
  colnames(p) <- paste0("p_", 1:m)

  row.names(w) <- row.names(y_model)
  colnames(w) <- paste0("p_", 1:j)

  res <- list(
    coefficients = beta_hat,
    residuals = aux_resid,
    fitted.values = y.fitted,
    nep = nep,
    nepk = nepk,
    vcov = var_beta,
    error.measure = accmeasure(y.fitted, y.values, which = errormeasure),
    support.stdUL ={if (length(support.signal) == 1) support.signal else NULL},
    support.matrix = zlim,
    p = p,
    w= w,
    lambda = lambda_hat,
    convergence = aux.convergence,
    p0 = p0,
    w0 = w0
    )

  class(res) <- "lmgce.fit"

  return(res)

}
