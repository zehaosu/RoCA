delta.mat.x <- function(bmat, vas, var = NULL, x) {
  p <- nrow(bmat)
  mat <- solve(diag(1, p) - bmat)
  if (is.null(var)) {
    smat <- mat %*% t(mat)
  } else {
    smat <- mat %*% diag(var) %*% t(mat)
  }
  delta <- NULL
  for (v in vas) {
    w <- rep(0, p)
    w[x] <- 1
    if (length(v) != 0) {
      coef <- solve(smat[v, v]) %*% smat[v, x]
      w[v] <- -coef
    }
    delta <- cbind(delta, t(mat) %*% w)
  }
  return(delta)
}

delta.mat.y <- function(bmat, vas, var = NULL, x, y) {
  p <- nrow(bmat)
  mat <- solve(diag(1, p) - bmat)
  if (is.null(var)) {
    smat <- mat %*% t(mat)
  } else {
    smat <- mat %*% diag(var) %*% t(mat)
  }
  delta <- NULL
  for (v in vas) {
    w <- rep(0, p)
    w[y] <- 1
    coef <- solve(smat[c(x, v), c(x, v)]) %*% smat[c(x, v), y]
    w[c(x, v)] <- -coef
    delta <- cbind(delta, t(mat) %*% w)
  }
  return(delta)
}

get_true_covmat <- function(bmat, vas, var = NULL, x, y) {
  delta.x <- delta.mat.x(bmat, vas, var, x)
  delta.y <- delta.mat.y(bmat, vas, var, x, y)
  k <- length(vas)
  covmat <- matrix(NA, nrow = k, ncol = k)
  if (is.null(var)) var <- 1
  for (i in 1:k) {
    for (j in i:k) {
      mat <- cbind(sqrt(var) * delta.x[, i],
                   sqrt(var) * delta.x[, j],
                   sqrt(var) * delta.y[, i],
                   sqrt(var) * delta.y[, j])
      u1 <- (mat[, 1] * mat[, 2])
      u2 <- (mat[, 1] * mat[, 3])
      u3 <- (mat[, 2] * mat[, 3])
      v1 <- (mat[, 3] * mat[, 4])
      v2 <- (mat[, 2] * mat[, 4])
      v3 <- (mat[, 1] * mat[, 4])
      tmp <- sum(apply(mat, 1, prod)) +
             sum(outer(u1, v1)) - sum(u1 * v1) +
             sum(outer(u2, v2)) - sum(u2 * v2) +
             sum(outer(u3, v3)) - sum(u3 * v3)
      covmat[i, j] <- tmp
    }
  }
  covmat[lower.tri(covmat)] <- t(covmat)[lower.tri(t(covmat))]
  if (length(vas) == 1) {
    scale.mat <- t(1 / apply(sqrt(var) * delta.x, 2, function(x) sum(x * x)))
  } else {
    scale.mat <- diag(1 / apply(sqrt(var) * delta.x, 2, function(x) sum(x * x)))
  }
  covmat <- scale.mat %*% covmat %*% scale.mat
  return(covmat)
}
