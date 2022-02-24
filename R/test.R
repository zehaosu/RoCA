covmat_boot <- function(data, index, vas = NULL, ns = NULL, x = NULL, y = NULL){
  l <- get_est_covmat(data = data[index, ], vas = vas, ns = ns,
                      x = x, y = y, method = "plugin")
  matrixcalc::vech(l$covmat)
}

get_covmat_covmat <- function(data, B = 100, vas, ns = NULL, x = NULL, y = NULL,
                              method = c("bootstrap", "plugin")) {
  if (is.null(ns)) ns <- colnames(data)
  method <- match.arg(method)
  
  if (method == "bootstrap") {
    C <- nrow(data) * cov(boot::boot(data = data, statistic = covmat_boot, R = B,
                            vas = vas, ns = ns, x = x, y = y)$t)
  }
  if (method == "plugin") {
    fml <- get_fml(vas = vas, ns = ns, x = x, y = y)
    res <- lapply(fml, get_resid, data = data)
    res.y <- sapply(res, function(r) r[["res.y"]])
    res.x <- sapply(res, function(r) r[["res.x"]])
    p <- length(fml)
    idx <- cbind(unlist(lapply(1:p, function(x) x:p)),
                 unlist(lapply(1:p, function(x) rep(x, p - x + 1))))
    C <- matrix(NA, nrow = nrow(idx), ncol = nrow(idx))
    for (i in 1:nrow(idx)) {
      for (j in i:nrow(idx)) {
        x1 <- mean(res.x[, idx[i, 1]] ^ 2)
        x2 <- mean(res.x[, idx[i, 2]] ^ 2)
        x3 <- mean(res.x[, idx[j, 1]] ^ 2)
        x4 <- mean(res.x[, idx[j, 2]] ^ 2)
        x5 <- mean(res.x[, idx[i, 1]] * res.y[, idx[i, 1]] * res.x[, idx[i, 2]] *
                     res.y[, idx[i, 2]])
        x6 <- mean(res.x[, idx[j, 1]] * res.y[, idx[j, 1]] * res.x[, idx[j, 2]] *
                     res.y[, idx[j, 2]])
        x7 <- cov(res.x[, idx[i, 2]] ^ 2, res.x[, idx[j, 2]] ^ 2)
        x8 <- cov(res.x[, idx[i, 2]] ^ 2, res.x[, idx[j, 1]] ^ 2)
        x9 <- cov(res.x[, idx[i, 1]] ^ 2, res.x[, idx[j, 2]] ^ 2)
        x10 <- cov(res.x[, idx[i, 1]] ^ 2, res.x[, idx[j, 1]] ^ 2)
        x11 <- cov(res.x[, idx[i, 1]] * res.y[, idx[i, 1]] * res.x[, idx[i, 2]] *
                     res.y[, idx[i, 2]], res.x[, idx[j, 2]] ^ 2)
        x12 <- cov(res.x[, idx[i, 1]] * res.y[, idx[i, 1]] * res.x[, idx[i, 2]] *
                     res.y[, idx[i, 2]], res.x[, idx[j, 1]] ^ 2)
        x13 <- cov(res.x[, idx[j, 1]] * res.y[, idx[j, 1]] * res.x[, idx[j, 2]] *
                     res.y[, idx[j, 2]], res.x[, idx[j, 2]] ^ 2)
        x14 <- cov(res.x[, idx[j, 1]] * res.y[, idx[j, 1]] * res.x[, idx[j, 2]] *
                     res.y[, idx[j, 2]], res.x[, idx[i, 2]] ^ 2)
        x15 <- cov(res.x[, idx[i, 1]] * res.y[, idx[i, 1]] * res.x[, idx[i, 2]] *
                     res.y[, idx[i, 2]], res.x[, idx[j, 1]] * res.y[, idx[j, 1]] * 
                     res.x[, idx[j, 2]] * res.y[, idx[j, 2]])
        
        num <- x1 * x2 * x3 * x4 * x15 -
          x1 * x2 * x3 * x6 * x11 -
          x1 * x2 * x4 * x6 * x12 -
          x1 * x3 * x4 * x5 * x13 -
          x2 * x3 * x4 * x5 * x14 +
          x5 * x6 * x1 * x3 * x7 +
          x5 * x6 * x1 * x4 * x8 +
          x5 * x6 * x2 * x3 * x9 +
          x5 * x6 * x2 * x4 * x10
        den <- (x1 * x2 * x3 * x4) ^ 2
        C[i, j] <- C[j, i] <- num / den
      }
    }
  }
  return(C)
}

rank_test <- function(sigma, n, C) {
  l <- svd(sigma)
  p <- nrow(sigma)
  pval <- numeric(p - 1)
  # need special modification for r = p - 1
  for(r in 1:(p - 1)){
    S <- expm::sqrtm(l$u[(r + 1):p, (r + 1):p] %*% t(l$u[(r + 1):p, (r + 1):p]))
    S.inv <- solve(S)
    A <- l$u[, (r + 1):p] %*% solve(l$u[(r + 1):p, (r + 1):p]) %*% S
    if (r == (p - 1)) {
      lambda <- S.inv ^ 2 * l$d[(r + 1):p] *
        (l$u[(r + 1):p, (r + 1):p]) ^ 2
      E <- as.matrix(p - r)
    } else {
      lambda <- S.inv %*% l$u[(r + 1):p, (r + 1):p] %*% diag(l$d[(r + 1):p]) %*%
        t(l$u[(r + 1):p, (r + 1):p]) %*% S.inv
      E <- matrixcalc::elimination.matrix(p - r)
    }
    D <- matrixcalc::duplication.matrix(p)
    G <- kronecker(t(A), t(A))
    omega <- E %*% G %*% D %*% C %*% t(D) %*% t(G) %*% t(E)
    chisq <- n * t(matrixcalc::vech(lambda)) %*% solve(omega) %*% matrixcalc::vech(lambda)
    pval[r] <- pchisq(chisq, df = (p - r) * (p - r + 1) / 2, lower.tail = FALSE)
  }
  return(pval)
}
