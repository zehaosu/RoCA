library(CVXR)

vech <- function(x, lower = TRUE) {
  if (all.equal(x, t(x)) != TRUE) stop("Not symmetric")
  if (lower) return(t(t(x[!upper.tri(x)])))
  else return(t(t(x[!lower.tri(x)])))
}

invvech <- function(x, lower = TRUE) {
  k <- (sqrt(1 + 8 * length(x)) - 1) / 2
  m <- matrix(0, nrow = k, ncol = k)
  if (lower) m[!upper.tri(m)] <- x
  else m[!lower.tri(m)] <- x
  d <- diag(m)
  return(m + t(m) - diag(d))
}

D.matrix <- function(k, lower = TRUE) {
  if (lower) return(matrixcalc::D.matrix(k))
  else {
    o <- matrix(0, nrow = k, ncol = k)
    o[!upper.tri(o)] <- 1:(k * (k + 1) / 2)
    o <- t(o)[!lower.tri(t(o))]
    return(matrixcalc::D.matrix(k)[, o])
  }
}

E.matrix <- function(k, lower = TRUE) {
  if (lower) return(matrixcalc::L.matrix(k))
  else {
    o <- matrix(0, nrow = k, ncol = k)
    o[!upper.tri(o)] <- 1:(k * (k + 1) / 2)
    o <- t(o)[!lower.tri(t(o))]
    return(matrixcalc::L.matrix(k)[o, ])
  }
}


irm <- function(sigma, C, r,
                w = 1, q = 1.1, t.max = 20, eps.1 = 1e-6, eps.2 = 1e-6,
                verbose = TRUE) {
  k <- nrow(sigma)
  stopifnot("Rank greater than matrix dimension" = (r <= k))
  
  C.inv <- solve(C)
  Q <- function(x){
    t(vech(sigma - x)) %*% C.inv %*% vech(sigma - x)
  }
  t <- 0
  # unconstrained problem
  S <- Variable(k, k, PSD = TRUE)
  I <- diag(1, nrow = k - r, ncol = k - r)
  obj <- quad_form(E.matrix(k) %*% (vec(sigma - S)), C.inv)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  S.best <- S.new <- result$getValue(S)
  lst <- eigen(S.new)
  u.old <- lst$values[r + 1]
  V <- lst$vectors[, (r + 1):k]
  Q.best <- result$value
  w <- w * q

  while (TRUE) {
    S <- Variable(k, k, PSD = TRUE)
    u <- Variable(1, nonneg = TRUE)
    A <- Variable(k - r, k - r, PSD = TRUE)
    I <- diag(1, nrow = k - r, ncol = k - r)
    obj <- quad_form(E.matrix(k) %*% (vec(sigma - S)), C.inv) + u * w
    constr <- list(u <= u.old,
                   A == u * I - t(V) %*% S %*% V)
    prob <- Problem(Minimize(obj), constr)
    result <- solve(prob)
    
    # update values
    S.old <- S.new
    S.new <- result$getValue(S)
    if (anyNA(S.new)) break
    V <- eigen(S.new)$vectors[, (r + 1):k]
    u.old <- result$getValue(u)
    w <- w * q
    t <- t + 1
    
    if (verbose) {
      cat(paste0("Iteration: ", t, "\n",
                 "Target function: ", round(Q(S.new), 6), "\n",
                 "Eigenvalue: ", u.old, "\n"))
    }
    if(t >= t.max | abs(u.old) < eps.1 | 
       abs((Q(S.new) - Q(S.old)) / Q(S.old)) < eps.2) {
      break
    }
  }
  if (anyNA(S.new)) return(Q(S.old))
  else return(Q(S.new))
  # return(S.best)
}
