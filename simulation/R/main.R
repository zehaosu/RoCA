require(pcalg)
require(igraph)

get_fml <- function(vas, ns, x = NULL, y = NULL){
  id <- match(c("X", "Y"), ns)
  if (is.na(id[1]) && is.null(x)) stop("\"X\" is not specified")
  if (is.na(id[2]) && is.null(y)) stop("\"Y\" is not specified")
  
  if (!anyNA(id)) {
    fml <- lapply(vas, function(s) as.formula(paste0(c('Y ~ X', ns[s]),
                                                     collapse = ' + ')))
  }
  else if (is.na(id[1]) && !is.na(id[2])) {
    fml <- lapply(vas, function(s) as.formula(paste("Y ~", 
                                              paste0(c(ns[x], ns[s]),
                                                     collapse = ' + '))))
  }
  else if (is.na(id[2]) && !is.na(id[1])) {
    fml <- lapply(vas, function(s) as.formula(paste(ns[y], "~ X +",
                                              paste0(ns[s], collapse = ' + '))))
  }
  else {
    fml <- lapply(vas, function(s) as.formula(paste(ns[y], "~",
                                                    paste0(c(ns[x], ns[s]),
                                                           collapse = ' + '))))
  }
  
  fml
}

# util function for bootstrap
get_est_boot <- function(data, index, fml){
  d <- data[index, ]
  est <- sapply(fml, function(f) coef(lm(f, data = d))[[2]])
  est
}

get_resid <- function(data, fml){
  fit1 <- lm(fml, data = data)
  # check for empty adjustment set case
  if(grepl("\\+", as.character(fml)[3])){
    fit2 <- lm(as.formula(sub("\\+", "~", as.character(fml)[3])), data = data)
    return(list(
      tau = coef(fit1)[[2]],
      res.y = unname(resid(fit1)),
      res.x = unname(resid(fit2))
    ))
  }
  else{
    return(list(
      tau = coef(fit1)[[2]],
      res.y = unname(resid(fit1)),
      res.x = data[, as.character(fml)[3]]
    ))
  }
}

get_est_covmat <- function(data, B = 100, vas, ns = NULL, x = NULL, y = NULL,
                           method = c("bootstrap", "plugin")){
  
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(ns)) ns <- colnames(data)
  fml <- get_fml(vas = vas, ns = ns, x = x, y = y)
  
  method <- match.arg(method)
  
  if (method == "bootstrap") {
    est_boot <- boot::boot(data = data, statistic = get_est_boot, R = B, fml = fml)
    est <- est_boot$t0
    covmat <- cov(est_boot$t) * nrow(data)
  }
  else {
    res <- lapply(fml, get_resid, data = data)
    est <- sapply(res, function(r) r[["tau"]])
    res.y <- sapply(res, function(r) r[["res.y"]])
    res.x <- sapply(res, function(r) r[["res.x"]])
    covmat <- matrix(NA, nrow = length(fml), ncol = length(fml))
    for (i in seq_along(fml)) {
      for (j in i:length(fml)) {
        num <- mean(res.y[, i] * res.x[, i] * res.y[, j] * res.x[, j])
        den <- mean(res.x[, i] ^ 2) * mean(res.x[, j] ^ 2)
        covmat[i, j] <- covmat[j, i] <- num / den
      }
    }
    covmat <- covmat
  }
  
  list(est = est,
       covmat = covmat)
}

# import matrixcalc::vech
ic_recon <- function(sigma, k, n, aic = FALSE){
  dec <- svd(sigma)
  if (k == 1) rec.sigma <- dec$d[1] * dec$u[, 1]  %*% t(dec$u[, 1])
  else rec.sigma <- dec$u[, 1:k] %*% diag(dec$d[1:k]) %*% t(dec$u[, 1:k])
  p <- nrow(sigma)
  if(aic) {
    return(n * sum(matrixcalc::vech(sigma - rec.sigma) ^ 2) + 
         2 * (p * k - k * (k - 1) / 2))
  }
  else {
    return(n * sum(matrixcalc::vech(sigma - rec.sigma) ^ 2) + 
             log(n) * (p * k - k * (k - 1) / 2))
  }
}

ic_recon_svd <- function(sigma, k, n, rho){
  p <- nrow(sigma)
  if (k == p) return(log(n) * (p * k - k * (k - 1) / 2))
  l <- svd(sigma)
  S <- expm::sqrtm(l$u[(k + 1):p, (k + 1):p] %*% t(l$u[(k + 1):p, (k + 1):p]))
  S.inv <- solve(S)
  A <- l$u[, (k + 1):p] %*% solve(l$u[(k + 1):p, (k + 1):p]) %*% S
  if (k == (p - 1)) {
    lambda <- S.inv ^ 2 * l$d[(k + 1):p] *
      (l$u[(k + 1):p, (k + 1):p]) ^ 2
    E <- as.matrix(p - k)
  } else {
    lambda <- S.inv %*% l$u[(k + 1):p, (k + 1):p] %*% diag(l$d[(k + 1):p]) %*%
      t(l$u[(k + 1):p, (k + 1):p]) %*% S.inv
    E <- matrixcalc::elimination.matrix(p - k)
  }
  D <- matrixcalc::duplication.matrix(p)
  G <- kronecker(t(A), t(A))
  omega <- E %*% G %*% D %*% rho %*% t(D) %*% t(G) %*% t(E)
  chisq <- n * t(matrixcalc::vech(lambda)) %*% solve(omega) %*% matrixcalc::vech(lambda)
  return(drop(chisq + log(n) * (p * k - k * (k - 1) / 2)))
}

simulate_amat <- function(n, prob){
  amat <- matrix(0, nrow = n, ncol = n)
  for(i in seq_len(n - 1)) {
    amat[i + 1, 1:i] <- rbinom(i, 1, prob)
  }
  amat
}

## import pcalg::possDe pcalg::adjustment
set_xy <- function(amat){
  ## x <- sample(which(colSums(amat) > 0), 1)
  p <- ncol(amat)
  nd <- sapply(1:p, function(i)
    length(pcalg::possDe(m = amat, x = i, type = "dag", possible = FALSE))) - 1
  x <- sample(1:p, size = 1, prob = nd / sum(nd))
  if (!is.null(names(x))) names(x) <- NULL
  xd <- setdiff(pcalg::possDe(m = amat, x = x, type = "dag", possible = FALSE), x)
  if (length(xd) == 1) y <- xd
  else y <- sample(xd, 1)
  ## vas <- pcalg::adjustment(amat = amat, amat.type = "dag",
  ##                          x = x, y = y, set.type = set.type)
  list(x = x, y = y)
}

tau <- function(bmat, x, y) {
  g <- igraph::graph_from_adjacency_matrix(t(bmat), mode = "directed",
                                           weighted = TRUE)
  # topological order of nodes
  topo.ord <- as.vector(igraph::topo_sort(g))
  bmat.ord <- bmat[topo.ord, topo.ord]
  x.ord <- which(topo.ord == x)
  y.ord <- which(topo.ord == y)
  p <- ncol(bmat)
  tau <- matrix(0, p, 1)
  tau[x.ord] <- 1
  if (y.ord - x.ord > 1) {
    for (i in (x.ord + 1):y.ord) {
      tau[i] <- bmat.ord[i, ] %*% tau
    }
    return(tau[y.ord])
  }
  else {
    return(bmat[y, x])
  }
}

beta <- function(cmat, x, y, z) {
  z_ <- c(x, z)
  beta_vec <- solve(cmat[z_, z_]) %*% cmat[z_, y]
  return(beta_vec[1])
}

get_weight <- function(m, min = 0.1, max = 1) {
  w <- runif(m, min = min, max = max)
  sgn <- sample(c(-1, 1), size = m, replace = TRUE)
  w * sgn
}

simulate_graph <- function(n, d, method = "er", min = 1, max = 2) {
  counter <- 1
  g <- randDAG(n, d, method = method, weighted = TRUE,
               wFUN = list(get_weight, min = min, max = max))
  w <- t(as(g, "matrix"))
  dimnames(w) <- NULL
  a <- w
  a[a != 0] <- 1
  ww <- t(as(dag2cpdag(g), "matrix"))
  while(TRUE) {
  l <- set_xy(a)
  nonforb <- setdiff(1:nrow(ww), c(l$x, pcalg:::forbiddenNodes(ww, l$x, l$y, "cpdag")))
  vas <- pcalg::adjustment(amat = ww, amat.type = "cpdag", x = l$x, y = l$y, set.type = "minimal")
  vas <- union(unname(vas), list(nonforb))
  if (length(vas) >= 2) break
  if (counter >= 100) stop("Tried too many combinations... no output given")
  counter <- counter + 1
  }
  return(c(list(amat = a, bmat = w), l))
}

name_to_id <- function(name, ns){
  if(!is.list(name)) name <- list(name)
  lapply(name, function(k) match(k, ns))
}

id_to_name <- function(id, ns){
  if(!is.list(id)) id <- list(id)
  lapply(id, function(k) ns[k])
}
