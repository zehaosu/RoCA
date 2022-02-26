getFml <- function(x, y, Z, V.name){
  fml <- lapply(Z, function(s) as.formula(paste(V.name[y], "~",
                                                    paste0(c(V.name[x], V.name[s]),
                                                           collapse = ' + '))))
  fml
}

# lm fitting and retrieving residuals
getResid <- function(data, fml){
  fit1 <- lm(fml, data = data)
  # check for empty adjustment set case
  if(grepl("\\+", as.character(fml)[3])){
    fit2 <- lm(as.formula(sub("\\+", "~", as.character(fml)[3])), data = data)
    return(list(
      tau = coef(fit1)[[2]],
      res.y = unname(resid(fit1)),
      res.x = unname(resid(fit2))
    ))
  } else {
    return(list(
      tau = coef(fit1)[[2]],
      res.y = unname(resid(fit1)),
      res.x = data[, as.character(fml)[3]]
    ))
  }
}

# computing TCE estimate and asymptotic covmat estimate
getEstCov <- function(data, x, y, Z, V.name = NULL){
  
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(V.name)) V.name <- colnames(data)
    fml <- getFml(Z = Z, V.name = V.name, x = x, y = y)
    
    res <- lapply(fml, getResid, data = data)
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
  
  list(est = est,
       covmat = covmat)
}

# getting rank estimation via BIC
ICRankEst <- function(sigma, n){
  p <- nrow(sigma)
  dec <- svd(sigma)
    ic <- sapply(1:p, function(k) {
      if (k == 1) rec.sigma <- dec$d[1] * dec$u[, 1]  %*% t(dec$u[, 1])
      else rec.sigma <- dec$u[, 1:k] %*% diag(dec$d[1:k]) %*% t(dec$u[, 1:k])
      sigma.diff <- sigma - rec.sigma
      return(n * crossprod(sigma.diff[lower.tri(sigma.diff, diag = TRUE)]) + 
               log(n) * (p * k - k * (k - 1) / 2))
    })
  return(which.min(ic))
}

# getting adjustment sets with a unique node each
getUniqueZ <- function(vas) {
  ns <- Reduce(union, vas)
  vas_unique <- list()
  for (i in seq_along(vas)) {
    ns_other <- Reduce(union, vas[-i])
    if (!all(vas[[i]] %in% ns_other)) {
      vas_unique <- c(vas_unique, vas[i])
    }
  }
  vas_rem <- setdiff(vas, vas_unique)
  ns_unique <- Reduce(union, vas_unique)
  while (length(ns_unique) < length(ns)) {
    vas_rm <- numeric(0L)
    for (i in seq_along(vas_rem)) {
      if (all(vas_rem[[i]] %in% ns_unique)) vas_rm <- c(vas_rm, i)
    }
    if (length(vas_rm) > 0) vas_rem <- vas_rem[-vas_rm]
    if (length(vas_rem) == 0) break
    v <- vas_rem[[1]]
    vas_rem <- vas_rem[-1]
    vas_rm <- numeric(0L)
    for (i in seq_along(vas_unique)) {
      if (length(vas_rm) > 0) {
        ns_other <- Reduce(union, c(vas_unique[-c(vas_rm, i)], list(v)))
      } else {
        ns_other <- Reduce(union, c(vas_unique[-i], list(v)))
      }
      if (all(vas_unique[[i]] %in% ns_other)) vas_rm <- c(vas_rm, i)
    }
    if (length(vas_rm) > 0) {
      vas_unique <- c(vas_unique[-vas_rm], list(v))
    } else {
      vas_unique <- c(vas_unique, list(v))
    }
    ns_unique <- Reduce(union, vas_unique)
  }
  return(vas_unique)
}

# main function for performing the test
waldTest <- function(data, x, y, Z, rank.est, contrast = NULL) {
  l <- getEstCov(data = data, x = x, y = y, Z = Z)
  n <- nrow(data)
  sigma <- l$covmat
  
  # no. of adjustment sets used
  k <- nrow(sigma)
  
  if (!is.null(contrast)) {
    if (any(dim(contrast) != c(k - 1, k))) {
      contrast <- NULL
      message(paste0("Wrong dimension of contrast matrix: expected (",
                     k - 1, ", ", k, ") and got (", dim(contrast)[1], ", ",
                     dim(contrast)[2], ")\nUsing default contrast matrix"))
    }
    if (!(all.equal(rowSums(contrast), rep(0, k - 1)) == TRUE)) {
      contrast <- NULL
      message("Input matrix is not a valid contrast matrix\nUsing default contrast matrix")
    }
  }
  
  if (is.null(contrast)) {
    # provide default contrast matrix
    C <- sapply(1:(k - 1),
                function(j) c(rep(0, j - 1), c(1, -1), rep(0, k - 1 - j)))
    C <- t(C)
  } else {
    C <- contrast
  }
  
  v <- C %*% l$est
  S <- C %*% sigma %*% t(C)
  
  if (rank.est) {
      r <- ICRankEst(sigma = sigma, n = n)
      if (r == 1) stop("Rank of Sigma estimated to be one; test not available!")
      r <- ICRankEst(sigma = S, n = n)
  } else {
      r <- k - 1
  }
  
  tmp <- svd(S)
  if (r > 1) {
    S.inv <- tmp$u[, 1:r] %*% diag(1/tmp$d[1:r]) %*% t(tmp$u[, 1:r])
  } else {
    S.inv <- 1/tmp$d[1] * tmp$u[, 1] %*% t(tmp$u[, 1])
  }
  csq <- n * t(v) %*% S.inv %*% v
  
  pval <- drop(pchisq(csq, df = r, lower.tail = FALSE))
  
  c(pval = pval, statistic = drop(csq), rank = as.integer(r))
}
