wald_test <- function(data, B = 100, vas, ns = NULL,
                      method = c("bootstrap", "plugin"),
                      amat = NULL, x = NULL, y = NULL,
                      rank.method = c("bic", "aic", "svd", "none"),
                      level = 0.05,
                      rank = FALSE, statistic = FALSE) {
    l <- get_est_covmat(data = data, B = B, vas = vas, ns = ns,
                        method = method, x = x, y = y)
    n <- nrow(data)
  
    sigma <- l$covmat
    # no. of valid adjustment sets used
    p <- nrow(sigma)
    if(p <= 1) stop(paste('Test not available for', p, 'adjustment set'))
    
    # contrast matrix
    C <- sapply(1:(p - 1),
                function(j) c(rep(0, j - 1), c(1, -1), rep(0, p - 1 - j)))
    C <- t(C)
    v <- C %*% l$est
    S <- C %*% sigma %*% t(C)
    
    rank.method <- match.arg(rank.method)
    # what if p = 2?
    if (rank.method == "bic") {
        r <- which.min(sapply(1:(p - 1),
                              function(k) ic_recon(S, k = k, n = n)))
    }
    if (rank.method == "aic") {
        r <- which.min(sapply(1:(p - 1),
                              function(k) ic_recon(S, k = k, n = n,
                                                   aic = TRUE)))
    }
    if (rank.method == "svd") {
        G <- get_covmat_covmat(data = data, vas = vas, ns = ns,
                               method = "plugin")
        P <- matrixcalc::elimination.matrix(p - 1) %*%
            kronecker(G, G) %*% matrixcalc::duplication.matrix(p)
        G <- P %*% G %*% t(P)
        r <- which.min(sapply(1:(p-1),
                              function(k) ic_recon_svd(S, k = k, n = n,
                                                       rho = C)))
    }
    if (rank.method == "none") {
        r <- p - 1
    }

    # need to screen for rank of sigma as well!
    if (rank.method == "bic") {
        r.max <- which.min(sapply(1:p,
                                  function(k) ic_recon(sigma, k = k,
                                                       n = n))) - 1
        if (r.max == 0) {
            message("Rank is zero; no comparison available")
            return(NA)
        }
    }
    
    tmp <- svd(S)
    
    if (r > 1) {
        S.inv <- tmp$u[, 1:r] %*% diag(1/tmp$d[1:r]) %*%
            t(tmp$u[, 1:r])
    } else if (r == 1) {
        S.inv <- 1/tmp$d[1] * tmp$u[, 1] %*% t(tmp$u[, 1])
    }
  
    csq <- n * t(v) %*% S.inv %*% v
  
    pval <- drop(pchisq(csq, df = r, lower.tail = FALSE))
  
    if (rank & statistic) c(pval = pval, rank = as.integer(r),
                            statistic = drop(csq))
    else if(rank) c(pval = pval, rank = as.integer(r))
    else if(statistic) c(pval = pval, statistic = drop(csq))
    else pval
}
