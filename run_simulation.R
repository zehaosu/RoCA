run_simulation_inner <- function(n, m, bmat, x, y, family, param, M,
                                 level, max.time) {
    # get unweighted adjacency matrix
    amat <- bmat
    amat[amat != 0] <- 1
  
    dat <- simulate_from_dag(amat = bmat, n = n,
                             family = family, param = param)
  
    if (family == "normal") {
        # tranform data to format required by ges
        score <- new("GaussL0penObsScore", dat)
        # estimating the CPDAG C
        cpdag <- ges(score)$essgraph
        amat.learn <- t(as(cpdag, "matrix"))
        g.type <- "cpdag"
    }
    else {
        g.learn <- pcalg::lingam(dat)
        amat.learn <- g.learn$Bpruned
        g.type <- "dag"
    }
    amat.learn[amat.learn != 0] <- 1
    
    if (!(y %in% possDe(amat.learn, x, type = g.type,
                        possible = FALSE))) {
        stop("Y is not a descendent of X")
    }
  
    vas.learn <- pcalg::adjustment(amat.learn, amat.type = g.type,
                                   x = x, y = y, set.type = "minimal")
    ## get non-forbidden nodes relative to (x, y)
    nonforb <- setdiff(1:nrow(amat.learn),
                       c(x, pcalg:::forbiddenNodes(amat.learn, x, y,
                                                   g.type)))
    # use minimal sets only with a unique node
    if (length(vas.learn) == 1) {
        vas.test.min.unique <- unname(vas.learn)
    } else {
        vas.test.min.unique <- get_unique_node_vas(unname(vas.learn))
    }
    
    vas.test.minprune <- union(vas.test.min.unique, list(nonforb))
    vas.test.all <- R.utils::withTimeout({
        pcalg::adjustment(amat.learn, amat.type = g.type,
                          x = x, y = y, set.type = "all")
    }, timeout = max.time, onTimeout = "warning")
    
    if (class(vas.test.all) == "character") vas.test.all <- NULL
    
    if (!is.null(vas.test.all) & (length(vas.test.all) < m)) {
        
        if (length(vas.test.all) == 1 ||
            length(vas.test.minprune) == 1) {
            stop("1 adjustment set available; not testable")
        }
  
        valid.all <- sapply(vas.test.all,
                            function(z) gac(amat = amat, x = x, y = y,
                                            z = z, type = "dag")$gac)

        valid.minprune <- sapply(vas.test.minprune,
                                 function(z) gac(amat = amat, x = x,
                                                 y = y, z = z,
                                                 type = "dag")$gac)
    
        pval <- replicate(M, {
            dat_ <-  simulate_from_dag(amat = bmat, n = m,
                                       family = family, param = param)
            c(pval.all = wald_test(data = dat_, vas = vas.test.all,
                                   x = x, y = y,
                                   method = "plugin",
                                   rank.method = "bic"),
              pval.minprune = wald_test(data = dat_,
                                        vas = vas.test.minprune,
                                        x = x, y = y,
                                        method = "plugin",
                                        rank.method = "none"))
        })
        
        if (family == "normal" || family == "t") {
            dmat <- diag(param ^ 2)
        }
        if (family == "uniform") dmat <- diag(param ^ 2 / 3)
        if (family == "logistic") dmat <- diag(pi ^ 2 / 3 * param ^ 2)
        
        if (all(valid.all)) {
            region.all <- 1
        } else if (sum(valid.all) > 0) {
            region.all <- 3
        } else {
            cmat <- solve(diag(ncol(amat)) - bmat)
            cmat <- cmat %*% dmat %*% t(cmat)
            betas <- sapply(vas.test.all, beta, cmat = cmat,
                            x = x, y = y)
            if (all(abs(betas - mean(betas)) <
                    .Machine$double.eps ^ 0.5)) {
                region.all <- 2
            } else {
                region.all <- 3
            }
        }
    
        if (all(valid.minprune)) {
            region.minprune <- 1
        } else if (sum(valid.minprune) > 0) {
            region.minprune <- 3
        } else {
            cmat <- solve(diag(ncol(amat)) - bmat)
            cmat <- cmat %*% dmat %*% t(cmat)
            betas <- sapply(vas.test.minprune, beta, cmat = cmat,
                            x = x, y = y)
            if (all(abs(betas - mean(betas)) <
                    .Machine$double.eps ^ 0.5)) {
                region.minprune <- 2
            } else {
                region.minprune <- 3
            }
        }
        
        area <- apply(pval, 1, function(y) {
            idx <- !is.na(y)
            y.ord <- sort(c(0, 1, y[idx]))
            MESS::auc(y.ord,
                      sapply(y.ord,
                             function(x) mean(y <= x, na.rm = TRUE)))
        })
        names(area) <- c("auc.all", "auc.minprune")
        count <- apply(pval, 1,
                       function(y) sum(y <= level, na.rm = TRUE))
        names(count) <- c("rej.all", "rej.minprune")
        count.n <- apply(pval, 1, function(y) sum(!is.na(y)))
        names(count.n) <- c("n.all", "n.minprune")
        return(c(area, count, count.n,
                 region.all = region.all,
                 region.minprune = region.minprune,
                 k.all = length(vas.test.all),
                 k.minprune = length(vas.test.minprune)))
    } else {
        
        if (length(vas.test.minprune) == 1) {
            stop("1 adjustment set available; not testable")
        }

        valid.minprune <- sapply(vas.test.minprune,
                                 function(z) gac(amat = amat, x = x,
                                                 y = y, z = z,
                                                 type = "dag")$gac)
        
        pval.minprune <- replicate(M, {
            dat_ <-  simulate_from_dag(amat = bmat, n = m,
                                       family = family, param = param)
            wald_test(data = dat_, vas = vas.test.minprune,
                      x = x, y = y,
                      method = "plugin", rank.method = "none")
        })
        
        if (family == "normal" || family == "t") {
            dmat <- diag(param ^ 2)
        }
        if (family == "uniform") dmat <- diag(param ^ 2 / 3)
        if (family == "logistic") dmat <- diag(pi ^ 2 / 3 * param ^ 2)

        if (all(valid.minprune)) {
            region.minprune <- 1
        } else if (sum(valid.minprune) > 0) {
            region.minprune <- 3
        } else {
            cmat <- solve(diag(ncol(amat)) - bmat)
            cmat <- cmat %*% dmat %*% t(cmat)
            betas <- sapply(vas.test.minprune, beta, cmat = cmat,
                            x = x, y = y)
            if (all(abs(betas - mean(betas)) <
                    .Machine$double.eps ^ 0.5)) {
                region.minprune <- 2
            } else {
                region.minprune <- 3
            }
        }

        idx <- !is.na(pval.minprune)
        pval.minprune.ord <- sort(c(0, 1, pval.minprune[idx]))
        auc.minprune <- MESS::auc(pval.minprune.ord,
                                  sapply(pval.minprune.ord,
                                         function(x) {
                                             mean(pval.minprune <= x,
                                                  na.rm = TRUE)
                                         }))
        rej.minprune <- sum(pval.minprune <= level, na.rm = TRUE)
        n.minprune <- sum(!is.na(pval.minprune))
        
        return(c(auc.all = NA, auc.minprune = auc.minprune,
                 rej.all = NA, rej.minprune = rej.minprune,
                 n.all = NA, n.minprune = n.minprune,
                 region.all = NA, region.minprune = region.minprune,
                 k.all = NA, k.minprune = length(vas.test.minprune)))
    }
}


run_simulation <- function(n, m, bmat, x, y, family = "normal",
                           param = NULL, R = 20, M = 100,
                           level = 0.05, max.time = 60, dir = NULL) {
    family <- match.arg(family,
                        c("normal", "uniform", "logistic", "t"))
    p <- nrow(bmat)
    if (is.null(param)) {
        if (family == "normal") {
            param <- sqrt(runif(p, min = 0.5, max = 1.5))
        }
        if (family == "uniform") {
            param <- runif(p, min = 1.2, max = 2.1)
        }
        if (family == "logistic") {
            param <- runif(p, min = 0.4, max = 0.7)
        }
        if (family == "t") {
            param <- sqrt(runif(p, min = 0.5, max = 1.5))
        }
    }
    out <- NULL
    for(i in 1:R) {
        res <- try(run_simulation_inner(n = n, m = m, bmat = bmat,
                                        x = x, y = y,
                                        family = family,
                                        param = param, M = M,
                                        level = level,
                                        max.time = max.time),
                   silent = TRUE)
        if (class(res) == "try-error") next
        out <- rbind(out, res)
    }

    return(as.data.frame(out))
}
