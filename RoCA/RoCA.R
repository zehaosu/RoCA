RoCA <- function(data, x, y, g = NULL,
                    strategy = c("all", "min+"),
                     contrast = NULL, time.max = 60)
{
    ## Purpose: Perform *Ro*bustness test of total (causal) effect estimates via *C*ovariate *A*djustment
    ## --------------------------------------------------------------------------------------------------
    ## Arguments:
    ## - data: data for testing
    ## - x: column number of exposure variable
    ## - y: column number of outcome variable
    ## - g: candidate causal DAG, an object of "matrix", "graphNEL" or "igraph"
    ## - strategy: one of "all" or "min+"
    ## - contrast: contrast matrix of appropriate dimensions
    ##             if NULL, the default contrast matrix (see paper) will be used
    ## - time.max: time limit on adjustment set extraction

    require(pcalg)
    require(igraph)
  
    # create adjancency matrix from graph object
    g.class <- class(g)[1]
    if (g.class == "graphNEL") {
        amat <- as(g, "matrix")
        dimnames(amat) <- NULL
    } else if (g.class == "igraph") {
        amat <- t(as_adjacency_matrix(g, names = FALSE, sparse = FALSE))
    } else if (g.class == "matrix") {
        amat <- g
        dimnames(amat) <- NULL
    } else {
        stop("Incompatible class for candidate graph")
    }
    
    # set time limit
    setTimeLimit(cpu = time.max, elapsed = time.max, transient = TRUE)
    on.exit({
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    })
    
    set.type <- ifelse(strategy == "all", "all", "minimal")
    
    Z <- tryCatch({
        adjustment(amat = amat, amat.type = "dag", x = x, y = y,
                   set.type = set.type)
    }, error = function(e) {
        if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
            ## when reaching timeout, apply "min+" as adjustment set method
            ## override orginal method.set
            message("Timeout when getting adjustment sets, using strategy 'min+' instead")
            strategy <<- "min+"
            ## get minimal sets first
            adjustment(amat = amat, amat.type = "dag", x = x, y = y,
                       set.type = "minimal")
        } else {
            # error not related to timeout
            stop(e)
        }
    })
    
    if (strategy == "min+") {
        ## use minimal sets only with a unique node
        if (length(Z) == 1) {
            Z.unique <- unname(Z)
        } else {
            Z.unique <- getUniqueZ(unname(Z))
        }
        ## append all non-forbidden nodes
        nonforb <- setdiff(1:nrow(amat),
                           c(x, pcalg:::forbiddenNodes(amat, x, y, "dag")))
        Z <- c(Z.unique, list(nonforb))
    }
    
    if (length(Z) == 1) {
        stop("1 adjustment set available; not testable")
    }
  
    ## call the function for test
    test.info <- waldTest(data = data, x = x, y = y, Z = Z, rank.est = (strategy == "all"), contrast = contrast)
    
    res <- split(unname(test.info), names(test.info))
    res$Z <- Z
    res$strategy <- strategy
    res$contrast <- contrast
    class(res) <- "RoCA"
    res
}

print.RoCA <- function(x, digits = 4) {
    cat("Null hypothesis: Total effect estimates identify the same quantity\n")
    z <- c("p-value" = x$pval, statistic = x$statistic, df = x$rank)
    print.default(c(format(z[1:2], digits = digits), format(z[3], digits = 0)),
                  print.gap = 2L, quote = FALSE)
    cat(paste0(length(x$Z), " adjustment sets used, ", x$rank,
               " over-identifying constraints estimated\n"))
}
