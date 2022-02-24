get_rand_error <- function(family, param, n) {
  switch(family,
         normal = rnorm(n, mean = 0, sd = param),
         uniform = runif(n, min = -param, max = param),
         logistic = rlogis(n, location = 0, scale = param),
         t = rt(n, df = 5) / sqrt(5 / 3) * param)
}

get_param <- function(family) {
  switch(family,
         normal = sqrt(runif(1, min = 0.5, max = 1.5)),
         uniform = runif(1, min = 1.2, max = 2.1),
         logistic = runif(1, min = 0.4, max = 0.7),
         t = sqrt(runif(1, min = 0.5, max = 1.5)))
}

# import igraph graph_from adjacency_matrix topo_sort
# family can be "normal" or "uniform" or "logistic" or  "t" or " mixed"
simulate_from_dag <- function(amat, n = 1000, ns = NULL, seed = NA,
                              family = "normal", param = NULL, 
                              perc = NULL, distr = NULL,
                              set.coef = FALSE){
  g <- igraph::graph_from_adjacency_matrix(t(amat), mode = "directed",
                                           weighted = TRUE)
  # topological order of nodes
  topo.ord <- as.vector(igraph::topo_sort(g))
  
  # check if weight matrix is unnamed
  # if (is.null(dimnames(amat)) && is.null(ns))
  #   message('Unnamed nodes in DAG')
  if (!is.null(ns)) dimnames(amat) <- list(ns, ns)
  
  # retrieve node names
  if (!is.null(dimnames(amat))) ns <- dimnames(amat)[[1]]
  
  p <- ncol(amat)
  data <- matrix(NA, nrow = n, ncol = p)
  
  lst <- list()
  
  # set seed if available
  if (!is.na(seed)) set.seed(seed)
  
  family <- match.arg(family, c("normal", "uniform", "logistic", "t", "mixed"))
  if (family == "normal") {
    if(is.null(param)) {
      param <- sqrt(runif(p, min = 0.5, max = 1.5))
      lst[["param"]] <- param
    }
    eps <- matrix(rnorm(n * p, mean = 0, sd = param), nrow = n, byrow = TRUE)
  }
  if (family == "uniform") {
    if(is.null(param)) {
      param <- runif(p, min = 1.2, max = 2.1)
      lst[["param"]] <- param
    }
    eps <- matrix(runif(n * p, min = -param, max = param), nrow = n, byrow = TRUE)
  }
  if (family == "logistic") {
    if(is.null(param)) {
      param <- runif(p, min = 0.4, max = 0.7)
      lst[["param"]] <- param
    }
    eps <- matrix(rlogis(n * p, location = 0, scale = param), nrow = n, byrow = TRUE)
  }
  if (family == "t") {
    if(is.null(param)) {
      param <- sqrt(runif(p, min = 0.5, max = 1.5))
      lst[["param"]] <- param
    }
    # variance = 5 / (5 - 2) = 5 / 3
    eps <- matrix(rt(n * p, df = 5) / sqrt(5 / 3) * param, nrow = n, byrow = TRUE)
  }
  if (family == "mixed") {
    if (is.null(distr)) {
      perc <- perc / sum(perc)
      # plus one to make from *1* to *4*
      distr <- findInterval(1:p, p * cumsum(perc),
                            rightmost.closed = TRUE, left.open = TRUE) + 1 
      # shuffle
      distr <- c("normal", "uniform", "logistic", "t")[sample(distr)]
      lst[["distr"]] <- distr
      }
    if (is.null(param)) {
      param <- sapply(distr, get_param, USE.NAMES = FALSE)
      lst[["param"]] <- param
    }
    eps <- mapply(get_rand_error, family = distr, param = param, n = n)
  }
  
  # get edge coefficients if set.coef is flagged true
  if (set.coef) {
    n.edge <- sum(amat != 0)
    amat[amat != 0] <- runif(n.edge, min = 0.1, max = 2) * sample(c(-1, 1), n.edge, replace = TRUE)
    lst[["amat"]] <- amat
  }
  
  for (i in topo.ord) {
    pa <- which(amat[i, ] != 0) # parent nodes
    if (length(pa) == 0) data[, i] <- eps[, i]
    else if (length(pa) == 1) data[, i] <- data[, pa] * amat[i, pa] + eps[, i]
    else data[, i] <- data[, pa] %*% amat[i, pa] + eps[, i]
  }
  if (!is.null(ns)) dimnames(data) <- list(NULL, ns)
  
  lst[["data"]] <- as.data.frame(data)
  
  if (length(lst) > 1) return(lst)
  else return(lst[[1]])
}
