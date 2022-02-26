# randomly choosing x and y given the adjacency matrix
selectVars <- function(amat){
  # x <- sample(which(colSums(amat) > 0), 1)
  p <- ncol(amat)
  nd <- sapply(1:p, function(i)
    length(pcalg::possDe(m = amat, x = i, type = "dag", possible = FALSE))) - 1
  x <- sample(1:p, size = 1, prob = nd / sum(nd))
  if (!is.null(names(x))) names(x) <- NULL
  xd <- setdiff(pcalg::possDe(m = amat, x = x, type = "dag", possible = FALSE), x)
  if (length(xd) == 1) y <- xd
  else y <- sample(xd, 1)
  list(x = x, y = y)
}

# uniformly assigning edge weights
getWgt <- function(m, min = 0.1, max = 1) {
  w <- runif(m, min = min, max = max)
  sgn <- sample(c(-1, 1), size = m, replace = TRUE)
  w * sgn
}

# sampling unweighted and weighted adj matrix and setting x, y
simDAG <- function(n, d, method = "er", min = 1, max = 2) {
  g <- randDAG(n, d, method = method, weighted = TRUE,
               wFUN = list(getWgt, min = min, max = max))
  w <- t(as(g, "matrix"))
  dimnames(w) <- NULL
  a <- w
  a[a != 0] <- 1
  l <- selectVars(a)
  return(c(list(amat = a, bmat = w), l))
}

# helper for simFromDAG
getRandErr <- function(family, param, n) {
  switch(family,
         normal = rnorm(n, mean = 0, sd = param),
         uniform = runif(n, min = -param, max = param),
         logistic = rlogis(n, location = 0, scale = param),
         t = rt(n, df = 5) / sqrt(5 / 3) * param)
}

# helper for simFromDAG
getParam <- function(family) {
  switch(family,
         normal = sqrt(runif(1, min = 0.5, max = 1.5)),
         uniform = runif(1, min = 1.2, max = 2.1),
         logistic = runif(1, min = 0.4, max = 0.7),
         t = sqrt(runif(1, min = 0.5, max = 1.5)))
}

# simulating data from a given DAG
# family can be "normal" or "uniform" or "logistic" or  "t" or " mixed"
simFromDAG <- function(amat, n = 1000, ns = NULL, seed = NA,
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
      param <- sapply(distr, getParam, USE.NAMES = FALSE)
      lst[["param"]] <- param
    }
    eps <- mapply(getRandErr, family = distr, param = param, n = n)
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
