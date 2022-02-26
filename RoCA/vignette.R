source("utils.R")
source("RoCA.R")
source("utils-vignette.R")

library(pcalg)
library(igraph)

set.seed(32)
# random ER graph of 20 nodes and expected neighbourhood size 6
l <- simDAG(20, 6)
n <- 1000
dat <- simFromDAG(amat = l$bmat, n = n, family = "normal", param = 1)
RoCA(data = dat, x = l$x, y = l$y,
         g = l$amat, strategy = "min+")

# random ER graph of 7 nodes and expected neighbourhood size 2
l <- simDAG(7, 2)
dat <- simFromDAG(amat = l$bmat, n = n, family = "normal", param = 1)
RoCA(data = dat, x = l$x, y = l$y,
         g = l$amat, strategy = "all")

# 8 adjustment sets used by "all"
# feeding a user-defined contrast matrix...
contr.alt <- t(sapply(1:7, function(j)
  c(rep(-1 / 7, j - 1), 1, rep(-1 / 7, 8 - j))))
RoCA(data = dat, x = l$x, y = l$y,
         g = l$amat, strategy = "all",
         contrast = contr.alt)
# slightly different result
