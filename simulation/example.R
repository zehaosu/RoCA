source("path/setup.R")

amat.0 <- matrix(
  c(0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  10, 10, byrow = TRUE)

amat.1 <- matrix(
  c(0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  10, 10, byrow = TRUE)

amat.2 <- matrix(
  c(0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  10, 10, byrow = TRUE)

amat.3 <- matrix(
  c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  10, 10, byrow = TRUE)

amat.3_ <- matrix(
  c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  10, 10, byrow = TRUE)

ns <- c("X", "Y", "A1", "A2", "B1", "B2", "V", "D", "R", "F")

vas.all.0 <- adjustment(amat.0, "dag", 1, 2, "all")
vas.all.1 <- adjustment(amat.1, "dag", 1, 2, "all")
vas.all.3 <- adjustment(amat.3, "dag", 1, 2, "all")
vas.all.3_ <- adjustment(amat.3_, "dag", 1, 2, "all")

vas.min.0 <- adjustment(amat.0, "dag", 1, 2, "minimal")
vas.min.1 <- adjustment(amat.1, "dag", 1, 2, "minimal")
vas.min.3 <- adjustment(amat.3, "dag", 1, 2, "minimal")
vas.min.3_ <- adjustment(amat.3_, "dag", 1, 2, "minimal")

nonforb.0 <- setdiff(1:nrow(amat.0),
                     c(1, pcalg:::forbiddenNodes(amat.0, 1, 2, "dag")))
nonforb.1 <- setdiff(1:nrow(amat.1),
                     c(1, pcalg:::forbiddenNodes(amat.1, 1, 2, "dag")))
nonforb.3 <- setdiff(1:nrow(amat.3),
                     c(1, pcalg:::forbiddenNodes(amat.3, 1, 2, "dag")))
nonforb.3_ <- setdiff(1:nrow(amat.3_),
                      c(1,
                        pcalg:::forbiddenNodes(amat.3_, 1, 2, "dag")))

vas.minpruned.0 <- c(get_unique_node_vas(vas.min.0), list(nonforb.0))
vas.minpruned.1 <- c(get_unique_node_vas(vas.min.1), list(nonforb.1))
vas.minpruned.3 <- c(get_unique_node_vas(vas.min.3), list(nonforb.3))
vas.minpruned.3_ <- c(vas.min.3_, list(nonforb.3_))

library(foreach)
set.seed(2022)
N <- c(25, 100, 400, 1600)
df <- foreach(n = N, .combine = "rbind") %do% {
  foreach(i = 1:100, .combine = "rbind") %do% {
      dat <- simulate_from_dag(amat = amat.0, n = n, family = "normal",
                               param = 1, ns = ns)
      c(
        pval.minpruned.0 = wald_test(data = dat, vas = vas.minpruned.0,
                                     x = 1, y = 2, method = "plugin",
                                     rank.method = "none"),
        pval.minpruned.1 = wald_test(data = dat, vas = vas.minpruned.1,
                                     x = 1, y = 2, method = "plugin",
                                     rank.method = "none"),
        pval.minpruned.3 = wald_test(data = dat, vas = vas.minpruned.3,
                                     x = 1, y = 2, method = "plugin",
                                     rank.method = "none"),
        pval.minpruned.3_ = wald_test(data = dat,
                                      vas = vas.minpruned.3_,
                                      x = 1, y = 2, method = "plugin",
                                      rank.method = "none"),
        pval.all.0 = wald_test(data = dat, vas = vas.all.0,
                               x = 1, y = 2, method = "plugin",
                               rank.method = "bic"),
        pval.all.1 = wald_test(data = dat, vas = vas.all.1,
                               x = 1, y = 2, method = "plugin",
                               rank.method = "bic"),
        pval.all.3 = wald_test(data = dat, vas = vas.all.3,
                               x = 1, y = 2, method = "plugin",
                               rank.method = "bic"),
        pval.all.3_ = wald_test(data = dat, vas = vas.all.3_,
                                x = 1, y = 2, method = "plugin",
                                rank.method = "bic"),
        n = n
      )
  }
}
df <- as.data.frame(df)

write.table(df, file = "path/cache/example.csv", sep = ",",
            col.names = TRUE, row.names = FALSE)
