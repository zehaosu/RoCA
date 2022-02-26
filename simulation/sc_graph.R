library(readxl)
library(dplyr)
library(ggplot2)
library(plotrix)

source("path/to/setup.R")

dat <- read_xls("path/to/data/1. cd3cd28.xls")
dat <- dat %>% log %>% scale(scale = FALSE)

amat <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
                   1,0,0,0,0,0,0,1,1,0,0,
                   0,0,0,0,1,0,0,0,0,0,0,
                   0,0,1,0,1,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,
                   0,1,0,0,0,0,0,1,0,0,0,
                   0,0,0,0,1,0,0,1,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,
                   0,0,1,1,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,1,1,0,0,
                   0,0,0,0,0,0,0,1,1,0,0), ncol=11))

amat.sachs <- t(matrix(c(0,0,0,0,0,0,0,1,1,0,0,
                1,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,1,0,0,0,0,0,0,
                0,0,1,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,1,0,0,0,
                0,0,0,0,0,1,0,1,0,0,0,
                0,0,0,0,0,0,0,0,1,0,0,
                0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,0,0,0,1,1,0,0), ncol=11))


names <- c("Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA",
           "PKC", "p38", "JNK")
colnames(dat) <- names
rownames(amat) <- names
colnames(amat) <- names
rownames(amat.sachs) <- names
colnames(amat.sachs) <- names

descendants <- lapply(1:11, possDe, m = amat, type = "dag")

df <- data.frame(x = NULL, y = NULL)

for (i in seq_along(descendants)) {
    if (length(descendants[[i]]) == 1) next
    de.rm <- setdiff(descendants[[i]], i)
    for (j in seq_along(de.rm)) {
        df <- rbind(df, c(i, de.rm[j]))
    }
}

res <- apply(df, 1, function(x) {
    vas <- adjustment(amat, "dag", x = x[1], y = x[2],
                      set.type = "all")
    print(length(vas))
    wald_test(dat, vas = vas, ns = names,
              method = "plugin", x = x[1], y = x[2],
              rank.method = "bic", rank = TRUE, statistic = TRUE)
})

res[1, ] * nrow(df) # Bonferroni adj pval

res.min <- apply(df, 1, function(x) {
    vas <- adjustment(amat, "dag", x = x[1], y = x[2],
                      set.type = "minimal")
    vas <- union(unname(vas), 
                 list(setdiff(1:nrow(amat),
                              c(x[1],
                                pcalg:::forbiddenNodes(amat, x[1],
                                                       x[2], "dag")))))
    wald_test(dat, vas = vas, ns = names,
              method = "plugin", x = x[1], y = x[2],
              rank.method = "none", rank = TRUE, statistic = TRUE)
})

res.min[1, ] * nrow(df) # Bonferroni adj pval
