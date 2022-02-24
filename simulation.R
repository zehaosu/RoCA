no.tasks <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
task.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
fname <- paste0("path/cache/graph-",
                sprintf("%04d", task.id), ".csv")

source("path/setup.R")
source("path/run_simulation.R")

set.seed(2022)
seeds <- sample(1:1e6, size = no.tasks, replace = FALSE)
set.seed(seeds[task.id])

g.size <- sample(c(10, 20), 1)
nbhd.size <- sample(c(2, 3, 4, 5), 1)
l <- simulate_graph(n = g.size, d = nbhd.size)
p <- nrow(l$bmat)
family <- sample(c("normal", "uniform", "logistic", "t"), 1)
if (family == "normal") param <- sqrt(runif(p, min = 0.5, max = 1.5))
if (family == "uniform") param <- runif(p, min = 1.2, max = 2.1)
if (family == "logistic") param <- runif(p, min = 0.4, max = 0.7)
if (family == "t") param <- sqrt(runif(p, min = 0.5, max = 1.5))

for (n in c(100, 400)){
  for (m in c(50, 100, 200, 400)) {
      df <- run_slow_simulation(n = n, m = m, bmat = l$bmat,
                                x = l$x, y = l$y,
                                family = family, param = param,
                                R = 20, max.time = 60)
    if(is.null(df) || nrow(df) == 0) next
    df$n <- n
    df$m <- m
    df$family <- family
    df$g.size <- g.size
    df$nbhd.size <- nbhd.size
    df$seed <- seeds[task.id]
    write.table(df, file = fname, append = TRUE, sep = ",",
                col.names = !file.exists(fname), row.names = FALSE)
  }
}
