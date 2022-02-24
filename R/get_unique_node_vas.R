get_unique_node_vas <- function(vas) {
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

# test <- function(vas) {
#   vas_unique <- get_unique_node_vas(vas)
#   ns <- Reduce(union, vas)
#   ns_unique <- Reduce(union, vas_unique)
#   has_unique <- sapply(seq_along(vas_unique), function(i) {
#     ns_other <- Reduce(union, vas_unique[-i])
#     !all(vas_unique[[i]] %in% ns_other)})
#   print(paste("Cover all?", length(ns_unique) == length(ns)))
#   print(paste("Has unique?", all(has_unique)))
#   print(paste("Sets reduced from", length(vas), "to", length(vas_unique)))
# }
# 
# vas <- replicate(50, {
#   len <- sample(1:40, 1)
#   sample(1:100, len)
# }, simplify = FALSE)
# test(vas)