get_modifiedDAGs <- function(DAG, TARGETS, PARENTS) {
  K <- length(PARENTS)
  DAGs <- lapply(1:K, function(k) modify_parents(DAG, TARGETS[k,], PARENTS[[k]]))
  return(DAGs)
}
