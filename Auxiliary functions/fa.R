fa <- function(node, DAG) {
  pa <- which(DAG[,node] != 0)
  fa <- c(node, pa)
  return(fa)
}
