as_graphNEL <- function(DAG) {
  q <- ncol(DAG)
  nodes <- as.character(1:q)
  ft <- which(DAG != 0, T)
  graphNEL <- graph::ftM2graphNEL(ft, V = nodes, edgemode = "directed")
  return(graphNEL)
}
