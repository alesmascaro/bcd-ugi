modify_parents <- function(DAG, kTARGETS, kPARENTS) {
  t <- which(kTARGETS != 0)
  if (length(t) != 0) {
    DAG[,t] <- kPARENTS[,t]
  }
  return(DAG)
}
