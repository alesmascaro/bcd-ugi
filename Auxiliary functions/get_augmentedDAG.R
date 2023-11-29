get_augmentedDAG <- function(DAG, TARGETS, PARENTS) {
  AugDAG <- DAG
  whcTARGETS <- which(TARGETS == 1)
  if (length(whcTARGETS) != 0) {
    AugDAG[,whcTARGETS] <- PARENTS[,whcTARGETS]
  }
  AugDAG <- cbind(rbind(AugDAG, TARGETS),0)
  return(AugDAG)
}
