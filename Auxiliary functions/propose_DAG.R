propose_DAG <- function(currentDAG, currentTARGETS, currentPARENTS, fast) {
  A <- currentDAG
  q <- ncol(A)
  K <- nrow(currentTARGETS)
  A_na <- A
  diag(A_na) <- NA

  # Define the set of possible operations
  # The cardinality of O will change depending on how many edges are present in the DAG

  id_set = c()
  dd_set = c()
  rd_set = c()



  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }

  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }

  ## set of nodes for rd
  set_rd = which(A_na == 1, TRUE)
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  # browser()
  if (fast == TRUE) {
    O = rbind(id_set, dd_set, rd_set)
    # browser()
    max.it <- 0
    repeat {
      # max.it <- max.it + 1
      m <- sample(1:nrow(O), 1)
      possibleDAG <- operation(O[m,1], A, O[m,2:3])
      modifiedDAGs <- get_modifiedDAGs(possibleDAG, currentTARGETS, currentPARENTS)
      acyclicity_ok <- prod(sapply(1:length(modifiedDAGs), function(k) gRbase::is.DAG(modifiedDAGs[[k]])))
      verify <- acyclicity_ok
      if (verify == TRUE) {
        break
      }
    }
    A_next <- possibleDAG
  } else {
    if (length(id_set) != 0) {
      is_id_valid <- vector("logical", nrow(id_set))
      for (i in 1:nrow(id_set)) {
        possibleDAG <- operation(id_set[i,1], A, id_set[i,2:3])
        modifiedDAGs <- get_modifiedDAGs(possibleDAG, currentTARGETS, currentPARENTS)
        is_id_valid[i] <- as.logical(prod(sapply(1:length(modifiedDAGs), function(k) gRbase::is.DAG(modifiedDAGs[[k]]))))
      }
      id_set <- id_set[is_id_valid,]
    }
    if (length(rd_set) != 0) {
      is_rd_valid <- vector("logical", nrow(rd_set))
      for (i in 1:nrow(rd_set)) {
        possibleDAG <- operation(rd_set[i,1], A, rd_set[i,2:3])
        modifiedDAGs <- get_modifiedDAGs(possibleDAG, currentTARGETS, currentPARENTS)
        is_rd_valid[i] <- as.logical(prod(sapply(1:length(modifiedDAGs), function(k) gRbase::is.DAG(modifiedDAGs[[k]]))))
        if (max_indegree != FALSE) {
          max_pa <- max(colSums(possibleDAG))
          indegree_ok <- max_pa <= max_indegree
        }
        is_rd_valid[i] <- as.logical(is_rd_valid[i]*indegree_ok)
      }
      rd_set <- rd_set[is_rd_valid,]
    }

    O <- rbind(id_set, dd_set, rd_set)
    m <- sample(1:nrow(O), 1)
    A_next <- operation(O[m,1], A, O[m,2:3])
  }

  op.type <- O[m,1]
  if (op.type == 3) {
    op.node <- O[m,-1]
  } else {
    op.node <- O[m,3]
  }
  proposed.opcard <- nrow(O)

  return(list(proposedDAG = A_next,
              op.type = op.type, op.node = op.node, op.card = proposed.opcard))

}
