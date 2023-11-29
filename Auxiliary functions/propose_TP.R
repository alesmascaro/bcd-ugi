propose_TP <- function(k, currentDAG, currentTARGETS, currentPARENTS, fast) {
  q <- ncol(currentDAG)
  modifiedDAG <- modify_parents(currentDAG, currentTARGETS[k,], currentPARENTS[[k]])
  AugmentedDAG <- get_augmentedDAG(modifiedDAG, currentTARGETS[k,], currentPARENTS[[k]])

  whccurrentTARGETS <- which(currentTARGETS[k,] == 1)
  whcnocurrentTARGETS <- setdiff(1:q, whccurrentTARGETS)

  # Define the set of possible "ordinary" operations
  A_na <- modifiedDAG
  diag(A_na) <- NA

  id_set = c()
  dd_set = c()
  rd_set = c()


  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  set_id = set_id[which(set_id[,2] %in% whccurrentTARGETS),]
  if(length(set_id) != 0){
    if (length(set_id) == 2) {
      id_set = cbind(1, t(as.matrix(set_id)))
    } else {
      id_set = cbind(1, set_id)
    }
  }

  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  set_dd = set_dd[which(set_dd[,2] %in% whccurrentTARGETS),]
  if(length(set_dd) != 0){
    if (length(set_dd) == 2) {
      dd_set = cbind(2, t(as.matrix(set_dd)))
    } else {
      dd_set = cbind(2, set_dd)
    }
  }

  ## set of nodes for rd
  set_rd = which(A_na == 1, TRUE)
  set_rd = set_rd[intersect(which(set_rd[,1] %in% whccurrentTARGETS), which(set_rd[,2] %in% whccurrentTARGETS)),]
  if(length(set_rd) != 0){
    if (length(set_rd) == 2) {
      rd_set = cbind(3, t(as.matrix(set_rd)))
    } else {
      rd_set = cbind(3, set_rd)
    }
  }

  ### Simple Operations
  zerostep_O <- c(0,0,0)
  onestep_list <- list(id_set, dd_set, rd_set)

  # browser()

  whc_onestep <- which(sapply(onestep_list, length) != 0)
  if (length(whc_onestep) != 0) {
    onestep_O = do.call(rbind, onestep_list[whc_onestep])
  } else {
    onestep_O = c()
  }

  # Add soft interventions and TARGET removal

  si_set <- c()
  sr_set <- c()

  set_si <- whcnocurrentTARGETS
  if (length(set_si) != 0) {
    si_set <- cbind(1, q+1, set_si)
    length_si <- nrow(si_set)
  }

  set_sr <- which(apply(currentDAG == modifiedDAG, 2, prod) == 1)
  set_sr <- intersect(set_sr, whccurrentTARGETS)
  if (length(set_sr) != 0) {
    sr_set <- cbind(2,q+1, set_sr)
    length_sr <- nrow(sr_set)
  }

  # browser()

  soft_O <- rbind(si_set, sr_set)

  # browser()
  if (fast == TRUE) {
    if (length(onestep_O) != 0) {
      O <- rbind(onestep_O, soft_O)
    } else {
      O <- soft_O
    }
    opcard <- nrow(O)
    repeat {
      m <- sample(1:nrow(O), 1)
      op.type <- O[m,1]
      possibleAugmentedDAG <- operation(O[m,1], AugmentedDAG, O[m,2:3])
      possibleDAG <- possibleAugmentedDAG[1:q, 1:q]
      acyclicity_ok <- gRbase::is.DAG(possibleDAG)
      verify <- acyclicity_ok
      if (verify == TRUE) {
        break
      }
    }
  } else {
    if (length(id_set) != 0) {
      is_id_valid <- vector("logical", nrow(id_set))
      for (i in 1:nrow(id_set)) {
        possibleAugmentedDAG <- operation(id_set[i,1], AugmentedDAG, id_set[i,2:3])
        possibleDAG <- possibleAugmentedDAG[1:q, 1:q]
        is_id_valid[i] <- as.logical(gRbase::is.DAG(possibleDAG))
      }
      id_set <- id_set[is_id_valid,]
    }
    if (length(rd_set) != 0) {
      is_rd_valid <- vector("logical", nrow(rd_set))
      for (i in 1:nrow(rd_set)) {
        possibleAugmentedDAG <- operation(rd_set[i,1], AugmentedDAG, rd_set[i,2:3])
        possibleDAG <- possibleAugmentedDAG[1:q, 1:q]
        is_rd_valid[i] <- as.logical(gRbase::is.DAG(possibleDAG))
      }
      rd_set <- rd_set[is_rd_valid,]
    }
    onestep_list <- list(id_set, dd_set, rd_set)
    whc_onestep <- which(sapply(onestep_list, length) != 0)
    if (length(whc_onestep) != 0) {
      onestep_O = do.call(rbind, onestep_list[whc_onestep])
    } else {
      onestep_O = c()
    }
    if (length(onestep_O) != 0) {
      O <- rbind(onestep_O, soft_O)
    } else {
      O <- soft_O
    }
    opcard <- nrow(O)

    m <- sample(1:nrow(O), 1)
    op.type <- O[m,1]
    possibleAugmentedDAG <- operation(O[m,1], AugmentedDAG, O[m,2:3])
    possibleDAG <- possibleAugmentedDAG[1:q, 1:q]
  }

  proposedTARGETS <- currentTARGETS
  proposedTARGETS[k,] <- possibleAugmentedDAG[q+1,1:q]
  whcproposedTARGETS <- which(proposedTARGETS[k,] == 1)
  proposed_kthPARENTS <- matrix(0,q,q)
  proposed_kthPARENTS[,whcproposedTARGETS] <- possibleDAG[,whcproposedTARGETS]

  whc_changedTARGETS <- which(proposedTARGETS[k,] != currentTARGETS[k,])
  whc_changedPARENTS <- which(colSums(proposed_kthPARENTS == currentPARENTS[[k]]) != q)
  op.node <- union(whc_changedTARGETS, whc_changedPARENTS)

  out <- list(proposedTARGETS = proposedTARGETS, proposed_kthPARENTS = proposed_kthPARENTS,
              op.node = op.node, op.card = opcard)

  return(out)
}
