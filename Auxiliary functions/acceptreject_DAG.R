acceptreject_DAG <- function(tXX_list, n_vec, currentTARGETS, currentPARENTS,
                             currentDAG, proposedDAG, node, op.type, op.card,
                             DW_a, DW_U, a_w, b_w,
                             fast, ordering) {
  q <- ncol(currentDAG)
  ## Logproposal ratio approximated to 1
  if (fast == TRUE) {
    logproposal.ratio <- log(1)
  } else {
    logproposal.ratio <- log(op.card/get_opDAGcard(proposedDAG, currentTARGETS, currentPARENTS))
  }

  ## Compute logprior ratio
  logDAGprior.ratio <- lbeta(a_w + sum(proposedDAG), b_w + q*(q-1)/2 - sum(proposedDAG)) -
    lbeta(a_w + sum(currentDAG), b_w + q*(q-1)/2 - sum(currentDAG))

  ## Compute log marginal likelihood ratio
  currentDAGs_list <- get_modifiedDAGs(currentDAG, currentTARGETS, currentPARENTS)
  proposedDAGs_list <- get_modifiedDAGs(proposedDAG, currentTARGETS, currentPARENTS)
  K <- length(currentDAGs_list)
  q <- ncol(currentDAG)

  if (length(node) == 1) { # I'm not reversing an edge: just compare the changed node

    ### Define families of changed node
    currentfamilies_list <- lapply(currentDAGs_list, function(k) fa(node, k))
    proposedfamilies_list <- lapply(proposedDAGs_list, function(k) fa(node, k))
    ### Define subsets of tXX based on the changed node and different DAGs (families) for each k
    current_fa.tXX <- lapply(1:K, function(k) tXX_list[[k]][currentfamilies_list[[k]], currentfamilies_list[[k]]])
    proposed_fa.tXX <- lapply(1:K, function(k) tXX_list[[k]][proposedfamilies_list[[k]], proposedfamilies_list[[k]]])

    ### Initialize lml
    currentDW_intblockslml <- 0
    proposedDW_intblockslml <- 0
    ### Compute lml for interventional blocks (if present)
    int_blocks <- which(currentTARGETS[,node] != 0)
    arethereint <- sum(int_blocks) != 0
    if (arethereint) {
      currentDW_intblockslml <- sapply(int_blocks, function(k) DW_blocklml(node, currentDAGs_list[[k]], current_fa.tXX[[k]], n_vec[k], DW_a, DW_U))
      proposedDW_intblockslml <- sapply(int_blocks, function(k) DW_blocklml(node, proposedDAGs_list[[k]], proposed_fa.tXX[[k]], n_vec[k], DW_a, DW_U))
      if (sum(currentDW_intblockslml) != sum(proposedDW_intblockslml)) print("Oh cazzo")
    }
    ### Compute lml for observational blocks (always present)
    obs_blocks <- setdiff(1:K, int_blocks)
    n_obs <- sum(n_vec[obs_blocks])
    #### Compute the corresponding statitic for the observational part from list of tXX
    current_fa.tXX.obs <- Reduce("+", current_fa.tXX[obs_blocks])
    proposed_fa.tXX.obs <- Reduce("+", proposed_fa.tXX[obs_blocks])

    currentDW_obsblockslml <- DW_blocklml(node, currentDAG, current_fa.tXX.obs, n_obs, DW_a, DW_U)
    proposedDW_obsblockslml <- DW_blocklml(node, proposedDAG, proposed_fa.tXX.obs, n_obs, DW_a, DW_U)


    ## Compute total current and proposed lml
    current_lml <- sum(c(as.numeric(currentDW_obsblockslml), as.numeric(currentDW_intblockslml)))
    proposed_lml <- sum(c(as.numeric(proposedDW_obsblockslml), as.numeric(proposedDW_intblockslml)))

  } else {

    currentfamilies1_list <- lapply(currentDAGs_list, function(k) fa(node[1], k))
    currentfamilies2_list <- lapply(currentDAGs_list, function(k) fa(node[2], k))
    proposedfamilies1_list <- lapply(proposedDAGs_list, function(k) fa(node[1], k))
    proposedfamilies2_list <- lapply(proposedDAGs_list, function(k) fa(node[2], k))

    current_fa1.tXX <- lapply(1:K, function(k) tXX_list[[k]][currentfamilies1_list[[k]], currentfamilies1_list[[k]]])
    current_fa2.tXX <- lapply(1:K, function(k) tXX_list[[k]][currentfamilies2_list[[k]], currentfamilies2_list[[k]]])
    proposed_fa1.tXX <- lapply(1:K, function(k) tXX_list[[k]][proposedfamilies1_list[[k]], proposedfamilies1_list[[k]]])
    proposed_fa2.tXX <- lapply(1:K, function(k) tXX_list[[k]][proposedfamilies2_list[[k]], proposedfamilies2_list[[k]]])


    current1DW_intblockslml <- 0
    current2DW_intblockslml <- 0
    proposed1DW_intblockslml <- 0
    proposed2DW_intblockslml <- 0

    ### First node
    int1_blocks <- which(currentTARGETS[,node[1]] != 0)
    arethereint1 <- sum(int1_blocks) != 0
    if (arethereint1) {
      current1DW_intblockslml <- sapply(int1_blocks, function(k) DW_blocklml(node[1], currentDAGs_list[[k]], current_fa1.tXX[[k]], n_vec[k], DW_a, DW_U))
      proposed1DW_intblockslml <- sapply(int1_blocks, function(k) DW_blocklml(node[1], proposedDAGs_list[[k]], proposed_fa1.tXX[[k]], n_vec[k], DW_a, DW_U))
    }
    obs1_blocks <- setdiff(1:K, int1_blocks)
    n_obs1 <- sum(n_vec[obs1_blocks])

    current_fa1.tXX.obs <- Reduce("+", current_fa1.tXX[obs1_blocks])
    proposed_fa1.tXX.obs <- Reduce("+", proposed_fa1.tXX[obs1_blocks])
    current1DW_obsblockslml <- DW_blocklml(node[1], currentDAG, current_fa1.tXX.obs, n_obs1, DW_a, DW_U)
    proposed1DW_obsblockslml <- DW_blocklml(node[1], proposedDAG, proposed_fa1.tXX.obs, n_obs1, DW_a, DW_U)

    ### Second node
    int2_blocks <- which(currentTARGETS[,node[2]] != 0)
    arethereint2 <- sum(int2_blocks) != 0
    if (arethereint2) {
      current2DW_intblockslml <- sapply(int2_blocks, function(k) DW_blocklml(node[2], currentDAGs_list[[k]], current_fa2.tXX[[k]], n_vec[k], DW_a, DW_U))
      proposed2DW_intblockslml <- sapply(int2_blocks, function(k) DW_blocklml(node[2], proposedDAGs_list[[k]], proposed_fa2.tXX[[k]], n_vec[k], DW_a, DW_U))
    }
    obs2_blocks <- setdiff(1:K, int2_blocks)
    n_obs2 <- sum(n_vec[obs2_blocks])

    current_fa2.tXX.obs <- Reduce("+", current_fa2.tXX[obs2_blocks])
    proposed_fa2.tXX.obs <- Reduce("+", proposed_fa2.tXX[obs2_blocks])
    current2DW_obsblockslml <- DW_blocklml(node[2], currentDAG, current_fa2.tXX.obs, n_obs2, DW_a, DW_U)
    proposed2DW_obsblockslml <- DW_blocklml(node[2], proposedDAG, proposed_fa2.tXX.obs, n_obs2, DW_a, DW_U)

    ## Get final values for current and proposed lml
    current_lml <- sum(c(current1DW_obsblockslml, current1DW_intblockslml, current2DW_obsblockslml, current2DW_intblockslml))
    proposed_lml <- sum(c(proposed1DW_obsblockslml, proposed1DW_intblockslml, proposed2DW_obsblockslml, proposed2DW_intblockslml))
  }

  acp.ratio <- min(0, proposed_lml - current_lml + logDAGprior.ratio +
                     logproposal.ratio)
  is.accepted <- log(stats::runif(1)) < acp.ratio

  return(is.accepted)
}
