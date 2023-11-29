acceptreject_kthTP <- function(k, currentDAG, currentTARGETS, currentPARENTS,
                               proposedTARGETS, proposed_kthPARENTS,
                               tXX_list, n_vec, DW_a, DW_U, alpha, beta, a, b,
                               node, op.card, fast) {

  # modifiedDAG <- modify_parents(currentDAG, currentTARGETS[k,], currentPARENTS[[k]])
  # if (gRbase::is.DAG(modifiedDAG) == F) {
  #   return(FALSE)
  # }

  q <- ncol(currentDAG)
  K <- length(tXX_list)

  ### Log-proposal ratio
  if (fast == TRUE) {
    logproposal.ratio <- log(1)
  } else if (fast == FALSE) {
    logproposal.ratio <- log(op.card/get_TPcard(k, currentDAG, proposedTARGETS, proposed_kthPARENTS))
  }

  ### Marginal prior ratio for TARGETS
  prop_margTARGETSprior <- lbeta(a + sum(proposedTARGETS[k,]), b + q - sum(proposedTARGETS[k,])) - lbeta(a,b)
  current_margTARGETprior <- lbeta(a + sum(currentTARGETS[k,]), b + q - sum(currentTARGETS[k,])) - lbeta(a,b)

  logTARGETSprior.ratio <- prop_margTARGETSprior - current_margTARGETprior

  ### Marginal prior probability for PARENTS
  prop_margPARENTSprior <- sum(lbeta(alpha + colSums(proposed_kthPARENTS),
                                     beta + q - colSums(proposed_kthPARENTS)) -
                                 lbeta(alpha, beta))
  current_margPARENTSprior <- sum(lbeta(alpha + colSums(currentPARENTS[[k]]),
                                        beta + q - colSums(currentPARENTS[[k]])) -
                                    lbeta(alpha, beta))
  logPARENTSprior.ratio <- prop_margPARENTSprior - current_margPARENTSprior

  ### Marginal likelihood ratio
  currentDAGs <- get_modifiedDAGs(currentDAG, currentTARGETS, currentPARENTS)
  proposedPARENTS <- currentPARENTS
  proposedPARENTS[[k]] <- proposed_kthPARENTS
  proposedDAGs <- get_modifiedDAGs(currentDAG, proposedTARGETS, proposedPARENTS)

  if (length(node) == 0) {
    return(TRUE)
  } else {
    current_nodeslml <- vector("double", length(node))
    proposed_nodeslml <- vector("double", length(node))

    for (j in 1:length(node)) {
      ### Define families of changed node
      currentfamilies <- lapply(currentDAGs, function(k) fa(node[j], k))
      proposedfamilies <- lapply(proposedDAGs, function(k) fa(node[j], k))
      ### Define subsets of tXX based on the changed node and different DAGs (families) for each context
      current_fa.tXX <- lapply(1:K, function(k) tXX_list[[k]][currentfamilies[[k]], currentfamilies[[k]]])
      proposed_fa.tXX <- lapply(1:K, function(k) tXX_list[[k]][proposedfamilies[[k]], proposedfamilies[[k]]])
      ### Initialize lml
      currentDW_intblockslml <- 0
      proposedDW_intblockslml <- 0
      ### Compute lml for interventional blocks (if present)
      current_intblocks <- which(currentTARGETS[,node[j]] != 0)
      proposed_intblocks <- which(proposedTARGETS[,node[j]] != 0)
      current_arethereint <- sum(current_intblocks) != 0
      proposed_arethereint <- sum(proposed_intblocks) != 0
      if (current_arethereint) {
        currentDW_intblockslml <- sapply(current_intblocks, function(k) DW_blocklml(node[j], currentDAGs[[k]], current_fa.tXX[[k]], n_vec[k], DW_a, DW_U))
      }
      if (proposed_arethereint) {
        proposedDW_intblockslml <- sapply(proposed_intblocks, function(k) DW_blocklml(node[j], proposedDAGs[[k]], proposed_fa.tXX[[k]], n_vec[k], DW_a, DW_U))
      }
      ### Compute lml for observational blocks (always present)
      current_obsblocks <- setdiff(1:K, current_intblocks)
      proposed_obsblocks <- setdiff(1:K, proposed_intblocks)
      current_nobs <- sum(n_vec[current_obsblocks])
      proposed_nobs <- sum(n_vec[proposed_obsblocks])
      #### Compute the corresponding statitic for the observational part from list of tXX
      # browser()
      current_fa.tXX.obs <- Reduce("+", current_fa.tXX[current_obsblocks])
      proposed_fa.tXX.obs <- Reduce("+", proposed_fa.tXX[proposed_obsblocks])

      currentDW_obsblockslml <- DW_blocklml(node[j], currentDAG, current_fa.tXX.obs, current_nobs, DW_a, DW_U)
      proposedDW_obsblockslml <- DW_blocklml(node[j], currentDAG, proposed_fa.tXX.obs, proposed_nobs, DW_a, DW_U)

      ## Compute total current and proposed lml
      current_nodeslml[j] <- sum(c(as.numeric(currentDW_obsblockslml), as.numeric(currentDW_intblockslml)))
      proposed_nodeslml[j] <- sum(c(as.numeric(proposedDW_obsblockslml), as.numeric(proposedDW_intblockslml)))
    }

    current_lml <- sum(current_nodeslml)
    proposed_lml <- sum(proposed_nodeslml)
  }

  acp.ratio <- min(0, proposed_lml - current_lml + logTARGETSprior.ratio + logPARENTSprior.ratio +
                     logproposal.ratio)
  is.accepted <- log(stats::runif(1)) < acp.ratio
  return(is.accepted)
}
