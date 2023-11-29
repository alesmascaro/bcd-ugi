learn_DTP <- function(S, burn, data, initialDAG,
                      DW_a, DW_U, alpha, beta, a, b, a_w, b_w,
                      fast = FALSE) {
  n.iter <- burn + S
  q <- ncol(data[[1]])

  tXX_list <- lapply(data, function(X) crossprod(X))
  n_vec <- sapply(data, nrow)
  K <- length(data)

  ## Initialize arrays and current values
  DAGs <- array(0, dim = c(q,q,n.iter))
  TARGETs <- array(0, dim = c(K,q,n.iter))
  PARENTs <- vector("list", K)
  for (i in (1:K)) {
    PARENTs[[i]] <- array(0, dim = c(q,q,n.iter))
  }

  prop.opcard_mat <- matrix(0, K, n.iter)

  # currentDAG <- matrix(0, ncol = q, nrow = q)
  currentDAG <- initialDAG
  currentTARGETS <- matrix(0, ncol = q, nrow = K)
  currentPARENTS <- vector("list", K)
  for (i in 1:K) {
    currentPARENTS[[i]] <- matrix(0, ncol = q, nrow = q)
  }

  pb <- utils::txtProgressBar(min = 2, max = n.iter, style = 3)

  for (i in 1:n.iter) {
    perm <- sample(1:K, K)

    for (k in perm) {
      if (k == 1) {
        propDAG <- propose_DAG(currentDAG, currentTARGETS, currentPARENTS, fast)
        prop.opcard_mat[k,i] <- propDAG$op.card
        is.DAG.accepted <- acceptreject_DAG(tXX_list, n_vec, currentTARGETS, currentPARENTS,
                                            currentDAG, propDAG$proposedDAG,
                                            node = propDAG$op.node, op.type = propDAG$op.type,
                                            op.card = propDAG$op.card,
                                            DW_a, DW_U, a_w, b_w,
                                            fast)

        # browser()

        if (is.DAG.accepted) {
          currentDAG <- propDAG$proposedDAG
        }

        DAGs[,,i] <- currentDAG
      } else {
        propkth_TP <- propose_TP(k,currentDAG, currentTARGETS, currentPARENTS, fast)
        prop.opcard_mat[k,i] <- propkth_TP$op.card
        is.kthTP.accepted <- acceptreject_kthTP(k, currentDAG, currentTARGETS, currentPARENTS,
                                                propkth_TP$proposedTARGETS, propkth_TP$proposed_kthPARENTS,
                                                tXX_list, n_vec, DW_a, DW_U, alpha, beta, a, b,
                                                node = propkth_TP$op.node, op.card = propkth_TP$op.card,
                                                fast)

        if (is.kthTP.accepted) {
          currentTARGETS[k,] <- propkth_TP$proposedTARGET[k,]
          currentPARENTS[[k]] <- propkth_TP$proposed_kthPARENTS
        }
        PARENTs[[k]][,,i] <- currentPARENTS[[k]]
      }
    }
    TARGETs[,,i] <- currentTARGETS
    # browser()
    utils::setTxtProgressBar(pb, i)
    close(pb)
  }

  ## Remove burnin
  DAGs <- DAGs[,,(burn+1):n.iter]
  TARGETs <- TARGETs[,,(burn+1):n.iter]
  for (k in 1:K) {
    PARENTs[[k]] <- PARENTs[[k]][,,(burn+1):n.iter]
  }
  prop.opcard_mat <- prop.opcard_mat[,(burn+1):n.iter]

  out <- list(DAGs = DAGs, TARGETs = TARGETs, PARENTs = PARENTs)
  return(out)
}
