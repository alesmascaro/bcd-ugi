DW_blocklml <- function(node, DAG, fa.tXX, n, a, U) {
  j <- node
  pa <- which(DAG[,node] != 0)
  fa <- c(j, pa)
  q <- ncol(DAG)

  pa_card <- length(pa)

  a_j <- (a+pa_card-q+1)


  if (pa_card == 0) {
    U_jj <- U[j,j]
    Upost_jj <- U_jj + fa.tXX

    prior.normcost <- -lgamma(a_j/2) + a_j/2*log(U_jj/2)
    post.normcost <- -lgamma(a_j/2 + n/2) + (a_j/2 + n/2)*log(Upost_jj/2)
  } else {
    Uprior_pa <- U[pa,pa]
    Uprior_fa <- U[fa,fa]
    Uprior_paj.j <- U[pa,j]
    if (all(Uprior_pa[!diag(pa_card)] == 0)) {
      Uprior_pa.inv <- 1/Uprior_pa
      Uprior_pa.det <- prod(diag(as.matrix(Uprior_pa)))
      Uprior_jj <- as.numeric(U[j,j])
    } else {
      Uprior_pa.inv <- chol2inv(chol(Uprior_pa))
      Uprior_pa.det <- det(as.matrix(Uprior_pa))
      Uprior_jj <- as.numeric(U[j,j] - t(Uprior_paj.j)%*%Uprior_pa.inv%*%Uprior_paj.j)
    }

    Upost_fa <- Uprior_fa + fa.tXX
    Upost_pa <- Upost_fa[-1, -1]
    Upost_paj.j <- Upost_fa[-1,1]
    Upost_jj <- as.numeric(Upost_fa[1,1] - t(Upost_paj.j)%*%chol2inv(chol(Upost_pa))%*%Upost_paj.j)

    prior.normcost <- -lgamma(a_j/2) + ((a_j)/2*(log(Uprior_jj))+(1/2)*log(Uprior_pa.det)) + (n/2)*log(2)
    post.normcost <- -lgamma(a_j/2 + n/2) + ((a_j+n)/2*(log(Upost_jj))+(1/2)*log(det(as.matrix(Upost_pa))))

  }

  blocklml <- -n/2*log(2*pi) + prior.normcost - post.normcost

  return(blocklml)
}
