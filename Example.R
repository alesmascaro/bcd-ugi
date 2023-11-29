q <- 10
w <- 0.2
DAG1 <- matrix(0,q,q)
set.seed(1)
DAG1[lower.tri(DAG1)] <- sample(c(0,1), q*(q-1)/2, T, prob = c(1-w, w))
gRbase::plot(as(DAG1, "graphNEL"))
set.seed(1)
L1 <- runif(q^2, 0.3, 1)*sample(c(-1,1), q^2, T)*DAG1
L1 <- L1 + diag(1,q)
D1 <- diag(1,q)
Sigma1 <- crossprod(solve(L1))

DAG2 <- DAG1
DAG2[2,4] <- DAG2[6,4] <- 1
gRbase::plot(as(DAG2, "graphNEL"))
L2 <- L1
set.seed(2)
L2[c(2,6,9),4] <- runif(3, 0.3, 1)*sample(c(-1,1), 3, T)
D2 <- D1
Sigma2 <- crossprod(solve(L2))

n <- 500
set.seed(1)
X1 <- mvtnorm::rmvnorm(n, sigma = Sigma1)
X2 <- mvtnorm::rmvnorm(n, sigma = Sigma2)

files.sources <- list.files(path=paste0(getwd(), "/Auxiliary functions") ,full.names = TRUE)
sapply(files.sources, source)
source(paste0(getwd(), "/learn_DTP.R"))

out <- learn_DTP(5000, 2000, list(X1, X2), matrix(0,q,q), q-1+1/2/n, diag(1,q), 1, 1, 1, 1, 1, 1,fast = TRUE)

round(apply(out$DAGs, c(1,2), mean))
DAG1
round(apply(out$TARGETs, c(1,2), mean))
round(apply(out$PARENTs[[2]], c(1,2), mean))
