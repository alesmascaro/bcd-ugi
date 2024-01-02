  ## Source necessary files and packages
files.sources <- list.files(path=paste0(getwd(), "/Auxiliary functions") ,full.names = TRUE)
sapply(files.sources, source)
source(paste0(getwd(), "/learn_DTP.R"))

# install.packages("gRbase")
# install.packages("BiocManager")
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")

# Generate DAG ------------------------------------------------------------

q <- 10
w <- 0.2
DAG1 <- matrix(0,q,q)
set.seed(1)
DAG1[lower.tri(DAG1)] <- sample(c(0,1), q*(q-1)/2, T, prob = c(1-w, w))

Rgraphviz::plot(as_graphNEL(DAG1))

# Generate DAG and DAG parameters -----------------------------------------

set.seed(1)
L1 <- runif(q^2, 0.3, 1)*sample(c(-1,1), q^2, T)*DAG1
L1 <- L1 + diag(1,q)
D1 <- diag(1,q)
Sigma1 <- crossprod(solve(L1))

DAG2 <- DAG1
DAG2[2,4] <- DAG2[6,4] <- 1
Rgraphviz::plot(as_graphNEL(DAG2))

L2 <- L1
set.seed(2)
L2[c(2,6,9),4] <- runif(3, 0.3, 1)*sample(c(-1,1), 3, T)
D2 <- D1
Sigma2 <- crossprod(solve(L2))


# Generate data -----------------------------------------------------------

n <- 500
set.seed(1)
X1 <- mvtnorm::rmvnorm(n, sigma = Sigma1)
X2 <- mvtnorm::rmvnorm(n, sigma = Sigma2)
data <- list(X1, X2)

# Use learn_DTP -----------------------------------------------------------

  ## Specify prior hyperparams
    ### Parameter prior
DW_a <- q-1+1/(2*n)
DW_U <- diag(1,q)
    ### Induced parents prior
a_phi = b_phi = 1
    ### Targets prior
a_eta = b_eta = 1
    ### DAG prior
a_pi = b_pi = 1

  ## Run MCMC
S = 30000
burn = 20000
initial_DAG <- matrix(0,q,q)
set.seed(1)
out <- learn_DTP(S, burn, data, initial_DAG,
                 DW_a, DW_U, a_phi, b_phi, a_eta, b_eta, a_pi, b_pi, fast = FALSE)


# Output ------------------------------------------------------------------

round(apply(out$DAGs, c(1,2), mean), 2)
round(apply(out$TARGETs, c(1,2), mean), 1)
round(apply(out$PARENTs[[2]], c(1,2), mean), 1)

MPM_DAG <- round(apply(out$DAGs, c(1,2), mean))
par(mfrow = c(1,2))
Rgraphviz::plot(as_graphNEL(MPM_DAG))
Rgraphviz::plot(as_graphNEL(DAG1))

