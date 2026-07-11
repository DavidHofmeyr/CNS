if(system.file(package = "FNN") == "") install.packages("FNN")
if(system.file(package = "Matrix") == "") install.packages("Matrix")
if(system.file(package = "matrixStats") == "") install.packages("matrixStats")
if(system.file(package = "igraph") == "") install.packages("igraph")
library(FNN)
library(Matrix)
library(matrixStats)
library(igraph)


### Function CNS implements the Custering by Non-parametric
### Smoothing algorithm, initially described in Hofmeyr, D.P. (2025) currently on Arxiv,
### and slightly modified based on comments from reviewers.
## Arguments:
## X = data matrix with observations row-wise
## kmax = maximum number of clusters
## nn = options for number of nearest neighbours
## lams = options for weight given to initial solution
## distance = metric for nearest neighbours; one of "Euclidean" and "cosine"
## iters = number of iterations for approximate matrix inversion, with larger
##          values being better approximation. Default is 50, following what is used in the paper


### Output:
## list with fields
## $probabilities = matrix of cluster membership probabilities
## $cluster = cluster allocations
## $vals = array storing the values of the model selection criterion for all settings considered

### All default values for settings are as they are in the revised paper, currently under submission


CNS <- function(X, kmin = 2, kmax = 50, nn = floor(log(nrow(X)))*1:4, lams = seq(0.005, 0.05, length = 10), distance = 'Euclidean', iters = 50, init = "graph", oq = TRUE){
  d <- ncol(X)
  n <- nrow(X)

  if(distance=='cosine'){
    X <- X - matrix(colMeans(X), n, d, byrow = TRUE)
    X <- X/sqrt(rowSums(X^2))
  }


  nns0 <- get.knn(X, max(nn)-1)
  nns0$nn.index <- cbind(1:n, nns0$nn.index)
  if(distance=="cosine") nns0$nn.dist <- nns0$nn.dist^2
  nns0$nn.dist <- nns0$nn.dist + 1e-10
  nns0$nn.dist <- cbind(0, nns0$nn.dist)
  nns <- nns0$nn.index


  ### vals stores the values of the model selection criterion used to select a final solution
  vals <- array(0, dim = c(kmax, length(nn), length(lams)))

  ### Qs stores the information needed to produce the cluster probability vectors for all solutions
  Qs <- array(dim = c(n, kmax, length(nn), length(lams)))

  if(iters==Inf){
    for(ni in 1:length(nn)){
      L <- sparseMatrix(i = rep(1:n, nn[ni]), j = nns[,1:nn[ni]], dims = c(n,n))
      G <- graph_from_adjacency_matrix(L)

      ### Find "local maxima" in the neighbour relations
      cs <- colSums(L)

      uix <- numeric(kmax)
      uix[1] <- which.max(cs)
      if(init=="graph"){
        ds <- apply(rbind(distances(G, v = uix[1], mode = "in"), distances(G, v = uix[1], mode = "out")), 2, min)
        ds[!is.finite(ds)] <- 1e100
      }
      else{
        if(distance=="Euclidean") ds <- sqrt(colSums((t(X)-X[uix[1],])^2))
        else if(distance=="cosine") ds <- colSums((t(X)-X[uix[1],])^2)
      }
      for(k in 2:kmax){
        uix[k] <- which.max(ds*cs)
        if(init=="graph"){
          dnew <- apply(rbind(distances(G, v = uix[k], mode = "in"), distances(G, v = uix[k], mode = "out")), 2, min)
          dnew[!is.finite(dnew)] <- 1e100
        }
        else{
          if(distance=="Euclidean") dnew <- sqrt(colSums((t(X)-X[uix[k],])^2))
          else if(distance=="cosine") dnew <- colSums((t(X)-X[uix[k],])^2)
        }
        ds <- ds*(ds<dnew) + dnew*(dnew<=ds)
      }


      ### pix stores the uix indices in matrix form
      pix <- sparseMatrix(i=uix, j=1:length(uix), dims = c(n, length(uix)))

      for(li in 1:length(lams)){
        lam <- lams[li]

        ### Compute information needed for solution with current settings of nn and lambda
        I <- sparseMatrix(1:n,1:n) - (1-lam)*L/nn[ni]

        if(oq) Qs[,1:min(kmax, length(uix)),ni,li] <- matrix(oQ(lam*solve(I, pix, sparse = TRUE), min(kmax, length(uix))))
        else Qs[,1:min(kmax, length(uix)),ni,li] <- matrix(lam*solve(I, pix, sparse = TRUE))


        ### store the model selection criterion values for all number of clusters up to kmax for these settings
        vals[,ni,li] <- c(rep(0, kmin-1), sapply(kmin:kmax, function(k){
          ref <- (1-lam)*(1/nn[ni] + 1/n - 2/sqrt(n*nn[ni]))
          (mean(rowMaxs(1/k+Qs[,1:k,ni,li]-rowSums(Qs[,1:k,ni,li])/k)) - (n-k+k^2)/n/k)/ref
        }))
      }
    }
  }
  else{
    for(ni in 1:length(nn)){
      L <- sparseMatrix(i = rep(1:n, nn[ni]), j = nns[,1:nn[ni]], dims = c(n,n))
      G <- graph_from_adjacency_matrix(L)
      ### Find "local maxima" in the neighbour relations
      cs <- colSums(L)

      uix <- numeric(kmax)
      uix[1] <- which.max(cs)
      if(init=="graph"){
        ds <- apply(rbind(distances(G, v = uix[1], mode = "in"), distances(G, v = uix[1], mode = "out")), 2, min)
        ds[!is.finite(ds)] <- 1e100
      }
      else{
        if(distance=="Euclidean") ds <- sqrt(colSums((t(X)-X[uix[1],])^2))
        else if(distance=="cosine") ds <- colSums((t(X)-X[uix[1],])^2)
      }
      for(k in 2:kmax){
        uix[k] <- which.max(ds*cs)
        if(init=="graph"){
          dnew <- apply(rbind(distances(G, v = uix[k], mode = "in"), distances(G, v = uix[k], mode = "out")), 2, min)
          dnew[!is.finite(dnew)] <- 1e100
        }
        else{
          if(distance=="Euclidean") dnew <- sqrt(colSums((t(X)-X[uix[k],])^2))
          else if(distance=="cosine") dnew <- colSums((t(X)-X[uix[k],])^2)
        }
        ds <- ds*(ds<dnew) + dnew*(dnew<=ds)
      }


      ### pix stores the uix indices in matrix form
      pix <- sparseMatrix(i=uix, j=1:length(uix), dims = c(n, length(uix)))

      qs <- appr_solv(L/nn[ni], pix, 1-lams, iter = iters)

      for(li in 1:length(lams)){

        if(oq) Qs[,1:min(kmax, length(uix)),ni,li] <- oQ(lams[li]*qs[,,li])
        else Qs[,1:min(kmax, length(uix)),ni,li] <- lams[li]*qs[,,li]

        ### store the model selection criterion values for all number of clusters up to kmax for these settings
        vals[,ni,li] <- c(rep(0, kmin-1), sapply(kmin:kmax, function(k){
          ref <- (1-lams[li])*(1/nn[ni] + 1/n - 2/sqrt(n*nn[ni]))
          (mean(rowMaxs(1/k+Qs[,1:k,ni,li]-rowSums(Qs[,1:k,ni,li])/k)) - (n-k+k^2)/n/k)/ref
        }))
      }
    }
  }


  ### if some numerical issues occurred, ensure maximum can be found without issue
  vals[is.na(vals)] <- -1


  ### find the best solution and compute the matrix of cluster probabilities
  parms <- which(vals==max(vals), arr.ind = TRUE)[1,]
  if(parms[1] > 1) prob <- 1/parms[1]+Qs[,1:parms[1],parms[2],parms[3]]-rowSums(Qs[,1:parms[1],parms[2],parms[3]])/parms[1]
  else prob <- matrix(1, n, 1)

  ### return
  list(probabilities = prob, clusters = apply(prob, 1, which.max), vals = vals, nn = nn, lams = lams, X = X)

}


oQ <- function(Q, q = ncol(Q)){
  cs <- colSums(Q)
  sim <- (t(Q)%*%Q)
  diag(sim) <- Inf
  ord <- numeric(q)
  ord[1] <- which.max(cs)
  ds <- sim[ord[1],]
  ds[ord[1]] <- Inf
  for(k in 2:q){
    ord[k] <- which.min(ds/cs^2)
    ds <- ds*(ds > sim[ord[k],]) + sim[ord[k],]*(ds <= sim[ord[k],])
  }
  Q[,ord]
}

appr_solv <- function(W, F0, lams, iter = 100){
  out <- array(0, dim = c(nrow(F0), ncol(F0), length(lams)))
  add <- F0
  for(it in 0:iter){
    for(j in 1:length(lams)) out[,,j] <- out[,,j] + as.matrix(add)*lams[j]^it
    add <- W%*%add
  }
  for(j in 1:length(lams)) out[,,j] <- out[,,j] + as.matrix(add)*(lams[j]^(it+1))/(1-lams[j])
  out
}


