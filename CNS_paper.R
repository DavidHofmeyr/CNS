if(system.file(package = "FNN") == "") install.packages("FNN")
if(system.file(package = "Matrix") == "") install.packages("Matrix")

### Function CNS implements the Custering by Nonparametric
### Smoothing algorithm in Hofmeyr, D.P. (2025), Arxiv:2503.09134
## Arguments:
## X = data matrix with observations row-wise
## kmax = maximum number of clusters
## nn = options for number of nearest neighbours
## lams = options for weight given to initial solution
## k0 = maximum number of potential points to use as "informative"
##      members of a cluster in the initial solution. This is purely
##      in place for computational reasons.
### Output:
## list with fields
## $probabilities = matrix of cluster membership probabilities
## $cluster = cluster allocations
## $vals = array storing the values of the model selection criterion for all settings considered

### All default values for settings are as they were used in the paper


CNS <- function(X, kmax = 30, nn = c(5, 7, 9, 11, 13, 15), lams = c(.01, .02), k0 = 300){
  d <- ncol(X)
  n <- nrow(X)
  
  ### Find neighbours up to largest number over which to fit models
  nns0 <- get.knn(X, max(nn))
  nns <- nns0$nn.index
  
  ### vals stores the values of the model selection criterion used to select a final solution
  vals <- array(0, dim = c(kmax, length(nn), length(lams)))
  
  ### Qs stores the information needed to produce the cluster probability vectors for all solutions
  Qs <- array(dim = c(n, kmax, length(nn), length(lams)))
  
  ### Loop over different options
  for(ni in 1:length(nn)){
    L <- sparseMatrix(i = rep(1:n, nn[ni]), j = nns[,1:nn[ni]], dims = c(n,n))
    
    ### Find "local maxima" in the neighbour relations
    cs <- colSums(L)
    
    mds <- which(cs >= apply(matrix(cs[nns[,1:nn[ni]]], n, nn[ni]), 1, max))
    
    ### uix stores the indices for the initialisation
    if(length(mds)==1) uix <- rep(mds, kmax)
    else if(length(mds) > k0){
      uix <- mds[order(get.knn(X[mds,], 1)$nn.dist*cs[mds], decreasing = TRUE)[1:k0]]
    }
    else uix <- mds
    
    ### pix stores the uix indices in matrix form
    pix <- sparseMatrix(i=uix, j=1:length(uix), dims = c(n, length(uix)))
    
    for(li in 1:length(lams)){
      lam <- lams[li]
      
      ### Compute information needed for solution with current settings of nn and lambda
      I <- sparseMatrix(1:n,1:n) - (1-lam)*L/nn[ni]
      
      Qs[,1:min(kmax, length(uix)),ni,li] <- matrix(oQ(lam*solve(I, pix, sparse = TRUE), min(kmax, length(uix))))

      ### kopt is the value for K which maximises the normalisation term for model selection
      kopt <- sqrt((n-(n-lam)*nn[ni]/(nn[ni]+1-lam))/(1-lam))
      
      ### store the model selection criterion values for all number of clusters up to kmax for these settings
      vals[,ni,li] <- c(0, sapply(2:kmax, function(k){
        ref <- (n-lam)*lam*(nn[ni]/kopt+1-lam)/n/(nn[ni]-(1-lam)*(nn[ni]-lam))+kopt/n*lam - (n-kopt+kopt^2)/n/kopt
        (mean(apply(1/k+Qs[,1:k,ni,li]-rowSums(Qs[,1:k,ni,li])/k, 1, max)) - (n-k+k^2)/n/k)/ref
      }))
      
    }
  }
  
  ### if some numerical issues occurred, ensure maximum can be found without issue
  vals[is.na(vals)] <- -1
  
  
  ### find the best solution and compute the matrix of cluster probabilities
  parms <- which(vals==max(vals), arr.ind = T)[1,]
  if(parms[1] > 1) prob <- 1/parms[1]+Qs[,1:parms[1],parms[2],parms[3]]-rowSums(Qs[,1:parms[1],parms[2],parms[3]])/parms[1]
  else prob <- matrix(1, n, 1)
  
  ### return
  list(probabilities = prob, clusters = apply(prob, 1, which.max), vals = vals)
  
}
