
library(snowfall)
library(parallel)
library(Matrix) 

# sfInit(parallel = TRUE, cpus = detectCores() - 6)
# sfLibrary(ADSIHT)
# sfLibrary(Matrix)

MT_DSIHT <- function(num, X, n, K, p, ic.coef, ic.scale,
                     coef1, coef2, kappa, eta, center , scale) {
  
  y_list <- lapply(X, function(z) z[, num])
  x_list <- lapply(X, function(z) z[, -num])
  
  fit <- ADSIHT.ML(x_list, y_list, method = "ols", 
                   ic.coef = ic.coef, ic.scale = ic.scale, 
                   coef1 = coef1, coef2 = coef2, 
                   kappa = kappa, eta=eta, L = 15 ,
                   center = center, scale= scale )
  beta <- fit[["beta"]][[which.min(fit[["ic"]])]]
  
  if( center == TRUE){
    temp <- unlist(y_list) - as.matrix(bdiag(x_list)) %*% beta - rep(fit[["intercept"]][[which.min(fit[["ic"]])]], times = n)
  }else{
    temp <- unlist(y_list) - as.matrix(bdiag(x_list)) %*% beta 
  }
  
  sigma <- rep(0, K)
  
  for (i in 1:K) {
    
    if (i == 1) {
      sigma[i] <- sum(temp[1:n[i]]^2) / n[i]
    } else {
      sigma[i] <- sum(temp[(sum(n[1:(i-1)]) + 1):sum(n[1:i])]^2) / n[i]
    }
  }
  omega <- lapply(1:K, function(i) {
    
    temp = -1/sigma[i]*beta[((i-1)*(p-1)+1):(i*(p-1))]
    return(append(temp, 1/sigma[i], after = num-1))
    
    })
  
  return(omega)
}


#### Main function of MIGHT ####
JGML <- function(X, ic.coef = 1, ic.scale = 2, coef1 = 1, coef2 = 0.1, 
                 kappa = 0.9, eta=0.9, center , scale){
  
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x)) #sample size of each task
  
  res <- sfLapply(1:p, function(num) MT_DSIHT(num, X = X, n = n, K = K, p = p, 
                                              ic.coef = ic.coef, ic.scale = ic.scale,
                                              coef1 = coef1, coef2 = coef2, 
                                              kappa = kappa, eta=eta, 
                                              center = center, scale= scale))
  
  res_Omega <- lapply(1:K, function(i) matrix(unlist(purrr::map(res, ~.x[[i]])), ncol = p))
  res_Omega <- lapply(res_Omega, function(A) {
    n <- nrow(A)  
    for (i in 2:n) {
      
      for (j in 1:(i-1)) {
        
        temp = ifelse(abs(A[i, j]) >= abs(A[j, i]), A[j, i], A[i, j]) 
        A[i, j] = A[j, i] = temp 
      }}
    
    return(A)
  })
  
  return(res_Omega)
}

# with a given s0
MT_DSIHT_s0 <- function(num, X, n, K, p, s0, ic.coef, ic.scale, coef1, coef2, kappa) {
  
  y_list <- lapply(X, function(z) z[, num])
  x_list <- lapply(X, function(z) z[, -num])
  
  fit <- ADSIHT.ML(x_list, y_list, method = "ols", ic.coef = ic.coef, ic.scale = ic.scale, coef1 = coef1, coef2 = coef2, kappa = kappa, s0 = s0)
  
  beta <- fit[["beta"]][[which.min(fit[["ic"]])]]
  temp <- unlist(y_list) - as.matrix(bdiag(x_list)) %*% beta-rep(fit[["intercept"]][[which.min(fit[["ic"]])]], times = n)

  sigma <- rep(0, K)
  for (i in 1:K) {
    
    if (i == 1) {
      
      sigma[i] <- sum(temp[1:n[i]]^2) / n[i]
    } else {
      
      sigma[i] <- sum(temp[(sum(n[1:(i-1)]) + 1):sum(n[1:i])]^2) / n[i]
    }
  }
  
  omega <- lapply(1:K, function(i) {
    
    temp = -1/sigma[i]*beta[((i-1)*(p-1)+1):(i*(p-1))]
    
    return(append(temp, 1/sigma[i], after = num-1))
  })
  
  return(omega)
}

JGML_s0 <- function(X, s0 = c(4, 5), ic.coef = 2, ic.scale = 1, coef1 = 1, coef2 = 0.1, kappa = 0.9) {
  
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x))
  
  res <- sfLapply(1:p, function(num) MT_DSIHT_s0(num, X = X, n = n, K = K, p = p, s0 = s0, 
                                                 ic.coef = ic.coef, ic.scale = ic.scale,
                                                 coef1 = coef1, coef2 = coef2, kappa = kappa))
  res_Omega <- lapply(1:K, function(i) matrix(unlist(purrr::map(res, ~.x[[i]])), ncol = p))
  
  res_Omega <- lapply(res_Omega, function(A) {
    #A = matrix(-5:10,4,4)
    n <- nrow(A)   
    for (i in 2:n) {
      
      for (j in 1:(i-1)) {
        
          temp = ifelse(abs(A[i, j]) >= abs(A[j, i]), A[j, i], A[i, j]) 
          A[i, j] = A[j, i] = temp 
      }}
    
    return(A)
    })
  
  return(res_Omega)
}
