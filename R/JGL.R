# EBIC
Calculate_ebic <- function(Theta, data_list, n, p, gamma) {
  
  K <- length(data_list)
  bic_values <- numeric(K)

  for (k in 1:K) {
   
    S <- cov(data_list[[k]])
    
    theta <- Theta[[k]]
    df <- sum(abs(theta[upper.tri(theta, diag = TRUE)]) > 1e-6)
    # df <- sum(abs(theta) > 1e-6)
    # 
    log_likelihood <- log(det(theta)) - sum(diag(S %*% theta))
    # BIC
    bic_values[k] <- -n[k] * log_likelihood + df * (log(n[k]) + gamma*log(p))
  }

  # return BIC
  return(sum(bic_values))
}

##### JEM 2011 #####
JointGrpah_BIC <- function(X, gamma){
  
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x))
  
  X_list <- as.matrix(do.call(rbind, X))
  Y <- rep(1:length(X), each = nrow(X[[1]]))
  
  lambda <- exp(seq(-2, 2, length.out = 5))
  bic_results <- rep(1e10, length(lambda))

  for (i in seq_along(lambda)) {
      # fit
      fit <- jointgraph_train(X_list, Y, lambda_value = lambda[i], 
                              adaptive_weight=array(1, c(length(unique(Y)), ncol(X_list), ncol(X_list))))
      res <- lapply(1:length(unique(Y)), function(i) fit$OMEGA[i,,])
      # BIC
      bic_results[i] <- Calculate_ebic(res, X, n, p, gamma)
  }

  # best lambda1 & lambda2
  optimal_indices <- which(bic_results == min(bic_results, na.rm = TRUE), arr.ind = TRUE)[1]
  optimal_lambda <- lambda[optimal_indices]

  # fit with optimal paras
  best_fit <- jointgraph_train(X_list, Y, lambda_value = optimal_lambda, adaptive_weight=array(1, c(length(unique(Y)), ncol(X_list), ncol(X_list))))
  return(lapply(1:length(unique(Y)), function(i) best_fit$OMEGA[i,,]))

}

##### GGL 2014 #####
JGL_EBIC <- function(X, penalty = "group", tol = 1e-3, maxiter = 300, gamma = 0) {
  
  # grid search lambda1 & lambda2
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x))
  
  lambda1_values <- exp(seq(log(1e-4), log(1e-1), length.out = 5))
  lambda2_values <- exp(seq(log(1e-4), log(1e-1), length.out = 5))
  # lambda1_values <- exp(seq(log(ifelse (nrow(X[[1]])>=ncol(X[[1]]), 1e-6, 1e-4)), log(0.5), length.out = 5))
  # lambda2_values <- exp(seq(log(ifelse (nrow(X[[1]])>=ncol(X[[1]]), 1e-6, 1e-4)), log(0.5), length.out = 10))
  bic_results <- matrix(NA, nrow = length(lambda1_values), ncol = length(lambda2_values))

  for (i in seq_along(lambda1_values)) {
    
    for (j in seq_along(lambda2_values)) {
      
      lambda1 <- lambda1_values[i]
      lambda2 <- lambda2_values[j]

      # fit
      fit <- JGL(Y = X, penalty = penalty, lambda1 = lambda1, lambda2 = lambda2,
                 return.whole.theta = T, tol = tol, maxiter = maxiter, truncate = 1e-5)

      # BIC
      bic_results[i, j] <- Calculate_ebic(fit$theta, X, n, p, gamma)
    }
  }

  # best lambda1 & lambda2
  optimal_indices <- which(bic_results == min(bic_results, na.rm = TRUE), arr.ind = TRUE)[1, ]
  optimal_lambda1 <- lambda1_values[optimal_indices[1]]
  optimal_lambda2 <- lambda2_values[optimal_indices[2]]

  # fit with optimal paras
  best_fit <- JGL(Y = X, penalty = "group", lambda1 = optimal_lambda1, lambda2 = optimal_lambda2,  return.whole.theta = T, tol = 1e-3)
  return( list(res=best_fit, lambda1= lambda1, lambda2=lambda2 ) )
}



jewel_BIC <- function(X, gamma = 0) {
  
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x))
  
  lambda1_values <- exp(seq(log(1e-3), log(1), length.out = 5))
  lambda2_values <- exp(seq(log(1e-3), log(1), length.out = 5))
  bic_results <- matrix(1e10, nrow = length(lambda1_values), ncol = length(lambda2_values))

  for (i in seq_along(lambda1_values)) {
    
    for (j in seq_along(lambda2_values)) {
      
      lambda1 <- lambda1_values[i]
      lambda2 <- lambda2_values[j]
      # fit
      fit <- fasjem(X, lambda = lambda1, epsilon = lambda2)

      # BIC
      bic_results[i, j] <- Calculate_ebic(fit, X, n, p, gamma)
    }
  }

  # best lambda1 & lambda2
  optimal_indices <- which(bic_results == min(bic_results, na.rm = TRUE), arr.ind = TRUE)
  optimal_lambda1 <- lambda1_values[optimal_indices[1]]
  optimal_lambda2 <- lambda2_values[optimal_indices[2]]
 
  # fit with optimal paras
  best_fit <- jewel(X, lambda1 = optimal_lambda1, lambda2 = optimal_lambda2, verbose = F)
  return(best_fit)
}


##### Fasjem 2017 #####
library(fasjem)

fasjem_EBIC <- function(X, lambda_grid = 10^( seq(-3, 0, length.out = 5) ),
                        epsilon_grid = 10^( seq(-2, 1, length.out = 5) ),
                        gamma=1 ) {
  
  bic_results <- matrix(NA, nrow = length(lambda_grid), ncol = length(epsilon_grid))
  best_lambda <- NULL
  best_epsilon <- NULL
  best_theta_list <- NULL
  p <- ncol(X[[1]])
  K <- length(X)
  n <- sapply(X, function(x) nrow(x))

  for (i in seq_along(lambda_grid)) {
    
    for (j in seq_along(epsilon_grid)) {
      
      lambda <- lambda_grid[i]
      epsilon <- epsilon_grid[j]
      result <- fasjem(X, lambda = lambda, epsilon = epsilon,
                       gamma= lambda, rho= 1)
      bic_results[i, j] <- Calculate_ebic(result, X, n, p, gamma)
    }
  }
  optimal_indices <- which(bic_results == min(bic_results, na.rm = TRUE), arr.ind = TRUE)

  if( all(is.na(optimal_indices)) ){
    
    optimal_lambda = lambda_grid[ length(lambda_grid)-1]
    optimal_epsilon = epsilon_grid[ length(epsilon_grid )-1]
  }else{
    
    optimal_lambda <- lambda_grid[optimal_indices[1]]
    optimal_epsilon <- epsilon_grid[optimal_indices[2]]
  }
  
  result <- fasjem(X, lambda = optimal_lambda, gamma= optimal_lambda,
                   epsilon = optimal_epsilon, rho= 1)
  return(result)
}
