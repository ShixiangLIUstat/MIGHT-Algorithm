library(glasso)
library(clime) 

# EBIC 
calculate_ebic <- function(theta, S, n, gamma = 0.5) {
  
  p <- ncol(S) #dimension
  log_likelihood <- n * (sum(diag(S %*% theta)) - log(det(theta)))
  non_zero <- sum(abs(theta[upper.tri(theta, diag = TRUE)]) > 1e-6) #non zero entry
  penalty <- non_zero * log(n) + gamma * non_zero * log(p)
  ebic <- log_likelihood + penalty
  
  return(ebic)
}

##### Glasso #####

# separate estimation by glasso
tune_glasso_single <- function(data, lambda_grid, gamma = 0.5) {
  
  n <- nrow(data)
  S <- cov(data)  
  best_ebic <- Inf
  best_lambda <- NULL
  best_theta <- NULL

  for (lambda in lambda_grid) {
    
    result <- glasso(S, rho = lambda)
    theta <- result$wi # inverse covariance matrix
    ebic <- calculate_ebic(theta, S, n, gamma)

    if (ebic < best_ebic) {
      
      best_ebic <- ebic
      best_lambda <- lambda
      best_theta <- theta
    }
  }

  return(list(best_lambda = best_lambda, best_theta = best_theta))
}

# multi glasso
tune_glasso_multiple <- function(data_list, lambda_grid, gamma = 0.5) {
  
  results <- list()

  for (i in seq_along(data_list)) {
    
    result <- tune_glasso_single(data_list[[i]], lambda_grid, gamma)
    results[[i]] <- result$best_theta
  }

  return(results)
}


##### Nodewise-HTP #####
tune_HTP_single <- function(X, ic.coef) {
  
  p <- ncol(X)
  n <- nrow(X)
  res <- lapply(1:p, function(i) {
    
      y <- X[, i]
      x <- X[, -i]
      fit <- AdaL0::bess3(x, y, threshold = 1000, ic.coef = ic.coef)
      beta <- fit$beta[, fit$best.size]
      intercept <- fit$intercept[fit$best.size]
      sigma <- sum((y - x%*%beta - intercept)^2)/n
      temp <- -1/(sigma+1e-2)*beta
      
      return(append(temp, 1/(sigma+1e-2), after = i-1))
  })
  
  return(matrix(unlist(res), ncol = p))
}

tune_HTP_multiple <- function(data_list, ic.coef = 1) {
  
  results <- list()

  for (k in seq_along(data_list)) {
    
    result <- tune_HTP_single(data_list[[k]], ic.coef = ic.coef)
    n <- nrow(result)  
    for (i in 2:n) {
      
      for (j in 1:(i-1)) {
        
        temp = ifelse(abs(result[i, j]) >= abs(result[j, i]), 
                      result[j, i], result[i, j]) 
        result[i, j] = result[j, i] = temp 
      }
    }
    
    results[[k]] <- result
  }

  return(results)
}


##### Compute entry-wise MCC #####

calculate_mcc <- function(mat1, mat2) {

  mat1_bin <- (abs(mat1) > 1e-6)
  mat2_bin <- (abs(mat2) > 1e-6)

  upper_tri_1 <- mat1_bin[upper.tri(mat1_bin)]
  upper_tri_2 <- mat2_bin[upper.tri(mat2_bin)]
   
  return( mccr(upper_tri_1, upper_tri_2 ) )
}


list_mcc <- function(list1, list2) {
  
  # same length
  if (length(list1) != length(list2)) {
    stop("The two lists must have the same length!")
  }

  # MCC
  mcc_values <- sapply(seq_along(list1), function(i) {
    calculate_mcc(list1[[i]], list2[[i]])
  })

  return(mcc_values)
}

##### Calute entry-wise TPR & TNR #####
calculate_tpr_tnr <- function(mat_est, mat_true) {
  
  mat_true_bin <- abs(mat_true) > 1e-6
  mat_est_bin <- abs(mat_est) > 1e-6

  true_upper <- mat_true_bin[upper.tri(mat_true_bin, diag = TRUE)]
  est_upper <- mat_est_bin[upper.tri(mat_est_bin, diag = TRUE)]

  # TP, TN, FP, FN
  tp <- sum(true_upper & est_upper) # True Positive
  fn <- sum(true_upper & !est_upper) # False Negative
  tn <- sum(!true_upper & !est_upper) # True Negative
  fp <- sum(!true_upper & est_upper) # False Positive

  # TPR, TNR
  tpr <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)  
  tnr <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)

  return(list(TPR = tpr, TNR = tnr,
              TP = tp, FN=fn,
              TN=tn, FP=fp ))
}


list_tpr_tnr <- function(list_est, list_true) {

  if (length(list_true) != length(list_est)) {
    stop("The two lists must have the same length!")
  }

  results <- lapply(seq_along(list_true), function(i) {
    calculate_tpr_tnr(list_est[[i]], list_true[[i]])
  })

  tpr_values <- mean(sapply(results, function(res) res$TPR))
  tnr_values <- mean(sapply(results, function(res) res$TNR))

  return(c(TPR = tpr_values, TNR = tnr_values))
}


##### Compute all data measure #####
# Conclude: Frobenius norm, max column-wise ell_2 norm, 
#           entry-wise tpr-tnr-mcc, covariate-wise tpr-tnr-mcc, 

my_measure = function(list_est, list_true){
  
  K = length(list_est);
  p = nrow(list_est[[1]])
  
  # 1. mean Frobenius norm 
  F_err = mean( sapply(1:K, function(i) 
                      norm(list_est[[i]]-list_true[[i]], type = "F") 
                ) )
  
  
  # 2. max column-wise ell_2 norm
  temp= matrix(0,p,p)
  
  for(k in 1:K){
    
    temp = temp + (list_est[[k]]-list_true[[k]])^2 
  }
  L2max = max( colSums(temp) ) 
  
  
  # 3. entry-wise (i.e., edge-wise) tpr-tnr 
  temp = list_tpr_tnr(list_est, list_true) 
  temp = unname(temp)
  
  entry_tpr = temp[1]
  entry_tnr = temp[2]
  
  
  # 4. entry-wise (i.e., edge-wise) mcc
  trueloc = unlist( lapply(1:K, function(k){
    abs( list_true[[k]][upper.tri(list_true[[k]])] ) >=1e-6 
  }) )
  
  estloc = unlist( lapply(1:K, function(k){
    abs( list_est[[k]][upper.tri(list_est[[k]])] ) >=1e-6 
  }) )
  
  entry_mcc = mccr(trueloc, estloc )
  
  
  # 5. covariate-wise tpr-tnr
  trueloc = estloc = matrix(0,p-1,p) # Consider nodewise regression...
  for(k in 1:K){                     # ...with the K-task learning...
    for( j in 1:p ){                 # ...w.r.t. each covariate X_j
      
      trueloc[,j] = trueloc [,j]+ list_true[[k]][-j,j]
      estloc[,j] = estloc[,j] + list_est[[k]][-j,j]
    }
  }
  
  trueloc = abs(trueloc) >=1e-6 
  estloc = abs(estloc) >=1e-6 
  table(est=estloc, real =trueloc)
  
  TP <- sum( trueloc &  estloc ) # True Positive
  FN <- sum( trueloc & !estloc) # False Negative
  TN <- sum(!trueloc & !estloc) # True Negative
  FP <- sum(!trueloc &  estloc) # False Positive
  
  cov_tpr = TP/ max(TP+FN,1)
  cov_tnr = TN/ max(TN+FP,1)
  
  
  # 6. covariate-wise mcc
  cov_mcc = mccr(trueloc, estloc )
  
  
  res = c(F_err, L2max, 
          entry_tpr, entry_tnr, entry_mcc, 
          cov_tpr, cov_tnr, cov_mcc)
  names(res) = c("Frobenius", "max_L2",
                 "entry_TPR", "entry_TNR", "entry_MCC",
                 "cov_TPR", "cov_TNR", "cov_MCC")
  
  return(res)
}
 