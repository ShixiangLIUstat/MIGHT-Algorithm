
##### generate_Gamma_L #####

# get symmetric element-wise sparse 0-1 matrix with diagnal entries always at 0
# generate support through binominal r.v. (Erdos-Renyi graph)
# generate around link*p*p support 

generate_Gamma_L <- function(p, link = 0.1,  max_attempts = 100) {
  
  #set.seed(123) # For reproducibility
  attempt <- 0
  Gamma <- matrix(0, nrow = p, ncol = p)
  min_nonzero_elements <- ceiling(0.5 * p) # Minimum number of nonzero elements
  
  repeat {
    
    attempt <- attempt + 1
    # Step 1: Generate symmetric Γ matrix
    Gamma <- matrix(0, nrow = p, ncol = p)
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) { 
        
        Gamma[j, i] = Gamma[i, j] = rbinom(1, 1, link) # Ensure symmetry
      }
    }
    
    # Check if Γ has at least min_nonzero_elements nonzero elements
    if (sum(Gamma != 0) >= min_nonzero_elements) break
    if (attempt >= max_attempts) {
      stop("Failed to generate a valid Γ matrix after max attempts.")
    }
  } 
  
  # Step 3: Extract non-zero elements in the upper triangular part (excluding diagonal)
  upper_triangle_nonzero <- Gamma[upper.tri(Gamma)]
  nonzero_values <- length(upper_triangle_nonzero[upper_triangle_nonzero != 0])
  
  return(list(Gamma = Gamma, nonzero_values = nonzero_values))
}


##### generate_Omega_L #####

# change part (rate0) support of Gamma to non-support 
# get the specific signal strength (lam) of Gamma
# ensure PD matrix
# lam    IT IS MEANINGLESS, DO NOT USE IT
# rate0  the rate represents the proportion of support positions removed from Gamma
# stren  the signal strength of variance, the bigger, the singal strength weaker


generate_Omega_L <- function(Gamma, lam=1, rate0 = 0.2, stren = 0){
  
  p <- nrow(Gamma)
  upper_triangle_indices <- which(upper.tri(Gamma), arr.ind = TRUE)
  
  non_zero_indices <- upper_triangle_indices[Gamma[upper_triangle_indices] != 0, ] #support index
  
  num_to_zero <- ceiling( rate0 * nrow(non_zero_indices))
  zero_indices <- sample(1:nrow(non_zero_indices), num_to_zero)
  
  for (i in 1:nrow(non_zero_indices) ) {
    
    row <- non_zero_indices[i, 1]
    col <- non_zero_indices[i, 2]
    
    if( i %in% zero_indices){
      Gamma[row, col] = Gamma[col, row] = 0
    }else{
      temp = runif(1, min=0.5, max=1) * ( 2*rbinom(1,1,0.5)-1 ) 
      Gamma[row, col] = Gamma[col, row] =  temp * lam
    }
  }
  
  eigenvalues <- eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eigenvalues)
  lambda <- max(0, -lambda_min)+0.1
  
  Omega <- Gamma + (lambda + stren) * diag(p)
  return(Omega)
}
