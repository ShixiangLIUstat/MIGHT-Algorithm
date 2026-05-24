
##### generate_Gamma_L #####

# get symmetric element-wise sparse 0-1 matrix with diagnal entries always at 0
# generate support through binominal r.v. (Erdos-Renyi graph, "ER") or band pattern ("band")
# generate around link*p*p support 

generate_Gamma_L <- function(p, link = 0.1, type= "ER" ,
                             max_attempts = 100) {
  
  #set.seed(123) # For reproducibility
  attempt <- 0 
  min_nonzero_elements <- ceiling(0.5 * p) # Minimum number of nonzero elements
  
  # Step 1: Generate symmetric Γ matrix
  repeat {
  
    attempt <- attempt + 1
    Gamma <- matrix(0, nrow = p, ncol = p)
    if(type == "ER"){   # Erdos-Renyi
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) { 
          Gamma[j, i] = Gamma[i, j] = rbinom(1, 1, link) # Ensure symmetry
        }
      }
    }else if(type=="band"){ # if band: |i-j|=1,2
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) { 
          Gamma[j, i] = Gamma[i, j] = ifelse( abs(i-j) <=2, 1, 0) # Band
        }
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
# lam     scaling parameter (IT IS MEANINGLESS, DO NOT USE IT)
# rate0   the rate represents the proportion of support positions removed from Gamma
# stren   the signal strength of variance, the bigger, the singal strength weaker


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



##### generate_Omega_L_Hub #####
# generate A number block-diagonal precision matrices with block-star connection
# change their upperleft and bottomright randomly and create most A^2 different pattern
# inspired by Ma and Michailidis 2016 JMLR
# known K first (therefore we could get K matrix randomly)
#
# p       dimension of matrix
# K       the requirement number of precision matrices
# A       the created base graphs
# bn      the number of blocks
# stren   the signal strength of variance, the bigger, the singal strength weaker
# seed    random seed
#
generate_Omega_L_Hub = function(p, K, A, bn, stren=0, seed){
  set.seed(seed)
  locp = p/bn
  
  # step 1: generate base graphs
  BaseList = list()
  for(ind in 1:A){
    BaseList[[ind]] = matrix(0, p, p)
  }
  
  for(blockind in 1:bn){    # for each block
    for(baseind in 1:A ){  # for each base graph
      center = sample(1:locp,1)
      block = matrix(0, locp, locp)
      block[center, ] = block[, center] = 
        runif(locp, min=0.5, max=1) * ( 2*rbinom(1,1,0.5)-1 )  
      block[center, center] =0
      
      myind = ( (blockind-1)*locp +1 ) : ( blockind*locp )
      
      BaseList[[baseind]][myind, myind] = block
    }
  }
  
  # step 2: generate different precision matrices
  res = list()
  for(ind in 1:K){
    temp = matrix(0, p, p)
    upleft = sample(1:A,1)
    bottomright = sample(1:A,1)
    
    temp[ 1:(p/2), 1:(p/2) ] =
      BaseList[[upleft]][ 1:(p/2), 1:(p/2) ]  # generate upperleft from one base graph
    temp[ (p/2+1):p, (p/2+1):p ] =
      BaseList[[bottomright]][ (p/2+1):p, (p/2+1):p ] # generate bottomright from another
    
    eigenvalues <- eigen(temp, symmetric = TRUE, only.values = TRUE)$values
    lambda_min <- min(eigenvalues)
    lambda <- max(0, -lambda_min)+0.1
  
    res[[ind]] =temp + (lambda + stren) * diag(p)
  }
  
  # return K precision matrices
  # thier pattern all come from A base graphs with random combination 
  return(res)
}



##### generate_Omega_Scale_free #####
# generate A number block-diagonal precision matrices with Scale_free network
# change their upperleft and bottomright randomly and create most A^2 different pattern
# inspired by Ma and Michailidis 2016 JMLR 
#
# p       dimension of matrix
# K       the requirement number of precision matrices
# A       the created base graphs
# bn      the number of blocks
# stren   the signal strength of variance, the bigger, the singal strength weaker
# seed    random seed
#
generate_Omega_Scale_free = function(p, K, A, bn=2, stren=0, seed){
  set.seed(seed)
  locp = p/bn
  
  # step 1: generate base graphs
  BaseList = list()
  for(ind in 1:A){
    BaseList[[ind]] = matrix(0, p, p)
  }
  
  for(blockind in 1:bn){    # for each block
    for(baseind in 1:A ){  # for each base graph
      
      block = matrix(0, locp, locp)
      
      # get PA graph for this block
      gtemp = sample_pa(n = locp, power = 1, m = 1, directed = FALSE)
      adj_matrix = as_adjacency_matrix(gtemp, sparse = FALSE)
      # apply(adj_matrix,2,sum)
      
      for( i in 1:(locp-1)){
        for(j in (i+1): locp ){
          if(adj_matrix[i,j]!= 0){
            block[i,j] = block[j,i] =
              runif(1, min=0.5, max=1) * ( 2*rbinom(1,1,0.5)-1 ) 
          } 
        }
      }
       
      myind = ( (blockind-1)*locp +1 ) : ( blockind*locp )
      BaseList[[baseind]][myind, myind] = block
    }
  }
  
  # step 2: generate different precision matrices
  res = list()
  for(ind in 1:K){
    temp = matrix(0, p, p)
    upleft = sample(1:A,1)
    bottomright = sample(1:A,1)
    
    temp[ 1:(p/2), 1:(p/2) ] =
      BaseList[[upleft]][ 1:(p/2), 1:(p/2) ]  # generate upperleft from one base graph
    temp[ (p/2+1):p, (p/2+1):p ] =
      BaseList[[bottomright]][ (p/2+1):p, (p/2+1):p ] # generate bottomright from another
    
    eigenvalues <- eigen(temp, symmetric = TRUE, only.values = TRUE)$values
    lambda_min <- min(eigenvalues)
    lambda <- max(0, -lambda_min)+0.1
    
    res[[ind]] =temp + (lambda + stren) * diag(p)
  }
  
  # return K precision matrices
  # thier pattern all come from A base graphs with random combination 
  return(res)
}
