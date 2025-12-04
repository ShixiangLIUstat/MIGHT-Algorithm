
library(JGL)
library(mvnfast)
library(ADSIHT)
library(snowfall)
library(Matrix)
library(fasjem) 
library(simule)
library(jointgraph)
library(mccr)
library(stringr)
library(glasso)
library(ggpubr)
library(latex2exp)

source('./generate.R')
source('./JGL MIGHT.R')
source('./JGL.R')
source('./separate.R')

# link    (0< <1) The probability of the main ER graph
# gamma   (>=0) The tunning parameter in EBIC
# lam     (DO NOT CHANGE) JUST TAKE IT AS 1
# rate0   rate of support each graph randomly drop from the main graph
# stren   the signal strength of variance, the bigger, the singal strength weaker

simu_L <- function(seed, n, K, p, link=0.1, gamma=1,
                     lam=1, rate0=0.5, stren=0, scale =13) {
  
  set.seed(seed)

  #generate support location
  data <- generate_Gamma_L(p, link=link )
  
  #generate the precsion matrix for each task 
  Precision <- lapply(1:K, function(x) 
    generate_Omega_L(data$Gamma, lam=lam, rate0 = rate0, stren=stren) )
  
  #get Covariance of each task
  Sigma <- lapply(Precision, function(x) solve(x))
  
  X <- lapply(Sigma, function(x) rmvn(n, rep(0, p), x))
  
  res = matrix(0,6, 8)
  colnames(res) = c("Frobenius", "max_L2",
                      "entry_TPR", "entry_TNR", "entry_MCC",
                      "cov_TPR", "cov_TNR", "cov_MCC")
  rownames(res) = c("MIGHT", "GGL", "Sep Classo",
                    "JEM", "Sep Node", "FJEM")
  
  # 1. MIGHT: our method 
  result <- JGML(X, ic.coef = 4.8, ic.scale = 1.2, coef1 = 1, coef2 = 1, 
                 kappa = 0.9, center = F, scale= scale)
  res[1,] = my_measure(result, Precision)
  
  
  # 2. GGL: Witten 2014 
  result2 = JGL_EBIC(X, penalty = "group", tol = 1e-3, 
                     maxiter = 200, gamma = gamma)
  res[2,]  = my_measure(result2$res$theta, Precision)
  
  
  # 3. Separate Glasso
  result3 <- tune_glasso_multiple(X, lambda_grid = seq(ifelse (nrow(X[[1]])>ncol(X[[1]]), 1e-4, 1e-2), 
                                                       1, length.out = 50), 
                                  gamma = gamma)
  res[3,] = my_measure(result3, Precision)
  
  
  # 4. JEM: Guo 2011
  result4 <- JointGrpah_BIC(X, gamma = gamma)
  res[4,] = my_measure(result4, Precision)
  
  
  # 5. Separate Nodewise
  result5 <- tune_HTP_multiple(X, ic.coef = 4.8)
  res[5,] = my_measure(result5, Precision)
  
  
  # 6. FJEM: Wang 2017 
  result6 = fasjem_EBIC(X)
  res[6,] = my_measure(result6, Precision)
  
  
  return( res )
}
