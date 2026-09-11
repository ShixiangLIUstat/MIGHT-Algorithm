
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
library(AdaL0)
library(igraph)

source('./generate.R')
source('./JGL MIGHT.R')
source('./JGL.R')
source('./separate.R')

# link    (0< <1) The probability of the main ER graph (only used in ER graph)
# gamma   (>=0) The tunning parameter in EBIC
# type    (default "ER") structure of graphs: ER, band, hub, and scale
# distr   (default "norm") distribution 
# mytime  (default False) record computational time
# lam     (DO NOT CHANGE) JUST TAKE IT AS 1
# rate0   rate of support each graph randomly drops from the main graph (used in ER and band graph)
# stren   the signal strength of variance, the bigger, the signal strength weaker
# scale   the scaling parameter used in MIGHT
#
#
simu_Ltime <- function(seed, n, K, p, link=0.05, gamma=1,
                       type="ER", distr = "norm", mytime = F,
                       lam=1, rate0=0.5, stren=0, scale =13) {
  
  set.seed(seed)
  
  #generate support location
  
  if(type %in% c("ER", "band") ){
    data <- generate_Gamma_L(p, link=link, type=type )
    
    #generate the precision matrix for each task 
    Precision <- lapply(1:K, function(x) 
      generate_Omega_L(data$Gamma, lam=lam, rate0 = rate0, stren=stren) )
  }else if(type=="hub"){
    
    # generate with randomly changed hub structure 
    Precision = generate_Omega_L_Hub(p=p, K=K, A=3, bn=40, stren=stren, seed=seed)
  }else if(type=="scale"){
    
    # generate with randomly changed PA structure 
    Precision = generate_Omega_Scale_free(p=p, K=K, A=3, bn=4, stren=stren, seed=seed) 
  }
  
  #get Covariance of each task
  Sigma <- lapply(Precision, function(x) solve(x))
  
  if(distr == "norm"){
    X <- lapply(Sigma, function(x) rmvn(n, rep(0, p), x))
  }else{ # generate t distribution
    X <- lapply(Sigma, function(x){
      # t distribution with degree distr and variance 1
      temp = matrix( rt(n*p,df=distr)/sqrt( distr/(distr-2) ), p, n ) 
      ev = eigen(x)
      Sighalf = ev$vectors %*% diag(sqrt(ev$values))  
      res = Sighalf %*% temp
      return( t(res) )
    } )
  }
  
  
  res = matrix(0,6, 9)
  colnames(res) = c("Frobenius", "max_L2",
                    "entry_TPR", "entry_TNR", "entry_MCC",
                    "cov_TPR", "cov_TNR", "cov_MCC", "time")
  rownames(res) = c("MIGHT", "GGL", "Sep Glasso",
                    "JEM", "Sep Node", "FJEM")
  
  # 1. MIGHT: our method 
  start_time = Sys.time()
  result <- JGML(X, ic.coef = 4.8, ic.scale = 1.2, coef1 = 1, coef2 = 1, 
                 kappa = 0.9, center = F, scale= scale)
  end_time = Sys.time()
  res[1,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[1,1:8] = my_measure(result, Precision)
  
  
  # 2. GGL: Witten 2014 
  start_time = Sys.time()
  result2 = JGL_EBIC(X, penalty = "group", tol = 1e-3, 
                     maxiter = 200, gamma = gamma)
  end_time = Sys.time()
  res[2,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[2,1:8]  = my_measure(result2$res$theta, Precision)
  
  
  # 3. Separate Glasso
  start_time = Sys.time()
  result3 <- tune_glasso_multiple(X, lambda_grid = 10^( seq(-3, 0, length.out = 5) ), 
                                  gamma = gamma)
  end_time = Sys.time()
  res[3,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[3,1:8] = my_measure(result3, Precision)
  
  
  # 4. JEM: Guo 2011
  start_time = Sys.time()
  result4 <- JointGrpah_BIC(X, gamma = gamma)
  end_time = Sys.time()
  res[4,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[4,1:8] = my_measure(result4, Precision)
  
  
  # 5. Separate Nodewise
  start_time = Sys.time()
  result5 <- tune_HTP_multiple(X, ic.coef = 4.8)
  end_time = Sys.time()
  res[5,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[5,1:8] = my_measure(result5, Precision)
  
  
  # 6. FJEM: Wang 2017 
  start_time = Sys.time()
  result6 = fasjem_EBIC(X)
  end_time = Sys.time()
  res[6,9] = as.numeric( difftime(end_time, start_time, units = "secs") )
  res[6,1:8] = my_measure(result6, Precision)
  
  if(mytime==F){
    return( res[,1:8] )
  }else{
    return( res )
  }
}

