rm(list = ls())
# install.packages("mvnfast")
# install.packages("remotes")
# remotes::install_github("Weilin23/jointgraph")

# install.packages("BH_1.90.0-1.tar.gz", repos = NULL, type = "source")
# install.packages("JGL_2.3.2.tar.gz", repos = NULL, type = "source")
# install.packages("ADSIHT_0.2.1.tar.gz", repos = NULL, type = "source")
# install.packages("fasjem_1.1.2.tar.gz", repos = NULL, type = "source")
# install.packages("RcppArmadillo_15.2.4-1.tar.gz", repos = NULL, type = "source")

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
library(latex2exp)

source('./simuL.R')

##### p=200 case #####

# consider performance of competing methods under larger p=200 
# and different (n,K) setting: (100,5), (100,10), (200,10)
# and alternative graph structures: ER; band; block hub; block scale free
# 
# Each setting follows MC times repetition

# Example
# seed=123; n=100; K=10; p=200; link=0.05; gamma=1; 
# lam=1; rate0 = 0.5; stren=0; scale=22; type="ER";
# t1= Sys.time()
# myhat = simu_L(seed=123, n=100, K=10, p=100,
#        link=0.05, gamma=1, lam=1, rate0= 0.5, stren=0, scale=17 )
# t2=Sys.time()
# t2-t1


 ##### 1. n=100; K=5 #####
  MC = 50; 
  n=100; K=5; p=200;
  scale = 17
  TRY1005 = list() 
  
  for(mc in 1:MC){
  ERtemp = simu_L(seed=mc, n=100, K=5, p=200, link=0.05, gamma=1, type="ER",
                     lam=1, rate0=0.5, stren=0, scale = scale)
  
  Bandtemp = simu_L(seed=mc, n=100, K=5, p=200, link=0.05, gamma=1, type="band",
                   lam=1, rate0=0.5, stren=0, scale = scale)
  
  Hubtemp = simu_L(seed=mc, n=100, K=5, p=200, link=0.05, gamma=1, type="hub",
                  lam=1, rate0=0.5, stren=0, scale = scale)

  PAtemp = simu_L(seed=mc, n=100, K=5, p=200, link=0.05, gamma=1, type="scale",
                          distr = "norm", lam=1, stren=0, scale = scale)

  TRY1005[[mc]] = list( ER   = ERtemp,
                        band = Bandtemp,
                        hub  = Hubtemp,
                        scalefree = PAtemp ) 
  
  save(TRY1005, file = "rev_n100K5.RData")
  }
  
  
  ##### 2. n=100; K=10 #####
  MC = 50; 
  n=100; K=10; p=200;
  
  TRY10010 = list() 
  
  for(mc in 1:MC){
  ERtemp = simu_L(seed=mc, n=100, K=10, p=200, link=0.05, gamma=1, type="ER",
                  lam=1, rate0=0.5, stren=0, scale = 17)
  
  Bandtemp = simu_L(seed=mc, n=100, K=10, p=200, link=0.05, gamma=1, type="band",
                    lam=1, rate0=0.5, stren=0, scale = 17)
  
  Hubtemp = simu_L(seed=mc, n=100, K=10, p=200, link=0.05, gamma=1, type="hub",
                   lam=1, rate0=0.5, stren=0, scale = 17)
   
  PAtemp = simu_L(seed=mc, n=100, K=10, p=200, link=0.05, gamma=1, type="scale",
                          distr = "norm", lam=1, stren=0, scale = 17)
   
  TRY10010[[mc]] = list( ER   = ERtemp,
                        band = Bandtemp,
                        hub  = Hubtemp,
                        scalefree = PAtemp ) 
  
  save(TRY10010, file = "rev_n100K10.RData")
  }



 ##### 3. n=200; K=10 #####
  MC = 50; 
  n=200; K=10; p=200;
  
  TRY20010 = list() 
  
  for(mc in 1:MC){
  ERtemp = simu_L(seed=mc, n=200, K=10, p=200, link=0.05, gamma=1, type="ER",
                  lam=1, rate0=0.5, stren=0, scale = 17)
  
  Bandtemp = simu_L(seed=mc, n=200, K=10, p=200, link=0.05, gamma=1, type="band",
                    lam=1, rate0=0.5, stren=0, scale = 17)
  
  Hubtemp = simu_L(seed=mc, n=200, K=10, p=200, link=0.05, gamma=1, type="hub",
                   lam=1, rate0=0.5, stren=0, scale = 17)
   
  PAtemp = simu_L(seed=mc, n=100, K=10, p=200, link=0.05, gamma=1, type="scale",
                          distr = "norm", lam=1, stren=0, scale = 17)
   
  TRY20010[[mc]] = list( ER   = ERtemp,
                        band = Bandtemp,
                        hub  = Hubtemp,
                        scalefree = PAtemp ) 
  
  save(TRY20010, file = "rev_n200K10.RData")
  }


 ##### 4. heavy tail t distribution #####
MC = 50
scale = 15
Heavy = list()

for( mc in 1:MC ){
  temp = list()
  
  t1= Sys.time()
  temp[[1]] = simu_L(seed=10*mc+1, n=100, K=10, p=100, type="ER", distr=3,
                     link=0.1, gamma=1, lam=1, rate0= 0.5, stren=0, scale= scale ) 
  temp[[1]]; t2=Sys.time() 
  
  t1= Sys.time()
  temp[[2]] = simu_L(seed=10*mc+2, n=100, K=10, p=100, type="ER", distr=6,
                     link=0.1, gamma=1, lam=1, rate0= 0.5, stren=0, scale= scale )
  temp[[2]]; t2=Sys.time() 
  
  t1= Sys.time()
  temp[[3]] = simu_L(seed=10*mc+3, n=100, K=10, p=100, type="ER", distr=9,
                     link=0.1, gamma=1, lam=1, rate0= 0.5, stren=0, scale= scale)
  temp[[3]]; t2=Sys.time() 

  Heavy[[mc]] = temp 
  save(Heavy, file = "rev_robust.RData")
}


