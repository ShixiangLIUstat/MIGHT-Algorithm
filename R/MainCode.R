rm(list = ls())

# install.packages("remotes")
# remotes::install_github("Weilin23/jointgraph")
# install.packages("JGL_2.3.2.tar.gz", repos = NULL, type = "source")
# install.packages("ADSIHT_0.2.1.tar.gz", repos = NULL, type = "source")
# install.packages("fasjem_1.1.2.tar.gz", repos = NULL, type = "source")


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

source('./simuL.R')
# seed=123; n=100; K=10; p=100;link=0.1; gamma=1; lam=1; rate0 = 0.2; stren=0; scale=22

##### 1. Fix n=100, p=100, K=10, rate0=0.2, 0.5, 0.8 ##### 

sfInit(parallel = TRUE, cpus = detectCores() - 4)  
sfLibrary(ADSIHT)
sfLibrary(Matrix)
sfExport("JGML")    
sfExportAll()

TRY = list()
ra = c(0.2, 0.5, 0.8); MC=100

Sys.time()
for( i in 1:3){
  
  ttt = list()
  for(mc in 1:MC){
    
    ttt[[mc]] = simu_L(seed=1000*i+ 2*mc, n=100, K=10, p=100,
                         link=0.1, gamma=1, lam=1, rate0= ra[i], stren=0, scale=22 )
    cat(mc,"/",MC,", ", i, "/", 3, "\r")
  }
  
  TRY[[i]] = ttt
  Sys.time()
}

# save(TRY, file = "100Final1.RData")
sfStop()

# load("100Final1.RData")

# gettable: get the average measure:
#           TPR TNR MCC *100
#           get the truncated average error norm of FastJEM
# order:    MIGHT, GGL, JEM, FJEM, GLASSO, NODEWISE
get_table = function(R1, mysd =T, rate=0.8){
  
  MC = length(R1)
  mymean = myse = matrix(0,6,8)
  myout = mymean
  
  for( i in 1:6){
    
    for(j in 1:8){
      
      temp = rep(0,MC)
      for( k in 1:MC){
        
        temp[k] = R1[[k]][i,j]
      }
      if(j==2){temp = temp/10}
      if(j>2 ){temp = temp*100}
      if( (i==6) & (j %in% c(1,2)) ){
        q_high <- quantile(temp, rate)
        temp = temp[temp<q_high]
      }
      
      mymean[i,j] = mean(temp)
      myse[i,j] = sd(temp)/sqrt(length(temp))
      
      if(mysd == T){
        
        myout[i,j] = str_c( sprintf("%.3f", mymean[i,j]), " (",
                          sprintf("%.3f", myse[i,j]), ")" )
      }else{
        myout[i,j] = mymean[i,j]
        
      }
    }
  }
  
  reper = matrix(0,6,8)
  reper[1:4,] = myout[c(1,2,4,6),]
  reper[5:6,] = myout[c(3,5),]
  return(reper)
  # return( list(mean = mymean, se=myse))
}

ex1 = rbind(get_table(TRY[[1]]),
            get_table(TRY[[2]]),
            get_table(TRY[[3]]) )
rownames(ex1) = rep( c("MIGHT", "GGL", "JEM", "FJEM",
                       "Sep Glasso", "Sep Node"), 3 )
colnames(ex1) = c( "Frobenius Norm", "Max L2 Norm",
                   "TPR-Edge", "TNR-Edge", "MCC-Edge",
                   "TPR-Ngbr", "TNR-Ngbr", "MCC-Ngbr" )
write.csv(ex1, "Maintable1.csv", row.names = T)



##### 2. Fix n=100, p=100, K=10, signal influence ##### 

sfInit(parallel = TRUE, cpus = detectCores() - 4) 
sfLibrary(ADSIHT)
sfLibrary(Matrix)
sfExport("JGML")    
sfExportAll()

Sys.time()
MC= 50
STR= 1/seq(1, 5, l=5)-0.1
R21=list()

for(l in 1:length(STR)){
  
  temp=list()
  for(mc in 1:MC){
    
    temp[[mc]] = simu_L(seed=2000*l+20*mc, n=100, K=10, p=100,
                          link=0.1, gamma=1, rate0= 0.5, stren = STR[l], scale=22)
    cat(mc,"/",MC,", ", l, "/", length(STR), "\r")
  }
  
  R21[[l]]=temp
}
Sys.time()
sfStop()

# load("100newR21310.RData")


R100 = R21

R21= list()
for( i in 1: length(R100) ){
  
  R21[[i]] = get_table(R100[[i]], mysd =F, rate=0.8) 
}

str2 = 1/seq(1, 5, l=5)
pn = length(str2)
plotr1 = data.frame("Signal"=rep(0,24*pn),        # minimal signal strength
                    "Methods" = rep("c" ,24*pn),  # estimate method
                    "Metric" = rep("c" , 24*pn),  # metric
                    "value" = rep(0,24*pn) )      # value

for( i in 1:pn ){ # for each signal point
  
  plotr1$Metric[ (i*24-23):(i*24)] = rep( c("Frobenius Norm", "Max L2 Norm",
                                            "MCC-Edge", "MCC-Ngbr"), 6 )
  
  plotr1$Methods[(i*24-23):(i*24)] = c( rep("MIGHT", 4),
                                        rep("GGL", 4), rep("JEM", 4), rep("FJEM", 4),
                                        rep("Sep Glasso", 4), rep("Sep Node", 4) )
  
  plotr1$Signal[(i*24-23):(i*24)] = rep( 1/str2[i], 24 )
  
  vtemp = sapply(R21, function(x) t(x[,c(1, 2, 5, 8)])  )
  
  if(i>0){
    ttemp = vtemp[,i]
    ttemp[13:14] = NaN
    plotr1$value[(i*24-23):(i*24)] =  ttemp
  }else{
    plotr1$value[(i*24-23):(i*24)] =  vtemp[,i]
  }
  
}


plotr1 = mutate( plotr1, Metric = factor(Metric, levels = c("Frobenius Norm", "Max L2 Norm",
                                                            "MCC-Edge", "MCC-Ngbr") ),
                 Methods = factor(Methods, levels = c("MIGHT", "GGL", "JEM", "FJEM",
                                                      "Sep Glasso", "Sep Node" )) )

p1 = ggplot(plotr1, aes(x = Signal, y = value, group = Methods,
                        color = Methods, linetype = Methods,
                        shape = Methods)) +
  geom_line(linewidth=0.8) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(15:18, 6, 3)) +
  scale_linetype_manual(values = rep(1, 6)) +
  facet_wrap(~Metric, scale = "free", ncol = 2) + 
  guides( color = guide_legend(nrow = 1),
          linetype = guide_legend(nrow = 1),
          shape = guide_legend(nrow = 1) ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",   
    legend.key.size = unit(30, "pt"),
    legend.box.spacing = unit(0, 'pt'),
    legend.text = element_text(size = 15),
    legend.box.margin = margin(5, 100, 0, 100),
    legend.spacing.x = unit(10, "pt"),
    axis.text = element_text(size = 15),   
    axis.title = element_text(size = 15),  
    strip.text = element_text(size = 15)  
    ) +
  ylab("") + xlab("Signal strength 1/r")
p1

ggsave(file = "./MainFig2.pdf", dpi = 1200, p1, height =3.5, width = 9)



##### 3. Fix p=50, n=100, asymptotic distribution ##### 

set.seed(33)
p=50; K=10; 
link=0.1; gamma=1; lam=1; rate0 = 0.5; stren=0

#Generate and Fix support location
data <- generate_Gamma_L(p, link=link )

#Generate and Fix the precsion matrix for each dataset
Precision <- lapply(1:K, function(x) 
  generate_Omega_L(data$Gamma, lam=lam, rate0 = rate0, stren=stren) )

#Generate and Fix Covariance of each dataset
Sigma <- lapply(Precision, function(x) solve(x))

#n=100
simuasy = function(seed, Sigma, n, gamma=1, scale =15 ) {
  
  set.seed(seed)
  K = length(Sigma)
  p = nrow(Sigma[[1]])
  
  X <- lapply(Sigma, function(x) rmvn(n, rep(0, p), x))
  
  res = list()
  
  # 1. MIGHT: our result 
  result <- JGML(X, ic.coef = 4.8, ic.scale = 1.2, coef1 = 1, coef2 = 1, 
                 kappa = 0.9, center = F, scale= scale )
  res[[1]] = list( result[[1]], result[[2]] )
  
  
  # 2. GGL: Witten 2014 
  result2 = JGL_EBIC(X, penalty = "group", tol = 1e-3,
                     maxiter = 200, gamma = gamma)$res$theta 
  res[[2]] = list( result2[[1]], result2[[2]] )
  
  
  # 3. JEM: Guo 2011 
  result3 <- JointGrpah_BIC(X, gamma = gamma)
  res[[3]] = list( result3[[1]], result3[[2]] )
  
  
  # 4. Separate nodewise: Shu 2024 
  result4 <- tune_HTP_multiple(X, ic.coef = 4.8 )
  res[[4]] = list( result4[[1]], result4[[2]] )
  
  temp = list( data = list(X[[1]], X[[2]]),
               res = res )
  
  return( temp )
}

sfInit(parallel = TRUE, cpus = detectCores() - 4)  
sfLibrary(ADSIHT) 
sfLibrary(Matrix)
sfExport("JGML")     
sfExportAll()        



Sys.time()
MC= 300
Res361 = list()

for(mc in 1:MC){
  
  Res361[[mc]] = simuasy(seed=3e4 + 2*mc, Sigma = Sigma,
                          n=100, gamma=1, scale =14 )
  
  cat(mc,"/",MC,"\r")
}

Sys.time()
sfStop()

save(Res361, file = "3R361.RData")

# load("3R361.RData")

#loc: 1:the precision's location;  2,3:the entry's location
asy = function(R,loc=c(2,5,6), truetheta=T2){
  
  MC = length(R)
  res = matrix(0,MC,4)
  colnames(res) = c("MIGHT", "GGL", "JEM", "IHT")
  myrow=loc[2]; mycol=loc[3]
  
  for( mc in 1:MC){
    
    data = R[[mc]][[1]][[loc[1]]]
    n= nrow(data); p = ncol(data)
    Sigmahat = 1/n * t(data) %*% data
    
    for(i in 1:4){ # 4 methods 
      
      Thetahat = R[[mc]][[2]][[i]][[loc[1]]]
      lochat = Thetahat[myrow, mycol]
      
      Shat = which(Thetahat[, mycol] !=0)
      ttmat = matrix(0,p,p)
      ttmat[Shat, Shat] = solve(Sigmahat[Shat, Shat])
      
      vu= ttmat[myrow, myrow] * ttmat[mycol, mycol] + ttmat[myrow, mycol]^2
      
      # ve=0
      # for(nn in 1:n){
      #   ve = ve + ( sum(ttmat[myrow,] * data[nn,]) * sum( data[nn,] * ttmat[,mycol]) )^2
      # }
      # ve = ve/n - lochat^2
      # 
      # vq=0
      # for(nn in 1:n){
      #   vq = vq + ( sum(Thetahat[myrow,] * data[nn,]) * sum( data[nn,] * Thetahat[,mycol]) )^2
      # }
      # vq = vq/n - lochat^2
      
      res[mc,i] = sqrt(n)*(lochat -truetheta )/sqrt(vu)
      
    }
  }
  
  return(res)
}

tran =function(res, ml){
  
  if( sum(res[,ml] %in% c(-Inf, Inf, NaN) )>0 ){
    
    ttt = which(res[,ml] %in% c(-Inf, Inf, NaN) )
    t1= res[-ttt,ml] 
  }else{
    
    t1 = res[,ml]  
  }
  
  return(t1)
}


mc = c("#D3D3D3", 
       #"#00008B",
       "black", "#FF8C00" ) 

## fig 1
pdf(file = "./1127_simu1.pdf",  width = 10, height = 7)

par(mfrow=c(3,4),  mai=c(0.6, 0.1, 0.5, 0.1), 
    cex.main=2, 
    cex.axis=1.5 )


mename = c("MIGHT", "GGL", "JEM", "Sep Node")

LOC = matrix(c(1,10,10,
               1,20,20,
               1,30,30), 3, 3, byrow=T)
for(mloc in 1:nrow(LOC)){
  Ti = Precision[[LOC[mloc,1]]][LOC[mloc,2],LOC[mloc,3]]
  Ti = round(Ti,3)
  res = asy(R=Res361, loc=LOC[mloc,], truetheta=Ti)
  temptex = sprintf("$\\Theta_{%d,%d}^{(%d)} = %f$", 
                    LOC[mloc,2], LOC[mloc,3], LOC[mloc,1], Ti)
  for(me in 1:4){
    xnorm = tran(res,me)
    if(abs(3-max(xnorm)) + abs(-3-min(xnorm)) < 2){
      xlim=c(-3,3)
    }else{ xlim = c(min(xnorm),max(xnorm)) }
    hist(xnorm, col = mc[1], freq = F, yaxt="n",
         breaks=seq(min(xnorm),max(xnorm),l=12),
         xlab="", xlim=xlim,
         main= mename[me], ylab="", ylim=c(0,0.6), mgp=c(2, 1, 0)
    )
    abline(v= mean(xnorm), col= mc[2],lwd = 2)
    abline(v= 0, col=mc[3],lwd = 2)
    truex = seq(xlim[1], xlim[2], l=5e2)
    lines(x=truex, y=dnorm(truex,sd=1 ),lwd = 2, col = mc[3])
    if(me == 4){
      center_x <- grconvertX(0.5, from="ndc", to="user")
      mtext(TeX(temptex), side=1, line=3.5, at=center_x, xpd=NA, cex=1.2)
    }
  }
}

dev.off()


#fig2
mc = c("#D3D3D3", 
       #"#00008B",
       "black", "#FF8C00" )

pdf(file = "./1127_simu2.pdf",  width = 10, height = 7)

par(mfrow=c(3,4),  mai=c(0.6, 0.1, 0.5, 0.1), 
    cex.main=2, 
    cex.axis=1.5 )


mename = c("MIGHT", "GGL", "JEM", "Sep Node")

LOC = matrix(c(2,12,42,
               2,24,5,
               2,36,5 ), 3, 3, byrow=T)
for(mloc in 1:nrow(LOC)){
  Ti = Precision[[LOC[mloc,1]]][LOC[mloc,2],LOC[mloc,3]]
  Ti = round(Ti,3)
  res = asy(R=Res361, loc=LOC[mloc,], truetheta=Ti)
  temptex = sprintf("$\\Theta_{%d,%d}^{(%d)} = %f$", 
                    LOC[mloc,2], LOC[mloc,3], LOC[mloc,1], Ti)
  for(me in 1:4){
    xnorm = tran(res,me)
    if(abs(3-max(xnorm)) + abs(-3-min(xnorm)) < 1.2){
      xlim=c(-3,3)
    }else{ xlim = c(min(xnorm),max(xnorm)) }
    hist(xnorm, col = mc[1], freq = F, yaxt="n",
         breaks=seq(min(xnorm),max(xnorm),l=12),
         xlab="", xlim=xlim,
         main= mename[me], ylab="", ylim=c(0,0.7), mgp=c(2, 1, 0)
    )
    abline(v= mean(xnorm), col= mc[2],lwd = 2)
    abline(v= 0, col=mc[3],lwd = 2)
    truex = seq(xlim[1], xlim[2], l=5e2)
    lines(x=truex, y=dnorm(truex,sd=1 ),lwd = 2, col = mc[3])
    if(me == 4){
      center_x <- grconvertX(0.5, from="ndc", to="user")
      mtext(TeX(temptex), side=1, line=3.5, at=center_x, xpd=NA, cex=1.2)
    }
  }
}

dev.off()
