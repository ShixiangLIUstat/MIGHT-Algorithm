rm(list = ls())
# package
library(mvnfast)
library(ADSIHT)
library(Matrix)

source('./simuL.R')


library(ADSIHT)
library(mvnfast)
set.seed(123)
n = 50; p = 10; K = 4
x_list <- lapply(1:K, function(x) rmvn(n, mu=rep(1, p),
                                       sigma = toeplitz( (x/2/K)^(1:p-1) ) ) )
fit = MIGHT(X=x_list, scale = 10, parallel=T, ncpus=1 )

# results
solve( toeplitz( 0.5^(1:p-1) ) )  # ground truth of the K-th graph
fit[[K]]                        # joint estimation of the K-th graph 

