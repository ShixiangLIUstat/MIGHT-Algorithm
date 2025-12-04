# MIGHT-Algorithm
Multi-task Iterative Graphical Hard Thresholding 

This repository contains R scripts (in the file "R") and a demo (in the file "demo") that demonstrate how to run the MIGHT method.

## Files

* `R/generate.R` — Generate multiple precision matrices with similar structure for simulation studies.
* `R/JGL.R` — Implementations of competing joint estimation methods: BIC-based GGL, JEM, and FJEM.
* `R/JGL MIGHT.R` — Implementation of the proposed MIGHT algorithm.
* `R/separate.R` — Implementations of competing separate estimation methods: BIC-based separate Glasso and separate nodewise regression.
* `R/simuL.R` — Code for a single simulation run (one repetition).
* `R/MainCode.R` — Main script to run MIGHT simulations: runs experiments, assembles result tables, and produces figures.
* `demo/MIGHT demo.R` — A runnable demo that illustrates the MIGHT algorithm on a concrete example.
* `package/ADSIHT_0.2.1.tar` — R package of the MIGHT method
* `package/AdaL0_1.0.1.tar` — R package of separate nodewise regression

## Requirements
R >= 4.0
Required packages: JGL, mvnfast, ADSIHT, snowfall, Matrix, fasjem, simule, jointgraph, AdaL0, glasso, mccr, stringr, ggpubr, latex2exp

## Run demo
In R:
```r
# package
library(mvnfast)
library(ADSIHT)
library(Matrix)

source('./simuL.R')
set.seed(123)
n = 50; p = 10; K = 4

# generate K precision matrix with the same structure (inverse of AR(1) Topeplitz)
x_list <- lapply(1:K, function(x) rmvn(n, mu=rep(1, p),
                                       sigma = toeplitz( (x/2/K)^(1:p-1) ) ) )

# joint estimation with MIGHT method
fit = MIGHT(X=x_list, scale = 10, parallel=T, ncpus=1 )

# results
solve( toeplitz( 0.5^(1:p-1) ) )  # ground truth of the K-th graph
fit[[K]]                        # joint estimation of the K-th graph 


