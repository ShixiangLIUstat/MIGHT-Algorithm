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
source("demo/demo.R")
