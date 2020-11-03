#!/usr/bin/Rscript --vanilla

library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(randomForest)

FitRandomForest <- function(x, y, ntrees=500, parallel_=TRUE){
  if (parallel_){
    ncores = getDoParWorkers()
    set.seed(98364)
    rf.model <- foreach(ntree = rep(floor(ntrees/ncores), ncores),
                        .combine = combine, 
                        .multicombine = TRUE, 
                        .packages = "randomForest") %dorng% {
                randomForest(x, y, importance = TRUE, ntree = ntree)  
                }

    return(rf.model)
  } else {
    set.seed(98364)
    rf.model <- randomForest(x, y, importance = TRUE, ntree = ntrees) 
        
    return(rf.model)
  }
}