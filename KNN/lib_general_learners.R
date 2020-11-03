#!/usr/bin/Rscript --vanilla

library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(kernlab)
library(e1071)
library(caret)
library(klaR)

FitKnn <- function(x, y, parallel_){
  ctrl <- trainControl(method = "cv", number = 10, allowParallel = parallel_)
  knngrid <- expand.grid(k = 1:30)
  set.seed(34535)
  knnregfit <- train(x = x, 
                     y = y, 
                     method = "knn",
                     trControl = ctrl,
                     tuneGrid = knngrid)
  knn.model <- knnregfit$finalModel
  
  return(knn.model)
}